/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
 * See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <cstring>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <unordered_map>
#include <functional>
#include <assert.h>

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_ACCESS
#include <unistd.h>  /* access() */
#endif

#include "ph5_concat.hpp"

/*----< Concatenator() >-----------------------------------------------------*/
Concatenator::Concatenator(int           nprocs,
                           int           rank,
                           MPI_Comm      comm,
                           MPI_Info      info,
                           size_t        num_input_files,
                           string const& output,
                           bool          append_mode,
                           bool          posix_open,
                           bool          in_memory_io,
                           bool          chunk_caching,
                           size_t        compress_threshold,
                           bool          one_process_create,
                           unsigned int  zip_level,
                           bool          enforce_contiguous,
                           size_t        buffer_size,
                           int           io_strategy,
                           string const& part_key_base) :
    comm(comm),
    info(info),
    nprocs(nprocs),
    rank(rank),
    num_input_files(num_input_files),
    io_buffer_size(buffer_size),
    zip(zip_level),
    enforce_contiguous(enforce_contiguous),
    posix_open(posix_open),
    in_memory_io(in_memory_io),
    chunk_caching(chunk_caching),
    compress_threshold(compress_threshold * 1048576),
    one_process_create(one_process_create),
    io_strategy(io_strategy),
    output_file_name(output),
    part_key_base(part_key_base)
{
    // chunk_size_threshold = 1*1024*1024; // Chunk size threshold (1 MiB or 1 M elements)
    chunk_size_threshold = 256*1024; // Chunk size threshold (256K elements)
    in_memory_cache_size = 512*1024*1024ull; // In-memory buffer increase (512 MiB)
    output_meta_cache_size = 128*1024*1024; // metadata cache size (128 MiB)
    raw_chunk_cache_size = 64*1024*1024; // raw chunk cache size (64 MiB)
    max_local_size_in_bytes = 0;

    original_num_groups = 0;
    num_groups_have_key = 0;
    num_groups = 0;

    original_total_num_datasets = 0;
    total_num_datasets = 0;

    /* whether to create a partition key dataset in each group */
    if (part_key_base.compare("") == 0)
        add_partition_key = false;
    else
        add_partition_key = true;

    /* spill's group index in groups[] array */
    spill_grp_no = -1;
    /* hash lookup tables, one per input file, built based on the user
     * indicated key base dataset in group /spill
     */
    lookup_table = new table[num_input_files];

    async_reqs = new MPI_Request[nprocs * 2];
    async_statuses = new MPI_Status[nprocs * 2];
    buffer = NULL;

    c_1d_2d = 0;
    o_1d = 0;
    r_1d = 0;
    w_1d = 0;
    o_2d = 0;
    r_2d = 0;
    w_2d = 0;
    o_f = 0;
    close_in_dsets = 0;
    close_out_dsets = 0;
    num_allreduce = 0;
    num_exscan = 0;
}

/*----< ~Concatenator() >----------------------------------------------------*/
Concatenator::~Concatenator()
{
    size_t ii, jj;

    for (ii=0; ii<num_groups; ii++) {
        for (jj=0; jj<groups[ii].num_dsets; jj++) {
            groups[ii].dsets[jj].in_dset_ids.clear();
            groups[ii].dsets[jj].in_dim0.clear();
        }
        delete [] groups[ii].dsets;
    }
    delete [] groups;
    delete [] lookup_table;

    delete []async_reqs;
    delete []async_statuses;

    if (buffer != NULL)
        delete []buffer;
}

/*----< set_rawdata_cache() >------------------------------------------------*/
/* This subroutine increases the raw data chunk cache size in hope to improve
 * decompression time. However, experiments show it is not as effective as
 * doing in-memory I/O. However, we observed no significant difference in
 * performance for reading 2579 small-sized datasets.
 *
 * Note from HDF5 document for H5Pset_cache(), it says the followings.
 * Raw dataset chunk caching is not currently supported when using the MPI I/O
 * and MPI POSIX file drivers in read/write mode; see H5Pset_fapl_mpio and
 * H5Pset_fapl_mpiposix, respectively. When using one of these file drivers,
 * all calls to H5Dread and H5Dwrite will access the disk directly, and
 * H5Pset_cache will have no effect on performance.
 */
static
int set_rawdata_cache(hid_t  fapl_id,
                      size_t rdcc_nslots,  /* Number of slots in hash table */
                      size_t rdcc_nbytes,  /* Size of chunk cache in bytes */
                      double w0)
{
    int err, err_exit=0;

    /* set the raw data chunk cache to improve decompression time */
    int mdc_nelmts;      /* Dummy parameter in API, no longer used by HDF5 */
    err = H5Pget_cache(fapl_id, &mdc_nelmts, &rdcc_nslots, &rdcc_nbytes, &w0);
    if (err < 0) HANDLE_ERROR("H5Pget_cache")

    /* increasing cache size for write seems no effect */

    /* rdcc_nslots should be a prime number and approximately 100 times number
     * of chunks that can fit in rdcc_nbytes. The default value used by HDF5 is
     * 521. Others can be 10007, 67231. However, experiments show that
     * increasing rdcc_nslots for this read operation actually performs worse.
     */

    err = H5Pset_cache(fapl_id, mdc_nelmts, rdcc_nslots, rdcc_nbytes, w0);
    if (err < 0) HANDLE_ERROR("H5Pset_cache")
fn_exit:
    return err_exit;
}

/*----< op_func() >----------------------------------------------------------*/
/* Operator function to be called by H5Ovisit. */
static
herr_t op_func(hid_t obj, const char *name, const H5O_info_t *info, void *me)
{
    return reinterpret_cast<Concatenator *>(me)->collect_metadata(obj, name, info);
}

/*----< construct_metadata() >-----------------------------------------------*/
/* Collect all metadata of groups and datasets from the assigned input files
 * and call MPI functions to calculated the aggregated sizes.
 */
int Concatenator::construct_metadata(vector<string> const &inputs)
{
    int err_exit=0;
    herr_t err;
    hid_t file_id, fapl_id;

    /* Go through all assigned input files and collect all metadata. */
    if (posix_open == true) {
        /* use POSIX I/O mode to open and read each file */
        if (in_memory_io == true) {
            fapl_id = H5Pcreate(H5P_FILE_ACCESS);
            if (fapl_id < 0) HANDLE_ERROR("H5Pcreate")

            /* Enabling HDF5 in-memory I/O results in faster reads. When
             * opening an existing file and in-memory-I/O is on, HDF5 reads the
             * entire file into an internal buffer, no matter how big is the
             * file. Thus, argument 'increment' of H5Pset_fapl_core take no
             * effect on reading an existing file.
             */
            err = H5Pset_fapl_core(fapl_id, 0, 0);
            if (err < 0) HANDLE_ERROR("H5Pset_fapl_core")

            if (chunk_caching) {
                err = set_rawdata_cache(fapl_id, 521, raw_chunk_cache_size, 1.0);
                if (err < 0) HANDLE_ERROR("set_rawdata_cache")
            }
        } else {
            fapl_id = H5P_DEFAULT;
        }
    } else {
        /* Use MPI-IO mode to open files. Note HDF5 in-memory I/O feature is
         * not available when using MPI-IO driver
         */
        fapl_id = H5Pcreate(H5P_FILE_ACCESS);
        if (fapl_id < 0) HANDLE_ERROR("H5Pcreate")

        /* use MPI_COMM_SELF to open files assigned */
        err = H5Pset_fapl_mpio(fapl_id, MPI_COMM_SELF, info);
        if (err < 0) HANDLE_ERROR("H5Pset_fapl_mpio")

        /* Sets metadata I/O mode for read operations to independent */
        err = H5Pset_all_coll_metadata_ops(fapl_id, false);
        if (err < 0) HANDLE_ERROR("H5Pset_all_coll_metadata_ops")
    }

    /* iterate all input files assigned to this process */
    for (auto it = inputs.begin(); it != inputs.end(); it++) {
        /* open the assigned file in read-only mode */
        file_id = H5Fopen(it->c_str(), H5F_ACC_RDONLY, fapl_id);
        if (file_id < 0) HANDLE_ERROR(string("H5Fopen ") + *it)

        /* iterate all data objects (groups and datasets in each group) to
         * collect their metadata
         */
#if defined HAS_H5OVISIT3 && HAS_H5OVISIT3
        err = H5Ovisit3(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, this, H5O_INFO_ALL);
        if (err < 0) HANDLE_ERROR(string("H5Ovisit3") + *it)
#else
        err = H5Ovisit(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, this);
        if (err < 0) HANDLE_ERROR(string("H5Ovisit") + *it)
#endif

        /* if partition key dataset is to be added and group /spill cannot be
         * found, then error out.
         */
        if (add_partition_key && spill_grp_no == -1) {
            fprintf(stderr, "Error: group /spill cannot be found in file %s\n",
                    it->c_str());
            err_exit = -1;
            goto fn_exit;
        }

        /* Instead of closing the input files, we keep the file handles
         * and re-use them later when concatenating datasets.
         */
        input_files.insert(make_pair(*it, file_id));
    }

    if (fapl_id != H5P_DEFAULT) {
        err = H5Pclose(fapl_id);
        if (err < 0) HANDLE_ERROR("H5Pclose")
    }

    /* check if partitioning key base dataset in each group is found */
    if (add_partition_key) {
        for (size_t ii=0; ii<num_groups; ii++) {
            if (groups[ii].key_base == NULL && rank == 0) {
                printf("[%d] Warning: partition key base '%s' cannot be found in group '%s'\n",
                       rank, part_key_base.c_str(), groups[ii].name.c_str());
            }
        }
        /* error out if the partition key base dataset is not found in group
         * /spill
         */
        if (groups[spill_grp_no].key_base == NULL) {
            fprintf(stderr, "Error: partiition key base dataset '/spill/%s' cannot be found\n", part_key_base.c_str());
            err_exit = -1;
            goto fn_exit;
        }
    }

    /* Sum the most significant dimension for all datasets. */
    accumulate_dimensions();

    /* Allocate I/O buffer space for read and write. */
    if (io_buffer_size == 0) io_buffer_size = max_local_size_in_bytes;

    buffer = new char[io_buffer_size];

fn_exit:
    return err_exit;
}

/*----< set_metadata_cache() >-----------------------------------------------*/
static
int set_metadata_cache(hid_t  file_id,
                       size_t cache_size,
                       double min_clean_fraction,
                       double dirty_bytes_threshold)
{
    int err, err_exit=0;
    H5AC_cache_config_t config;

    /* Set the metadata cache size which may improve data object open time */
    config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    err = H5Fget_mdc_config(file_id, &config);
    if (err < 0) HANDLE_ERROR("H5Fget_mdc_config")

    config.max_size = cache_size;
    config.min_size = cache_size;
    config.initial_size = cache_size;
    config.min_clean_fraction = min_clean_fraction;
    config.dirty_bytes_threshold = dirty_bytes_threshold * cache_size * min_clean_fraction;
    config.decr_mode = H5C_decr__off;
    err = H5Fset_mdc_config(file_id, &config);
    if (err < 0) HANDLE_ERROR("H5Fset_mdc_config")
fn_exit:
    return err_exit;
}

/*----< create_partition_key() >---------------------------------------------*/
/* Create a new partition key dataset in a group by inheriting most of the
 * metadata of the base dataset.
 */
int Concatenator::create_partition_key(GrpInfo &grp)
{
    DSInfo_t &seq = grp.dsets[grp.num_dsets++];

    /* copy contents of base dataset over to key dataset, seq */
    seq = *grp.key_base;

    seq.is_key_base = false;
    seq.is_key_seq  = true;
    seq.name        = part_key_base + ".seq";
    seq.type_id     = H5T_STD_I64LE;
    seq.type_size   = H5Tget_size(H5T_STD_I64LE);
    if (seq.global_dims[0] == 0)
        seq.layout = H5D_COMPACT;
    else
        seq.layout = (enforce_contiguous == true) ? H5D_CONTIGUOUS : H5D_CHUNKED;

    herr_t err = create_dataset(grp.id, seq, false);
    if (err < 0) RETURN_ERROR("create_dataset", seq.name.c_str())
    grp.seq_dset = &seq;
    total_num_datasets++;

    return 0;
}

/*----< file_create() >------------------------------------------------------*/
int Concatenator::file_create()
{
    int err_exit=0, file_exist=1;
    size_t ii, jj, kk;
    herr_t err;
    hid_t group_id=0, fapl_id;

#if defined PROFILE_HDF5_INTERNAL && PROFILE_HDF5_INTERNAL
    MPI_Barrier(comm);
    H5Vreset();
    H5Venable();
#endif

    /* Create the output file using MPI-IO driver */
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id < 0) HANDLE_ERROR("H5Pcreate")

    err = H5Pset_fapl_mpio(fapl_id, comm, info);
    if (err < 0) HANDLE_ERROR("H5Pset_fapl_mpio")

    /* set collective mode for metadata reads and writes */
    err = H5Pset_all_coll_metadata_ops(fapl_id, true);
    if (err < 0) HANDLE_ERROR("H5Pset_all_coll_metadata_ops")

    err = H5Pset_coll_metadata_write(fapl_id, true);
    if (err < 0) HANDLE_ERROR("H5Pset_coll_metadata_write")

    /* use a large metadata block size */
    err = H5Pset_meta_block_size(fapl_id, 4194304);
    if (err < 0) HANDLE_ERROR("H5Pset_meta_block_size")

#ifdef HAVE_ACCESS
    /* if access() is available, use it to check whether file already exists
     * rank 0 calls access() and broadcasts file_exist */
    if (rank == 0 && access(output_file_name.c_str(), F_OK) == -1)
        file_exist = 0;

    MPI_Bcast(&file_exist, 1, MPI_INT, 0, comm);
    if (file_exist) {
        cout<<output_file_name.c_str()<<" already exists." <<endl;
        return -1;
    }
#endif

    if (one_process_create == true) {
        /* Root process creates the file first, followed by all other processes
         * open the file
         */
        if (rank == 0) {
            hid_t one_fapl_id = H5Pcreate(H5P_FILE_ACCESS);
            if (one_fapl_id < 0) HANDLE_ERROR("H5Pcreate")

            /* use a large metadata block size */
            err = H5Pset_meta_block_size(one_fapl_id, 4194304);
            if (err < 0) HANDLE_ERROR("H5Pset_meta_block_size")

            /* create a new file */
            output_file_id = H5Fcreate(output_file_name.c_str(), H5F_ACC_EXCL,
                                       H5P_DEFAULT, one_fapl_id);
            if (output_file_id < 0)
                HANDLE_ERROR(string("H5Fcreate in exclusive mode ") + output_file_name.c_str())

            err = H5Pclose(one_fapl_id);
            if (err < 0) HANDLE_ERROR("H5Pclose")

            /* rank 0 creates all the groups and datasets */
            for (ii=0; ii<num_groups; ii++) {
                /* create a new group */
                group_id = H5Gcreate2(output_file_id, groups[ii].name.c_str(),
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (group_id < 0) HANDLE_ERROR("H5Gcreate2")
                groups[ii].id = group_id;

                for (jj=0; jj<groups[ii].num_dsets; jj++) {
                    /* Create a new dataset in the output file */
                    bool toFill = (groups[ii].dsets[jj].type_class == H5T_STRING);
                    err = create_dataset(group_id, groups[ii].dsets[jj], toFill);
                    if (err < 0) HANDLE_ERROR("create_dataset")

                    /* zero-sized datasets have been closed in create_dataset().
                     * Now, close non-zero sized datasets.
                     */
                    if (groups[ii].dsets[jj].global_dims[0] > 0) {
                        err = H5Dclose(groups[ii].dsets[jj].out_dset_id);
                        if (err < 0) HANDLE_ERROR("H5Dclose")
                    }
                }
                err = H5Gclose(group_id);
                if (err < 0) HANDLE_ERROR("H5Gclose")
            }

            /* must flush as the file may still be cached locally */
            err = H5Fflush(output_file_id, H5F_SCOPE_LOCAL);
            if (err < 0) HANDLE_ERROR("H5Fflush")

            err = H5Fclose(output_file_id);
            if (err < 0) HANDLE_ERROR("H5Fclose")
        }

        /* all processes wait until root completes the file create */
        MPI_Barrier(comm);

        /* all ranks collectively open the newly created file */
        output_file_id = H5Fopen(output_file_name.c_str(), H5F_ACC_RDWR,
                                 fapl_id);
        if (output_file_id < 0)
            HANDLE_ERROR(string("H5Fopen ") + output_file_name)

        err = set_metadata_cache(output_file_id, output_meta_cache_size,
                                 0.3, 0.45);
        if (err < 0) HANDLE_ERROR("set_metadata_cache")

        /* all ranks collectively open all non-zero datasets */
        for (ii=0; ii<num_groups; ii++) {
            /* create a new group */
            group_id = H5Gopen(output_file_id, groups[ii].name.c_str(),
                               H5P_DEFAULT);
            if (group_id < 0) HANDLE_ERROR("H5Gopen")
            groups[ii].id = group_id;

            for (jj=0; jj<groups[ii].num_dsets; jj++) {
                /* open only globally non-zero datasets */
                if (groups[ii].dsets[jj].global_dims[0] > 0) {
                    hid_t id;
                    id = H5Dopen(group_id, groups[ii].dsets[jj].name.c_str(),
                                 H5P_DEFAULT);
                    if (id < 0) HANDLE_ERROR(string("H5Dopen ")+groups[ii].dsets[jj].name)
                    groups[ii].dsets[jj].out_dset_id = id;
                    calculate_chunk_size(groups[ii].dsets[jj]);
                }
            }

            /* After all the datasets have been created, create a new partition
             * key dataset in each group. Note the partition key base dataset
             * may be missing in a group. In this case, the key seq dataset
             * will not be created for that group.
             */
            if (add_partition_key && groups[ii].key_base != NULL) {
                err = create_partition_key(groups[ii]);
                if (err < 0) HANDLE_ERROR("create_partition_key")
            }

            err = H5Gclose(group_id);
            if (err < 0) HANDLE_ERROR("H5Gclose")
        }
    }
    else { /* all processes collectively create file and datasets */
        if (chunk_caching) {
            err = set_rawdata_cache(fapl_id, 67231, raw_chunk_cache_size, 1.0);
            if (err < 0) HANDLE_ERROR("set_rawdata_cache")
        }

        /* collectively create a new output file */
        output_file_id = H5Fcreate(output_file_name.c_str(), H5F_ACC_EXCL,
                                   H5P_DEFAULT, fapl_id);
        if (output_file_id < 0)
            HANDLE_ERROR(string("H5Fcreate in exclusive mode ") + output_file_name.c_str())

        err = set_metadata_cache(output_file_id, output_meta_cache_size,
                                 0.3, 0.45);
        if (err < 0) HANDLE_ERROR("set_metadata_cache")

        /* Collectively create all the groups and datasets */
        for (ii=0; ii<num_groups; ii++) {
            group_id = H5Gcreate2(output_file_id, groups[ii].name.c_str(),
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (group_id < 0) HANDLE_ERROR("H5Gcreate2")
            groups[ii].id = group_id;

            for (jj=0; jj<groups[ii].num_dsets; jj++) {
                /* Create a new dataset in the output file. Keep the handle
                 * opened so that it can be used later when concatenating
                 * datasets.
                 */
                bool toFill = (groups[ii].dsets[jj].type_class == H5T_STRING);
                err = create_dataset(group_id, groups[ii].dsets[jj], toFill);
                if (err < 0) HANDLE_ERROR("create_dataset")
            }

            /* After all the datasets have been created, create a new partition
             * key dataset in each group. Note the partition key base dataset
             * may be missing in a group. In this case, the key seq dataset
             * will not be created for that group.
             */
            if (add_partition_key && groups[ii].key_base != NULL) {
                err = create_partition_key(groups[ii]);
                if (err < 0) HANDLE_ERROR("create_partition_key")
            }

            err = H5Gclose(group_id);
            if (err < 0) HANDLE_ERROR("H5Gclose")
        }
    }
    err = H5Pclose(fapl_id);
    if (err < 0) HANDLE_ERROR("H5Pclose")

    /* save the original total number of datasets */
    original_total_num_datasets = total_num_datasets;

    /* remove zero-sized datasets from each group from groups[].dsets[] */
    for (ii=0; ii<num_groups; ii++) {
        kk = -1;
        for (jj=0; jj<groups[ii].num_dsets; jj++) {
            if (groups[ii].dsets[jj].global_dims[0] == 0) {
                groups[ii].dsets[jj].in_dset_ids.clear();
                total_num_datasets--;
            }
            else if (++kk < jj) {
                groups[ii].dsets[kk] = groups[ii].dsets[jj];
                if (groups[ii].dsets[kk].is_key_base)
                    groups[ii].key_base = groups[ii].dsets + kk;
                else if (groups[ii].dsets[kk].is_key_seq)
                    groups[ii].seq_dset = groups[ii].dsets + kk;
            }
        }
        groups[ii].num_dsets = kk + 1;
    }
    /* remove empty groups from groups[] */
    num_groups_have_key = 0;
    kk = -1;
    for (ii=0; ii<num_groups; ii++) {
        if (groups[ii].num_dsets == 0) {
            delete [] groups[ii].dsets;
            continue;
        }
        else if (++kk < ii)
            groups[kk] = groups[ii];

        if (groups[kk].key_base != NULL)
            num_groups_have_key++;
    }
    num_groups = kk + 1;

#if defined DEBUG && DEBUG
    if (add_partition_key) {
        for (ii=0; ii<num_groups; ii++) {
            if (groups[ii].key_base != NULL) {
                assert(groups[ii].seq_dset->is_key_seq);
            }
        }
    }
#endif

    /* Calculate the number of collective writes for each dataset based on the
     * I/O buffer size.
     */
    calculate_num_writes();

    /* use MPI collective I/O */
    dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    if (dxpl_id < 0) HANDLE_ERROR("H5Pcreate")
    err = H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    if (err < 0) HANDLE_ERROR("H5Pset_dxpl_mpio")

#if defined PROFILE_HDF5_INTERNAL && PROFILE_HDF5_INTERNAL
    MPI_Barrier(comm);
    H5Vprint();
    H5Vdisable();
#endif
fn_exit:
    return err_exit;
}

/*----< file_open() >--------------------------------------------------------*/
int Concatenator::file_open()
{
    int err_exit=0, file_exist=1;
    size_t ii, jj, kk;
    herr_t err;
    hid_t group_id=0, fapl_id;

#if defined PROFILE_HDF5_INTERNAL && PROFILE_HDF5_INTERNAL
    MPI_Barrier(comm);
    H5Vreset();
    H5Venable();
#endif

#ifdef HAVE_ACCESS
    /* if access() is available, use it to check whether file already exists
     * rank 0 calls access() and broadcasts file_exist */
    if (rank == 0 && access(output_file_name.c_str(), F_OK) == -1)
        file_exist = 0;

    MPI_Bcast(&file_exist, 1, MPI_INT, 0, comm);
    if (!file_exist) {
        cout<<output_file_name.c_str()<<" does not exists." <<endl;
        return -1;
    }
#endif

    original_total_num_datasets = 0;
    total_num_datasets = 0;

    /* remove zero-sized datasets from each group from groups[].dsets[] */
    for (ii=0; ii<num_groups; ii++) {
        /* save the original total number of datasets */
        original_total_num_datasets += groups[ii].num_dsets;
        kk = -1;
        for (jj=0; jj<groups[ii].num_dsets; jj++) {
            if (groups[ii].dsets[jj].global_dims[0] == 0) {
                groups[ii].dsets[jj].in_dset_ids.clear();
            }
            else if (++kk < jj) {
                groups[ii].dsets[kk] = groups[ii].dsets[jj];
                if (groups[ii].dsets[kk].is_key_base)
                    groups[ii].key_base = groups[ii].dsets + kk;
                else if (groups[ii].dsets[kk].is_key_seq)
                    groups[ii].seq_dset = groups[ii].dsets + kk;
            }
        }
        groups[ii].num_dsets = kk + 1;
        total_num_datasets += groups[ii].num_dsets;
    }
    /* remove empty groups from groups[] */
    num_groups_have_key = 0;
    kk = -1;
    for (ii=0; ii<num_groups; ii++) {
        if (groups[ii].num_dsets == 0) {
            delete [] groups[ii].dsets;
            continue;
        }
        else if (++kk < ii)
            groups[kk] = groups[ii];

        if (groups[kk].key_base != NULL)
            num_groups_have_key++;
    }
    num_groups = kk + 1;

    /* Open the output file using MPI-IO driver */
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id < 0) HANDLE_ERROR("H5Pcreate")

    err = H5Pset_fapl_mpio(fapl_id, comm, info);
    if (err < 0) HANDLE_ERROR("H5Pset_fapl_mpio")

    /* set collective mode for metadata reads and writes */
    err = H5Pset_all_coll_metadata_ops(fapl_id, true);
    if (err < 0) HANDLE_ERROR("H5Pset_all_coll_metadata_ops")

    err = H5Pset_coll_metadata_write(fapl_id, true);
    if (err < 0) HANDLE_ERROR("H5Pset_coll_metadata_write")

    if (chunk_caching) {
        err = set_rawdata_cache(fapl_id, 67231, raw_chunk_cache_size, 1.0);
        if (err < 0) HANDLE_ERROR("set_rawdata_cache")
    }

    /* collectively open the output file */
    output_file_id = H5Fopen(output_file_name.c_str(), H5F_ACC_RDWR, fapl_id);
    if (output_file_id < 0)
        HANDLE_ERROR(string("H5Fopen ") + output_file_name)

    err = set_metadata_cache(output_file_id, output_meta_cache_size, 0.3, 0.45);
    if (err < 0) HANDLE_ERROR("set_metadata_cache")

    /* Collectively open all the groups and datasets */
    for (ii=0; ii<num_groups; ii++) {

        group_id = H5Gopen(output_file_id, groups[ii].name.c_str(), H5P_DEFAULT);
        if (group_id < 0) HANDLE_ERROR("H5Gopen")
        groups[ii].id = group_id;

        for (jj=0; jj<groups[ii].num_dsets; jj++) {
            /* Open dataset jj in group ii. Keep the dataset opened, so that it
             * can be used later when appending.
             */
            err = open_dataset(group_id, groups[ii].dsets[jj]);
            if (err < 0) HANDLE_ERROR("open_dataset")
        }

        /* check if this group contains a key dataset */
        if (add_partition_key && groups[ii].key_base != NULL) {
            /* key dataset is the last one in dsets[] */

            DSInfo_t &seq = groups[ii].dsets[groups[ii].num_dsets];
            /* copy contents of key base dataset over to key dataset seq */
            seq = *groups[ii].key_base;
            seq.name = part_key_base + ".seq";

            /* zero-sized dataset has been removed earlier */
            assert(seq.global_dims[0] > 0);

	    /* Open key dataset. In append mode, the partition key dataset must
             * have already existed
             */
	    seq.out_dset_id = H5Dopen(group_id, seq.name.c_str(), H5P_DEFAULT);
            if (seq.out_dset_id < 0)
                HANDLE_ERROR(string("H5Dopen on partition key dataset ")+seq.name)
            /* seq.global_dims has been increased in key_base.global_dims */
            err = H5Dset_extent(seq.out_dset_id, seq.global_dims);
            if (err < 0) RETURN_ERROR("H5Dset_extent",seq.name.c_str())

            seq.is_key_base = false;
            seq.is_key_seq  = true;
            seq.type_id     = H5T_STD_I64LE;
            seq.type_size   = H5Tget_size(H5T_STD_I64LE);
            groups[ii].seq_dset = &seq;
            groups[ii].num_dsets++;
            total_num_datasets++;
        }

        err = H5Gclose(group_id);
        if (err < 0) HANDLE_ERROR("H5Gclose")
    }
    err = H5Pclose(fapl_id);
    if (err < 0) HANDLE_ERROR("H5Pclose")

    /* Calculate the number of collective writes for each dataset based on the
     * I/O buffer size.
     */
    calculate_num_writes();

    /* use MPI collective I/O */
    dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    if (dxpl_id < 0) HANDLE_ERROR("H5Pcreate")
    err = H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
    if (err < 0) HANDLE_ERROR("H5Pset_dxpl_mpio")

#if defined PROFILE_HDF5_INTERNAL && PROFILE_HDF5_INTERNAL
    MPI_Barrier(comm);
    H5Vprint();
    H5Vdisable();
#endif
fn_exit:
    return err_exit;
}

/*----< collect_metadata() >-------------------------------------------------*/
int Concatenator::collect_metadata(hid_t             obj_id,
                                   char const       *name,
                                   H5O_info_t const *obj_info)
{
    static int grp_no=-1, file_no=-1;
    herr_t err;
    unordered_map<string, DSInfo_t> grp;

    if (obj_info->type == H5O_TYPE_GROUP) {
        /* This object is a group. Retrieve metadata about this group */
        H5G_info_t grp_info;
        err = H5Gget_info_by_name(obj_id, name, &grp_info, H5P_DEFAULT);
        if (err < 0) RETURN_ERROR("H5Gget_info_by_name", name)

        /* Root group is '.' also first object read by a new input file */
        if (!strcmp(name, ".")) {
            if (file_no == -1) { /* set only when reading the 1st file */
                num_groups = grp_info.nlinks;
                original_num_groups = num_groups;
                groups = new GrpInfo[num_groups];
            }
            file_no++;   /* number of input file */
            grp_no = -1; /* restart from -1 for a new input file */
            return 0;
        }
        else grp_no++;

        if (file_no == 0) { /* set only when reading the 1st file */
            groups[grp_no].name.assign(name);
            /* allocate space for dataset metadata in this group.
             * +1 is for an additional partitioning key dataset */
            groups[grp_no].dsets = new DSInfo_t[grp_info.nlinks + 1];
            groups[grp_no].key_base = NULL;
            groups[grp_no].seq_dset = NULL;

            /* save the groups[] array index of group /spill */
            if (!strcmp(name, "spill")) spill_grp_no = grp_no;
        }
        groups[grp_no].num_dsets = 0; /* to be increased in later visit */
    }
    else if (obj_info->type == H5O_TYPE_DATASET) {
        /* This object is a dataset */
        DSInfo_t &dset = groups[grp_no].dsets[groups[grp_no].num_dsets];
        groups[grp_no].num_dsets++;

        hid_t dset_id;

        /* retrieve group name and dataset name */
        char *name_copy = strdup(name);
        char *group_name = strtok(name_copy, "/");
        char *dataset_name = strtok(NULL, "/");

        /* group should have already been registered, because groups are
         * visited first */
        if (groups[grp_no].name.compare(group_name))  {
            printf("Expecting group name %s of dataset %s but got %s\n",
                   groups[grp_no].name.c_str(), dataset_name, group_name);
            return -1;
        }

        /* Open the dataset. Note that obj_id is not the dataset/group ID. */
        dset_id = H5Dopen(obj_id, name, H5P_DEFAULT);
        if (dset_id < 0) RETURN_ERROR("H5Dopen",name)

        /* obtain dataset's dimension sizes into dset_dims[] */
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) RETURN_ERROR("H5Dget_space",name)
        hsize_t dset_dims[2];
        int ndims = H5Sget_simple_extent_dims(space_id, dset_dims, NULL);
        if (ndims < 0) RETURN_ERROR("H5Sget_simple_extent_dims",name)
        err = H5Sclose(space_id);
        if (err < 0) RETURN_ERROR("H5Sclose",name)

        if (file_no == 0) { /* This is the first input file */
            dset.name.assign(dataset_name);
            dset.local_dims[0] = dset_dims[0];
            dset.local_dims[1] = dset_dims[1];
            dset.type_id = H5Dget_type(dset_id);
            if (dset.type_id < 0) RETURN_ERROR("H5Dget_type",name)
            dset.type_class = H5Tget_class(dset.type_id);
            if (dset.type_class < 0) RETURN_ERROR("H5Tget_class", name);

            if (dset.type_class == H5T_STRING) {
                /* special treatment for string datasets of variable length */
                dset.type_id = H5Tcopy(H5T_C_S1);
                if (dset.type_id < 0) RETURN_ERROR("H5Tcopy",name)
                err = H5Tset_size(dset.type_id, MAX_STR_LEN);
                if (err < 0) RETURN_ERROR("H5Tset_size",name)
                dset.type_size = MAX_STR_LEN;
            }
            else {
                dset.type_size = H5Tget_size(dset.type_id);
                if (dset.type_size  < 0) RETURN_ERROR("H5Tget_size",name)
            }
            dset.layout = (enforce_contiguous == true) ? H5D_CONTIGUOUS
                                                       : H5D_CHUNKED;
            dset.is_key_seq  = false;
            dset.is_key_base = false;

            if (add_partition_key) {
                string key_name = part_key_base + ".seq";
                /* check if key dataset 'seq' has already existed */
                if (key_name.compare(dataset_name) == 0) {
                    printf("Error: dataset '%s' already existed in group %s\n",
                           dataset_name,group_name);
                    return -1;
                }
                /* check if this is the partition key base dataset */
                if (part_key_base.compare(dataset_name) == 0) {
                    assert(groups[grp_no].key_base == NULL);
                    dset.is_key_base = true;
                    groups[grp_no].key_base = &dset;
                }
            }
            total_num_datasets++;
        }
        else {
            /* accumulate this dataset's dimension sizes across all assigned
             * input files */
            dset.local_dims[0] += dset_dims[0];

            assert(dset_dims[1] == dset.local_dims[1]);
            assert(dset.name.compare(dataset_name) == 0);
        }
        /* Note all datasets in the same group share the size of dimension 0 */
        groups[grp_no].shared_dim0 = dset.local_dims[0];

        /* a dataset may be of zero size */
        if (dset_dims[0] > 0) {
            dset.in_dset_ids.push_back(dset_id);  /* save dataset ID */
            dset.in_dim0.push_back(dset_dims[0]); /* save dataset dim0 size */
            hsize_t dset_size = dset_dims[0] * dset_dims[1] * dset.type_size;
            /* find the max size among all datasets */
            if (max_local_size_in_bytes < dset_size)
                max_local_size_in_bytes = dset_size;
        }
        else {
            /* This dataset is of zero size in this input file. Keep the vector
             * size the same as num_input_files, in case this dataset is of
             * size zero but non-zero sized in another file.
             */
            dset.in_dset_ids.push_back(-1);
            dset.in_dim0.push_back(0);
            /* No need to keep zero-size dataset open */
            err = H5Dclose(dset_id);
            if (err < 0) RETURN_ERROR("H5Dclose",name)
        }

        free(name_copy);
    }

    return 0;
}

/*----< calculate_num_writes() >---------------------------------------------*/
void Concatenator::calculate_num_writes()
{
    size_t ii, jj, kk;
    size_t *numWrites;

    numWrites = new size_t[total_num_datasets];

    /* calculate number of rounds of collective writes per dataset, given the
     * I/O buffer size
     */
    kk = 0;
    for (ii=0; ii<num_groups; ii++) {
        for (jj=0; jj<groups[ii].num_dsets; jj++) {
            DSInfo_t &dset = groups[ii].dsets[jj];
            assert(dset.global_dims[0] > 0);
            size_t nrows_in_buffer = io_buffer_size
                                   / (dset.local_dims[1] * dset.type_size);
            numWrites[kk++] = (dset.local_dims[0] % nrows_in_buffer == 0) ?
                              (dset.local_dims[0] / nrows_in_buffer) :
                              (dset.local_dims[0] / nrows_in_buffer) + 1;
        }
    }
    assert(kk == total_num_datasets);

    /* find the maximumsr among all processes */
    MPI_Allreduce(MPI_IN_PLACE, numWrites, total_num_datasets, MPI_SIZE_T,
                  MPI_MAX, comm);
    num_allreduce++;

    /* copy to local objects */
    kk = 0;
    for (ii=0; ii<num_groups; ii++)
        for (jj=0; jj<groups[ii].num_dsets; jj++)
            groups[ii].dsets[jj].num_writes = numWrites[kk++];

    delete []numWrites;
}

#if defined DEBUG && DEBUG
/*----< check_h5_objects() >-------------------------------------------------*/
static
int check_h5_objects(const char *filename, hid_t fid)
{
    ssize_t num_objs;
    size_t ii, howmany;
    H5I_type_t ot;
    hid_t *objs;
    char obj_name[1024];

    num_objs = H5Fget_obj_count(fid, H5F_OBJ_ALL);
    if (num_objs <= 0) return num_objs;
    /* ignore FILE object */
    if (num_objs > 1) printf("%zd object(s) open\n", num_objs);

    objs = (hid_t*) malloc(num_objs * sizeof(hid_t));
    howmany = H5Fget_obj_ids(fid, H5F_OBJ_ALL, num_objs, objs);
    if (howmany > 1) printf("open objects:\n");

    for (ii=0; ii<howmany; ii++) {
        string type_name="";
        ot = H5Iget_type(objs[ii]);
             if (ot == H5I_FILE)      continue; /*  type_name = "H5I_FILE"; */
        else if (ot == H5I_GROUP)     type_name = "H5I_GROUP";
        else if (ot == H5I_DATATYPE)  type_name = "H5I_DATATYPE";
        else if (ot == H5I_DATASPACE) type_name = "H5I_DATASPACE";
        else if (ot == H5I_DATASET)   type_name = "H5I_DATASET";
        else if (ot == H5I_ATTR)      type_name = "H5I_ATTR";
        H5Iget_name(objs[ii], obj_name, 1024);
        printf("Still opened in file %s %4zd: type %s, name %s\n",
               filename, ii, type_name.c_str(), obj_name);
    }
    free(objs);
    return howmany;
}
#endif

/*----< close_output_file() >------------------------------------------------*/
int Concatenator::close_output_file()
{
    int err_exit=0;

#if defined DEBUG && DEBUG
    check_h5_objects(output_file_name.c_str(), output_file_id);
#endif

    herr_t err = H5Fclose(output_file_id);
    if (err < 0) HANDLE_ERROR("H5Fclose")

    err = H5Pclose(dxpl_id);
    if (err < 0) HANDLE_ERROR("H5Pclose")

fn_exit:
    return err_exit;
}

/*----< close_input_files() >------------------------------------------------*/
int Concatenator::close_input_files()
{
    int err_exit=0;

    /* Close all the input files. */
    while (input_files.size() > 0) {
        auto infile = input_files.begin();
#if defined DEBUG && DEBUG
        check_h5_objects(infile->first.c_str(), infile->second);
#endif
        herr_t err = H5Fclose(infile->second);
        if (err < 0) HANDLE_ERROR("H5Fclose")
        input_files.erase(infile);
    }

fn_exit:
    return err_exit;
}

/*----< accumulate_dimensions() >--------------------------------------------*/
void Concatenator::accumulate_dimensions()
{
    size_t ii, jj;
    hsize_t *my_lens, *my_offs;

    my_lens = new hsize_t[num_groups];
    my_offs = new hsize_t[num_groups];

    /* copy local aggregated sizes to a contiguous array */
    for (ii=0; ii<num_groups; ii++) {
        my_lens[ii] = groups[ii].shared_dim0;
        my_offs[ii] = 0;
    }

    /* prefix scan to get the starting offset for each group */
    MPI_Exscan(my_lens, my_offs, num_groups, MPI_HSIZE_T, MPI_SUM, comm);
    num_exscan++;

    /* Sum the dataset 1st dimension sizes among all processes. */
    MPI_Allreduce(MPI_IN_PLACE, my_lens, num_groups, MPI_HSIZE_T,
                  MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, &max_local_size_in_bytes, 1, MPI_HSIZE_T,
                  MPI_MAX, comm);
    num_allreduce += 2;

    /* calculate global dataset metadata */
    for (ii=0; ii<num_groups; ii++) {
        /* all datasets in the same group share the same 1st dimension size */
        groups[ii].shared_dim0 = my_lens[ii]; /* update to global_dims[0] */

        DSInfo_t *dsets = groups[ii].dsets;
        for (jj=0; jj<groups[ii].num_dsets; jj++) {
            dsets[jj].global_off     = my_offs[ii];
            dsets[jj].global_dims[0] = my_lens[ii];
            dsets[jj].global_dims[1] = dsets[jj].local_dims[1];
            dsets[jj].size_in_bytes  = dsets[jj].global_dims[0] *
                                       dsets[jj].global_dims[1] *
                                       dsets[jj].type_size;

            /* If the aggregated size is 0, do not chunk the dataset. */
            if (dsets[jj].global_dims[0] == 0)
                dsets[jj].layout = H5D_COMPACT;
            else if (dsets[jj].size_in_bytes < compress_threshold)
                /* Disable compression for datasets of size smaller than
                 * MIN_DATASET_SIZE. This can reduce the number of chunks, thus
                 * the size of metadata, and cost of metadata operations.
                 */
                dsets[jj].layout = H5D_CONTIGUOUS;
        }
    }

    delete []my_offs;
    delete []my_lens;
}

/*----< calculate_chunk_size() >---------------------------------------------*/
void Concatenator::calculate_chunk_size(DSInfo_t &dset)
{
    size_t chunk_len;

    /* no zero-sized dataset should call this subroutine */
    assert(dset.global_dims[0] > 0);

    /* no chunking along the 2nd dimension */
    dset.chunk_dims[1] = dset.global_dims[1];

#define USE_NELEM_BASED_CHUNKING
#ifdef USE_NELEM_BASED_CHUNKING
    chunk_len = chunk_size_threshold;
#else
    chunk_len = chunk_size_threshold / dset.type_size;
#endif

#if USE_MIN_DSET_SIZE
    /* set chunk size to global dim size if smaller than the threshold */
    if (dset.global_dims[1] == 1) /* 1D dataset */
        /* use whichever is smaller */
        dset.chunk_dims[0] = (dset.global_dims[0] < chunk_len) ?
                              dset.global_dims[0] : chunk_len;
    else /* 2D dataset */
        dset.chunk_dims[0] = (dset.global_dims[0] < 128) ?
                              dset.global_dims[0] : 128;
#else
    /* no matter the dim size if larger or smaller than the threshold */
    if (dset.global_dims[1] == 1) /* 1D dataset */
        dset.chunk_dims[0] = chunk_len;
    else { /* 2D dataset */
        if (dset.global_dims[1] > 1024)
            dset.chunk_dims[0] = 128;
        else
            dset.chunk_dims[0] = chunk_len;
    }
#endif

    dset.chunk_size = dset.chunk_dims[0] * dset.chunk_dims[1]
                    * dset.type_size;
}

/*----< create_dataset() >---------------------------------------------------*/
int Concatenator::create_dataset(hid_t     group_id,
                                 DSInfo_t& dset,
                                 bool      toFill)
{
    int err_exit=0;
    herr_t err;
    hid_t space_id, dcpl_id;
#if defined PROFILE && PROFILE
    double start;
    start = MPI_Wtime();
#endif

    /* Create a space to define dataset's dimensions. */
    if (dset.layout == H5D_CHUNKED) {
        hsize_t maxdims[2] = {H5S_UNLIMITED, dset.global_dims[1]};
        space_id = H5Screate_simple(2, dset.global_dims, maxdims);
    }
    else /* fixed-size dataset (H5D_COMPACT or H5D_CONTIGUOUS)*/
        space_id = H5Screate_simple(2, dset.global_dims, NULL);

    if (space_id < 0) HANDLE_ERROR("H5Screate_simple")

    /* Create chunking & compression setting. */
    if (dset.layout == H5D_CHUNKED) {
        dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (dcpl_id < 0) HANDLE_ERROR("H5Pcreate")

        calculate_chunk_size(dset);

        /* H5Pset_chunk changes/sets the layout to H5D_CHUNKED */
        err = H5Pset_chunk(dcpl_id, 2, dset.chunk_dims);
        if (err < 0) {
            cout<<"size: "<<dset.global_dims[0]<<" x "<<dset.global_dims[1]
                <<" chunks: "<<dset.chunk_dims[0]<<" x "<<dset.chunk_dims[1]
                <<endl;
            HANDLE_ERROR("H5Pcreate")
        }

#if 0
        /* NOTE: A workaround for avoiding 'monotonically non-decreasing
         * offset' issue. Once a file is collectively opened (created), when
         * many datasets with different data types are created (written), an
         * assertion failure occurs from H5Dchunk.c. By setting the compression
         * property as follows, we can avoid the issue. NOTE: 10/26/2019, this
         * setting turns out to cause wrong data in the partial chunk. For now,
         * we do not know why but the 'monotonically non-decreasing offset'
         * issue is not reproduced even without
         * H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS setting. Temporarily, the
         * below setting is commented out. We will find out the reason asap.
         */
        unsigned int chunk_info;

        err = H5Pget_chunk_opts(dcpl_id, &chunk_info);
        if (err < 0) HANDLE_ERROR("H5Pget_chunk_opts")

        chunk_info |= H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS;
        err = H5Pset_chunk_opts(dcpl_id, chunk_info);
        if (err < 0) HANDLE_ERROR("H5Pset_chunk_opts")
#endif

        if (zip > 0) {
            err = H5Pset_deflate(dcpl_id, zip);
            if (err < 0) HANDLE_ERROR("H5Pset_deflate")
        }

        if (toFill) {
            /* HDF5 default fill value is 0, skip call H5Pset_fill_value() */
            err = H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_ALLOC);
            if (err < 0) HANDLE_ERROR("H5Pset_fill_time")
        }
        else {
            /* Never fill in the chunk in advance. */
            err = H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);
            if (err < 0) HANDLE_ERROR("H5Pset_fill_time")
        }

        /* setting H5D_ALLOC_TIME_LATE takes no effect for parallel I/O on
         * compressed datasets, which HDF5 always uses H5D_ALLOC_TIME_EARLY.
         */
        err = H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_LATE);
        if (err < 0) HANDLE_ERROR("H5Pset_alloc_time")
    }
    else if (dset.layout == H5D_COMPACT) {
        dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (dcpl_id < 0) HANDLE_ERROR("H5Pcreate")

        err = H5Pset_layout(dcpl_id, dset.layout);
        if (err < 0) HANDLE_ERROR("H5Pset_layout")
    }
    else if (dset.layout == H5D_CONTIGUOUS) {
        dcpl_id = H5P_DEFAULT;
    }
    else {
        cout<<"dset.layout should've set to CHUNKED, COMPACT, or CONTIGUOUS..."<<endl;
        return -1;
    }

    /* Create a dataset in the output file. */
    dset.out_dset_id = H5Dcreate2(group_id, dset.name.c_str(), dset.type_id,
                                  space_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    if (dset.out_dset_id < 0) HANDLE_ERROR(string("H5Dcreate2 ")+dset.name)

    if (dcpl_id != H5P_DEFAULT) {
        err = H5Pclose(dcpl_id);
        if (err < 0) HANDLE_ERROR("H5Pclose")
    }
    err = H5Sclose(space_id);
    if (err < 0) HANDLE_ERROR("H5Sclose")

    /* no need to keep zero-sized datasets open */
    if (dset.global_dims[0] == 0) {
        if (dset.type_id != H5T_STD_I64LE) {
            /* the partitioning key dataset is created using the predefined
             * datatype H5T_STD_I64LE. It is illegal to close predefined type.
             */
            err = H5Tclose(dset.type_id);
            if (err < 0) HANDLE_ERROR("H5Tclose")
        }
        err = H5Dclose(dset.out_dset_id);
        if (err < 0) HANDLE_ERROR("H5Dclose")
    }
fn_exit:
#if defined PROFILE && PROFILE
    c_1d_2d += MPI_Wtime() - start;
#endif
    return err_exit;
}

/*----< open_dataset() >-----------------------------------------------------*/
int Concatenator::open_dataset(hid_t     group_id,
                               DSInfo_t& dset)
{
    int err_exit=0, ndims;
    herr_t err;
    hid_t space_id;
    hsize_t dset_dims[2];
#if defined PROFILE && PROFILE
    double start;
    start = MPI_Wtime();
#endif

    /* only non-zero sized datasets will reach here */
    assert(dset.global_dims[0] > 0);

    /* Open an existing dataset in the output file. */
    dset.out_dset_id = H5Dopen(group_id, dset.name.c_str(), H5P_DEFAULT);
    if (dset.out_dset_id < 0) HANDLE_ERROR(string("H5Dopen ")+dset.name)

    /* collect the current dataset size */
    space_id = H5Dget_space(dset.out_dset_id);
    if (space_id < 0) HANDLE_ERROR("H5Dget_space")
    ndims = H5Sget_simple_extent_dims(space_id, dset_dims, NULL);
    if (ndims < 0) RETURN_ERROR("H5Sget_simple_extent_dims",dset.name.c_str())
    err = H5Sclose(space_id);
    if (err < 0) RETURN_ERROR("H5Sclose",dset.name.c_str())

    /* extend the dim[0] size by dset_dims[0] */
    if (dset_dims[0] > 0) {
        dset.global_dims[0] += dset_dims[0];
        err = H5Dset_extent(dset.out_dset_id, dset.global_dims);
        if (err < 0) RETURN_ERROR("H5Dset_extent",dset.name.c_str())

        /* update the starting offset of the appending */
        dset.global_off += dset_dims[0];
    }

fn_exit:
#if defined PROFILE && PROFILE
    c_1d_2d += MPI_Wtime() - start;
#endif
    return err_exit;
}

