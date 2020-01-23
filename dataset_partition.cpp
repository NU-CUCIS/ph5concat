/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
 * See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <assert.h>
#include "ph5_concat.hpp"

int Concatenator::concat_small_datasets(std::vector<std::string> const &inputs)
{
    int err_exit=0;
    size_t ii, jj, kk;
    herr_t err;

    /* concatenate only 1D datasets */
    err = concat_datasets(false);
    if (err < 0) HANDLE_ERROR("H5Dclose")

    /* close large 2D datasets, as they will be collective opened later */
    for (ii=0; ii<num_groups; ii++) {
        for (jj=0; jj<groups[ii].num_dsets; jj++) {
            DSInfo_t &dset = groups[ii].dsets[jj];
            if (dset.global_dims[1] > 1) {
                for (kk=0; kk<dset.in_dset_ids.size(); kk++) {
                    err = H5Dclose(dset.in_dset_ids[kk]);
                    if (err < 0) HANDLE_ERROR("H5Dclose")
                }
            }
        }
    }

    /* Close all the input files */
    close_input_files();

fn_exit:
    return err_exit;
}

int Concatenator::open_input_files(std::vector<std::string> files,
                                   bool collective_io)
{
    int err_exit=0;
    herr_t err;
    hid_t file_id;
    hid_t plist_id;
#if defined PROFILE && PROFILE
    double start;
#endif

#if defined PROFILE && PROFILE
    start = MPI_Wtime();
#endif

    /* remove all elements previous added */
    input_files.clear();

    /* Open all the input files */
    for (auto file = files.begin(); file < files.end(); file++) {
        auto entry = input_files.find(*file);
        if (entry == input_files.end()) { /* new file */
            if (collective_io == true) {
                plist_id = H5Pcreate(H5P_FILE_ACCESS);
                if (plist_id < 0) HANDLE_ERROR("H5Pcreate")
                err = H5Pset_fapl_mpio(plist_id, comm, info);
                if (err < 0) HANDLE_ERROR("H5Pset_fapl_mpio")
                err = H5Pset_all_coll_metadata_ops(plist_id, true);
                if (err < 0) HANDLE_ERROR("H5Pset_all_coll_metadata_ops")
                file_id = H5Fopen(file->c_str(), H5F_ACC_RDONLY, plist_id);
                err = H5Pclose(plist_id);
                if (err < 0) HANDLE_ERROR("H5Pclose")
            }
            else {
                file_id = H5Fopen(file->c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            }
            if (file_id < 0) HANDLE_ERROR("H5Fopen")
            input_files.insert(std::make_pair(*file, file_id));
        }
        else {
            /* error: repeated file detected */
            std::cout<<"["<<__FILE__<<"]["<<__FUNCTION__<<"]["<<__LINE__<<"] "<<*file<<std::endl;
            err_exit = -1;
            goto fn_exit;
        }
    }
#if defined PROFILE && PROFILE
    o_f += MPI_Wtime() - start;
#endif

fn_exit:
    return err_exit;
}

int Concatenator::concat_large_datasets(std::vector<std::string> const &inputs)
{
    int err_exit=0;
    size_t ii, jj, kk;
    herr_t err;
    hid_t file_id;

#if defined PROFILE_HDF5_INTERNAL && PROFILE_HDF5_INTERNAL
    MPI_Barrier(comm);
    H5Vreset();
    H5Venable();
#endif

    /* All processes open all input files in MPI-I/O mode, and read large
     * datasets collectively */
    open_input_files(inputs, true);

    /* Retrieve the file handle of output file */
    file_id = output_file_id;

    for (ii=0; ii<num_groups; ii++) {
        for (jj=0; jj<groups[ii].num_dsets; jj++) {
            DSInfo_t &dset = groups[ii].dsets[jj];

            if (dset.global_dims[1] == 1) /* skip 1D datasets */
                continue;

            assert(dset.global_dims[0] > 0);
            dset.cur_offset = 0;
            dset.cur_chunk_offset = 0;

            for (kk=0; kk<inputs.size(); kk++) {
                hsize_t counts[2], offsets[2], dset_size[2];
                auto ff = input_files.find(inputs[kk]);
                /* Open a dataset from one input file and read the entire local
                 * dataset. We assume any local dataset fits into the memory
                 * space of the node. */
                std::string path = "/" + groups[ii].name + "/" + dset.name;
                hid_t id;
                id = H5Dopen(ff->second, path.c_str(), H5P_DEFAULT);
                if (id < 0) HANDLE_ERROR("H5Dopen")

                err = read_2d_dataset(id, dset, counts, offsets, dset_size);
                if (err < 0) {
                    std::cout<<"read_2d_dataset failed."<<std::endl;
                    return -1;
                }

                /* Close the input dataset. */
                err = H5Dclose(id);
                if (err < 0) HANDLE_ERROR("H5Dclose")

                err = write_2d_dataset(dset, counts, offsets, dset_size);
                if (err < 0) {
                    std::cout<<"write_2d_dataset failed."<<std::endl;
                    return -1;
                }
            }

            /* Close the output dataset. */
            err = H5Dclose(dset.out_dset_id);
            dset.out_dset_id = -1;
            if (err < 0) HANDLE_ERROR("H5Dclose")
        }
    }

#if defined PROFILE_HDF5_INTERNAL && PROFILE_HDF5_INTERNAL
    H5Vprint();
    H5Vdisable();
#endif

fn_exit:
    return err_exit;
}

int Concatenator::read_dataset(int       file_no, /* input file number */
                               DSInfo_t &dset)
{
    int err_exit=0;
    herr_t err;
    hid_t dset_id;
#if defined PROFILE && PROFILE
    double start;

    start = MPI_Wtime();
#endif

    dset_id = dset.in_dset_ids[file_no];

    /* Read the whole dataset. */
    err = H5Dread(dset_id, dset.type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  &buffer[dset.cur_offset]);
    if (err < 0) HANDLE_ERROR("H5Dread")

#if defined PROFILE && PROFILE
    if (dset.global_dims[1] == 1)
        r_1d += MPI_Wtime() - start;
    else
        r_2d += MPI_Wtime() - start;
#endif

    /* Close the input dataset. */
    err = H5Dclose(dset_id);
    if (err < 0) HANDLE_ERROR("H5Dclose")

    dset.cur_offset += dset.in_dim0[file_no] * dset.local_dims[1] * dset.type_size;

fn_exit:
    return err_exit;
}

int Concatenator::write_dataset(DSInfo_t &dset,
                                hsize_t  *offs,
                                hsize_t  *lens)
{
    int err_exit=0;
    herr_t err;
    hid_t space_id=-1, memspace_id=-1;
#if defined PROFILE && PROFILE
    double start;

    MPI_Barrier(comm);
    start = MPI_Wtime();
#endif

#if defined DEBUG && DEBUG
    std::cout<<"R"<<rank<<" "<<dset.name.c_str()<<" type: "<<dset.type_id
             <<" type_size: "<<dset.type_size<<" num_rows x cols: "
             <<dset.global_dims[0]<<" x "<<dset.global_dims[1]<<" write off: "
             <<offs[0]<<" x "<<offs[1]<<" lens: "<<lens[0]<<" x "<<lens[1]
             <<std::endl;
#endif

    /* Setup the memory space and the file space for the wrtie. */
    memspace_id = H5Screate_simple(2, lens, NULL);
    if (memspace_id < 0) HANDLE_ERROR("H5Screate_simple")
    space_id = H5Dget_space(dset.out_dset_id);
    if (space_id < 0) HANDLE_ERROR("H5Dget_space")
    err = H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offs, NULL, lens, NULL);
    if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab")

    err = H5Dwrite(dset.out_dset_id, dset.type_id, memspace_id, space_id,
                   dxpl_id, buffer);
    if (err < 0) HANDLE_ERROR("H5Dwrite")

#if defined PROFILE && PROFILE
    if (dset.global_dims[1] == 1)
        w_1d += MPI_Wtime() - start;
    else
        w_2d += MPI_Wtime() - start;
#endif

fn_exit:
    if (memspace_id >= 0)
        H5Sclose(memspace_id);
    if (space_id >= 0)
        H5Sclose(space_id);

    return err_exit;
}

int Concatenator::read_2d_dataset(hid_t     dset_id,
                                  DSInfo_t &dset,
                                  hsize_t  *counts,
                                  hsize_t  *offsets,
                                  hsize_t  *dset_size)
{
    int err_exit=0;
    herr_t err;
    hid_t space_id;
    hid_t memspace_id;
    hsize_t read_offsets[2], new_size;
#if defined PROFILE && PROFILE
    double start;
    start = MPI_Wtime();
#endif

    /* First we initialize the counts to 0 here.
     * If something goes wrong, we will not try to write anything. */
    counts[0] = 0;
    counts[1] = 0;
    offsets[0] = 0;
    offsets[1] = 0;

    space_id = H5Dget_space(dset_id);
    if (space_id < 0) HANDLE_ERROR("H5Dget_space")
    err = H5Sget_simple_extent_dims(space_id, dset_size, NULL);
    if (err < 0) HANDLE_ERROR("H5Sget_simple_extent_dims")

    /* Compute the read length. Note that the offsets are computed for write. */
    numerology(dset, dset_size, offsets, counts);

    /* Calculate read offsets. */
    read_offsets[0] = offsets[0] - (dset.cur_offset / (dset.global_dims[1] * dset.type_size));
    read_offsets[1] = 0;

    err = H5Sselect_hyperslab(space_id, H5S_SELECT_SET, read_offsets, NULL, counts, NULL);
    if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab")

    /* Create a memspac for collective read. */
    memspace_id = H5Screate_simple(2, counts, NULL);
    if (memspace_id < 0) HANDLE_ERROR("H5Screate_simple")

    new_size = counts[0] * counts[1] * dset.type_size;
    if (io_buffer_size < new_size) {
        io_buffer_size = new_size;
        delete [] buffer;
        buffer = new char[io_buffer_size];
    }

    err = H5Dread(dset_id, dset.type_id, memspace_id, space_id, dxpl_id, buffer);
    if (err < 0) HANDLE_ERROR("H5Dread")

    err = H5Sclose(memspace_id);
    if (err < 0) HANDLE_ERROR("H5Sclose")
    err = H5Sclose(space_id);
    if (err < 0) HANDLE_ERROR("H5Sclose")

#if defined PROFILE && PROFILE
    r_2d += MPI_Wtime() - start;
#endif

fn_exit:
    return err_exit;
}

int Concatenator::write_2d_dataset(DSInfo_t &dset,
                                   hsize_t  *counts,
                                   hsize_t  *offsets,
                                   hsize_t  *dset_size)
{
    int err_exit=0;
    herr_t err;
    hid_t space_id;
    hid_t memspace_id;
    hsize_t dset_size_in_bytes;
#if defined PROFILE && PROFILE
    double start;
#endif

    assert(dset.global_dims[0] > 0);
#if defined DEBUG && DEBUG
    std::cout<<"R"<<rank<<" "<<dset.name.c_str()<<" type: "<<dset.type_id
             <<" type_size: "<<dset.type_size<<" num_rows x cols: "
             <<dset.global_dims[0]<<" x "<<dset.global_dims[1]
             <<" chunk_size: "<<dset.chunk_size<<std::endl;
#endif
    /* Setup the memory space and the file space for the wrtie. */
    memspace_id = H5Screate_simple(2, counts, NULL);
    if (memspace_id < 0) HANDLE_ERROR("H5Screate_simple")
    space_id = H5Dget_space(dset.out_dset_id);
    if (space_id < 0) HANDLE_ERROR("H5Dget_space")
    err = H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offsets, NULL, counts, NULL);
    if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab")

    /* Write the data. */
#if defined PROFILE && PROFILE
    start = MPI_Wtime();
#endif
    err = H5Dwrite(dset.out_dset_id, dset.type_id, memspace_id, space_id, dxpl_id, buffer);
    if (err < 0) HANDLE_ERROR("H5Dwrite")
#if defined PROFILE && PROFILE
    w_2d += MPI_Wtime() - start;
#endif

    err = H5Sclose(memspace_id);
    if (err < 0) HANDLE_ERROR("H5Sclose")
    err = H5Sclose(space_id);
    if (err < 0) HANDLE_ERROR("H5Sclose")

    /* Update the current offset in the output file.
     * Since we go through all the input file one after another,
     * We need to keep the offset for the next write. */
    dset_size_in_bytes = dset_size[0] * dset_size[1] * dset.type_size;
    dset.cur_offset += dset_size_in_bytes;
    dset.cur_chunk_offset = dset.cur_offset - (dset.cur_offset % dset.chunk_size);
fn_exit:
    return err_exit;
}

int Concatenator::numerology(DSInfo_t &dset,
                             hsize_t  *dset_size, /* dataset size in a file */
                             hsize_t  *offsets,
                             hsize_t  *counts)
{
    size_t num_chunks;
    size_t my_starting_chunk;
    size_t num_my_chunks;
    size_t evenly_partitioned;
    size_t start_chunk;
    size_t extra_chunks;
    hsize_t dset_size_in_bytes;
    hsize_t my_offset, my_length;

    /* 0. Assertions. */
    if (dset.global_dims[0] == 0) {
        counts[0] = 0;
        counts[1] = 0;
        offsets[0] = 0;
        offsets[1] = 0;
        return 0;
    }

    if (dset.chunk_size <= 0) {
        std::cout<<"/"<<dset.name.c_str()<<"/"<<dset.name.c_str()<<
                   " has an invalid chunk_size."<<std::endl;
        return -1;
    }

    /* 1. Calculate the number of chunks that will be processed by this rank. */
    dset_size_in_bytes = dset_size[0] * dset_size[1] * dset.type_size;
    dset_size_in_bytes += (dset.cur_offset - dset.cur_chunk_offset);

    num_chunks = dset_size_in_bytes / dset.chunk_size;
    if (dset_size_in_bytes % dset.chunk_size > 0)
        num_chunks++;
    evenly_partitioned = num_chunks / nprocs;
    extra_chunks = num_chunks % nprocs;
    if (extra_chunks > 0)
        num_my_chunks = ((extra_chunks > static_cast<unsigned int>(rank)) ? evenly_partitioned + 1 : evenly_partitioned);
    else
        num_my_chunks = evenly_partitioned;

    /* 2. Calculate my_offset and my_length for a collective write. */
    if (num_my_chunks > 0 ) {
        start_chunk = dset.cur_chunk_offset / dset.chunk_size;
        my_starting_chunk = start_chunk + (evenly_partitioned * rank);
        my_starting_chunk += (extra_chunks > static_cast<unsigned int>(rank)) ? rank : extra_chunks;

        my_offset = (dset.cur_offset > my_starting_chunk * dset.chunk_size) ?
                    dset.cur_offset :
                    my_starting_chunk * dset.chunk_size;

        my_length = num_my_chunks * dset.chunk_size;

        if (my_offset == dset.cur_offset) // first chunk
            my_length -= (dset.cur_offset - dset.cur_chunk_offset);

        dset_size_in_bytes = dset_size[0] * dset_size[1] * dset.type_size;
        if ((my_offset + my_length) > (dset.cur_offset + dset_size_in_bytes)) // last chunk
            my_length -= ((my_offset + my_length) - (dset.cur_offset + dset_size_in_bytes));

        offsets[0] = my_offset / (dset.global_dims[1] * dset.type_size);
        offsets[1] = 0;

        counts[0] = my_length / (dset.global_dims[1] * dset.type_size);
        counts[1] = dset.global_dims[1];
    }
    else {
        offsets[0] = 0;
        offsets[1] = 0;

        counts[0] = 0;
        counts[1] = 0;
    }

    return 0;
}
