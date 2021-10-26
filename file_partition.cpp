/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
 * See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <assert.h>
#include <string.h> /* strncpy() */

#include "ph5_concat.hpp"

#if defined PROFILE && PROFILE
    #define SET_TIMER(ts) { ts = MPI_Wtime(); }
    #define GET_TIMER(ts, t) { t += MPI_Wtime() - ts; }
#else
    #define SET_TIMER(ts)
    #define GET_TIMER(ts, t)
#endif

/*----< write_partition_key_dataset() >--------------------------------------*/
/* For each group, read the key base dataset, generate key seq dataset, write
 * both key base and key seq datasets to the output file.
 */
int Concatenator::write_partition_key_dataset()
{
    int err_exit=0;
    herr_t err;
    hsize_t ii, jj, kk, start[2];

    /* num_groups now is the number of non-zero size groups */
    for (ii=0; ii<num_groups; ii++) {
        assert(groups[ii].shared_dim0 > 0);

        DSInfo_t     *base = groups[ii].key_base;
        DSInfo_t     *seq = groups[ii].seq_dset;
        unsigned int *base_buf, *base_ptr;
        int64_t      *seq_buf, *seq_ptr;

        /* allocate buffer to read the key base dataset. Note local_dims[0] is
         * the aggregated dim[0] size of assigned input files.
         */
        base_buf = new unsigned int [base->local_dims[0]];
        base_ptr = base_buf;
        seq_buf = new int64_t [base->local_dims[0]];
        seq_ptr = seq_buf;

        /* iterate all input files assigned to this process */
        for (jj=0; jj<num_input_files; jj++) {
            if (base->in_dim0[jj] == 0) continue;

            hid_t in_dset_id = base->in_dset_ids[jj];

            /* read the entire partition base dataset */
            err = H5Dread(in_dset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, base_ptr);
            if (err < 0) HANDLE_ERROR("H5Dread")

            /* close input key base dataset, as it is no longer used */
            err = H5Dclose(in_dset_id);
            if (err < 0) HANDLE_ERROR("H5Dclose");

            /* populate contents of key seq dataset from lookup table jj */
            for (kk=0; kk<base->in_dim0[jj]; kk++)
                seq_ptr[kk] = lookup_table[jj][base_ptr[kk]];

            base_ptr += base->in_dim0[jj];
            seq_ptr  += base->in_dim0[jj];
        }

        /* write partition key base dataset */
        start[0] = base->global_off;
        start[1] = 0;
        err = write_dataset_2D(*base, start, base->local_dims, base_buf);
        if (err < 0) HANDLE_ERROR("write_dataset_2D")

        /* close partition key base dataset */
        err = H5Tclose(base->type_id);
        if (err < 0) HANDLE_ERROR("H5Tclose")
        err = H5Dclose(base->out_dset_id);
        if (err < 0) HANDLE_ERROR("H5Dclose");
        base->out_dset_id = -1;

        /* write partition key seq dataset */
        err = write_dataset_2D(*seq, start, base->local_dims, seq_buf);
        if (err < 0) HANDLE_ERROR("write_dataset_2D")

        /* close partition key seq dataset */
        err = H5Dclose(seq->out_dset_id);
        if (err < 0) HANDLE_ERROR("H5Dclose");
        seq->out_dset_id = -1;

        delete [] base_buf;
        delete [] seq_buf;
    }

fn_exit:
    return err_exit;
}

/*----< generate_partition_key() >-------------------------------------------*/
/* For each input file, read the entire dataset /spill/key_base into memory,
 * and use its contents to build the key lookup table. The contents of
 * /spill/key_base must be a list of unique integers. Lookup table can hash
 * values in key_base to global (concatenated) unique keys. It will later
 * be used to generate partition key seq dataset in each group.
 */
int Concatenator::generate_partition_key(GrpInfo &grp)  /* group /spill */
{
    int err_exit=0;
    size_t ii, jj;
    herr_t err;
    int64_t seq_off;
    unsigned int *base_buf, *buf_ptr;

    assert(grp.key_base != NULL);
    assert(grp.key_base->global_dims[0] > 0);

    /* allocate buffer to read the partition key base dataset. Note
     * local_dims[0] is the aggregated dim[0] size of assigned input files.
     */
    base_buf = new unsigned int [grp.key_base->local_dims[0]];
    buf_ptr = base_buf;
    seq_off = grp.key_base->global_off;

    /* iterate all input files assigned to this process */
    for (ii=0; ii<num_input_files; ii++) {
        if (grp.key_base->in_dim0[ii] == 0) continue;

        hid_t in_dset_id = grp.key_base->in_dset_ids[ii];

        /* read the entire partition key base dataset */
        err = H5Dread(in_dset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, buf_ptr);
        if (err < 0) HANDLE_ERROR("H5Dread")

        /* use contents of key base dataset to construct a lookup hash table */
        for (jj=0; jj<grp.key_base->in_dim0[ii]; jj++)
            lookup_table[ii][buf_ptr[jj]] = seq_off + jj;

        buf_ptr += grp.key_base->in_dim0[ii];
        seq_off += grp.key_base->in_dim0[ii];
    }

fn_exit:
    delete [] base_buf;
    return err_exit;
}

/*----< concat_datasets() >--------------------------------------------------*/
int Concatenator::concat_datasets(bool process_large_dsets)
{
    int err_exit=0;
    herr_t err;
    hsize_t ii, jj, kk;
#if defined PROFILE && PROFILE
    double ts;
#endif
#if defined PROFILE_HDF5_INTERNAL && PROFILE_HDF5_INTERNAL
    MPI_Barrier(comm);
    H5Vreset();
    H5Venable();
#endif

    if (add_partition_key && !process_large_dsets)  {
        /* Read partition key base dataset from group /spill and generate the
         * lookup table.
         */
        err = generate_partition_key(groups[spill_grp_no]);
        if (err < 0) HANDLE_ERROR("generate_partition_key")
    }

    /* Traverse over all the groups and datasets. */
    for (ii=0; ii<num_groups; ii++) {
        size_t num_dsets = groups[ii].num_dsets;

        if (add_partition_key && !process_large_dsets)  {
            /* Skip writing the partition key seq dataset. -1 is because key
             * seq is the last dataset added.
             */
            num_dsets--;
        }

        for (jj=0; jj<num_dsets; jj++) {
            /* Skip writing the key base dataset */
            if (groups[ii].dsets[jj].is_key_base) continue;

            DSInfo_t &dset = groups[ii].dsets[jj];
            assert(dset.global_dims[0] > 0);

            if (process_large_dsets) {
                if (dset.global_dims[1] == 1)
                    continue;
            }
            else if (dset.global_dims[1] > 1)
                continue;

            hsize_t row_size = dset.global_dims[1] * dset.type_size;
            hsize_t nrows_in_buf = io_buffer_size / row_size;
            hsize_t nrows_written = 0;
            hsize_t out_offset = 0;
            size_t num_writes = 0;

            /* iterate input files assigned to this process */
            for (kk=0; kk<num_input_files; kk++) {
                if (dset.in_dim0[kk] == 0) continue;
                hid_t in_dset_id = dset.in_dset_ids[kk];
                hsize_t in_dset_nrows = dset.in_dim0[kk];
                hsize_t offset_in_dset = 0;

                /* do multiple rounds of reads and writes in case the I/O
                 * buffer size is smaller than the dataset size
                 */
                while (offset_in_dset < in_dset_nrows) {
                    hsize_t rows_to_read;
                    if (nrows_in_buf - nrows_written < in_dset_nrows - offset_in_dset)
                        rows_to_read = nrows_in_buf - nrows_written;
                    else
                        rows_to_read = in_dset_nrows - offset_in_dset;

                    /* Read from file kk into I/O buffer */
                    err = read_dataset2(dset, kk, in_dset_nrows,
                                        offset_in_dset, rows_to_read,
                                        nrows_written * row_size);
                    if (err < 0) HANDLE_ERROR("read_dataset2 /"+groups[ii].name+"/"+dset.name)

                    /* Update the offsets. */
                    nrows_written += rows_to_read;
                    offset_in_dset += rows_to_read;

                    /* Write the buffer out once full or done with read */
                    if (nrows_written == nrows_in_buf) {
                        /* Calculate the write hyperslab offset. */
                        hsize_t offs[2], lens[2];
                        offs[0] = dset.global_off + out_offset;
                        offs[1] = 0;
                        lens[0] = nrows_written;
                        lens[1] = dset.global_dims[1];

                        /* Write it collectively. */
                        err = write_dataset_2D(dset, offs, lens, buffer);
                        if (err < 0) HANDLE_ERROR("write_dataset_2D")

                        out_offset += nrows_written;
                        nrows_written = 0;
                        num_writes++;
                    }
                }

                /* Done with this input dataset */
                SET_TIMER(ts)
                err = H5Dclose(in_dset_id);
                if (err < 0) HANDLE_ERROR("H5Dclose")
                GET_TIMER(ts, close_in_dsets)
            } /* end of file loop kk */

            /* Flush out the leftover in the I/O buffer */
            if (nrows_written > 0) {
                /* Calculate the write hyperslab */
                hsize_t offs[2], lens[2];
                offs[0] = dset.global_off + out_offset;
                offs[1] = 0;
                lens[0] = nrows_written;
                lens[1] = dset.global_dims[1];

                /* Write it collectively. */
                err = write_dataset_2D(dset, offs, lens, buffer);
                if (err < 0) HANDLE_ERROR("write_dataset_2D")

                out_offset += nrows_written;
                nrows_written = 0;
                num_writes++;
            }

            /* Participate collective writes with zero-sized amount */
            while (num_writes < dset.num_writes) {
                hsize_t offs[2]={0, 0}, lens[2]={0, 0};
                err = write_dataset_2D(dset, offs, lens, buffer);
                if (err < 0) HANDLE_ERROR("write_dataset_2D")
                num_writes++;
            }

            /* Close the output dataset. */
            SET_TIMER(ts)
            err = H5Tclose(dset.type_id);
            if (err < 0) HANDLE_ERROR("H5Tclose")
            err = H5Dclose(dset.out_dset_id);
            if (err < 0) HANDLE_ERROR("H5Dclose");
            dset.out_dset_id = -1;
            GET_TIMER(ts, close_out_dsets)
        } /* end of dataset loop */
    } /* end of group loop */

#if defined PROFILE_HDF5_INTERNAL && PROFILE_HDF5_INTERNAL
    MPI_Barrier(comm);
    H5Vprint();
    H5Vdisable();
#endif

fn_exit:
    return err_exit;
}

/*----< read_dataset2() >----------------------------------------------------*/
int Concatenator::read_dataset2(DSInfo_t &dset,
                                size_t    file_no,
                                hsize_t   all_rows,
                                hsize_t   round_off,
                                hsize_t   round_len,
                                hsize_t   mem_off)
{
    int err_exit=0;
    herr_t err;
    hid_t dset_id, space_id=H5S_ALL, memspace_id=H5S_ALL;
#if defined PROFILE && PROFILE
    double ts, timing=0.0;
#endif
    SET_TIMER(ts)

    dset_id = dset.in_dset_ids[file_no];

    if (all_rows != round_len) { /* read a subset of dataset */
        hsize_t one[2]={1,1}, offs[2], lens[2];

        space_id = H5Dget_space(dset_id);
        if (space_id < 0) HANDLE_ERROR("H5Dget_space")

        offs[0] = round_off;
        offs[1] = 0;
        lens[0] = round_len;
        lens[1] = dset.global_dims[1];
        err = H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offs, NULL,
                                  one, lens);
        if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab");

        memspace_id = H5Screate_simple(2, lens, NULL);
        if (memspace_id < 0) HANDLE_ERROR("H5Screate_simple")
    }

    /* Read dataset. */
    if (dset.type_class == H5T_STRING) {
        /* Below assumes the string datasets stored in the input files are
         * defined of variable length. However, HDF5 1.10.6 and 1.12.0 have
         * not supported parallel I/O for variable-length datasets. We read
         * the string datasets one file per process and then write the
         * concatenated datasets as fixed-length. The string contents are
         * copied over to a fixed-length array.
         */
        htri_t isVLEN = H5Tis_variable_str(dset.type_id);
        char **rdata = (char **) malloc (round_len * sizeof (char *));
        if (space_id == H5S_ALL) {
            space_id = H5Dget_space(dset_id);
            if (space_id < 0) HANDLE_ERROR("H5Dget_space")
        }

        if (isVLEN) {
            /* prepare reading the variable-length dataset */
            hid_t memtype = H5Tcopy(H5T_C_S1);
            if (memtype < 0) HANDLE_ERROR("H5Tcopy");
            err = H5Tset_size(memtype, H5T_VARIABLE);
            if (err < 0) HANDLE_ERROR("H5Tset_size");
            err = H5Dread(dset_id, memtype, memspace_id, space_id, H5P_DEFAULT, rdata);
            if (err < 0) HANDLE_ERROR("H5Dread "+dset.name);
            err = H5Tclose(memtype);
            if (err < 0) HANDLE_ERROR("H5Tclose");

            /* copy rdata[*] to buffer+mem_off */
            hsize_t ii;
            char *ptr = buffer+mem_off;
            for (ii=0; ii<round_len; ii++) {
                strncpy(ptr, rdata[ii], MAX_STR_LEN);
                ptr += MAX_STR_LEN;
            }

            /* Close and release resources internally allocated by HDF5 for
             * variable-length strings.
             */
            err = H5Dvlen_reclaim (dset.type_id, space_id, H5P_DEFAULT, rdata);
            if (err < 0) HANDLE_ERROR("H5Dvlen_reclaim");
        }
        else {
            hsize_t ii;
            hsize_t dims[2];
            int ndims = H5Sget_simple_extent_dims(space_id, dims, NULL);
            assert(ndims == 2);

            rdata[0] = (char*) malloc(round_len * MAX_STR_LEN);
            for (ii=1; ii<round_len; ii++)
                rdata[ii] = rdata[ii-1] + MAX_STR_LEN;

            /* prepare reading the fixed-length dataset */
            hid_t memtype = H5Tcopy(H5T_C_S1);
            if (memtype < 0) HANDLE_ERROR("H5Tcopy");
            err = H5Tset_size(memtype, MAX_STR_LEN);
            if (err < 0) HANDLE_ERROR("H5Tset_size");
            err = H5Dread(dset_id, memtype, memspace_id, space_id, H5P_DEFAULT,
                          rdata[0]);
            if (err < 0) HANDLE_ERROR("H5Dread "+dset.name);
            err = H5Tclose(memtype);
            if (err < 0) HANDLE_ERROR("H5Tclose");

            /* copy rdata[*] to buffer+mem_off */
            char *ptr = buffer+mem_off;
            for (ii=0; ii<round_len; ii++) {
                strncpy(ptr, rdata[ii], MAX_STR_LEN);
                ptr += MAX_STR_LEN;
            }
            free(rdata[0]);
        }
        free (rdata);
    }
    else {
        err = H5Dread(dset_id, dset.type_id, memspace_id, space_id, H5P_DEFAULT,
                      buffer+mem_off);
        if (err < 0) HANDLE_ERROR("H5Dread "+dset.name)
    }

fn_exit:
    if (memspace_id != H5S_ALL)
        H5Sclose(memspace_id);
    if (space_id != H5S_ALL)
        H5Sclose(space_id);

    GET_TIMER(ts, timing)
#if defined PROFILE && PROFILE
    if (dset.global_dims[1] == 1)
        r_1d += timing;
    else
        r_2d += timing;
#endif

    return err_exit;
}

/*----< write_dataset_2D() >-------------------------------------------------*/
int Concatenator::write_dataset_2D(DSInfo_t &dset,
                                   hsize_t  *offs,
                                   hsize_t  *lens,
                                   void     *wbuf)
{
    int err_exit=0;
    herr_t err;
    hid_t space_id=-1, memspace_id=-1;
    hsize_t one[2]={1,1};
#if defined PROFILE && PROFILE
    double ts, timing=0.0;
#endif
    MPI_Barrier(comm);
    SET_TIMER(ts)

    /* Setup memory space */
    memspace_id = H5Screate_simple(2, lens, NULL);
    if (memspace_id < 0) HANDLE_ERROR("H5Screate_simple")

    /* Setup hyperslab file space */
    space_id = H5Dget_space(dset.out_dset_id);
    if (space_id < 0) HANDLE_ERROR("H5Dget_space")
    err = H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offs, NULL, one, lens);
    if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab")

    if (dset.type_class == H5T_STRING) {
        /* fixed-length string dataset */
        hid_t memtype = H5Tcopy(H5T_C_S1);
        if (memtype < 0) HANDLE_ERROR("H5Tcopy");
        err = H5Tset_size(memtype, MAX_STR_LEN);
        if (err < 0) HANDLE_ERROR("H5Tset_size");
        /* Write the data. */
        err = H5Dwrite(dset.out_dset_id, memtype, memspace_id, space_id,
                       dxpl_id, wbuf);
        err = H5Tclose(memtype);
        if (err < 0) HANDLE_ERROR("H5Tclose");
    }
    else {
        /* Write the data. */
        err = H5Dwrite(dset.out_dset_id, dset.type_id, memspace_id, space_id,
                       dxpl_id, wbuf);
    }
    if (err < 0) HANDLE_DSET_ERROR("H5Dwrite", dset.name)

fn_exit:
    if (memspace_id >= 0)
        H5Sclose(memspace_id);
    if (space_id >= 0)
        H5Sclose(space_id);

    GET_TIMER(ts, timing)
#if defined PROFILE && PROFILE
    if (dset.global_dims[1] == 1)
        w_1d += timing;
    else
        w_2d += timing;
#endif

    return err_exit;
}
