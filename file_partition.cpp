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

#if defined PROFILE && PROFILE
    #define SET_TIMER(ts) { ts = MPI_Wtime(); }
    #define GET_TIMER(ts, t) { t += MPI_Wtime() - ts; }
#else
    #define SET_TIMER(ts)
    #define GET_TIMER(ts, t)
#endif

/*----< write_partition_key_datasets() >-------------------------------------*/
int Concatenator::write_partition_key_datasets()
{
    int err_exit=0;
    herr_t err;
    hsize_t ii, jj, kk, start[2]={0,0}, count[2]={1,1};
    hsize_t *seq_offs, *seq_lens, *cnt_offs, *cnt_lens, *extents;
#if defined PROFILE && PROFILE
    double ts;
#endif

    /* some groups may not contain key-base datasets */
    cnt_lens = new hsize_t[num_groups_have_key * 5];
    seq_lens = cnt_lens + num_groups_have_key;
    cnt_offs = seq_lens + num_groups_have_key;
    seq_offs = cnt_offs + num_groups_have_key;
    extents  = seq_offs + num_groups_have_key;

    /* prepare the buffers for MPI communication */
    kk = 0;
    for (ii=0; ii<num_groups; ii++) {
        if (groups[ii].key_base == NULL) continue;
        extents[kk]  = groups[ii].cnt_dset->cnt_len;
        cnt_lens[kk] = groups[ii].cnt_dset->cnt_len;
        cnt_offs[kk] = 0;
        seq_lens[kk] = groups[ii].seq_dset->seq_len;
        seq_offs[kk] = 0;
        kk++;
    }

    /* get starting offsets for writing key.cnt and key.seq in parallel */
    MPI_Exscan(cnt_lens, cnt_offs, num_groups_have_key*2, MPI_HSIZE_T, MPI_SUM,
               comm);
    num_exscan++;

#ifdef DELAY_KEY_DSET_CREATION
    /* get the globally concatenated size of dataset key.cnt */
    MPI_Allreduce(MPI_IN_PLACE, extents, num_groups_have_key, MPI_HSIZE_T,
                  MPI_SUM, comm);
    num_allreduce++;

    kk = 0;
    for (ii=0; ii<num_groups; ii++) {
        if (groups[ii].key_base == NULL) continue;
        assert(groups[ii].seq_dset != NULL);
        assert(groups[ii].cnt_dset != NULL);

        hid_t grp_id;
        grp_id = H5Gopen(output_file_id, groups[ii].name.c_str(), H5P_DEFAULT);
        if (grp_id < 0) HANDLE_ERROR("H5Gopen")

        DSInfo_t *seq_dset = groups[ii].seq_dset;
        DSInfo_t *cnt_dset = groups[ii].cnt_dset;

        /* add offset to local sequence IDs to make them global unique */
        for (jj=0; jj<seq_dset->local_dims[0]; jj++)
            seq_dset->seq_buf[jj] += seq_offs[kk];

        /* create dataset key.seq */
        err = create_dataset(grp_id, *seq_dset, false);
        if (err < 0) HANDLE_ERROR("create_dataset")

        /* write dataset key.seq */
        start[0] = groups[ii].seq_dset->global_off;
        err = write_dataset_2D(*seq_dset, start, seq_dset->local_dims,
                               seq_dset->seq_buf);
        if (err < 0) HANDLE_ERROR("write_dataset_2D")
        delete [] seq_dset->seq_buf;

        err = H5Dclose(seq_dset->out_dset_id);
        if (err < 0) HANDLE_ERROR("H5Dclose");

        /* create dataset key.cnt */
        cnt_dset->global_dims[0] = extents[kk];
        err = create_dataset(grp_id, *cnt_dset, false);
        if (err < 0) HANDLE_ERROR("create_dataset")

        /* write dataset key.cnt */
        start[0] = cnt_offs[kk];
        count[0] = cnt_lens[kk];
        err = write_dataset_2D(*cnt_dset, start, count, cnt_dset->cnt_buf);
        if (err < 0) HANDLE_ERROR("write_dataset_2D")
        delete [] cnt_dset->cnt_buf;

        err = H5Gclose(grp_id);
        if (err < 0) HANDLE_ERROR("H5Gclose")
        kk++;
    }
#else
    if (set_extent) {
        /* get the globally concatenated size of dataset key.cnt */
        MPI_Allreduce(MPI_IN_PLACE, extents, num_groups_have_key, MPI_HSIZE_T,
                      MPI_SUM, comm);
        num_allreduce++;
    }

    kk = 0;
    for (ii=0; ii<num_groups; ii++) {
        if (groups[ii].key_base == NULL) continue;
        assert(groups[ii].seq_dset != NULL);
        assert(groups[ii].cnt_dset != NULL);

        DSInfo_t *seq_dset = groups[ii].seq_dset;
        DSInfo_t *cnt_dset = groups[ii].cnt_dset;

        /* add offset to local sequence IDs to make them global unique */
        for (jj=0; jj<seq_dset->local_dims[0]; jj++)
            seq_dset->seq_buf[jj] += seq_offs[kk];

        /* write dataset key.seq */
        start[0] = groups[ii].seq_dset->global_off;
        err = write_dataset_2D(*seq_dset, start, seq_dset->local_dims,
                               seq_dset->seq_buf);
        if (err < 0) HANDLE_ERROR("write_dataset_2D")
        delete [] seq_dset->seq_buf;

        err = H5Dclose(seq_dset->out_dset_id);
        if (err < 0) HANDLE_ERROR("H5Dclose");

        /* write dataset key.cnt */
        start[0] = cnt_offs[kk];
        count[0] = cnt_lens[kk];
        err = write_dataset_2D(*cnt_dset, start, count, cnt_dset->cnt_buf);
        if (err < 0) HANDLE_ERROR("write_dataset_2D")
        delete [] cnt_dset->cnt_buf;

#ifdef SET_EXTET_INTERLEAVE_WRITE
/* Interleaving H5Dset_extent() with H5Dwrite() in the same loop has shown
 * very expensive on Cori @NERSC. Below we extract calls to H5Dset_extent()
 * in a separate loop.
 */
        /* set dataset key.cnt to its true size. Warning: very expensive !!!
         * For now, we rely on the default fill value of 0 in key.cnt to tell
         * its true size. Luckily, value 0 will not affect the calculation of
         * data partitioning when using key.cnt (which stores counts of
         * consecutive elements sharing the same key.
         */
        if (set_extent) {
            SET_TIMER(ts)
            count[0] = extents[kk];
            err = H5Dset_extent(cnt_dset->out_dset_id, count);
            if (err < 0) HANDLE_ERROR("H5Dset_extent")
            GET_TIMER(ts, t_set_extent)
        }
        err = H5Dclose(cnt_dset->out_dset_id);
        if (err < 0) HANDLE_ERROR("H5Dclose");
#else
        if (!set_extent) {
            err = H5Dclose(cnt_dset->out_dset_id);
            if (err < 0) HANDLE_ERROR("H5Dclose");
        }
#endif
        kk++;
    }

#ifndef SET_EXTET_INTERLEAVE_WRITE
    if (set_extent) {
        /* Set dataset key.cnt to its true size.
         * Note calling H5Dset_extent() is expensive. Calling H5Dset_extent()
         * inside the loop with H5Dwrite() called interleaved is even more
         * expensive (like 10 times slower) !!!
         * Below, we use a separate loop to call H5Dset_extent() altogether
         * after all calls to H5Dwrite() are complete.
         */
        SET_TIMER(ts)
        kk = 0;
        for (ii=0; ii<num_groups; ii++) {
            if (groups[ii].key_base == NULL) continue;
            count[0] = extents[kk];
            err = H5Dset_extent(groups[ii].cnt_dset->out_dset_id, count);
            if (err < 0) HANDLE_ERROR("H5Dset_extent")

            err = H5Dclose(groups[ii].cnt_dset->out_dset_id);
            if (err < 0) HANDLE_ERROR("H5Dclose");
            kk++;
        }
        GET_TIMER(ts, t_set_extent)
    }
#endif
#endif

fn_exit:
    delete []cnt_lens;
    return err_exit;
}

/*----< generate_partition_keys() >------------------------------------------*/
int Concatenator::generate_partition_keys(GrpInfo  &grp,
                                          DSInfo_t &base)
{
    int err_exit=0;
    herr_t err;
    size_t ii, jj, qq, nn;
    void *base_buf;
    unsigned int   *base_ibuf=NULL, *base_iptr=NULL;
    unsigned short *base_sbuf=NULL, *base_sptr=NULL;
    hsize_t start[2];

    assert(grp.key_base != NULL);
    assert(base.global_dims[0] > 0);

    /* allocate buffers for datasets base, key.seq, and key.cnt */
    grp.seq_dset->seq_buf = new long long [base.local_dims[0]];
    grp.cnt_dset->cnt_buf = new long long [base.local_dims[0]];

    if (base.type_size == 4) {
        base_ibuf = new unsigned int [base.local_dims[0]];
        base_iptr = base_ibuf;
        base_buf  = base_ibuf;
    }
    else if (base.type_size == 2) {
        base_sbuf = new unsigned short [base.local_dims[0]];
        base_sptr = base_sbuf;
        base_buf  = base_sbuf;
    }
    else {
        delete [] grp.seq_dset->seq_buf;
        delete [] grp.cnt_dset->cnt_buf;
        printf("Error in %s line %d: unsupported datatype size %zd for key base dataset\n",
               __func__,__LINE__,base.type_size);
        return -1;
    }

    qq = 0;
    nn = 0;
    grp.seq_dset->seq_len = 0;
    for (ii=0; ii<base.in_dset_ids.size(); ii++) { /* iterate all input files */
        hid_t in_dset_id = base.in_dset_ids[ii];
        hsize_t in_dset_nrows = base.in_dim0[ii];

        /* read the whole dataset 'base' */
        if (base.type_size == 4)
            err = H5Dread(in_dset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, base_iptr);
        else
            err = H5Dread(in_dset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, base_sptr);
        if (err < 0) HANDLE_ERROR("H5Dread")

        /* close input dataset base, as it has been entirely read */
        err = H5Dclose(in_dset_id);
        if (err < 0) HANDLE_ERROR("H5Dclose");

        /* use contents of 'base' to calculate contents of dataset 'key.seq'
         * dataset 'key.cnt' */
        grp.seq_dset->seq_buf[qq++] = grp.seq_dset->seq_len;
        grp.cnt_dset->cnt_buf[nn] = 1;
        if (base.type_size == 4) {
            for (jj=1; jj<in_dset_nrows; jj++) {
                if (base_iptr[jj] == base_iptr[jj-1]) /* repeated ID */
                    grp.cnt_dset->cnt_buf[nn]++;      /* increment count */
                else {
                    grp.seq_dset->seq_len++; /* increment no. of unique IDs */
                    grp.cnt_dset->cnt_buf[++nn] = 1; /* a new unique ID */
                }
                /* assign the unique ID */
                grp.seq_dset->seq_buf[qq++] = grp.seq_dset->seq_len;
            }
            base_iptr += in_dset_nrows;
        }
        else if (base.type_size == 2) {
            for (jj=1; jj<in_dset_nrows; jj++) {
                if (base_sptr[jj] == base_sptr[jj-1]) /* repeated ID */
                    grp.cnt_dset->cnt_buf[nn]++;      /* increment count */
                else {
                    grp.seq_dset->seq_len++; /* increment no. of unique IDs */
                    grp.cnt_dset->cnt_buf[++nn] = 1; /* a new unique ID */
                }
                /* assign the unique ID */
                grp.seq_dset->seq_buf[qq++] = grp.seq_dset->seq_len;
            }
            base_sptr += in_dset_nrows;
        }
        nn++;                    /* increment one for next file */
        grp.seq_dset->seq_len++; /* increment one for next file */
    }
    grp.cnt_dset->cnt_len = nn; /* number of elements in cnt_buf */
    /* Now, seq_len is the number of unique IDs in seq_buf */

    /* write dataset base */
    start[0] = base.global_off;
    start[1] = 0;
    err = write_dataset_2D(base, start, base.local_dims, base_buf);
    if (err < 0) HANDLE_ERROR("write_dataset_2D")

    /* close datasets key base */
    err = H5Tclose(base.type_id);
    if (err < 0) HANDLE_ERROR("H5Tclose")
    err = H5Dclose(base.out_dset_id);
    if (err < 0) HANDLE_ERROR("H5Dclose");
    base.out_dset_id = -1;

fn_exit:
    if (base.type_size == 4)
        delete []base_ibuf;
    else if (base.type_size == 2)
        delete []base_sbuf;
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

    /* Traverse over all the groups and datasets. */
    for (ii=0; ii<num_groups; ii++) {
        size_t num_dsets = groups[ii].num_dsets;
        bool check_key_base = false;
        if (groups[ii].key_base != NULL) {
            check_key_base = !process_large_dsets;
            num_dsets -= 2;
            /* key.seq and key.cnt are the last 2 added. In this loop, their
             * contents are generated from the key base dataset and cached in
             * memory.  Later, they will be written to output file later in
             * write_partition_key_datasets().
             */
        }

        for (jj=0; jj<num_dsets; jj++) {
            DSInfo_t &dset = groups[ii].dsets[jj];
            assert(dset.global_dims[0] > 0);

            if (process_large_dsets) {
                if (dset.global_dims[1] == 1)
                    continue;
            }
            else if (dset.global_dims[1] > 1)
                continue;

            if (check_key_base && dset.is_key_base) {
                /* Use 'part_key_base' to generate the partitioning key
                 * datasets and write base dataset to output file.
                 */
                err = generate_partition_keys(groups[ii], dset);
                if (err < 0) HANDLE_ERROR("generate_partition_keys")
                check_key_base = false;
                continue;
            }

            hsize_t row_size = dset.global_dims[1] * dset.type_size;
            hsize_t nrows_in_buf = io_buffer_size / row_size;
            hsize_t nrows_written = 0;
            hsize_t out_offset = 0;
            size_t num_writes = 0;

            /* iterate input files assigned to this process */
            for (kk=0; kk<dset.in_dset_ids.size(); kk++) {
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
                    if (err < 0) HANDLE_ERROR("read_dataset2")

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
        hsize_t offs[2], lens[2];

        space_id = H5Dget_space(dset_id);
        if (space_id < 0) HANDLE_ERROR("H5Dget_space")

        offs[0] = round_off;
        offs[1] = 0;
        lens[0] = round_len;
        lens[1] = dset.global_dims[1];
        err = H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offs, NULL, lens,
                                  NULL);
        if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab");

        memspace_id = H5Screate_simple(2, lens, NULL);
        if (memspace_id < 0) HANDLE_ERROR("H5Screate_simple")
    }

    /* Read dataset. */
    err = H5Dread(dset_id, dset.type_id, memspace_id, space_id, H5P_DEFAULT,
                  buffer+mem_off);
    if (err < 0) HANDLE_ERROR("H5Dread")

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
    err = H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offs, NULL, lens, NULL);
    if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab")

    /* Write the data. */
    err = H5Dwrite(dset.out_dset_id, dset.type_id, memspace_id,
                   space_id, dxpl_id, wbuf);
    if (err < 0) HANDLE_ERROR("H5Dwrite")

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
