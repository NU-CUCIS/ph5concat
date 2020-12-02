/*
 * Copyright (C) 2020, Northwestern University and Fermi National Accelerator Laboratory
 * See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h> /* getopt() */
#include <limits.h> /* INT_MAX */
#include <assert.h>
#include <mpi.h>
#include "hdf5.h"
#include "zlib.h"

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

#define LINE_SIZE 64

static int verbose;

typedef struct {
    int      nDatasets;  /* number of datasets in this group */
    int     *nChunks;    /* number of chunks in each dataset */
    char   **dset_names; /* string names of datasets */
    void   **buf;        /* read buffers */
    double  *read_t;     /* read timings */
} NOvA_group;

/*----< read_dataset_names() >-----------------------------------------------*/
static
int read_dataset_names(int          rank,
                       int          seq_opt,
                       int          dset_opt,
                       char        *listfile,
                       NOvA_group **gList,           /* OUT */
                       int         *nTotalDatasets)  /* OUT */
{
    FILE *fptr;
    int d, j, g, len, nGroups, nDatasets, nDsetGrp, maxDsetGrp;
    char line[LINE_SIZE], gname[LINE_SIZE], *cur_gname;

    fptr = fopen(listfile, "r");
    if (fptr == NULL) {
        fprintf(stderr,"%d: Error: fail to open file %s (%s)\n",
               rank,listfile,strerror(errno));
        assert(0);
    }

    /* check number of datasets listed in the file */
    nGroups = 0;     /* number of groups */
    maxDsetGrp = 0;  /* max number of datasets among all groups */
    nDsetGrp = 0;    /* number of dataset in the current group */
    gname[0] = '\0';
    while (fgets(line, LINE_SIZE, fptr) != NULL) {
        if (line[0] == '\n') continue;   /* skip empty line */
        if (line[0] == '#') continue;    /* skip comment line */
        len = strlen(line);
        if (line[len-1] == '\n') line[--len] = '\0';
        while (len > 0 && line[len-1] == ' ') line[--len] = '\0';
        if (len == 0) continue; /* skip blank line */

        /* retrieve group name */
        cur_gname = strtok(line, "/");
        if (strcmp(gname, cur_gname)) { /* entering a new group */
            strcpy(gname, cur_gname);
            if (nDsetGrp > maxDsetGrp) maxDsetGrp = nDsetGrp;
            nDsetGrp = 0;
            nGroups++;
        }
        nDsetGrp++;
    }

    /* allocate an array of group objects */
    *gList = (NOvA_group*) malloc(nGroups * sizeof(NOvA_group));

    /* read dataset names */
    rewind(fptr);

    nGroups = -1;
    nDatasets = 0;
    gname[0] = '\0';
    while (fgets(line, LINE_SIZE, fptr) != NULL) {
        char name[LINE_SIZE];

        if (line[0] == '\n') continue;   /* skip empty line */
        if (line[0] == '#') continue;    /* skip comment line */
        len = strlen(line);
        if (line[len-1] == '\n') line[--len] = '\0';
        while (len > 0 && line[len-1] == ' ') line[--len] = '\0';
        if (len == 0) continue; /* skip blank line */

        /* now len is the true string length of line */

        /* retrieve group name */
        strcpy(name, line);
        cur_gname = strtok(name, "/");
        if (strcmp(gname, cur_gname)) { /* entering a new group */
            strcpy(gname, cur_gname);
            if (nGroups >= 0) (*gList)[nGroups].nDatasets = nDatasets;
            nGroups++;
            nDatasets = 0;
            /* allocate space to store names of datasets in this group */
            (*gList)[nGroups].dset_names = (char**) malloc(maxDsetGrp * sizeof(char*));
        }

        if (!strcmp(strtok(NULL, "/"), "evt.seq")) {
            /* add dataset evt.seq to first in the group */
            for (j=nDatasets; j>0; j--)
                (*gList)[nGroups].dset_names[j] = (*gList)[nGroups].dset_names[j-1];
            (*gList)[nGroups].dset_names[0] = (char*) malloc(len + 1);
            strcpy((*gList)[nGroups].dset_names[0], line);
        }
        else {
            (*gList)[nGroups].dset_names[nDatasets] = (char*) malloc(len + 1);
            strcpy((*gList)[nGroups].dset_names[nDatasets], line);
        }
        nDatasets++;
    }
    fclose(fptr);
    (*gList)[nGroups].nDatasets = nDatasets;
    nGroups++;

    *nTotalDatasets = 0;
    for (g=0; g<nGroups; g++) {
        *nTotalDatasets += (*gList)[g].nDatasets;
        /* allocate space to store number of chunks */
        (*gList)[g].nChunks = (int*) calloc((*gList)[g].nDatasets, sizeof(int));
        /* allocate array of read buffer pointers */
        (*gList)[g].buf = (void**) malloc((*gList)[g].nDatasets * sizeof(void*));
        /* allocate space to store read timings */
        (*gList)[g].read_t = (double*) malloc((*gList)[g].nDatasets * sizeof(double));
    }

    if (verbose && rank == 0) {
        for (g=0; g<nGroups; g++)
            for (d=0; d<(*gList)[g].nDatasets; d++)
                printf("group[%d].dataset[%d] %s\n", g, d, (*gList)[g].dset_names[d]);
    }

    if (rank == 0) {
        printf("Number of datasets to read = %d\n", *nTotalDatasets);
        printf("Number of groups = %d\n",nGroups);
        printf("Maximum number of datasets among groups = %d\n",maxDsetGrp);
        if (seq_opt == 0)
            printf("Read evt.seq method: root process reads and broadcasts\n");
        else if (seq_opt == 1)
            printf("Read evt.seq method: all processes read the entire evt.seq collectively\n");
        else
            printf("Read evt.seq method: root process reads evt.seq and scatters boundaries\n");

        if (dset_opt == 0)
            printf("Read datasets method: H5Dread\n");
        else
            printf("Read datasets method: MPI_File_read_all\n");
        printf("----------------------------------------------------\n");
    }
    fflush(stdout);

    return nGroups;
}

/*----< calculate_starts_ends() >--------------------------------------------*/
static
int calculate_starts_ends(hid_t    fd,
                          int      nprocs,
                          int      rank,
                          hsize_t *starts,  /* OUT */
                          hsize_t *ends)    /* OUT */
{
    herr_t err;
    hid_t seq, fspace;
    hsize_t j, ndims, dims[5], seq_len, my_count;

    /* open dataset '/spill/evt.seq' */
    seq = H5Dopen2(fd, "/spill/evt.seq", H5P_DEFAULT); assert(seq >= 0);

    /* inquire the size of '/spill/evt.seq' */
    fspace = H5Dget_space(seq); assert(fspace >= 0);
    ndims = H5Sget_simple_extent_ndims(fspace); assert(ndims == 2);
    err = H5Sget_simple_extent_dims(fspace, dims, NULL); assert(err >= 0);
    /* 2nd dimension of /spill/evt.seq is always 1 */
    assert(dims[1] == 1);

    seq_len = dims[0];
    if (verbose && rank == 0)
        printf("Size of /spill/evt.seq is %zd\n", (size_t)seq_len);

    /* calculate the range of event IDs assigned to this process. We name it
     * 'partition domain'. The domain ranges from starts[rank] to ends[rank].
     */
    my_count = seq_len / nprocs;
    for (j=0; j<seq_len % nprocs; j++) {
        starts[j] = my_count * j + j;
        ends[j] = starts[j] + my_count;
    }
    for (; j<nprocs; j++) {
        starts[j] = my_count * j + seq_len % nprocs;
        ends[j] = starts[j] + my_count - 1;
    }

    /* Note contents of /spill/evt.seq start with 0 and increment by 1, so
     * there is no need to read the contents of it. Check the contents only
     * in verbose mode.
     */
    if (verbose) {
        int64_t v, *seq_buf;
        printf("%d: /spill/evt.seq len=%zd starts[rank]=%zd ends[rank]=%zd\n",
               rank, (size_t)seq_len, (size_t)starts[rank], (size_t)ends[rank]);
        seq_buf = (int64_t*) malloc(seq_len * sizeof(int64_t));

        err = H5Dread(seq, H5T_STD_I64LE, fspace, fspace, H5P_DEFAULT, seq_buf);
        assert(err >= 0);

        for (v=0; v<seq_len; v++)
            if (seq_buf[v] != v) {
                printf("Error: /spill/evt.seq[%"PRId64"] expect %"PRId64" but got %"PRId64"\n",
                       v, v, seq_buf[v]);
                break;
            }
        free(seq_buf);
    }

    err = H5Sclose(fspace); assert(err >= 0);
    err = H5Dclose(seq); assert(err >= 0);
    return 1;
}

/* return the smallest index i, such that base[i] >= key */
static
size_t binary_search_min(int64_t   key,
                         int64_t *base,
                         size_t   nmemb)
{
    size_t low=0, high=nmemb;
    while (low != high) {
        size_t mid = (low + high) / 2;
        if (base[mid] < key)
            low = mid + 1;
        else
            high = mid;
    }
    return low;
}

/* return the largest index i, such that base[i] <= key */
static
size_t binary_search_max(int64_t  key,
                         int64_t *base,
                         size_t   nmemb)
{
    size_t low=0, high=nmemb;
    while (low != high) {
        size_t mid = (low + high) / 2;
        if (base[mid] <= key)
            low = mid + 1;
        else
            high = mid;
    }
    return (low-1);
}

/*----< read_evt_seq() >-----------------------------------------------------*/
static
int read_evt_seq(hid_t          fd,
                 NOvA_group    *group,
                 int            nprocs,
                 int            rank,
                 int            profile,
                 int            seq_opt,
                 hid_t          xfer_plist,
                 const hsize_t *starts,  /* IN:  [nprocs] */
                 const hsize_t *ends,    /* IN:  [nprocs] */
                 size_t        *low,     /* OUT */
                 size_t        *high)    /* OUT */
{
    herr_t err;
    hsize_t ndims, dims[5];
    int64_t *seq_buf=NULL;

    /* open dataset 'evt.seq', first in the group */
    hid_t seq = H5Dopen2(fd, group->dset_names[0], H5P_DEFAULT); assert(seq >= 0);

    /* inquire dimension sizes of 'evt.seq' in this group */
    hid_t fspace = H5Dget_space(seq); assert(fspace >= 0);
    ndims = H5Sget_simple_extent_ndims(fspace); assert(ndims == 2);
    err = H5Sget_simple_extent_dims(fspace, dims, NULL); assert(err>=0);
    assert(dims[1] == 1);
    err = H5Sclose(fspace); assert(err >= 0);

    /* data type of evt.seq is 64-bit integer */
    hid_t dtype = H5Dget_type(seq); assert(dtype >= 0);
    size_t dtype_size = H5Tget_size(dtype); assert(dtype_size > 0);
    err = H5Tclose(dtype); assert(err >= 0);

    if (seq_opt == 2) {
        /* rank 0 reads evt.seq, calculates lows and highs for all processes
         * and scatters
         */
        long long low_high[2], *bounds;
        if (rank == 0) {
            int j, k;
            seq_buf = (int64_t*) malloc(dims[0] * dims[1] * dtype_size);
            err = H5Dread(seq, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, seq_buf);
            assert(err >= 0);
            bounds = (long long*) malloc(nprocs * 2 * sizeof(long long));
            for (j=0, k=0; j<nprocs; j++) {
                bounds[k++] = binary_search_min(starts[j], seq_buf, dims[0]);
                bounds[k++] = binary_search_max(ends[j],   seq_buf, dims[0]);
            }
        }
        MPI_Scatter(bounds, 2, MPI_LONG_LONG, low_high, 2, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        *low  = low_high[0];
        *high = low_high[1];
        if (rank == 0) free(bounds);
    }
    else {
        if (seq_opt == 0) {
            /* rank 0 reads evt.seq and broadcasts it */
            seq_buf = (int64_t*) malloc(dims[0] * dims[1] * dtype_size);
            if (rank == 0) {
                err = H5Dread(seq, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                              seq_buf);
                assert(err >= 0);
            }
            MPI_Bcast(seq_buf, dims[0], MPI_LONG_LONG, 0, MPI_COMM_WORLD);
        }
        else {
            /* collective-read-the-whole-dataset is bad */
            seq_buf = (int64_t*) malloc(dims[0] * dims[1] * dtype_size);
            err = H5Dread(seq, H5T_STD_I64LE, H5S_ALL, H5S_ALL, xfer_plist, seq_buf);
            /* all-processes-read-the-whole-dataset is even worse */
            // err = H5Dread(seq, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, seq_buf);
            assert(err >= 0);
        }
        /* find the array index range from 'low' to 'high' that falls into
         * this process's partition domain.
         */
        *low  = binary_search_min(starts[rank], seq_buf, dims[0]);
        *high = binary_search_max(ends[rank],   seq_buf, dims[0]);
    }

    err = H5Dclose(seq); assert(err >= 0);

    if (verbose) {
        if (rank == 0)
            printf("%s[0]=%"PRId64" seq[len-1=%zd]=%"PRId64"\n",
                   group->dset_names[0],seq_buf[0],(size_t)dims[0]-1,seq_buf[dims[0]-1]);

        printf("%d: %s starts[rank]=%zd ends[rank]=%zd low=%zd high=%zd\n",
               rank, group->dset_names[0], (size_t)starts[rank],(size_t)ends[rank], *low, *high);
    }

    if (seq_buf != NULL) free(seq_buf);

    return 1;
}

/*----< read_hdf5() >--------------------------------------------------------*/
/* call H5Dread to read datasets one at a time */
static
int read_hdf5(hid_t       fd,
              int         rank,
              NOvA_group *group,
              hid_t       xfer_plist,
              size_t      low,
              size_t      high,
              int         profile)
{
    int d;
    herr_t  err;
    hid_t   dset, fspace, mspace, dtype;
    hsize_t ndims, dims[5], start[2], count[2];
    size_t read_len;

    /* iterate all the remaining datasets in this group */
    for (d=1; d<group->nDatasets; d++) {
        hsize_t one[2]={1,1};

        if (verbose && rank == 0) printf("dataset %s\n", group->dset_names[d]);

        /* open dataset */
        dset = H5Dopen2(fd, group->dset_names[d], H5P_DEFAULT); assert(dset >= 0);

        /* inquire dimension sizes of dset */
        fspace = H5Dget_space(dset); assert(fspace >= 0);
        ndims = H5Sget_simple_extent_ndims(fspace); assert(ndims == 2);
        err = H5Sget_simple_extent_dims(fspace, dims, NULL); assert(err>=0);

        /* set subarray/hyperslab access */
        start[0] = low;
        start[1] = 0;
        count[0] = high - low + 1;
        count[1] = dims[1];
        err = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, one,
                                  count); assert(err>=0);
        mspace = H5Screate_simple(2, count, NULL); assert(mspace>=0);

        /* get data type and size */
        dtype = H5Dget_type(dset); assert(dtype >= 0);
        size_t dtype_size = H5Tget_size(dtype); assert(dtype_size > 0);

        /* calculate read size in bytes and allocate read buffers */
        read_len = count[0] * count[1] * dtype_size;
        group->buf[d] = (void*) malloc(read_len);

        /* collectively read dataset d's contents */
        if (verbose)
            printf("%d: READ read_len=%zd dataset %s\n",rank,read_len,group->dset_names[d]);

        if (profile == 2) group->read_t[d]=MPI_Wtime();
        err = H5Dread(dset, dtype, mspace, fspace, xfer_plist, group->buf[d]);
        assert(err >= 0);

        err = H5Tclose(dtype);  assert(err >= 0);
        err = H5Sclose(mspace); assert(err >= 0);
        err = H5Sclose(fspace); assert(err >= 0);
        err = H5Dclose(dset);   assert(err >= 0);

        if (profile == 2) group->read_t[d]=MPI_Wtime()-group->read_t[d];
    }

    return 1;
}

typedef struct {
    MPI_Aint block_dsp;
    int      block_len;
    int      block_idx;
} off_len;

static int
off_compare(const void *a, const void *b)
{
    if (((off_len*)a)->block_dsp > ((off_len*)b)->block_dsp) return  1;
    if (((off_len*)a)->block_dsp < ((off_len*)b)->block_dsp) return -1;
    return 0;
}

/*----< read_mpio() >--------------------------------------------------------*/
/* call MPI_File_read_all to read raw chunks of each dataset ans decompress */
static
int read_mpio(hid_t       fd,
              int         rank,
              NOvA_group *group,
              MPI_File    fh,
              size_t      low,
              size_t      high,
              int         profile)
{
    int j, d, mpi_err;
    herr_t  err;
    hid_t   dset;
    size_t read_len;
    unsigned char *whole_chunk;
    MPI_Status status;

    for (d=1; d<group->nDatasets; d++) {
        hsize_t dims[5], chunk_dims[2];
        /* open dataset */
        dset = H5Dopen2(fd, group->dset_names[d], H5P_DEFAULT); assert(dset >= 0);

        /* find metadata of all the chunks of this dataset */

        /* inquire dimension sizes of dset */
        hid_t fspace = H5Dget_space(dset); assert(fspace >= 0);
        err = H5Sget_simple_extent_dims(fspace, dims, NULL); assert(err>=0);
        err = H5Sclose(fspace); assert(err >= 0);

        /* get data type and size */
        hid_t dtype = H5Dget_type(dset); assert(dtype >= 0);
        size_t dtype_size = H5Tget_size(dtype); assert(dtype_size > 0);
        err = H5Tclose(dtype); assert(err >= 0);

        /* calculate buffer read size in bytes and allocate read buffer */
        read_len = (high - low + 1) * dims[1] * dtype_size;
        group->buf[d] = (void*) malloc(read_len);

        /* inquire chunk size along each dimension */
        hid_t chunk_plist = H5Dget_create_plist(dset); assert(chunk_plist >= 0);
        int chunk_ndims = H5Pget_chunk(chunk_plist, 2, chunk_dims);
        assert(chunk_ndims == 2);
        assert(chunk_dims[1] == 1);
        err = H5Pclose(chunk_plist); assert(err>=0);

        hsize_t nChunks, offset[2];

        /* find the number of chunks to be read by this process */
        nChunks = (high / chunk_dims[0]) - (low / chunk_dims[0]) + 1;

        /* Note file offsets of chunks may not follow the increasing order of
         * chunk IDs read by this process. We must sort the offsets before
         * creating a file type. Construct an array of off-len-indx for such
         * sorting.
         */
        off_len *disp_indx = (off_len*) malloc(nChunks * sizeof(off_len));

        /* calculate the logical position of chunk's first element. See
         * https://hdf5.io/develop/group___h5_d.html#ga408a49c6ec59c5b65ce4c791f8d26cb0
         */
        offset[0] = (low / chunk_dims[0]) * chunk_dims[0];
        offset[1] = 0;
        read_len = 0;
        for (j=0; j<nChunks; j++) {
            hsize_t size;
            haddr_t addr;
            err = H5Dget_chunk_info_by_coord(dset, offset, NULL, &addr, &size); assert(err>=0);
            disp_indx[j].block_dsp = (MPI_Aint)addr;  /* chunk's file offset */
            disp_indx[j].block_len = (int)size;       /* compressed chunk size */
            disp_indx[j].block_idx = j;               /* chunk ID to be read by this process */
            read_len += size;
            offset[0] += chunk_dims[0];
        }
        err = H5Dclose(dset); assert(err >= 0);

        /* sort chunk offsets into an increasing order, as MPI-IO requires that
         * for file view
         */
        qsort(disp_indx, nChunks, sizeof(off_len), off_compare);

        /* allocate array_of_blocklengths[] and array_of_displacements[] */
        MPI_Aint *chunk_dsps = (MPI_Aint*) malloc(nChunks * sizeof(MPI_Aint));
        int *chunk_lens = (int*) malloc(nChunks * sizeof(int));
        for (j=0; j<nChunks; j++) {
            chunk_lens[j] = disp_indx[j].block_len;
            chunk_dsps[j] = disp_indx[j].block_dsp;
        }

        /* create the filetype */
        MPI_Datatype ftype;
        mpi_err = MPI_Type_create_hindexed(nChunks, chunk_lens, chunk_dsps, MPI_BYTE, &ftype);
        assert(mpi_err == MPI_SUCCESS);

        mpi_err = MPI_Type_commit(&ftype);
        assert(mpi_err == MPI_SUCCESS);
        free(chunk_lens);
        free(chunk_dsps);

        mpi_err = MPI_File_set_view(fh, 0, MPI_BYTE, ftype, "native", MPI_INFO_NULL);
        assert(mpi_err == MPI_SUCCESS);

        mpi_err = MPI_Type_free(&ftype);
        assert(mpi_err == MPI_SUCCESS);

        /* allocate buffer for reading the compressed chunks all at once */
        unsigned char *chunk_buf, *chunk_buf_ptr;
        chunk_buf = (unsigned char*) malloc(read_len);
        chunk_buf_ptr = chunk_buf;

        /* collective read */
        mpi_err = MPI_File_read_all(fh, chunk_buf, read_len, MPI_BYTE, &status);
        assert(mpi_err == MPI_SUCCESS);

        /* decompress each chunk into group->buf[d] */
        size_t whole_chunk_size = chunk_dims[0] * chunk_dims[1] * dtype_size;
        whole_chunk = (unsigned char*) malloc(whole_chunk_size);
        for (j=0; j<nChunks; j++) {
            int ret;
            z_stream z_strm;
            z_strm.zalloc    = Z_NULL;
            z_strm.zfree     = Z_NULL;
            z_strm.opaque    = Z_NULL;
            z_strm.avail_in  = disp_indx[j].block_len;
            z_strm.next_in   = chunk_buf_ptr;
            z_strm.avail_out = whole_chunk_size;
            z_strm.next_out  = whole_chunk;
            ret = inflateInit(&z_strm); assert(ret == Z_OK);
            ret = inflate(&z_strm, Z_SYNC_FLUSH); assert(ret == Z_OK || Z_STREAM_END);
            ret = inflateEnd(&z_strm); assert(ret == Z_OK);

            /* copy requested data to user buffer */
            if (disp_indx[j].block_idx == 0) { /* first chunk */
                size_t off = low % chunk_dims[0];
                size_t len;
                if (nChunks == 1)
                    len = high - low + 1;
                else
                    len = chunk_dims[0] - off;
                len *= chunk_dims[1] * dtype_size;
                off *= chunk_dims[1] * dtype_size;
                memcpy(group->buf[d], whole_chunk + off, len);
            }
            else if (disp_indx[j].block_idx == nChunks - 1) { /* last chunk */
                size_t len = (high+1) % chunk_dims[0];
                if (len == 0) len = chunk_dims[0];
                size_t off = high + 1 - len - low;
                off *= chunk_dims[1] * dtype_size;
                len *= chunk_dims[1] * dtype_size;
                memcpy(group->buf[d] + off, whole_chunk, len);
            }
            else { /* middle chunk, copy the full chunk */
                size_t off = chunk_dims[0] - low % chunk_dims[0];
                off += (disp_indx[j].block_idx - 1) * chunk_dims[0];
                off *= chunk_dims[1] * dtype_size;
                memcpy(group->buf[d] + off, whole_chunk, whole_chunk_size);
            }
            chunk_buf_ptr += disp_indx[j].block_len;
        }

        free(whole_chunk);
        free(chunk_buf);
        free(disp_indx);
    }
    return 1;
}

/*----< chunk_statistics() >-------------------------------------------------*/
static
int chunk_statistics(MPI_Comm    comm,
                     const char *infile,
                     NOvA_group *groups,
                     int         nGroups,
                     hsize_t    *starts,
                     hsize_t    *ends,
                     int         seq_opt,
                     int         profile)
{
    herr_t  err;
    hid_t   fd, dset;
    hsize_t ndims, dims[2];

    int g, d, j, k, nprocs, rank, nDatasets;
    int nchunks_shared=0, max_shared_chunks=0;
    int max_nchunks_read=0, min_nchunks_read=INT_MAX, total_nchunks=0;
    long long aggr_nchunks_read, my_nchunks_read=0;
    long long all_dset_size, all_evt_seq_size=0;
    long long *bounds, maxRead=0, minRead=LONG_MAX;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (rank == 0)
        bounds = (long long*) malloc(nprocs * 2 * sizeof(long long));

    /* collect statistics describing chunk contention */
    fd = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
    all_dset_size = 0;
    nDatasets = 0;
    for (g=0; g<nGroups; g++) {
        size_t low=0, high=0;

        nDatasets += groups[g].nDatasets;
        if (!strcmp(groups[g].dset_names[0], "/spill/evt.seq")) {
            groups[g].nChunks[0] = 0;
            if (rank == 0) {
                for (j=0, k=0; j<nprocs; j++) {
                    bounds[k++] = starts[j];
                    bounds[k++] = ends[j];
                }
            }
            low  = starts[rank];
            high = ends[rank];
        }
        else {
            /* rank 0 reads evt.seq, calculates lows and highs for all processes
             * and scatters them.
             */
            long long low_high[2];

            /* open dataset 'evt.seq', first in the group */
            hid_t seq = H5Dopen2(fd, groups[g].dset_names[0], H5P_DEFAULT); assert(seq >= 0);

            /* inquire the size of 'evt.seq' in this group */
            hid_t fspace = H5Dget_space(seq); assert(fspace >= 0);
            ndims = H5Sget_simple_extent_ndims(fspace); assert(ndims == 2);
            err = H5Sget_simple_extent_dims(fspace, dims, NULL); assert(err>=0);
            err = H5Sclose(fspace); assert(err >= 0);

            /* inquire chunk size along each dimension */
            hid_t chunk_plist;
            hsize_t chunk_dims[2];
            chunk_plist = H5Dget_create_plist(seq); assert(chunk_plist >= 0);
            int chunk_ndims = H5Pget_chunk(chunk_plist, 2, chunk_dims);
            assert(chunk_ndims == 2);
            assert(chunk_dims[1] == 1);
            err = H5Pclose(chunk_plist); assert(err>=0);
            groups[g].nChunks[0] = dims[0] / chunk_dims[0];
            if (dims[0] % chunk_dims[0]) groups[g].nChunks[0]++;

            if (rank == 0) {
                int64_t *seq_buf=NULL;

                my_nchunks_read += groups[g].nChunks[0];
                total_nchunks += groups[g].nChunks[0];

                /* data type of evt.seq is 64-bit integer */
                hid_t dtype = H5Dget_type(seq); assert(dtype >= 0);
                size_t dtype_size = H5Tget_size(dtype); assert(dtype_size > 0);
                all_dset_size    += dims[0] * dims[1] * dtype_size;
                all_evt_seq_size += dims[0] * dims[1] * dtype_size;
                err = H5Tclose(dtype); assert(err >= 0);

                seq_buf = (int64_t*) malloc(dims[0] * dims[1] * dtype_size);
                err = H5Dread(seq, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, seq_buf); assert(err >= 0);

                for (j=0, k=0; j<nprocs; j++) {
                    bounds[k++] = binary_search_min(starts[j], seq_buf, dims[0]);
                    bounds[k++] = binary_search_max(ends[j],   seq_buf, dims[0]);
                }
                free(seq_buf);

                if (seq_opt == 1) nchunks_shared += groups[g].nChunks[0];
            }
            else if (seq_opt == 1)
                my_nchunks_read += groups[g].nChunks[0];

            err = H5Dclose(seq); assert(err >= 0);

            MPI_Scatter(bounds, 2, MPI_LONG_LONG, low_high, 2, MPI_LONG_LONG, 0, comm);
            low  = low_high[0];
            high = low_high[1];
        }

        for (d=1; d<groups[g].nDatasets; d++) {
            hsize_t chunk_dims[2];

            /* open dataset */
            dset = H5Dopen2(fd, groups[g].dset_names[d], H5P_DEFAULT); assert(dset >= 0);

            /* inquire dimension sizes of dset */
            hid_t fspace = H5Dget_space(dset); assert(fspace >= 0);
            ndims = H5Sget_simple_extent_ndims(fspace); assert(ndims == 2);
            err = H5Sget_simple_extent_dims(fspace, dims, NULL); assert(err>=0);

            hid_t dtype = H5Dget_type(dset); assert(dtype >= 0);
            size_t dtype_size = H5Tget_size(dtype); assert(dtype_size > 0);
            err = H5Tclose(dtype); assert(err >= 0);

            /* inquire chunk size along each dimension */
            hid_t chunk_plist;
            chunk_plist = H5Dget_create_plist(dset); assert(chunk_plist >= 0);
            int chunk_ndims = H5Pget_chunk(chunk_plist, 2, chunk_dims);
            assert(chunk_ndims == 2);
            assert(chunk_dims[1] == 1);
            err = H5Pclose(chunk_plist); assert(err>=0);
            err = H5Dclose(dset); assert(err >= 0);

            /* calculate number of chunks of this dataset */
            groups[g].nChunks[d] = dims[0] / chunk_dims[0];
            if (dims[0] % chunk_dims[0]) groups[g].nChunks[d]++;
            total_nchunks += groups[g].nChunks[d];

            /* calculate number of chunks read by this process */
            int nchunks = (high / chunk_dims[0]) - (low / chunk_dims[0]) + 1;
            my_nchunks_read += nchunks;
            max_nchunks_read = MAX(max_nchunks_read, nchunks);
            min_nchunks_read = MIN(min_nchunks_read, nchunks);

            /* calculate nchunks_shared and max_shared_chunks */
            if (rank == 0) {
                /* bounds[] is populated on root only */
                int prev_chunk_id=-1, max_chunks=1;
                int k=0, chunk_id, prev_chunk=-1, prev_prev_chunk=-1;
                for (j=0; j<nprocs; j++) {
                    size_t read_len = (bounds[k+1] - bounds[k] + 1) * dims[1] * dtype_size;
                    maxRead = MAX(maxRead, read_len);
                    minRead = MIN(minRead, read_len);
                    chunk_id = bounds[k++] / chunk_dims[0];
                    if (chunk_id == prev_chunk) {
                        if (chunk_id > prev_prev_chunk) /* count only once */
                            nchunks_shared++;
                        prev_prev_chunk = prev_chunk;
                        prev_chunk = chunk_id;
                    }
                    else
                        prev_prev_chunk = -1;
                    prev_chunk = bounds[k++] / chunk_dims[0];

                    if (chunk_id == prev_chunk_id)
                        max_chunks++;
                    else {
                        max_shared_chunks = MAX(max_shared_chunks, max_chunks);
                        max_chunks = 1;
                        prev_chunk_id = prev_chunk;
                    }
                }
                max_shared_chunks = MAX(max_shared_chunks, max_chunks);
                all_dset_size += dims[0] * dims[1] * dtype_size;
            }
        }
    }
    H5Fclose(fd);
    if (rank == 0) free(bounds);

    MPI_Allreduce(&my_nchunks_read, &aggr_nchunks_read, 1, MPI_LONG_LONG, MPI_SUM, comm);

    if (rank == 0) {
        printf("Read amount MAX=%.2f MiB MIN=%.2f MiB (among all datasets, all processes)\n",
               (float)maxRead/1048576.0,(float)minRead/1048576.0);
//        printf("Read amount MAX=%lld B MIN=%lld B (among all datasets, all processes)\n", maxRead,minRead);
        printf("Total amount of datasets %.2f MiB (%.2f GiB)\n",
               (float)all_dset_size/1048576.0, (float)all_dset_size/1073741824.0);
        printf("Total amount of evt.seq  %.2f MiB (%.2f GiB)\n",
               (float)all_evt_seq_size/1048576.0, (float)all_evt_seq_size/1073741824.0);
        if (profile >= 1) {
            printf("total number of chunks in all %d datasets (exclude /spill/evt.seq): %d\n",
                   nDatasets, total_nchunks);
            printf("Aggregate number of chunks read by all processes: %lld\n",
                   aggr_nchunks_read);
            printf("        averaged among processes: %.2f\n", (float)aggr_nchunks_read/nprocs);
            printf("        averaged among processes among datasets: %.2f\n", (float)aggr_nchunks_read/nprocs/nDatasets);
            printf("Out of %d chunks, number of chunks read by two or more processes: %d\n",
                   total_nchunks,nchunks_shared);
            printf("Out of %d chunks, most shared chunk is read by number of processes: %d\n",
                   total_nchunks,max_shared_chunks);
            printf("----------------------------------------------------\n");
        }
        if (profile == 2) {
            printf("----------------------------------------------------\n");
            for (j=0, g=0; g<nGroups; g++)
                for (d=0; d<groups[g].nDatasets; d++)
                    printf("dataset[%3d] nchunks: %4d read time: %.4f sec. (%s)\n",
                           j++, groups[g].nChunks[d], groups[g].read_t[d], groups[g].dset_names[d]);
        }
    }
    fflush(stdout);
    MPI_Barrier(comm);

    if (profile >= 1)
        printf("rank %3d: number of chunks read=%lld (max=%d min=%d avg=%.2f among %d datasets, exclude evt.seq)\n",
               rank, my_nchunks_read, max_nchunks_read, min_nchunks_read,
               (float)my_nchunks_read/(float)nDatasets, nDatasets-nGroups);

    return 1;
}

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  [-h]           print this command usage message\n\
  [-v]           verbose mode (default: off)\n\
  [-p number]    performance profiling method (0, 1, or 2)\n\
                 0: report file open, close, read timings (default)\n\
                 1: report various chunk numebrs read information\n\
                 2: report read times for individual datasets\n\
  [-s number]    read method for evt.seq (0, 1, or 2)\n\
                 0: root process reads evt.seq and broadcasts (default)\n\
                 1: all processes read the entire evt.seq collectively\n\
                 2: root process reads evt.seq and scatters boundaries\n\
  [-m number]    read method for other datasets (0 or 1)\n\
                 0: use H5Dread (default)\n\
                 1: use MPI_file_read_all\n\
  [-l file_name] name of file containing dataset names to be read\n\
  [-i file_name] name of input HDF5 file\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v] [-p number] [-s number] [-m number] [-l file_name] [-i file_name]\n%s\n",
           progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) {
    hid_t   fd, fapl_id, xfer_plist;
    hsize_t *starts=NULL, *ends;
    herr_t  err;

    int seq_opt=0, dset_opt=0, profile=0;
    int c, d, g, nprocs, rank, nDatasets, nGroups, mpi_err;
    char *listfile=NULL, *infile=NULL;
    size_t low=0, high=0;
    double open_t, read_t, close_t;
    double max_open_t, max_read_t, max_close_t;
    double min_open_t, min_read_t, min_close_t;
    MPI_File fh;
    MPI_Info info = MPI_INFO_NULL;
    NOvA_group *groups=NULL;

    verbose = 0;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvp:s:m:l:i:")) != -1)
        switch(c) {
            case 'h': if (rank  == 0) usage(argv[0]);
                      goto fn_exit;
            case 'v': verbose = 1;
                      break;
            case 'p': profile = atoi(optarg);
                      break;
            case 's': seq_opt = atoi(optarg);
                      break;
            case 'm': dset_opt = atoi(optarg);
                      break;
            case 'l': listfile = strdup(optarg);
                      break;
            case 'i': infile = strdup(optarg);
                      break;
            default: break;
        }

    if (listfile == NULL) { /* list file name is mandatory */
        if (rank  == 0) {
            printf("list file is missing\n");
            usage(argv[0]);
        }
        goto fn_exit;
    }
    if (infile == NULL) { /* input file name is mandatory */
        if (rank  == 0) {
            printf("input file is missing\n");
            usage(argv[0]);
        }
        goto fn_exit;
    }

    if (rank == 0) {
        printf("Number of MPI processes = %d\n", nprocs);
        printf("Input dataset name file '%s'\n", listfile);
        printf("Input concatenated HDF5 file '%s'\n", infile);
    }

    /* read dataset names and get number of datasets, number of groups, maximum
     * number of datasets among groups
     */
    nGroups = read_dataset_names(rank, seq_opt, dset_opt, listfile,
                                 &groups, &nDatasets);

    starts = (hsize_t*) malloc(nprocs * 2 * sizeof(hsize_t));
    ends = starts + nprocs;

    MPI_Barrier(MPI_COMM_WORLD);
    open_t = MPI_Wtime();

    /* create file access property list and add MPI communicator */
    fapl_id = H5Pcreate(H5P_FILE_ACCESS); assert(fapl_id >= 0);
    err = H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL); assert(err >= 0);

    /* collectively open input file for reading */
    fd = H5Fopen(infile, H5F_ACC_RDONLY, fapl_id);
    if (fd < 0) {
        fprintf(stderr,"%d: Error: fail to open file %s (%s)\n",
                rank,  infile, strerror(errno));
        fflush(stderr);
        assert(0);
    }
    err = H5Pclose(fapl_id); assert(err >= 0);

    /* set MPI-IO hints */
    if (dset_opt > 0) {
        mpi_err = MPI_Info_create(&info); assert(mpi_err == MPI_SUCCESS);
        mpi_err = MPI_Info_set(info, "romio_cb_read", "enable"); assert(mpi_err == MPI_SUCCESS);
        mpi_err = MPI_File_open(MPI_COMM_WORLD, infile, MPI_MODE_RDONLY, info, &fh);
        assert(mpi_err == MPI_SUCCESS);
        mpi_err = MPI_Info_free(&info); assert(mpi_err == MPI_SUCCESS);
    }
    open_t = MPI_Wtime() - open_t;

    MPI_Barrier(MPI_COMM_WORLD);
    read_t = MPI_Wtime();

    /* calculate the range of event IDs responsible by all process */
    calculate_starts_ends(fd, nprocs, rank, starts, ends);

    /* set MPI-IO collective transfer mode */
    xfer_plist = H5Pcreate(H5P_DATASET_XFER); assert(xfer_plist>=0);
    err = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE); assert(err>=0);

    /* iterate all groups */
    for (g=0; g<nGroups; g++) {
        /* First, read dataset evt.seq, the first dataset in the group, and calculate
         * the event range (low and high indices) responsible by this process
         */
        if (profile == 2) groups[g].read_t[0]=MPI_Wtime();
        if (!strcmp(groups[g].dset_names[0], "/spill/evt.seq")) {
            /* no need to read /spill/evt.seq */
            low  = starts[rank];
            high = ends[rank];
        }
        else {
	    read_evt_seq(fd, groups+g, nprocs, rank, profile, seq_opt,
                         xfer_plist, starts, ends, &low, &high);
        }
        if (profile == 2) groups[g].read_t[0]=MPI_Wtime()-groups[g].read_t[0];

        if (dset_opt == 0)
            /* read the remaining datasets using H5Dread() */
            err = read_hdf5(fd, rank, groups+g, xfer_plist, low, high, profile);
        else if (dset_opt == 1)
            /* read the remaining datasets using MPI-IO, one dataset at a time */
            err = read_mpio(fd, rank, groups+g, fh, low, high, profile);

        /* This is where PandAna performs computation to identify events of
         * interest from the read buffers
         */

        /* free read allocated buffers all at once */
        for (d=1; d<groups[g].nDatasets; d++)
            free(groups[g].buf[d]);
    }
    err = H5Pclose(xfer_plist); assert(err>=0);

    read_t = MPI_Wtime() - read_t;

    /* close input file */
    MPI_Barrier(MPI_COMM_WORLD);
    close_t = MPI_Wtime();
    err = H5Fclose(fd); assert(err >= 0);

    if (dset_opt > 0) {
        mpi_err = MPI_File_close(&fh);
        assert(mpi_err == MPI_SUCCESS);
    }
    close_t = MPI_Wtime() - close_t;

    /* find the max/min timings among all processes */
    MPI_Allreduce(&open_t,  &max_open_t, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&read_t,  &max_read_t, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&close_t, &max_close_t,1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&open_t,  &min_open_t, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&read_t,  &min_read_t, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&close_t, &min_close_t,1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    if (profile == 2) {
        for (g=0; g<nGroups; g++)
            MPI_Allreduce(MPI_IN_PLACE, groups[g].read_t, groups[g].nDatasets,
                          MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        printf("MAX open_time=%.2f read_time=%.2f close_time=%.2f\n",
               max_open_t,max_read_t,max_close_t);
        printf("MIN open_time=%.2f read_time=%.2f close_time=%.2f\n",
               min_open_t,min_read_t,min_close_t);
        printf("----------------------------------------------------\n");
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    chunk_statistics(MPI_COMM_WORLD, infile, groups, nGroups, starts, ends, seq_opt, profile);

fn_exit:
    if (starts != NULL) free(starts);

    if (groups != NULL) {
        for (g=0; g<nGroups; g++) {
            for (d=0; d<groups[g].nDatasets; d++)
                free(groups[g].dset_names[d]);
            free(groups[g].dset_names);
            free(groups[g].nChunks);
            free(groups[g].buf);
            free(groups[g].read_t);
        }
        free(groups);
    }
    if (listfile != NULL) free(listfile);
    if (infile != NULL) free(infile);

    MPI_Finalize();
    return 0;
}
