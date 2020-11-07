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

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

#define LINE_SIZE 64

static
int read_dataset_names(int     rank,
                       char   *filename,
                       int    *nDatasets,  /* total number of datasets */
                       int    *nGroups,    /* total number of groups */
                       int    *maxDsetGrp, /* max number of datasets among groups */
                       char ***name_list)
{
    FILE *fptr;
    int i, nDsetGrp;
    char line[LINE_SIZE], gname[LINE_SIZE], *cur_gname;

    fptr = fopen(filename, "r");
    if (fptr == NULL) {
        fprintf(stderr,"%d: Error: fail to open file %s (%s)\n",
               rank,filename,strerror(errno));
        assert(0);
    }

    /* check number of datasets listed in the file */
    *nDatasets = 0;
    *nGroups = 0;
    *maxDsetGrp = 0;
    nDsetGrp = 0;
    gname[0] = '\0';
    while (fgets(line, LINE_SIZE, fptr) != NULL) {
        if (strlen(line) == 0) continue;
        if (line[0] == '\n') continue;
        if (line[0] == '#') continue;
        (*nDatasets)++;

        /* retrieve group name */
        cur_gname = strtok(line, "/");
        if (strcmp(gname, cur_gname)) { /* entering a new group */
            strcpy(gname, cur_gname);
            if (nDsetGrp > *maxDsetGrp) *maxDsetGrp = nDsetGrp;
            nDsetGrp = 0;
            (*nGroups)++;
        }
        nDsetGrp++;
    }

    *name_list = (char**) malloc(*nDatasets * sizeof(char*));
    (*name_list)[0] = (char*)  malloc(*nDatasets * LINE_SIZE * sizeof(char));
    for (i=1; i<*nDatasets; i++)
        (*name_list)[i] = (*name_list)[i-1] + LINE_SIZE;

   /* read dataset names */
    rewind(fptr);
    i = 0;
    while (fgets(line, LINE_SIZE, fptr) != NULL) {
        if (strlen(line) == 0) continue;
        if (line[0] == '\n') continue;
        if (line[0] == '#') continue;
        strcpy((*name_list)[i], line);
        (*name_list)[i][strlen(line)-1] = '\0';
        i++;
    }
    assert(i == *nDatasets);
    fclose(fptr);
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
  [-r number]    read method (0, 1, or 2)\n\
                 0: root process reads evt.seq and broadcasts (default)\n\
                 1: all processes read the entire evt.seq collectively\n\
                 2: root process reads evt.seq and scatters boundaries\n\
  [-l file_name] name of file containing dataset names to be read\n\
  [-i file_name] name of input HDF5 file\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v] [-p number] [-r number] [-l file_name] [-i file_name]\n%s\n",
           progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) {
    hid_t   fd, fapl_id, seq, dset, fspace, mspace, dtype, xfer_plist;
    hid_t   chunk_plist;
    hsize_t ndims, dims[5], seq_len, chunk_dims[2];
    hsize_t my_start, my_end, my_count, start[2], count[2];
    hsize_t *starts, *ends;
    herr_t  err;

    int r_opt=0, verbose=0, profile=0, nchunks_shared=0, max_shared_chunks=0;
    int first_chunk_id, last_chunk_id, nchunks, nchunks_read=0;
    int max_nchunks_read=0, min_nchunks_read=INT_MAX, total_nchunks=0;
    int chunk_ndims, nchunks_dim0, nchunks_dim1, *nchunks_dset;
    int c, d, j, nprocs, rank, nDatasets, nDsetGrp, nGroups, maxDsetGrp;
    char *listfile=NULL, *infile=NULL;
    char **name_list, gname[LINE_SIZE];
    size_t low, high, buf_len;
    long long *bounds, low_high[2], total_nchunks_read;
    double open_t, read_t, close_t;
    double max_open_t, max_read_t, max_close_t, *rt;
    double min_open_t, min_read_t, min_close_t;
    long long maxRead, minRead, maxRead_all, minRead_all, totalRead;
    MPI_Info info = MPI_INFO_NULL;

    /* read buffer */
    void **buf;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvp:r:l:i:")) != -1)
        switch(c) {
            case 'h': if (rank  == 0) usage(argv[0]);
                      goto fn_exit;
            case 'v': verbose = 1;
                      break;
            case 'p': profile = atoi(optarg);
                      break;
            case 'r': r_opt = atoi(optarg);
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

    /* read dataset names and get number of datasets, number of groups, maximum
     * number of datasets among groups
     */
    err = read_dataset_names(rank, listfile, &nDatasets, &nGroups, &maxDsetGrp,
                             &name_list);
    assert(err == 1);

    if (rank == 0) {
        printf("Number of MPI processes = %d\n", nprocs);
        printf("Input dataset name file '%s'\n", listfile);
        printf("Input concatenated HDF5 file '%s'\n", infile);
        printf("Number of datasets to read = %d\n", nDatasets);
        printf("Number of groups = %d\n",nGroups);
        printf("Maximum number of datasets among groups = %d\n",maxDsetGrp);
        if (r_opt == 0)
            printf("Read method: root process reads evt.seq and broadcasts\n");
        else if (r_opt == 1)
            printf("Read method: all processes read the entire evt.seq collectively\n");
        else
            printf("Read method: root process reads evt.seq and scatters boundaries\n");
        printf("----------------------------------------------------\n");
    }
    fflush(stdout);

    if (profile == 2) { /* collect read time of individual datasets */
        rt = (double*)malloc(nDatasets*sizeof(double));
        if (rank == 0) nchunks_dset = (int*)malloc(nDatasets*sizeof(int));
    }

    /* allocate array of read buffer pointers */
    buf = (void**) calloc(maxDsetGrp, sizeof(void*));

    if (r_opt == 2) {
        starts = (hsize_t*) malloc(nprocs * 2 * sizeof(hsize_t));
        ends = starts + nprocs;
        bounds = (long long*) malloc(nprocs * 2 * sizeof(long long));
    }

    MPI_Barrier(MPI_COMM_WORLD);
    open_t = MPI_Wtime();

    /* create file access property list */
    fapl_id = H5Pcreate(H5P_FILE_ACCESS); assert(fapl_id >= 0);
    /* add MPI communicator to the access property */
    err = H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info); assert(err >= 0);

    /* open input file for reading */
    fd = H5Fopen(infile, H5F_ACC_RDONLY, fapl_id);
    if (fd < 0) {
        fprintf(stderr,"%d: Error: fail to open file %s (%s)\n",
                rank,  infile, strerror(errno));
        fflush(stderr);
        assert(0);
    }
    err = H5Pclose(fapl_id); assert(err >= 0);
    open_t = MPI_Wtime() - open_t;

    MPI_Barrier(MPI_COMM_WORLD);
    read_t = MPI_Wtime();

    /* open dataset '/spill/evt.seq' */
    seq = H5Dopen2(fd, "/spill/evt.seq", H5P_DEFAULT); assert(seq >= 0);

    /* inquire the size of '/spill/evt.seq' */
    fspace = H5Dget_space(seq); assert(fspace >= 0);
    ndims = H5Sget_simple_extent_ndims(fspace); assert(ndims >= 0);
    err = H5Sget_simple_extent_dims(fspace, dims, NULL); assert(err >= 0);

    /* 2nd dimension of /spill/evt.seq is always 1 */
    assert(dims[1] == 1);

    seq_len = dims[0];
    if (verbose && rank == 0)
        printf("Size of /spill/evt is %zd\n", (size_t)seq_len);

    /* calculate the range of event IDs assigned to this process. We name it
     * 'partition domain'. The domain ranges from my_start to my_end.
     */
    my_count = seq_len / nprocs;
    my_start = my_count * rank;
    if (r_opt == 2) {
        for (j=0; j<seq_len % nprocs; j++) {
            starts[j] = my_count * j + j;
            ends[j] = starts[j] + my_count;
        }
        for (; j<nprocs; j++) {
            starts[j] = my_count * j + seq_len % nprocs;
            ends[j] = starts[j] + my_count - 1;
        }
    }
    if (rank < seq_len % nprocs) {
        my_start += rank;
        my_count++;
    }
    else {
        my_start += seq_len % nprocs;
    }
    my_end = my_start + my_count - 1;

    /* Note contents of /spill/evt.seq start with 0 and increment by 1, so
     * there is no need to read the contents of it. Check the contents only
     * in verbose mode.
     */
    if (verbose) {
        int64_t v, *seq_buf;
        printf("%d: /spill/evt.seq len=%zd my_start=%zd my_end=%zd\n",
               rank, (size_t)seq_len, (size_t)my_start, (size_t)my_end);
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

    /* set MPI-IO collective transfer mode */
    xfer_plist = H5Pcreate(H5P_DATASET_XFER); assert(xfer_plist>=0);
    err = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE); assert(err>=0);

    /* iterate all datasets */
    totalRead = 0;
    nDsetGrp = 0;
    gname[0]='\0';
    for (d=0; d<nDatasets; d++) {
        hsize_t one[2]={1,1};
        char *cur_gname, name[LINE_SIZE];

        /* retrieve group name */
        strcpy(name, name_list[d]);
        cur_gname = strtok(name, "/");

        if (strcmp(gname, cur_gname)) { /* entering a new group */
            int64_t *seq_buf;
            char seq_name[LINE_SIZE];

            if (verbose && rank == 0) printf("group[%d]: %s\n", d, cur_gname);

            strcpy(gname, cur_gname);
            sprintf(seq_name, "/%s/evt.seq", gname);

            /* open dataset 'evt.seq' */
            seq = H5Dopen2(fd, seq_name, H5P_DEFAULT); assert(seq >= 0);

            /* inquire the size of 'evt.seq' in this group */
            fspace = H5Dget_space(seq); assert(fspace >= 0);
            ndims = H5Sget_simple_extent_ndims(fspace); assert(ndims >= 0);
            err = H5Sget_simple_extent_dims(fspace, dims, NULL); assert(err>=0);

            buf_len = dims[0] * sizeof(int64_t);
            seq_buf = (int64_t*) malloc(buf_len);

            if (rank == 0 && profile == 2) {
                /* calculate the number of chunks of this dataset */
                chunk_plist = H5Dget_create_plist(seq); assert(chunk_plist >= 0);
                chunk_ndims = H5Pget_chunk(chunk_plist, 2, chunk_dims);
                err = H5Pclose(chunk_plist); assert(err>=0);
                nchunks_dset[d] = dims[0] / chunk_dims[0];
                if (dims[0] % chunk_dims[0]) nchunks_dset[d]++;
            }

            if (profile == 2) rt[d]=MPI_Wtime();
            if (r_opt == 2) {
                /* rank 0 reads evt.seq, calculates lows and highs for all
                 * processes and scatters
                 */
                if (rank == 0) {
                    err = H5Dread(seq, H5T_STD_I64LE, fspace, fspace, H5P_DEFAULT, seq_buf);
                    assert(err >= 0);
                    int k=0;
                    for (j=0; j<nprocs; j++) {
                        bounds[k++] = binary_search_min(starts[j], seq_buf, dims[0]);
                        bounds[k++] = binary_search_max(ends[j],   seq_buf, dims[0]);
                    }
                }
                MPI_Scatter(bounds, 2, MPI_LONG_LONG, low_high, 2, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
                low  = low_high[0];
                high = low_high[1];
            }
            else {
                if (r_opt == 0) {
                    /* rank 0 reads evt.seq and broadcasts it */
                    if (rank == 0) {
                        err = H5Dread(seq, H5T_STD_I64LE, fspace, fspace, H5P_DEFAULT,
                                      seq_buf);
                        assert(err >= 0);
                    }
                    MPI_Bcast(seq_buf, dims[0], MPI_LONG_LONG, 0, MPI_COMM_WORLD);
                }
                else {
                    /* collective-read-the-whole-dataset is bad */
                    err = H5Dread(seq, H5T_STD_I64LE, fspace, fspace, xfer_plist, seq_buf);
                    /* all-processes-read-the-whole-dataset is even worse */
                    // err = H5Dread(seq, H5T_STD_I64LE, fspace, fspace, H5P_DEFAULT, seq_buf);
                    assert(err >= 0);
                }
	        /* find the array index range from 'low' to 'high' that falls into
                 * this process's partition domain.
                 */
	        low = binary_search_min(my_start, seq_buf, dims[0]);
                high = binary_search_max(my_end, seq_buf, dims[0]);
            }
            if (profile == 2) rt[d]=MPI_Wtime()-rt[d];

            err = H5Sclose(fspace); assert(err >= 0);
            err = H5Dclose(seq); assert(err >= 0);
            if (verbose && rank == 0)
                printf("/%s/evt.seq[0]=%"PRId64" seq[len-1=%zd]=%"PRId64"\n",
                       gname,seq_buf[0],(size_t)dims[0]-1,seq_buf[dims[0]-1]);

            if (verbose)
                printf("%d: %s my_start=%zd my_end=%zd low=%zd high=%zd\n",
                       rank, gname, (size_t)my_start,(size_t)my_end, low,high);

            free(seq_buf);

            if (d == 0) /* 1st group */
                maxRead = minRead = buf_len;
            else {
                /* free read buffers allocated in previous group */
                for (j=0; j<nDsetGrp; j++) free(buf[j]);
                nDsetGrp = 0;
                maxRead = MAX(maxRead, buf_len);
                minRead = MIN(minRead, buf_len);
            }
            totalRead += buf_len;
        }
        /* skip dataset 'evt.seq', as it has already been read */
        if (0 == strcmp("evt.seq", strtok(NULL, "/"))) continue;

        /* open dataset */
        dset = H5Dopen2(fd, name_list[d], H5P_DEFAULT); assert(dset >= 0);

        /* inquire the size of dset */
        fspace = H5Dget_space(dset); assert(fspace >= 0);
        ndims = H5Sget_simple_extent_ndims(fspace); assert(ndims >= 0);
        err = H5Sget_simple_extent_dims(fspace, dims, NULL); assert(err>=0);

        /* set subarray/hyperslab access */
        start[0] = low;
        start[1] = 0;
        count[0] = high - low + 1;
        count[1] = dims[1];
        err = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, NULL, one,
                                  count); assert(err>=0);
        mspace = H5Screate_simple(2, count, NULL); assert(mspace>=0);

        /* check data type */
        dtype = H5Dget_type(dset); assert(dtype>=0);

        /* allocate read buffers */
        buf_len = count[0] * count[1];

        if (H5Tequal(dtype, H5T_STD_U8LE)) {
            buf_len *= sizeof(unsigned char);
            totalRead += dims[0] * dims[1] * sizeof(unsigned char);
        } else if (H5Tequal(dtype, H5T_STD_I16LE)) {
            buf_len *= sizeof(short);
            totalRead += dims[0] * dims[1] * sizeof(short);
        } else if (H5Tequal(dtype, H5T_STD_U16LE)) {
            buf_len *= sizeof(unsigned short);
            totalRead += dims[0] * dims[1] * sizeof(unsigned short);
        } else if (H5Tequal(dtype, H5T_STD_I32LE)) {
            buf_len *= sizeof(int);
            totalRead += dims[0] * dims[1] * sizeof(int);
        } else if (H5Tequal(dtype, H5T_STD_U32LE)) {
            buf_len *= sizeof(unsigned int);
            totalRead += dims[0] * dims[1] * sizeof(unsigned int);
        } else if (H5Tequal(dtype, H5T_STD_I64LE)) {
            buf_len *= sizeof(long long);
            totalRead += dims[0] * dims[1] * sizeof(long long);
        } else if (H5Tequal(dtype, H5T_STD_U64LE)) {
            buf_len *= sizeof(unsigned long long);
            totalRead += dims[0] * dims[1] * sizeof(unsigned long long);
        } else if (H5Tequal(dtype, H5T_IEEE_F32LE)) {
            buf_len *= sizeof(float);
            totalRead += dims[0] * dims[1] * sizeof(float);
        } else if (H5Tequal(dtype, H5T_IEEE_F64LE)) {
            buf_len *= sizeof(double);
            totalRead += dims[0] * dims[1] * sizeof(double);
        } else {
            printf("Error: invalid data type %s\n",  name_list[d]);
            break;
        }
        buf[nDsetGrp] = (void*) malloc(buf_len);

        maxRead = MAX(maxRead, buf_len);
        minRead = MIN(minRead, buf_len);

        /* collectively read its contents */
        if (verbose)
            printf("%d: READ buf_len=%zd dataset %s\n",rank,buf_len,name_list[d]);
        if (profile == 2) rt[d]=MPI_Wtime();
        err = H5Dread(dset, dtype, mspace, fspace, xfer_plist, buf[nDsetGrp]);
        if (profile == 2) rt[d]=MPI_Wtime()-rt[d];
        assert(err >= 0);

        if (profile >= 1) {
            /* inquire chunk size along each dimension */
            chunk_plist = H5Dget_create_plist(dset); assert(chunk_plist >= 0);
            chunk_ndims = H5Pget_chunk(chunk_plist, 2, chunk_dims);
            assert(chunk_ndims == 2);
            assert(chunk_dims[1] == 1);
            err = H5Pclose(chunk_plist); assert(err>=0);

            if (rank == 0 && profile == 2) {
                nchunks_dset[d] = dims[0] / chunk_dims[0];
                if (dims[0] % chunk_dims[0]) nchunks_dset[d]++;
            }

            if (rank == 0 && r_opt == 2) {
                int prev_chunk_id=-1, max_chunks=1;
                int k=0, chunk_id, prev_chunk=-1, prev_prev_chunk=-1;
                for (j=0; j<nprocs; j++) {
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
            }
            /* calculate number of chunks read */
            first_chunk_id = low  / chunk_dims[0];
            last_chunk_id  = high / chunk_dims[0];
            nchunks = last_chunk_id - first_chunk_id + 1;
            nchunks_read += nchunks;
            max_nchunks_read = MAX(max_nchunks_read, nchunks);
            min_nchunks_read = MIN(min_nchunks_read, nchunks);
            nchunks_dim0 = dims[0] / chunk_dims[0];
            if (dims[0] % chunk_dims[0]) nchunks_dim0++;
            nchunks_dim1 = dims[1] / chunk_dims[1];
            if (dims[1] % chunk_dims[1]) nchunks_dim1++;
            total_nchunks += nchunks_dim0 * nchunks_dim1;
        }

        err = H5Tclose(dtype); assert(err >= 0);
        err = H5Sclose(fspace); assert(err >= 0);
        err = H5Sclose(mspace); assert(err >= 0);
        err = H5Dclose(dset); assert(err >= 0);

        nDsetGrp++;
    }
    /* free read buffers allocated at the last group */
    for (j=0; j<nDsetGrp; j++) free(buf[j]);
    free(buf);

    if (r_opt == 2) {
        free(starts);
        free(bounds);
    }
    err = H5Pclose(xfer_plist); assert(err>=0);

    read_t = MPI_Wtime() - read_t;

    /* close input file */
    MPI_Barrier(MPI_COMM_WORLD);
    close_t = MPI_Wtime();
    err = H5Fclose(fd); assert(err >= 0);
    close_t = MPI_Wtime() - close_t;

    /* find the max/min timings among all processes */
    MPI_Allreduce(&open_t,  &max_open_t, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&read_t,  &max_read_t, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&close_t, &max_close_t,1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&open_t,  &min_open_t, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&read_t,  &min_read_t, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&close_t, &min_close_t,1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&maxRead, &maxRead_all,1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minRead, &minRead_all,1, MPI_LONG_LONG, MPI_MIN, MPI_COMM_WORLD);
    if (profile == 2)
        MPI_Allreduce(MPI_IN_PLACE, rt, nDatasets, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    total_nchunks_read = nchunks_read;
    MPI_Allreduce(MPI_IN_PLACE, &total_nchunks_read, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) {
        float maxR, minR;
        maxR = (float)maxRead_all / 1048576;
        minR = (float)minRead_all / 1048576;
        printf("Total read amount %.2f MiB (%.2f GiB)\n",
               (float)totalRead/1048576, (float)totalRead/1073741824);
        printf("Read amount MAX=%.2f MiB MIN=%.2f MiB\n", maxR,minR);
        printf("MAX open_time=%.2f read_time=%.2f close_time=%.2f\n",
               max_open_t,max_read_t,max_close_t);
        printf("MIN open_time=%.2f read_time=%.2f close_time=%.2f\n",
               min_open_t,min_read_t,min_close_t);
        if (profile >= 1) {
            printf("----------------------------------------------------\n");
            printf("total number of chunks in all %d datasets: %d\n",
                   nDatasets, total_nchunks);
            printf("Aggregate number of chunks read by all processes: %lld\n",
                   total_nchunks_read);
            if (r_opt == 2) {
                printf("Out of %d chunks, number of chunks read by two or more processes: %d\n",
                       total_nchunks,nchunks_shared);
                printf("Out of %d chunks, most shared chunk is read by number of processes: %d\n",
                       total_nchunks,max_shared_chunks);
            }
        }
        printf("----------------------------------------------------\n");
    }
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (profile >= 1)
        printf("rank %3d: number of chunks read=%d (max=%d min=%d avg=%.2f among %d datasets)\n",
               rank, nchunks_read, max_nchunks_read, min_nchunks_read,
               (float)nchunks_read/(float)nDatasets, nDatasets);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0 && profile == 2) {
        printf("----------------------------------------------------\n");
        for (d=0; d<nDatasets; d++)
            printf("dataset[%3d] nchunks: %4d read time: %.4f sec. (%s)\n",
                   d, nchunks_dset[d], rt[d], name_list[d]);
        free(nchunks_dset);
    }

    if (profile == 2) free(rt);
    free(name_list[0]);
    free(name_list);
    if (listfile != NULL) free(listfile);
    if (infile != NULL) free(infile);

fn_exit:
    MPI_Finalize();
    return 0;
}
