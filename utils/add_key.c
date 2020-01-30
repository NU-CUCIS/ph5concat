/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
 * See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strdup(), strlen() */
#include <unistd.h> /* getopt() */
#include <libgen.h> /* dirname() */
#include <assert.h> /* assert() */

#include <hdf5.h>

static int verbose;

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

#if defined PROFILE && PROFILE
#define SET_TIMER(ts) { ts = wtime(); }
#define GET_TIMER(ts, te, t) { te = wtime(); t += te - ts; ts = te; }

static double dopen_time, dcreate_time, dread_time, dwrite_time, dclose_time;

#include <sys/time.h>
static
double wtime(void) 
{
    double          now_time;
    struct timeval  etstart;
    struct timezone tzp;

    if (gettimeofday(&etstart, &tzp) == -1) {
        fprintf(stderr, "Error: calling gettimeofday() not successful.\n");
        return 0.0;
    }

    now_time = ((double)etstart.tv_sec) +              /* in seconds */
               ((double)etstart.tv_usec) / 1000000.0;  /* in microseconds */
    return now_time;
}
#else
#define SET_TIMER(ts)
#define GET_TIMER(ts, te, t)
#endif

#define RETURN_ERROR(func_name, obj_name) { \
    printf("Error in %s line %d: calling %s for object %s\n",__FILE__,__LINE__,func_name,obj_name); \
    return -1; \
}

#define HANDLE_ERROR(msg) { \
    printf("Error at line %d: %s\n",__LINE__, msg); \
    err_exit = -1; \
    goto fn_exit; \
}

/* parameters to be passed to the call back function */
struct op_data {
    hid_t    file_id;       /* file descriptor */
    char    *part_key_base; /* name of partition base dataset */
    size_t   num_groups;
    char   **key_base;      /* full path name of partition base dataset */
    herr_t   err;
};

/*----< collect_key_name() >-------------------------------------------------*/
/* call back function used in H5Ovisit() */
static
herr_t collect_key_name(hid_t             loc_id,   /* object ID */
                        const char       *name,     /* object name */
                        const H5O_info_t *info,     /* object metadata */
                        void             *operator) /* data passed from caller */
{
    struct op_data *it_op = (struct op_data*)operator;
    it_op->err = 0;

    if (info->type == H5O_TYPE_GROUP) {
        if (verbose) printf("group name=%s\n",name);
        return 0;
    }
    else if (info->type == H5O_TYPE_DATASET) {
        char *name_copy, *dataset_name, *key_name;
        name_copy = strdup(name);
        strtok(name_copy, "/");
        dataset_name = strtok(NULL, "/");
        if (strcmp(it_op->part_key_base, dataset_name)) {
            key_name = (char*) malloc(strlen(it_op->part_key_base)+10);
            sprintf(key_name, "%s.key.seq", it_op->part_key_base);
            if (strcmp(dataset_name, key_name) == 0) {
                printf("Error: partition key dataset %s already exists\n",
                       name);
                it_op->err = -1;
                free(key_name);
                free(name_copy);
                return 1;
		/* return a positive value causes the visit iterator to
		 * immediately return that positive value, indicating
		 * short-circuit success.
		 */
	    }
            sprintf(key_name, "%s.key.cnt", it_op->part_key_base);
            if (strcmp(dataset_name, key_name) == 0) {
                printf("Error: partition key dataset %s already exists\n",
                       name);
                it_op->err = -1;
                free(key_name);
                free(name_copy);
                return 1;
            }
            free(key_name);
            free(name_copy);
            return 0;
        }
        if (verbose) printf("\tdataset name=%s\n",name);
        it_op->key_base[it_op->num_groups++] = strdup(name);
        free(name_copy);
    }

    return 0;
}

/*----< incr_raw_data_chunk_cache() >----------------------------------------*/
static
int incr_raw_data_chunk_cache(hid_t fapl_id)
{
    int err_exit=0;
    herr_t err;

    /* set the raw data chunk cache to improve (de)compression time */
    int mdc_nelmts;      /* Dummy parameter in API, no longer used by HDF5 */
    size_t rdcc_nslots;  /* Number of slots in the hash table (default: 521)*/
    size_t rdcc_nbytes;  /* Size of chunk cache in bytes (default: 1 MiB)*/
    double w0;           /* policy to flush chunks (default: 0.75)
                          * 0: fully read/written chunks are treated no
                          *    differently than other chunks
                          * 1: fully read/written chunks are always preempted
                          *    before other chunks
                          */

    err = H5Pget_cache(fapl_id, &mdc_nelmts, &rdcc_nslots, &rdcc_nbytes, &w0);
    if (err < 0) HANDLE_ERROR("H5Pget_cache")

    /* increase cache size to 64 MiB */
    rdcc_nbytes = 67108864;

    /* set to 1 suggested by HDF5 for when only reads or writes chunks once */
    w0 = 1.;

    /* rdcc_nslots should be a prime number and approximately 100 times number
     * of chunks that can fit in rdcc_nbytes. The default value used by HDF5 is
     * 521. Bigger prime numbers are 10007, 67231, etc.
     */
    rdcc_nslots = 67231;

    err = H5Pset_cache(fapl_id, mdc_nelmts, rdcc_nslots, rdcc_nbytes, w0);
    if (err < 0) HANDLE_ERROR("H5Pset_cache")
fn_exit:
    return err_exit;
}

/*----< add_key() >----------------------------------------------------------*/
static
herr_t add_key(hid_t    file_id,    /* file ID */
               int      dry_run,    /* dry run for testing */
               size_t   num_groups, /* number of groups */
               char   **key_base,   /* full path name of key base dataset */
               size_t  *num_nonzero_groups, /* number of groups of size > 0 */
               size_t  *max_seq,    /* max key.seq size */
               size_t  *min_seq,    /* min key.seq size */
               size_t  *avg_seq,    /* avg key.seq size */
               size_t  *max_cnt,    /* max key.cnt size */
               size_t  *min_cnt,    /* min key.cnt size */
               size_t  *avg_cnt,    /* avg key.cnt size */
               double  *timing)
{
    int ii, err_exit=0;
    herr_t err=0;
#if defined PROFILE && PROFILE
    double ts=0.0, te=0.0;
#endif

    *avg_seq = 0;
    *avg_cnt = 0;
    *num_nonzero_groups = 0;
    for (ii=0; ii<num_groups; ii++) {
        size_t jj, nn, cnt_len;
        char *grp_name, *dset_name, *seq_name, *cnt_name;
        int ndims;
        hid_t run_id, subrun_id;
        hid_t base_id, type_id, space_id, dcpl_id, seq_id, cnt_id;
        hsize_t type_size, dset_dims[2], maxdims[2];
        unsigned int *run_buf, *subrun_buf, *base_ibuf;
        unsigned short *base_sbuf;
        unsigned long long *seq_buf, *cnt_buf;

        /* open the partition key base dataset */
        base_id = H5Dopen(file_id, key_base[ii], H5P_DEFAULT);
        if (base_id < 0) RETURN_ERROR("H5Dopen",key_base[ii])

        /* Get dimension sizes of base dataset */
        space_id = H5Dget_space(base_id);
        if (space_id < 0) RETURN_ERROR("H5Dget_space",key_base[ii])
        ndims = H5Sget_simple_extent_dims(space_id, dset_dims, maxdims);
        if (ndims < 0) RETURN_ERROR("H5Sget_simple_extent_dims",key_base[ii])
        err = H5Sclose(space_id);
        if (err < 0) RETURN_ERROR("H5Sclose",key_base[ii])

        /* key base dataset can only be 1D array */
        if (dset_dims[1] > 1) {
            H5Dclose(base_id);
            printf("Error: partition key base dataset can only be 1D  array\n");
            goto fn_exit;
        }

        seq_name = (char*) malloc(strlen(key_base[ii]) + 10);
        sprintf(seq_name, "/%s.key.seq",key_base[ii]);
        cnt_name = (char*) malloc(strlen(key_base[ii]) + 10);
        sprintf(cnt_name, "/%s.key.cnt",key_base[ii]);

        SET_TIMER(ts)
        if (dset_dims[0] == 0) { /* zero-size datasets */
            dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
            if (dcpl_id < 0) RETURN_ERROR("H5Pcreate", key_base[ii])

            /* use compact layout */
            err = H5Pset_layout(dcpl_id, H5D_COMPACT);
            if (err < 0) RETURN_ERROR("H5Pset_layout", key_base[ii])

            space_id = H5Screate_simple(2, dset_dims, NULL);
            if (err < 0) RETURN_ERROR("H5Screate_simple",key_base[ii])

            if (!dry_run) {
                /* create datasets key.seq */
                seq_id = H5Dcreate2(file_id, seq_name, H5T_STD_I64LE, space_id,
                                    H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
                if (seq_id < 0) RETURN_ERROR("H5Dcreate2", seq_name)
                err = H5Dclose(seq_id);
                if (err < 0) RETURN_ERROR("H5Dclose", seq_name)

                /* create datasets key.cnt */
                cnt_id = H5Dcreate2(file_id, cnt_name, H5T_STD_I64LE, space_id,
                                    H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
                if (seq_id < 0) RETURN_ERROR("H5Dcreate2", cnt_name)
                err = H5Dclose(cnt_id);
                if (err < 0) RETURN_ERROR("H5Dclose", cnt_name)
            }
            err = H5Pclose(dcpl_id);
            if (err < 0) RETURN_ERROR("H5Pclose", key_base[ii])
            err = H5Sclose(space_id);
            if (err < 0) RETURN_ERROR("H5Sclose",key_base[ii])
            err = H5Dclose(base_id);
            if (err < 0) RETURN_ERROR("H5Dclose", key_base[ii])

            free(seq_name);
            free(cnt_name);

            GET_TIMER(ts, te, timing[4])
            continue;
        }
        /* now the base dataset is non-zero sized */
        (*num_nonzero_groups)++;

        /* read dataset 'run' */
        char *name_copy = strdup(key_base[ii]);
        grp_name = dirname(name_copy);
        dset_name = (char*) malloc(strlen(grp_name)+10);

        sprintf(dset_name, "/%s/run", grp_name);
        run_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
        if (run_id < 0) RETURN_ERROR("H5Dopen", dset_name)

        run_buf = (unsigned int*) malloc(dset_dims[0] * sizeof(unsigned int));
        err = H5Dread(run_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      run_buf);
        if (err < 0) RETURN_ERROR("H5Dread", dset_name)
        err = H5Dclose(run_id);
        if (err < 0) RETURN_ERROR("H5Dclose", dset_name);

        /* read dataset 'subrun' */
        sprintf(dset_name, "/%s/subrun", grp_name);
        subrun_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
        if (subrun_id < 0) RETURN_ERROR("H5Dopen", dset_name)

        subrun_buf = (unsigned int*) malloc(dset_dims[0] * sizeof(unsigned int));
        err = H5Dread(subrun_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      subrun_buf);
        if (err < 0) RETURN_ERROR("H5Dread", dset_name)
        err = H5Dclose(subrun_id);
        if (err < 0) RETURN_ERROR("H5Dclose", dset_name);

        free(name_copy);
        free(dset_name);

        /* read base dataset */
        type_id = H5Dget_type(base_id);
        if (type_id < 0) RETURN_ERROR("H5Dget_type",key_base[ii])

        type_size = H5Tget_size(type_id);
        if (type_size < 0) RETURN_ERROR("H5Tget_size",key_base[ii])

        err = H5Tclose(type_id);
        if (err < 0) RETURN_ERROR("H5Tclose", key_base[ii])

        /* allocate read buffers and read the entire base dataset */
        if (type_size == 4) { /* base dataset is unsigned int */
            base_ibuf = (unsigned int*) malloc(dset_dims[0] *
                                               sizeof(unsigned int));
            err = H5Dread(base_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, base_ibuf);
        }
        else if (type_size == 2) { /* base dataset is unsigned short */
            base_sbuf = (unsigned short*) malloc(dset_dims[0] *
                                                 sizeof(unsigned short));
            err = H5Dread(base_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, base_sbuf);
        }
        else {
            printf("Error in %s line %d: unsupported datatype size %llu for key base dataset\n",
                   __func__,__LINE__,type_size);
            err = H5Sclose(space_id);
            if (err < 0) RETURN_ERROR("H5Sclose",key_base[ii])
            err = H5Dclose(base_id);
            if (err < 0) RETURN_ERROR("H5Dclose", key_base[ii])
            goto fn_exit;
        }
        if (err < 0) RETURN_ERROR("H5Dread", key_base[ii])

        /* close base dataset, as it has been entirely read */
        err = H5Dclose(base_id);
        if (err < 0) RETURN_ERROR("H5Dclose", key_base[ii]);
        GET_TIMER(ts, te, timing[5])

        /* allocate write buffers for key.seq and key.cnt */
        seq_buf = (unsigned long long*) malloc(dset_dims[0] *
                                               sizeof(unsigned long long));
        cnt_buf = (unsigned long long*) malloc(dset_dims[0] *
                                               sizeof(unsigned long long));

        /* use contents of 'base' to populate contents of dataset 'key.seq'
         * and 'key.cnt'
         */
        nn = 0;
        seq_buf[0] = 0;
        cnt_buf[0] = 1;
        if (type_size == 4) {
            for (jj=1; jj<dset_dims[0]; jj++) {
                if (   run_buf[jj] ==    run_buf[jj-1] &&
                    subrun_buf[jj] == subrun_buf[jj-1] &&
                     base_ibuf[jj] ==  base_ibuf[jj-1]) {
                    cnt_buf[nn]++; /* repeated ID, increment count */
                    seq_buf[jj] = seq_buf[jj-1];
                }
                else {
                    cnt_buf[++nn] = 1; /* a new unique ID */
                    seq_buf[jj] = seq_buf[jj-1] + 1;
                }
            }
        }
        else if (type_size == 2) {
            for (jj=1; jj<dset_dims[0]; jj++) {
                if (   run_buf[jj] ==    run_buf[jj-1] &&
                    subrun_buf[jj] == subrun_buf[jj-1] &&
                     base_sbuf[jj] ==  base_sbuf[jj-1]) {
                    cnt_buf[nn]++; /* repeated ID, increment count */
                    seq_buf[jj] = seq_buf[jj-1];
                }
                else {
                    cnt_buf[++nn] = 1; /* a new unique ID */
                    seq_buf[jj] = seq_buf[jj-1] + 1;
                }
            }
        }
        cnt_len = nn + 1; /* number of elements in cnt_buf */

        free(run_buf);
        free(subrun_buf);
        if (type_size == 4)
            free(base_ibuf);
        else if (type_size == 2)
            free(base_sbuf);
        GET_TIMER(ts, te, timing[6])

        /* Now the key datasets have been populated, the next step is to
         * create them and write them to file.
         */

        /* set chunk layout for two new key datasets */
        dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (dcpl_id < 0) RETURN_ERROR("H5Pcreate", key_base[ii])

        /* set compression level for two new key datasets */
        unsigned int zip = 6;
        err = H5Pset_deflate(dcpl_id, zip);
        if (err < 0) RETURN_ERROR("H5Pset_deflate", key_base[ii])

        /* no need to fill, as the entire datasets are written */
        err = H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);
        if (err < 0) RETURN_ERROR("H5Pset_fill_time", key_base[ii])

        /* make chunk size 1048576 bytes */
        hsize_t chunk_dims[2];;
        chunk_dims[1] = 1;

        /* create key.seq */
        chunk_dims[0] = (dset_dims[0] < 131072) ? dset_dims[0] : 131072;

        err = H5Pset_chunk(dcpl_id, 2, chunk_dims);
        if (err < 0) RETURN_ERROR("H5Pset_chunk", key_base[ii])

        space_id = H5Screate_simple(2, dset_dims, maxdims);
        if (space_id < 0) RETURN_ERROR("H5Screate_simple", seq_name)

        if (!dry_run) {
            seq_id = H5Dcreate2(file_id, seq_name, H5T_STD_I64LE, space_id,
                                H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
            if (seq_id < 0) RETURN_ERROR("H5Dcreate2", seq_name)
            GET_TIMER(ts, te, timing[7])

            /* write key.seq (entire dataset) */
            err = H5Dwrite(seq_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, seq_buf);
            if (err < 0) RETURN_ERROR("H5Dwrite", seq_name)
            GET_TIMER(ts, te, timing[8])

            err = H5Dclose(seq_id);
            if (err < 0) RETURN_ERROR("H5Dclose", seq_name)
            GET_TIMER(ts, te, timing[9])
        }
        err = H5Sclose(space_id);
        if (err < 0) RETURN_ERROR("H5Sclose",seq_name)

        free(seq_buf);

        if (*num_nonzero_groups == 1) {
            *min_seq = dset_dims[0];
            *min_cnt = cnt_len;
        }
        else {
            *min_seq = MIN(*min_seq, dset_dims[0]);
            *min_cnt = MIN(*min_cnt, cnt_len);
        }
        *max_seq = MAX(*max_seq, dset_dims[0]);
        *max_cnt = MAX(*max_cnt, cnt_len);
        *avg_seq += dset_dims[0];
        *avg_cnt += cnt_len;

        /* create key.cnt */
        dset_dims[0] = cnt_len; /* size of key.cnt may be smaller */
        chunk_dims[0] = (dset_dims[0] < 131072) ? dset_dims[0] : 131072;

        err = H5Pset_chunk(dcpl_id, 2, chunk_dims);
        if (err < 0) RETURN_ERROR("H5Pset_chunk", key_base[ii])

        space_id = H5Screate_simple(2, dset_dims, maxdims);
        if (space_id < 0) RETURN_ERROR("H5Screate_simple", cnt_name)

        if (!dry_run) {
            cnt_id = H5Dcreate2(file_id, cnt_name, H5T_STD_I64LE, space_id,
                                H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
            if (cnt_id < 0) RETURN_ERROR("H5Dcreate2", cnt_name)
            GET_TIMER(ts, te, timing[10])

            /* write key.cnt (entire dataset) */
            err = H5Dwrite(cnt_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, cnt_buf);
            if (err < 0) RETURN_ERROR("H5Dwrite", cnt_name)
            GET_TIMER(ts, te, timing[11])

            err = H5Dclose(cnt_id);
            if (err < 0) RETURN_ERROR("H5Dclose", cnt_name)
            GET_TIMER(ts, te, timing[12])
        }
        err = H5Sclose(space_id);
        if (err < 0) RETURN_ERROR("H5Sclose",cnt_name)

        free(cnt_buf);

        /* free up allocated object IDs */
        err = H5Pclose(dcpl_id);
        if (err < 0) RETURN_ERROR("H5Pclose", key_base[ii])

        free(seq_name);
        free(cnt_name);
        free(key_base[ii]);
    }

    if (*num_nonzero_groups > 0) {
        *avg_seq /= *num_nonzero_groups;
        *avg_cnt /= *num_nonzero_groups;
    }

fn_exit:
    return err_exit;
}

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
         char *type_name="";
         ot = H5Iget_type(objs[ii]);
              if (ot == H5I_FILE)      continue; /*  type_name = "H5I_FILE"; */
         else if (ot == H5I_GROUP)     type_name = "H5I_GROUP";
         else if (ot == H5I_DATATYPE)  type_name = "H5I_DATATYPE";
         else if (ot == H5I_DATASPACE) type_name = "H5I_DATASPACE";
         else if (ot == H5I_DATASET)   type_name = "H5I_DATASET";
         else if (ot == H5I_ATTR)      type_name = "H5I_ATTR";
         H5Iget_name(objs[ii], obj_name, 1024);
         printf("%s %4zd: type %s, name %s\n",
                filename, ii, type_name, obj_name);
    }
    free(objs);
    return howmany;
}

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  [-h]            print this command usage message\n\
  [-v]            verbose mode (default: off)\n\
  [-n]            dry run without creating key datasets (default: disabled)\n\
  -k datset_name  name of dataset used to generate partitioning keys (required)\n\
  file _name      HDF5 file name (required)\n\n\
  This utility program adds partitioning key datasets, key.seq and key.cnt,\n\
  to an HDF5 file.\n\n\
  Requirements for the HDF5 file:\n\
    1. contains multiple groups only at root level\n\
    2. each group may  contain multiple 2D datasets\n\
  Requirements for the partitioning base daatset:\n\
    1. the second dimension size of base dataset must be 1\n\
    2. if base dataset is missing in a group, the key datasets will not be\n\
       generated for that group\n\
    3. currently supports base dataset of type either 4-byte unsigned int or\n\
       2-byte short int\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v] -k dataset_name file_name\n%s\n", progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int c, err_exit=0, dry_run=0;
    char msg[1024], *infile=NULL, *part_key_base=NULL;
    herr_t err;
    hid_t file_id=-1, fapl_id=-1;
    H5G_info_t grp_info;
    size_t num_orig_groups, num_nonzero_groups;
    size_t max_seq=0, max_cnt=0, min_seq, min_cnt, avg_seq, avg_cnt;
    struct op_data it_op;
    double timing[13];

#if defined PROFILE && PROFILE
    double ts=0.0, te=0.0;
    dopen_time = dcreate_time = dread_time = dwrite_time = dclose_time = 0.0;
#endif
    for (c=0; c<13; c++) timing[c] = 0.0;

    verbose = 0; /* default is quiet */

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvnk:")) != -1)
        switch(c) {
            case 'h': usage(argv[0]);
                      err_exit = -1;
                      goto fn_exit;
            case 'v': verbose = 1;
                      break;
            case 'n': dry_run = 1;
                      break;
            case 'k': part_key_base = strdup(optarg);
                      break;
            default: break;
        }

    if (part_key_base == NULL) { /* input file name is mandatory */
        printf("Error: partition key base dataset name is missing.\n\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    it_op.part_key_base = part_key_base;

    if (argv[optind] == NULL) { /* input file name is mandatory */
        printf("input file is missing\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    infile = strdup(argv[optind]);

    if (verbose) printf("input file: %s\n", infile);

    SET_TIMER(ts)

    /* create file access property for read and write */
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id < 0) HANDLE_ERROR("H5Pcreate")

    err = incr_raw_data_chunk_cache(fapl_id);
    if (err < 0) HANDLE_ERROR("incr_raw_data_chunk_cache")

    /* open input file */
    if (dry_run)
        /* read-only mode  */
        file_id = H5Fopen(infile, H5F_ACC_RDONLY, fapl_id);
    else
        /* read-and-write mode  */
        file_id = H5Fopen(infile, H5F_ACC_RDWR, fapl_id);
    if (file_id < 0) {
        sprintf(msg, "Can't open input file %s\n", infile);
        HANDLE_ERROR(msg)
    }
    it_op.file_id = file_id;

    /* retrieve metadata about root group */
    err = H5Gget_info_by_name(file_id, "/", &grp_info, H5P_DEFAULT);
    if (err < 0) RETURN_ERROR("H5Gget_info_by_name", "root group")
    num_orig_groups = grp_info.nlinks;
    it_op.key_base = (char**) malloc(num_orig_groups * sizeof(char*));
    it_op.num_groups = 0;

    GET_TIMER(ts, te, timing[0])

    /* Iterate all objects and collect key base full path names */
    err = H5Ovisit(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, collect_key_name,
                   &it_op);
    if (err < 0) HANDLE_ERROR("H5Ovisit")
    if (it_op.err < 0) HANDLE_ERROR("H5Ovisit")

    assert(num_orig_groups >= it_op.num_groups);

    GET_TIMER(ts, te, timing[1])

    err = add_key(file_id, dry_run, it_op.num_groups, it_op.key_base,
                  &num_nonzero_groups, &max_seq, &min_seq, &avg_seq,
                  &max_cnt, &min_cnt, &avg_cnt, timing);
    if (err < 0) HANDLE_ERROR("add_key")

    GET_TIMER(ts, te, timing[2])
    
fn_exit:
    free(it_op.key_base);

    if (file_id >= 0) {
        check_h5_objects(infile, file_id);
        err = H5Fclose(file_id);
        if (err < 0) printf("Error at line %d: H5Fclose\n",__LINE__);
    }
    if (fapl_id >= 0) {
        err = H5Pclose(fapl_id);
        if (err < 0) printf("Error at line %d: H5Pclose\n",__LINE__);
    }

    GET_TIMER(ts, te, timing[3])

#if defined PROFILE && PROFILE
    if (err_exit == 0) {
        printf("-------------------------------------------------------\n");
        printf("Dry-run mode                      = %s\n", dry_run?"YES":"NO");
        printf("Input file name                   = %s\n", infile);
        printf("Partition key base dataset name   = %s\n", part_key_base);
        printf("number of groups in the file      = %zd\n", num_orig_groups);
        printf("number of groups contain key base = %zd\n", it_op.num_groups);
        printf("number of non-zero key bases      = %zd\n", num_nonzero_groups);
        printf("Partition key %s.key.seq :\n",part_key_base);
        printf("  max length among all groups     = %zd\n", max_seq);
        printf("  min length among all groups     = %zd\n", min_seq);
        printf("  avg length among all groups     = %zd\n", avg_seq);
        printf("Partition key %s.key.cnt :\n",part_key_base);
        printf("  max length among all groups     = %zd\n", max_cnt);
        printf("  min length among all groups     = %zd\n", min_cnt);
        printf("  avg length among all groups     = %zd\n", avg_cnt);
        printf("-------------------------------------------------------\n");
        printf("Open input file                   = %7.2f sec\n", timing[0]);
        printf("Collect partition key base names  = %7.2f sec\n", timing[1]);
        printf("Add partition key datasets        = %7.2f sec\n", timing[2]);
        printf("  Create zero-size key datasets   = %7.2f sec\n", timing[4]);
        printf("  Read run, subrun, base datasets = %7.2f sec\n", timing[5]);
        printf("  Generate key.seq and key.cnt    = %7.2f sec\n", timing[6]);
        printf("  Create key.seq                  = %7.2f sec\n", timing[7]);
        printf("  Write  key.seq                  = %7.2f sec\n", timing[8]);
        printf("  Close  key.seq                  = %7.2f sec\n", timing[9]);
        printf("  Create key.cnt                  = %7.2f sec\n", timing[10]);
        printf("  Write  key.cnt                  = %7.2f sec\n", timing[11]);
        printf("  Close  key.cnt                  = %7.2f sec\n", timing[12]);
        printf("File close                        = %7.2f sec\n", timing[3]);
        printf("-------------------------------------------------------\n");
        printf("End-to-end                        = %7.2f sec\n",
               timing[0] + timing[1] + timing[2] + timing[3]);
    }
#endif
    if (infile != NULL) free(infile);
    if (part_key_base != NULL) free(part_key_base);

    return (err_exit != 0);
}
