/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator
 * Laboratory
 * See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcmp(), strdup() */
#include <unistd.h> /* getopt() */
#include <assert.h> /* assert() */

#include <string>
#include <unordered_map>
#include <tuple>
using namespace std;

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

#define CALLBACK_ERROR(func_name, dset_name) { \
    printf("Error in %s line %d: calling %s for group %s dataset %s\n",__FILE__,__LINE__,func_name, grp_name, dset_name); \
    return -1; \
}

/* lookup table based on 3-tuple datasets */
typedef tuple<unsigned int, unsigned int, unsigned int> key;
struct hash3 : public unary_function<key, size_t> {
    size_t operator()(const key& k) const {
        return get<0>(k) ^ get<1>(k) ^ get<2>(k);
    }
};
typedef unordered_map<key, int64_t, hash3> table;

static
table build_lookup_table(hsize_t       len,
                         unsigned int *run,
                         unsigned int *subrun,
                         unsigned int *base)
{
    hsize_t i;
    table ret;

    for (i=0; i<len; i++)
        ret[key(run[i], subrun[i], base[i])] = i;

    return ret;
}

/* parameters to be passed to the call back function */
struct op_data {
    size_t   num_groups;
    char   **grp_names;     /* [num_groups] group names */
    herr_t   err;
};

/*----< gather_grp_names() >-------------------------------------------------*/
/* call back function used in H5Ovisit() */
static
herr_t gather_grp_names(hid_t             loc_id,/* object ID */
                        const char       *name,  /* object name */
                        const H5O_info_t *info,  /* object metadata */
                        void             *op)    /* data passed from caller */
{
    struct op_data *it_op = (struct op_data*)op;
    it_op->err = 0;

    /* skip if this is a dataset */
    if (info->type == H5O_TYPE_DATASET) return 0;

    /* info->type == H5O_TYPE_GROUP */
    if (verbose) printf("group name=%s\n",name);

    /* skip root group '.' the first object */
    if (!strcmp(name, ".")) return 0;

    it_op->grp_names[it_op->num_groups++] = strdup(name);
    return 0;
}

/*----< add_seq() >---------------------------------------------------------*/
static
herr_t add_seq(bool        dry_run,
               hid_t       file_id,
               const char *grp_name,
               const char *base_name,
               table       lookup_table,
               size_t     *num_nonzero_groups,
               size_t     *max_len,    /* max seq size */
               size_t     *min_len,    /* min seq size */
               size_t     *avg_len,    /* avg seq size */
               double     *timing)
{
    int ndims;
    size_t i;
    herr_t err;
    hid_t grp_id, run_id, srun_id, base_id, dcpl_id, space_id, seq_id;
    hsize_t dset_dims[2], maxdims[2];
    unsigned int *run_buf, *srun_buf, *base_buf;
    int64_t *seq_buf;
#if defined PROFILE && PROFILE
    double ts=0.0, te=0.0;
#endif

    SET_TIMER(ts)

    /* open group */
    grp_id = H5Gopen(file_id, grp_name, H5P_DEFAULT);
    if (grp_id < 0) CALLBACK_ERROR("H5Gopen", grp_name)

    /* open base dataset in this group */
    base_id = H5Dopen(grp_id, base_name, H5P_DEFAULT);
    if (base_id < 0) CALLBACK_ERROR("H5Dopen", base_name)

    /* Get dimension sizes of base dataset */
    space_id = H5Dget_space(base_id);
    if (space_id < 0) CALLBACK_ERROR("H5Dget_space", base_name)
    ndims = H5Sget_simple_extent_dims(space_id, dset_dims, maxdims);
    if (ndims != 2) CALLBACK_ERROR("H5Sget_simple_extent_dims", base_name)

    /* Get create property of base dataset when it was created */
    dcpl_id = H5Dget_create_plist(base_id);
    if (dcpl_id < 0) CALLBACK_ERROR("H5Dget_create_plist", base_name)

    /* if base dataset is zero-sized, just create dataset key_name */
    if (dset_dims[0] == 0) {
        err = H5Dclose(base_id);
        if (err < 0) CALLBACK_ERROR("H5Dclose", base_name)
        goto seq_create;
    }

    /* non-zero size group */
    (*num_nonzero_groups)++;

    if (*num_nonzero_groups == 1) {
        *max_len = dset_dims[0];
        *min_len = dset_dims[0];
        *avg_len = dset_dims[0];
    }
    else {
        *max_len = MAX(*max_len, dset_dims[0]);
        *min_len = MIN(*min_len, dset_dims[0]);
        *avg_len += dset_dims[0];
    }

    /* allocate buffers for run, subrun, base */
    run_buf  = (unsigned int*) malloc(dset_dims[0] * 3 * sizeof(unsigned int));
    srun_buf = run_buf  + dset_dims[0];
    base_buf = srun_buf + dset_dims[0];

    /* read the entire base dataset */
    err = H5Dread(base_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  base_buf);
    if (err < 0) CALLBACK_ERROR("H5Dread", base_name)

    if ((err = H5Dclose(base_id))  < 0) CALLBACK_ERROR("H5Dclose", base_name)

    /* open dataset run in this group */
    run_id = H5Dopen(grp_id, "run", H5P_DEFAULT);
    if (run_id < 0) CALLBACK_ERROR("H5Dopen", "run")

    /* read the entire dataset run */
    err = H5Dread(run_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  run_buf);
    if (err < 0) CALLBACK_ERROR("H5Dread", "run")

    if ((err = H5Dclose(run_id))  < 0) CALLBACK_ERROR("H5Dclose", "run")

    /* open datasets subrun in this group */
    srun_id = H5Dopen(grp_id, "subrun", H5P_DEFAULT);
    if (srun_id < 0) CALLBACK_ERROR("H5Dopen", "subrun")

    /* read the entire dataset subrun */
    err = H5Dread(srun_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  srun_buf);
    if (err < 0) CALLBACK_ERROR("H5Dread", "subrun")

    if ((err = H5Dclose(srun_id)) < 0) CALLBACK_ERROR("H5Dclose", "subrun")

    GET_TIMER(ts, te, timing[0])

    /* allocate buffer for dataset seq */
    seq_buf = (int64_t*) malloc(dset_dims[0] * sizeof(int64_t));

    /* table look up the seq values */
    for (i=0; i<dset_dims[0]; i++)
        seq_buf[i] = lookup_table[key(run_buf[i], srun_buf[i], base_buf[i])];

    free(run_buf);

    GET_TIMER(ts, te, timing[1])

seq_create:
    if (!dry_run) {
        char key_name[1024];
        sprintf(key_name, "%s.seq", base_name);

        /* create a new dataset */
        seq_id = H5Dcreate2(grp_id, key_name, H5T_STD_I64LE, space_id,
                            H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        if (seq_id < 0) CALLBACK_ERROR("H5Dcreate2", key_name)

        /* write seq (entire dataset) */
        if (dset_dims[0] > 0) {
            err = H5Dwrite(seq_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, seq_buf);
            if (err < 0) CALLBACK_ERROR("H5Dwrite", key_name)
        }

        GET_TIMER(ts, te, timing[2])

        /* close dataset */
        err = H5Dclose(seq_id);
        if (err < 0) CALLBACK_ERROR("H5Dclose", key_name)
    }
    if (dset_dims[0] > 0) free(seq_buf);

    /* close create property */
    err = H5Pclose(dcpl_id);
    if (err < 0) RETURN_ERROR("H5Pclose", base_name)
    /* close space */
    err = H5Sclose(space_id);
    if (err < 0) CALLBACK_ERROR("H5Sclose", base_name)
    /* close group */
    err = H5Gclose(grp_id);
    if (err < 0) CALLBACK_ERROR("H5Gclose", grp_name)

    GET_TIMER(ts, te, timing[3])

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
         printf("%s %4zd: type %s, name %s\n",
                filename, ii, type_name.c_str(), obj_name);
    }
    free(objs);
    return howmany;
}

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  [-h]          print this command usage message\n\
  [-v]          verbose mode (default: off)\n\
  [-n]          dry run without creating key datasets (default: disabled)\n\
  -k base_name  dataset name in group /spill to generate partitioning keys (required)\n\
  file_name     input/output HDF5 file name (required)\n\n\
  This utility program adds a new dataset in each group of the input file.\n\
  The new dataset, referred as the partition key dataset and to be named as\n\
  'base_name.seq', can be used for data partitioning in parallel read\n\
  operations. Its contents are generated based on the dataset 'base_name' in\n\
  group '/spill'. This base dataset must contain a list of unique integer\n\
  values, stored in an increasing order. An example is the dataset\n\
  '/spill/evt'. The data partitioning strategy for parallel reads is to\n\
  assign the dataset elements with the same 3-tuple of 'run', 'subrun', and\n\
  the base dataset to the same MPI process. Thus the partition key dataset\n\
  created in the output file stores a list of unique IDs corresponding to\n\
  the unique 3-tuples. The unique IDs are consistent among datasets across\n\
  all groups. Requirements for the HDF5 file:\n\
    1. must contain datasets /spill/run and /spill/subrun\n\
    2. contains multiple groups at root level\n\
    3. each group may contain multiple 2D datasets\n\
    4. all datasets in the same group share the 1st dimension size\n\
    5. each group must contain datasets run, subrun, and 'base_name'\n\
    6. the second dimension size of the 3 datasets must be 1\n\
    7. data type of the 3 datasets must be H5T_STD_U32LE\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v] -k base_name file_name\n%s\n", progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    bool dry_run=false;
    int c, err_exit=0, ndims;
    char msg[1024], *infile=NULL, dset_name[1024], *part_key_base=NULL;
    herr_t err;
    hid_t file_id=-1, fapl_id=-1, run_id, srun_id, base_id, space_id;
    hsize_t dset_dims[2], maxdims[2];
    unsigned int *run_buf=NULL, *srun_buf, *base_buf;
    H5G_info_t grp_info;
    size_t i, num_orig_groups=0, num_nonzero_groups=0;
    size_t max_seq=0, min_seq=0, avg_seq=0;
    struct op_data it_op;
    double timing[6], stime[6];
    table lookup_table;

#if defined PROFILE && PROFILE
    double ts=0.0, te=0.0;
#endif
    for (c=0; c<6; c++) timing[c] = 0.0;
    for (c=0; c<6; c++) stime[c] = 0.0;
    it_op.num_groups = 0;
    it_op.grp_names  = NULL;

    verbose = 0; /* default is quiet */

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvnk:")) != -1)
        switch(c) {
            case 'h': usage(argv[0]);
                      err_exit = -1;
                      goto fn_exit;
            case 'v': verbose = 1;
                      break;
            case 'k': part_key_base = strdup(optarg);
                      break;
            case 'n': dry_run = true;
                      break;
            default: break;
        }

    if (part_key_base == NULL) { /* key base dataset name is mandatory */
        printf("Error: partition key base dataset name is missing.\n\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }

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

    /* increase cache size for raw data */
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

    GET_TIMER(ts, te, timing[0])
#if defined PROFILE && PROFILE
    printf("Open input file                   = %7.2f sec\n", timing[0]);
#endif

    /* retrieve the number of groups */
    err = H5Gget_info_by_name(file_id, "/", &grp_info, H5P_DEFAULT);
    if (err < 0) RETURN_ERROR("H5Gget_info_by_name", "root group")
    num_orig_groups = grp_info.nlinks;
    it_op.grp_names = (char**) malloc(num_orig_groups * sizeof(char*));

    /* Iterate all objects to collect names of all group */
#if defined HAS_H5OVISIT3 && HAS_H5OVISIT3
    err = H5Ovisit3(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, gather_grp_names,
                    &it_op, H5O_INFO_ALL);
    if (err < 0) HANDLE_ERROR("H5Ovisit3")
#else
    err = H5Ovisit(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, gather_grp_names,
                   &it_op);
    if (err < 0) HANDLE_ERROR("H5Ovisit")
#endif
    if (it_op.err < 0) HANDLE_ERROR("H5Ovisit")

    GET_TIMER(ts, te, timing[1])
#if defined PROFILE && PROFILE
    printf("Collect group names               = %7.2f sec\n", timing[1]);
#endif

    /* open dataset /spill/run, /spill/subrun, /spill/base */
    run_id = H5Dopen(file_id, "/spill/run", H5P_DEFAULT);
    if (run_id < 0) RETURN_ERROR("H5Dopen", "/spill/run")
    srun_id = H5Dopen(file_id, "/spill/subrun", H5P_DEFAULT);
    if (srun_id < 0) RETURN_ERROR("H5Dopen", "/spill/subrun")
    sprintf(dset_name, "/spill/%s", part_key_base);
    base_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    if (base_id < 0) RETURN_ERROR("H5Dopen", dset_name)

    /* collect dimension sizes of /spill/base */
    space_id = H5Dget_space(base_id);
    if (space_id < 0) RETURN_ERROR("H5Dget_space", dset_name)
    ndims = H5Sget_simple_extent_dims(space_id, dset_dims, maxdims);
    if (ndims < 0) RETURN_ERROR("H5Sget_simple_extent_dims", dset_name)
    err = H5Sclose(space_id);
    if (err < 0) RETURN_ERROR("H5Sclose", dset_name)

    /* allocate buffers for /spill/run, /spill/subrun, /spill/base
     * All 3 are of the same size, as they belong to the same group.
     */
    run_buf  = (unsigned int*) malloc(dset_dims[0] * 3 * sizeof(unsigned int));
    srun_buf = run_buf  + dset_dims[0];
    base_buf = srun_buf + dset_dims[0];

    /* read the entire dataset /spill/run */
    err = H5Dread(run_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  run_buf);
    if (err < 0) RETURN_ERROR("H5Dread", "/spill/run")
    /* read the entire dataset /spill/subrun */
    err = H5Dread(srun_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  srun_buf);
    if (err < 0) RETURN_ERROR("H5Dread", "/spill/subrun")
    /* read the entire dataset /spill/base */
    err = H5Dread(base_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                  base_buf);
    if (err < 0) RETURN_ERROR("H5Dread", dset_name)

    /* close datasets */
    if ((err = H5Dclose(run_id))  < 0) RETURN_ERROR("H5Dclose", "/spill/run")
    if ((err = H5Dclose(srun_id)) < 0) RETURN_ERROR("H5Dclose", "/spill/subrun")
    if ((err = H5Dclose(base_id)) < 0) RETURN_ERROR("H5Dclose", dset_name)

    GET_TIMER(ts, te, timing[2])
#if defined PROFILE && PROFILE
    printf("Read 'spill' 3-tuples             = %7.2f sec\n", timing[2]);
#endif

    /* build a lookup table based on datasets run, subrun, and base */
    lookup_table = build_lookup_table(dset_dims[0], run_buf, srun_buf,
                                      base_buf);

    if (run_buf != NULL) free(run_buf);

    GET_TIMER(ts, te, timing[3])
#if defined PROFILE && PROFILE
    printf("Construct lookup table            = %7.2f sec\n", timing[3]);
#endif

    /* Iterate all groups and add a new dataset named 'part_key_base'.seq to
     * each group.
     */
    for (i=0; i<num_orig_groups; i++) {
        err = add_seq(dry_run, file_id, it_op.grp_names[i], part_key_base,
                      lookup_table, &num_nonzero_groups, &max_seq, &min_seq,
                      &avg_seq, stime);
        if (err < 0) RETURN_ERROR("add_seq", it_op.grp_names[i])
    }
    if (num_nonzero_groups > 0)
        avg_seq /= num_nonzero_groups;

    assert(num_orig_groups >= num_nonzero_groups);

    GET_TIMER(ts, te, timing[4])

fn_exit:
    if (it_op.grp_names != NULL) {
        for (i=0; i<num_orig_groups; i++)
            free(it_op.grp_names[i]);
        free(it_op.grp_names);
    }
    if (file_id >= 0) {
        check_h5_objects(infile, file_id);
        err = H5Fclose(file_id);
        if (err < 0) printf("Error at line %d: H5Fclose\n",__LINE__);
    }
    if (fapl_id >= 0) {
        err = H5Pclose(fapl_id);
        if (err < 0) printf("Error at line %d: H5Pclose\n",__LINE__);
    }

    GET_TIMER(ts, te, timing[5])

#if defined PROFILE && PROFILE
    if (err_exit == 0) {
        printf("-------------------------------------------------------\n");
        printf("Dry-run mode                      = %s\n", dry_run?"YES":"NO");
        printf("Input file name                   = %s\n", infile);
        printf("Partition key base dataset name   = %s\n", part_key_base);
        printf("number of groups in the file      = %zd\n", num_orig_groups);
        printf("number of non-zero size groups    = %zd\n", num_nonzero_groups);
        printf("Partition key name                = %s.seq\n",part_key_base);
        printf("  max length among all groups     = %zd\n", max_seq);
        printf("  min length among all groups     = %zd\n", min_seq);
        printf("  avg length among all groups     = %zd\n", avg_seq);
        printf("-------------------------------------------------------\n");
        printf("Open input file                   = %7.2f sec\n", timing[0]);
        printf("Collect group names               = %7.2f sec\n", timing[1]);
        printf("Read 'spill' 3-tuples             = %7.2f sec\n", timing[2]);
        printf("Construct lookup table            = %7.2f sec\n", timing[3]);
        printf("Add partition key datasets        = %7.2f sec\n", timing[4]);
        printf("  Read run, subrun, base datasets = %7.2f sec\n", stime[0]);
        printf("  Hash table lookup               = %7.2f sec\n", stime[1]);
        printf("  Write partition key seq         = %7.2f sec\n", stime[2]);
        printf("  Close partition key seq         = %7.2f sec\n", stime[3]);
        printf("File close                        = %7.2f sec\n", timing[5]);
        printf("-------------------------------------------------------\n");
        printf("End-to-end                        = %7.2f sec\n",
               timing[0] + timing[1] + timing[2] + timing[3] +
               timing[4] + timing[5]);
    }
#endif
    if (infile != NULL) free(infile);
    if (part_key_base != NULL) free(part_key_base);

    return (err_exit != 0);
}
