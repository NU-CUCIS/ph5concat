/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator
 * Laboratory
 * See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strdup() */
#include <unistd.h> /* getopt() */
#include <assert.h> /* assert() */

#include <string>
#include <vector>

#include <hdf5.h>

using namespace std;

static int verbose;

#define RETURN_ERROR(func_name, obj_name) {   \
    cout<<"Error at line "<<__LINE__<<" calling "<<func_name<<" for object "<<obj_name<<endl; \
    return -1;                                \
}

#define HANDLE_ERROR(msg) {                                   \
    cout<<"Error at line "<<__LINE__<<" failed: "<<msg<<endl; \
    err_exit = -1;                                            \
    goto fn_exit;                                             \
}

/* string-based basename */
string basename(string full_path)
{
    auto pos = full_path.find_last_of("/\\");
    if(pos != string::npos)
      return full_path.substr(pos+1);
    else
      return full_path;
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
    if (err < 0) HANDLE_ERROR("H5Pget_cache");

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
    if (err < 0) HANDLE_ERROR("H5Pset_cache");
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
        printf("Still opened in file %s %4zd: type %s, name %s\n",
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
  [-n]          dry-run mode (default: off)\n\
  -s src_path   full path of dataset who's first element's value will be used\n\
                to populate the new dataset in group '/spill' (required)\n\
  file_name     input/output HDF5 file name (required)\n\n\
  This utility program adds a new dataset in group '/spill' of the input\n\
  file. Argument 'src_path' should be in the form of '/group/dset'. The new\n\
  dataset to be created will be '/spill/dset', whose contents will be\n\
  single-valued, populated from the first element of '/group/dset'. This new\n\
  dataset is intended to be used as an additional index dataset for\n\
  generating partition key datasets. Requirements for the HDF5 file:\n\
    1. must contain group '/spill'\n\
    2. contains multiple groups at root level\n\
    3. each group may contain multiple 2D datasets\n\
    4. all datasets in the same group share the 1st dimension size\n\
    5. dataset src_path must exist, except for the one in group '/spill'\n\
    6. the second dimension size of the dataset 'dset' must be 1\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v|-n] -s src_path file_name\n%s\n", progname, USAGE);
}

int main(int argc, char **argv)
{
    char *infile=NULL, *src_path=NULL, *newsrc_name;
    int c, err_exit=0, ndims;
    double inVal;
    hid_t file_id=-1, fapl_id=-1, dcpl_id, dset, mspace, fspace;
    hid_t spill_id=-1, space_id=-1, run_id, newsrc_id, dtype;
    htri_t src_exist;
    herr_t err;
    hsize_t one[2]={1,1}, offs[2]={0,0}, lens[2]={1,1};
    hsize_t dset_dims[2], maxdims[2];
    bool dry_run=false;

    verbose = 0; /* default is quiet */

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvns:")) != -1) {
        switch(c) {
            case 'n': dry_run = true;
                      break;
            case 's': src_path = strdup(optarg);
                      break;
            case 'v': verbose = 1;
                      break;
            case 'h': usage(argv[0]);
                      return 0;
            default:  usage(argv[0]);
                      return 1;
        }
    }
    if (argv[optind] == NULL) { /* input file name is mandatory */
        fprintf(stderr,"Error: input file name is missing\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    if (src_path == NULL) { /* src_path dataset is required */
        fprintf(stderr,"Error: src_path dataset is required\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }

    infile = strdup(argv[optind]);
    if (verbose) printf("input file: %s\n", infile);

    /* create file access property for read and write */
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id < 0) HANDLE_ERROR("H5Pcreate");

    /* increase cache size for raw data */
    err = incr_raw_data_chunk_cache(fapl_id);
    if (err < 0) HANDLE_ERROR("incr_raw_data_chunk_cache");

    if (dry_run)
        /* read-only mode */
        file_id = H5Fopen(infile, H5F_ACC_RDONLY, fapl_id);
    else
        /* read-and-write mode  */
        file_id = H5Fopen(infile, H5F_ACC_RDWR, fapl_id);
    if (file_id < 0)
        HANDLE_ERROR(string("H5Fopen ")+infile);

    /* check if dataset src_path exists */
    src_exist = H5Lexists(file_id, src_path, H5P_DEFAULT);
    if (src_exist < 0) RETURN_ERROR("H5Lexists", src_path)
    if (src_exist == 0) {
        fprintf(stderr,"Error: src_path dataset '%s' does not exist\n", src_path);
        err_exit = -1;
        goto fn_exit;
    }

    /* create memspace */
    mspace = H5Screate_simple(2, lens, NULL);
    if (mspace < 0) HANDLE_ERROR("H5Screate_simple");

    /* open the src_path dataset */
    dset = H5Dopen(file_id, src_path, H5P_DEFAULT);
    if (dset < 0) RETURN_ERROR("H5Dopen", src_path);

    /* get the data type, as the new dataset will be of the same type */
    dtype = H5Dget_type(dset);
    if (dtype < 0) HANDLE_ERROR("H5Dget_type")

    /* check if this dataset is of zero size */
    fspace = H5Dget_space(dset);
    if (fspace < 0) RETURN_ERROR("H5Dget_space", src_path);
    ndims = H5Sget_simple_extent_dims(fspace, dset_dims, NULL);
    if (ndims < 0)
        HANDLE_ERROR(string("H5Sget_simple_extent_dims ")+src_path);
    if (dset_dims[0] == 0) {
        cerr << "Error at line "<<__LINE__<<": zero-sized dataset '"<<src_path<<"'"<<endl;
        err_exit = -1;
        goto fn_exit;
    }

    err = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offs, NULL, one, lens);
    if (err < 0) RETURN_ERROR("H5Sselect_hyperslab", src_path);

    /* read into an 8-type buffer (of type double) large enough for any type */
    err = H5Dread(dset, H5T_NATIVE_DOUBLE, mspace, fspace, H5P_DEFAULT, &inVal);
    if (err < 0) RETURN_ERROR("H5Dread", src_path);

    err = H5Dclose(dset);
    if (err < 0) RETURN_ERROR("H5Dclose", src_path);

    err = H5Sclose(mspace);
    if (err < 0) HANDLE_ERROR("H5Sclose")

    err = H5Sclose(fspace);
    if (err < 0) HANDLE_ERROR("H5Sclose")

    if (verbose) /* print values we retrieved from the src_path datasets */
        printf("src_path dataset %s = %e\n", src_path, inVal);

    /* open spill group */
    spill_id = H5Gopen(file_id, "/spill", H5P_DEFAULT);
    if (spill_id < 0) RETURN_ERROR("H5Gopen", "/spill");

    /* open run dataset in spill group to determine dataset size */
    run_id = H5Dopen(spill_id, "run", H5P_DEFAULT);
    if (run_id < 0) RETURN_ERROR("H5Dopen", "run");

    /* Get dimension sizes of dataset /spill/run */
    space_id = H5Dget_space(run_id);
    if (space_id < 0) RETURN_ERROR("H5Dget_space", "run");
    ndims = H5Sget_simple_extent_dims(space_id, dset_dims, maxdims);
    if (ndims != 2) RETURN_ERROR("H5Sget_simple_extent_dims", "run");

    /* Get create property of run dataset when it was created */
    dcpl_id = H5Dget_create_plist(run_id);
    if (dcpl_id < 0) RETURN_ERROR("H5Dget_create_plist", "run");

    /* close /spill/run */
    err = H5Dclose(run_id);
    if (err < 0) RETURN_ERROR("H5Dclose", "run");

    if (verbose) printf("/spill/run dims (%llu, 1)\n", dset_dims[0]);

    /* /spill/run should not be of zero-sized */
    if (dset_dims[0] == 0) RETURN_ERROR("Zero-sized group '/spill'", "run");

    newsrc_name = basename(src_path);
    src_exist = H5Lexists(spill_id, newsrc_name, H5P_DEFAULT);
    if (src_exist < 0) RETURN_ERROR("H5Lexists", newsrc_name)
    if (src_exist > 0) {
        fprintf(stderr,"Error: dataset '%s' already exists in group '/spill'\n",
                newsrc_name);
        err_exit = -1;
        goto fn_exit;
    }

    if (!dry_run) {
        /* Fill the entire new dataset with the same value. Data type casting
         * will be internally triggered by HDF5, if newsrc_name is not of
         * type double.
         */
        err = H5Pset_fill_value(dcpl_id, H5T_NATIVE_DOUBLE, &inVal);
        if (err < 0) RETURN_ERROR("H5Pset_fill_value", newsrc_name);

        err = H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_ALLOC);
        if (err < 0) RETURN_ERROR("H5Pset_fill_time", newsrc_name);

        /* create new dataset */
        newsrc_id = H5Dcreate2(spill_id, newsrc_name, dtype, space_id,
                               H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        if (newsrc_id < 0) RETURN_ERROR("H5Dcreate2", newsrc_name);

        /* close dataset */
        err = H5Dclose(newsrc_id);
        if (err < 0) RETURN_ERROR("H5Dclose", newsrc_name);
    }

    err = H5Pclose(dcpl_id);
    if (err < 0) HANDLE_ERROR("H5Pclose");
    err = H5Tclose(dtype);
    if (err < 0) HANDLE_ERROR("H5Tclose");

fn_exit:
    if (space_id != -1) {
        err = H5Sclose(space_id);
        if (err < 0) HANDLE_ERROR("H5Sclose")
    }
    if (spill_id != -1) {
        /* close spill group */
        err = H5Gclose(spill_id);
        if (err < 0) RETURN_ERROR("H5Gclose", "/spill");
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
    if (src_path != NULL) free(src_path);
    if (infile   != NULL) free(infile);

    return (err_exit == -1);
}

