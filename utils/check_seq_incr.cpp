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
#include <libgen.h> /* dirname(), basename() */
#include <assert.h> /* assert() */

#include <iostream>
#include <fstream>
#include <cerrno>
#include <string>
using namespace std;

#include <hdf5.h>

static int verbose;

#define HANDLE_ERROR(func_name, dset_name) { \
    printf("Error in %s line %d: calling %s for dataset %s\n",__FILE__,__LINE__,func_name, dset_name); \
    err_exit = -1; \
    goto fn_exit; \
}

/* parameters to be passed to the call back function */
struct op_data {
    char         *dname;
    unsigned int  max_id;
    herr_t        err;
};

/*----< check_order() >------------------------------------------------------*/
/* call back function used in H5Ovisit() */
static
herr_t check_order(hid_t             loc_id,/* object ID */
                   const char       *name,  /* object name */
                   const H5O_info_t *info,  /* object metadata */
                   void             *op)    /* data passed from caller */
{
    char *path=NULL, *dset_name, *grp_name, *base_grp, *dname;
    herr_t err=0;
    hid_t dset_id, space_id;
    hsize_t dset_dims[2];
    int ndims, err_exit=0;
    size_t i;
    struct op_data *it_op = (struct op_data*)op;
    unsigned int *buf=NULL;

    it_op->err = 0;

    /* skip if this is not a dataset */
    if (info->type != H5O_TYPE_DATASET) return 0;

    /* strip group name */
    path = strdup(name);
    dset_name = basename(path);
    grp_name = dirname(path);

    base_grp = dirname(it_op->dname);
    dname    = basename(it_op->dname);

    if (strcmp(grp_name, base_grp) == 0)
        goto fn_exit; /* check for /spill has been done in main() */

    if (strcmp(dset_name, dname) != 0)
        goto fn_exit; /* skip dataset that is not dname */

    /* Open the dataset. Note that loc_id is not the dataset/group ID. */
    dset_id = H5Dopen(loc_id, name, H5P_DEFAULT);
    if (dset_id < 0) HANDLE_ERROR("H5Dopen",name)

    /* Get dimension sizes */
    space_id = H5Dget_space(dset_id);
    if (space_id < 0) HANDLE_ERROR("H5Dget_space",name)
    ndims = H5Sget_simple_extent_dims(space_id, dset_dims, NULL);
    if (ndims < 0) HANDLE_ERROR("H5Sget_simple_extent_dims",name)
    err = H5Sclose(space_id);
    if (err < 0) HANDLE_ERROR("H5Sclose",name)

    /* skip zero-size dataset */
    if (dset_dims[0] == 0) {
        err = H5Dclose(dset_id);
        if (err < 0) HANDLE_ERROR("H5Dclose",name)
        goto fn_exit;
    }

    if (verbose)
        printf("Checking non-zero dataset %s in group %s\n",dset_name,grp_name);

    /* allocate read buffer */
    buf = new unsigned int [dset_dims[0]];

    /* read the entire dataset */
    err = H5Dread(dset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    if (err < 0) HANDLE_ERROR("H5Dread",name)
    err = H5Dclose(dset_id);
    if (err < 0) HANDLE_ERROR("H5Dclose",name)

    /*  check first element */
    if (buf[0] < 0 || buf[0] > it_op->max_id) {
        printf("Error: /%s/%s[0]=%u expected to be >=0 and <= %u\n",
               path,dset_name,buf[0], it_op->max_id);
        err_exit = -1;
        goto fn_exit;
    }
    for (i=1; i<dset_dims[0]; i++) {
        /* check if smaller than max_id */
        if (buf[i] > it_op->max_id) {
            printf("Error: /%s/%s[%zd]=%u expected to be <= %u\n",
                   path,dset_name,i,buf[i], it_op->max_id);
            err_exit = -1;
            goto fn_exit;
        }
        /* check if monotonically nondecreasing */
        if (buf[i] < buf[i-1]) {
            printf("Error: /%s/%s[%zd]=%u expects >= /%s/%s[%zd]=%u\n",
                   path,dset_name,i,buf[i], path,dset_name,i-1,buf[i-1]);
            err_exit = -1;
            goto fn_exit;
        }
    }

fn_exit:
    if (buf != NULL) delete [] buf;
    if (path != NULL) free(path);

    it_op->err = err_exit;
    return (err_exit == 0) ? 0 : 1;
}

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  [-h]       print this command usage message\n\
  [-v]       verbose mode (default: off)\n\
  [-d name]  partition key dataset name (default: '/spill/evt.seq')\n\
  infile     name of input HDF5 file (required)\n\n\
  This utility program checks the contents of dataset 'name' in each group of\n\
  the input file for whether the values are in a monotonically nondecreasing\n\
  order. In particular, the partition dataset specified at command-line option\n\
  '-d' is checked for whether the values start from 0 and increment by 1.\n\
  The same datasets in all other groups are checked only for a monotonically\n\
  nondecreasing order. Requirements for the input HDF5 file:\n\
    1. contains multiple groups at root level\n\
    2. each group may contain multiple 2D datasets\n\
    3. the 1st  dimension of all datasets in the same group share same size\n\
    4. each group must contain a partition key dataset named 'name'\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v|-d name] infile\n%s\n", progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int c, err_exit=0, ndims;
    char *infile=NULL, *dname=NULL;
    herr_t err;
    hid_t file_id=-1, dset_id, space_id;
    hsize_t i, dset_dims[2];
    unsigned int *buf=NULL;
    struct op_data it_op;

    verbose = 0; /* default is quiet */

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvd:")) != -1)
        switch(c) {
            case 'v': verbose = 1;
                      break;
            case 'd': dname = strdup(optarg);
                      break;
            case 'h': usage(argv[0]);
                      return 0;
            default:  usage(argv[0]);
                      return 1;
        }

    if (argv[optind] == NULL) { /* input file name is mandatory */
        printf("input file name is missing\n");
        usage(argv[0]);
        return 1;
    }
    infile = strdup(argv[optind]);
    if (verbose) printf("input list file: %s\n", infile);

    if (dname == NULL) dname = strdup("/spill/evt.seq");
    if (verbose) printf("dataset name: %s\n", dname);

    /* open file in read-only mode  */
    file_id = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {
        printf("Error at line %d: H5Fopen %s\n",__LINE__,infile);
        err_exit = -1;
        goto fn_exit;
    }

    /* Open the dataset */
    dset_id = H5Dopen(file_id, dname, H5P_DEFAULT);
    if (dset_id < 0) HANDLE_ERROR("H5Dopen",dname)

    /* Get dimension sizes */
    space_id = H5Dget_space(dset_id);
    if (space_id < 0) HANDLE_ERROR("H5Dget_space",dname)
    ndims = H5Sget_simple_extent_dims(space_id, dset_dims, NULL);
    if (ndims < 0) HANDLE_ERROR("H5Sget_simple_extent_dims",dname)
    err = H5Sclose(space_id);
    if (err < 0) HANDLE_ERROR("H5Sclose",dname)

    if (dset_dims[1] > 1) {
        printf("Error at line %d: H5Fopen %s\n",__LINE__,infile);
        printf("\tkey dataset %s 2nd dimension is expected of size 1 but got %llu\n",
               dname, dset_dims[1]);
        err_exit = -1;
        goto fn_exit;
    }

    if (verbose) printf("Checking dataset %s\n",dname);

    /* allocate read buffer */
    buf = new unsigned int [dset_dims[0]];

    /* read the entire dataset */
    err = H5Dread(dset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    if (err < 0) HANDLE_ERROR("H5Dread",dname)
    err = H5Dclose(dset_id);
    if (err < 0) HANDLE_ERROR("H5Dclose",dname)

    /* check if starts with 0 */
    if (buf[0] != 0) {
        printf("Error: %s[0] expects 0 but got %u\n",dname,buf[0]);
        err_exit = -1;
        goto fn_exit;
    }
    /* check if increments by 1 */
    for (i=1; i<dset_dims[0]; i++) {
        if (buf[i] != i) {
            printf("Error: %s[%lld] expects %lld but got %u\n",
                   dname,i,i,buf[i]);
            err_exit = -1;
            goto fn_exit;
        }
    }
    it_op.max_id = dset_dims[0] - 1;
    delete [] buf;
    buf = NULL;

    /* Iterate all objects in the input file */
    it_op.dname = dname;
#if defined HAS_H5OVISIT3 && HAS_H5OVISIT3
    err = H5Ovisit3(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, check_order, &it_op, H5O_INFO_ALL);
    if (err < 0) HANDLE_ERROR("H5Ovisit3", dname)
#else
    err = H5Ovisit(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, check_order, &it_op);
    if (err < 0) HANDLE_ERROR("H5Ovisit", dname)
#endif
    if (it_op.err < 0) HANDLE_ERROR("H5Ovisit", dname)

    err = H5Fclose(file_id);
    if (err < 0)
        printf("Error at line %d: H5Fclose %s\n",__LINE__, infile);

fn_exit:
    if (buf    != NULL) delete [] buf;
    if (infile != NULL) free(infile);
    if (dname  != NULL) free(dname);

    if (err_exit == 0)
        printf("All datasets pass the check\n");
    else
        printf("Check fails\n");

    return (err_exit != 0);
}
