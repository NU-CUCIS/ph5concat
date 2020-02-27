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
#include <libgen.h> /* basename() */
#include <assert.h> /* assert() */

#include <iostream>
#include <fstream>
#include <cerrno>
#include <vector>
#include <string>
#include <tuple>
#include <map>
using namespace std;

#include <hdf5.h>

static int verbose, debug;

#define HANDLE_ERROR(msg) { \
    printf("Error at line %d: %s\n",__LINE__, msg); \
    err_exit = -1; \
    goto fn_exit; \
}

#define CALLBACK_ERROR(func_name, dset_name) { \
    printf("Error in %s line %d: calling %s for dataset %s\n",__FILE__,__LINE__,func_name, dset_name); \
    err_exit = -1; \
    goto fn_exit; \
}

/* compare for 2-tuple (run, subrun) */
typedef tuple<unsigned int, unsigned int> key;
struct comp2 : public unary_function<key, size_t> {
    bool operator()(const key& lhs, const key& rhs) const {
        if (get<0>(lhs) < get<0>(rhs)) return true;
        if (get<0>(lhs) > get<0>(rhs)) return false;
        if (get<1>(lhs) < get<1>(rhs)) return true;
        return false;
    }
};
typedef map<key, string, comp2> tuple_map;

/* parameters to be passed to the call back function */
struct op_data {
    bool is_set_run;
    bool is_set_subrun;
    unsigned int run;
    unsigned int subrun;
    herr_t   err;
};

/*----< get_IDs() >----------------------------------------------------------*/
/* call back function used in H5Ovisit() */
static
herr_t get_IDs(hid_t             loc_id,/* object ID */
               const char       *name,  /* object name */
               const H5O_info_t *info,  /* object metadata */
               void             *op)    /* data passed from caller */
{
    herr_t err=0;
    hid_t dset_id, space_id;
    hsize_t dset_dims[2];
    int ndims, err_exit=0;
    size_t i;
    struct op_data *it_op = (struct op_data*)op;
    unsigned int *buf;

    it_op->err = 0;

    if (info->type != H5O_TYPE_DATASET) return 0;

    /* remove group name */
    char *path = strdup(name);
    char *dset_name = basename(path);
    if (strcmp(dset_name, "run") != 0 && strcmp(dset_name, "subrun") != 0)
        goto fn_exit; /* skip dataset that is not 'run' or 'subrun' */

    /* Open the dataset. Note that loc_id is not the dataset/group ID. */
    dset_id = H5Dopen(loc_id, name, H5P_DEFAULT);
    if (dset_id < 0) CALLBACK_ERROR("H5Dopen",name)

    /* Get dimension sizes */
    space_id = H5Dget_space(dset_id);
    if (space_id < 0) CALLBACK_ERROR("H5Dget_space",name)
    ndims = H5Sget_simple_extent_dims(space_id, dset_dims, NULL);
    if (ndims < 0) CALLBACK_ERROR("H5Sget_simple_extent_dims",name)
    err = H5Sclose(space_id);
    if (err < 0) CALLBACK_ERROR("H5Sclose",name)

    /* skip zero-size dataset */
    if (dset_dims[0] == 0) {
        err = H5Dclose(dset_id);
        if (err < 0) CALLBACK_ERROR("H5Dclose",name)
        goto fn_exit;
    }

    /* allocate read buffer */
    buf = new unsigned int [dset_dims[0]];

    /* read the entire dataset */
    err = H5Dread(dset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    if (err < 0) CALLBACK_ERROR("H5Dread",name)
    err = H5Dclose(dset_id);
    if (err < 0) CALLBACK_ERROR("H5Dclose",name)

    if (strcmp(dset_name, "run") == 0) {
        /* get the first value of run[0] */
        if (! it_op->is_set_run) {
            it_op->is_set_run = true;
            it_op->run = buf[0];
        }
        /* check run[] in all groups for consistency */
        for (i=0; i<dset_dims[0]; i++) {
            if (buf[i] != it_op->run) {
                printf("Error: inconsistent run ID %u, expecting %u\n",
                       buf[i], it_op->run);
                err_exit = -1;
                goto fn_exit;
            }
        }
    }
    else { /* dataset "subrun" */
        if (! it_op->is_set_subrun) {
            it_op->is_set_subrun = true;
            it_op->subrun = buf[0];
        }
        /* check subrun[] in all groups for consistency */
        for (i=0; i<dset_dims[0]; i++) {
            if (buf[i] != it_op->subrun) {
                printf("Error: inconsistent subrun ID %u, expecting %u\n",
                       buf[i], it_op->subrun);
                err_exit = -1;
                goto fn_exit;
            }
        }
    }
    delete [] buf;

fn_exit:
    free(path);

    it_op->err = err_exit;
    return (err_exit == 0) ? 0 : 1;
}

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  [-h]          print this command usage message\n\
  [-v]          verbose mode (default: off)\n\
  [-d]          debug mode (default: off)\n\
  [-o outfile]  output file name (default: 'out_list.txt')\n\
  infile        input file name contains a list of HDF5 file names (required)\n\n\
  This utility program re-order the files in infile into a sorted list based\n\
  on the increasing order of 'run' and 'subrun' IDs. Requirements for the\n\
  input HDF5 files:\n\
    1. must contain datasets '/spill/run' and '/spill/subrun'\n\
    2. may contain multiple groups at root level\n\
    3. each group may contain multiple 2D datasets\n\
    4. 1st dimension of all datasets in the same group share the same size\n\
    5. each group must contain datasets 'run' and 'subrun'\n\
    6. data type of datasets 'run' and 'subrun' must be H5T_STD_U32LE\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v|-d|-o outfile] infile\n%s\n", progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int c, err_exit=0;
    char *infile=NULL, *outfile=NULL;
    herr_t err;
    hid_t file_id=-1;
    size_t i;
    struct op_data it_op;
    ifstream in_fd;
    ofstream out_fd;
    string line;
    vector<string> in_list;
    tuple_map file_list;

    verbose = 0; /* default is quiet */
    debug   = 0; /* default is no */

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvdo:")) != -1)
        switch(c) {
            case 'h': usage(argv[0]);
                      err_exit = -1;
                      goto fn_exit;
            case 'v': verbose = 1;
                      break;
            case 'd': debug = 1;
                      break;
            case 'o': outfile = strdup(optarg);
                      break;
            default: break;
        }

    if (argv[optind] == NULL) { /* input file name is mandatory */
        printf("input file name is missing\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    infile = strdup(argv[optind]);
    if (outfile == NULL) outfile = strdup("out_list.txt");

    if (verbose) printf("input list file: %s\n", infile);

    /* open input file and catch error */
    try {
        in_fd.open(infile);
        if (!in_fd)
            throw ios_base::failure(strerror(errno));
    }
    catch (ifstream::failure& e) {
        cerr << "Error: opening file \""<<infile<<"\" (" << e.what() 
             << ")" << endl;
        err_exit = -1;
        goto fn_exit;
    }

    /* read input file contents */
    while (getline(in_fd, line)) {
        if (line.length() == 0)
            continue; /* skip empty lines */
        if (line.at(0) == '#')
            continue; /* skip comment line (start with #) */
        /* save file name to in_list */
        in_list.push_back(line);
    }
    in_fd.close();

    if (verbose)
        for (i=0; i<in_list.size(); i++)
            cout << "input HDF5 file: "<< in_list[i] << '\n'; 

    for (i=0; i<in_list.size(); i++) {
        /* open file in read-only mode */
        file_id = H5Fopen(in_list[i].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {
            printf("Error at line %d: H5Fopen %s\n",__LINE__,
                   in_list[i].c_str());
            continue;
        }

        if (debug) {
            /* Iterate all objects to collect run and subrun IDs */
            it_op.is_set_run = false;
            it_op.is_set_subrun = false;
            err = H5Ovisit(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, get_IDs,
                           &it_op);
            if (err < 0) HANDLE_ERROR("H5Ovisit")
            if (it_op.err < 0) HANDLE_ERROR("H5Ovisit")

            if (verbose)
                printf("File %zd: run ID = %d subrun ID = %u\n",
                       i, it_op.run,it_op.subrun);

            /* use sorted map in an increasing order of run and subrun */
            file_list[key(it_op.run, it_op.subrun)] = in_list[i];
        }
        else {
            unsigned int run, subrun;
            hsize_t offs[2]={0,0}, lens[2]={1,1};
            hid_t dset_id, mem_space_id, file_space_id;

            mem_space_id = H5Screate_simple(2, lens, NULL);
            if (mem_space_id < 0) HANDLE_ERROR("H5Screate_simple")

            /* open /spill/run */
            dset_id = H5Dopen(file_id, "/spill/run", H5P_DEFAULT);
            if (dset_id < 0) HANDLE_ERROR("H5Dopen")

            file_space_id = H5Dget_space(dset_id);
            if (file_space_id < 0) HANDLE_ERROR("H5Dget_space")
            err = H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offs, NULL,
                                      lens, NULL);
            if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab");

            err = H5Dread(dset_id, H5T_NATIVE_UINT, mem_space_id, file_space_id,
                          H5P_DEFAULT, &run);
            if (err < 0) HANDLE_ERROR("H5Dread")

            err = H5Dclose(dset_id);
            if (err < 0) HANDLE_ERROR("H5Dclose")
            H5Sclose(file_space_id);

            /* open /spill/subrun */
            dset_id = H5Dopen(file_id, "/spill/subrun", H5P_DEFAULT);
            if (dset_id < 0) HANDLE_ERROR("H5Dopen")

            file_space_id = H5Dget_space(dset_id);
            if (file_space_id < 0) HANDLE_ERROR("H5Dget_space")
            err = H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offs, NULL,
                                      lens, NULL);
            if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab");

            err = H5Dread(dset_id, H5T_NATIVE_UINT, mem_space_id, file_space_id,
                          H5P_DEFAULT, &subrun);
            if (err < 0) HANDLE_ERROR("H5Dread")

            err = H5Dclose(dset_id);
            if (err < 0) HANDLE_ERROR("H5Dclose")

            H5Sclose(file_space_id);
            H5Sclose(mem_space_id);

            if (verbose)
                printf("File %zd: run ID = %d subrun ID = %u\n",i,run,subrun);

            /* use sorted map in an increasing order of run and subrun */
            file_list[key(run, subrun)] = in_list[i];
        }
        err = H5Fclose(file_id);
        if (err < 0)
            printf("Error at line %d: H5Fclose %s\n",__LINE__,
                   in_list[i].c_str());
    }

    /* create the output file */
    try {
        out_fd.open(outfile, ios::out);
        if (!out_fd)
            throw ios_base::failure(strerror(errno));
    }
    catch (ifstream::failure& e) {
        cerr << "Error: creating file \""<<outfile<<"\" (" << e.what() 
             << ")" << endl;
        err_exit = -1;
        goto fn_exit;
    }

    if (verbose)
        printf("Number of files in output file list=%zd\n",file_list.size());

    /* print the sorted file names to the output file */
    for (auto f = file_list.begin(); f != file_list.end(); f++)
        out_fd << f->second << '\n'; 

    out_fd.close();

fn_exit:
    if (outfile != NULL) free(outfile);
    if (infile  != NULL) free(infile);

    return (err_exit != 0);
}
