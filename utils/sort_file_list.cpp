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
    cout<<"Error at line "<<__LINE__<<" for file "<<in_list[i]<<": "<<msg<<" failed."<<endl; \
    err_exit = -1; \
    goto fn_exit; \
}

#define CALLBACK_ERROR(func_name, dset_name) { \
    printf("Error in %s line %d: calling %s for dataset %s\n",__FILE__,__LINE__,func_name, dset_name); \
    err_exit = -1; \
    goto fn_exit; \
}

struct compv {
    bool operator()(const vector<long long>& lhs, const vector<long long>& rhs) const
    {
        assert(lhs.size() == rhs.size());

        for (auto i = 0u; i < lhs.size(); i++) {
            if (lhs[i] < rhs[i]) return true;
            if (lhs[i] > rhs[i]) return false;
        }
        return false;
    }
};

typedef map<vector<long long>, string, compv> tuple_map;

/* parameters to be passed to the call back function */
struct op_data {
    vector<string> names;
    vector<bool> is_set;
    vector<long long> event_index;
    herr_t   err;
    op_data() {}
    op_data(vector<string> _names)
      : names(_names)
    {
      is_set      = vector<bool     >(names.size(), false);
      event_index = vector<long long>(names.size(), 0    );
    }
    void add_level(string name)
    {
      names      .push_back(name);
      is_set     .push_back(false);
      event_index.push_back(0);
    }
    void reset()
    {
      for(auto iset = 0u; iset < is_set.size(); iset++)
        is_set[iset] = false;
    }
    int index(string name) const
    {
      for(size_t idx = 0; idx < names.size(); idx++)
        if(strcmp(name.c_str(), names[idx].c_str()) == 0) return idx;
      return -1;
    }
};

/* string-based basename */
string basename(string full_path)
{
    auto pos = full_path.find_last_of("/\\");
    if(pos != string::npos)
      return full_path.substr(pos+1);
    else
      return full_path;
}

/*----< get_IDs() >----------------------------------------------------------*/
/* call back function used in H5Ovisit() */
static
herr_t get_IDs(hid_t             loc_id,/* object ID */
               const char       *name,  /* object name */
               const H5O_info_t *info,  /* object metadata */
               void             *op)    /* data passed from caller */
{
    herr_t err=0;
    hid_t dset, fspace;
    hsize_t dset_dims[2];
    int ndims, err_exit=0;
    size_t i;
    struct op_data *it_op = (struct op_data*)op;
    long long *buf;

    it_op->err = 0;

    /* skip objects that are not HDF5 datasets */
    if (info->type != H5O_TYPE_DATASET) return 0;

    /* remove group name */
    char *path = strdup(name);
    char *dset_name = basename(path);

    /* skip dataset that is not in list of index names
       Note: dset_name is compared to basename of the level index
             full path to ensure consistency within file.
             Importantly, the current group does not contain
             a matching dataset, this group will be skipped.
     */
    int level_idx = -1;
    for (i=0; i<it_op->names.size(); i++) {
        if ((string) dset_name == basename(it_op->names[i])) {
            level_idx = i;
            break;
        }
    }
    if (level_idx < 0) goto fn_exit;

    /* Open the dataset. Note that loc_id is not the dataset/group ID. */
    dset = H5Dopen(loc_id, name, H5P_DEFAULT);
    if (dset < 0) CALLBACK_ERROR("H5Dopen",name);

    /* Get dimension sizes */
    fspace = H5Dget_space(dset);
    if (fspace < 0) CALLBACK_ERROR("H5Dget_space",name);
    ndims = H5Sget_simple_extent_dims(fspace, dset_dims, NULL);
    if (ndims < 0) CALLBACK_ERROR("H5Sget_simple_extent_dims",name);
    err = H5Sclose(fspace);
    if (err < 0) CALLBACK_ERROR("H5Sclose",name);

    /* skip zero-size dataset */
    if (dset_dims[0] == 0) {
        err = H5Dclose(dset);
        if (err < 0) CALLBACK_ERROR("H5Dclose",name);
        goto fn_exit;
    }

    /* allocate read buffer */
    buf = new long long [dset_dims[0]];

    /* read the entire dataset */
    err = H5Dread(dset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    if (err < 0) CALLBACK_ERROR("H5Dread",name);
    err = H5Dclose(dset);
    if (err < 0) CALLBACK_ERROR("H5Dclose",name);

    /* get the first value of level */
    if (! it_op->is_set[level_idx]) {
        it_op->is_set     [level_idx] = true;
        it_op->event_index[level_idx] = buf[0];
    }

    /* check level in all groups for consistency */
    for (i=0; i<dset_dims[0]; i++) {
        if (buf[i] != it_op->event_index[level_idx]) {
            printf("Warn at line %d: inconsistent '%s' ID %lld, expecting %lld\n",
                   __LINE__, name, buf[i], it_op->event_index[level_idx]);
            break;
        }
    }

    delete [] buf;

fn_exit:
    free(path);

    it_op->err = err_exit;
    return (err_exit == 0) ? 0 : 1;
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
  [-d]          debug mode (default: off)\n\
  [-k paths]    full paths to datasets, separated by comma, to be used for\n\
                multi-index sorting. Only the first elements of the datasets\n\
                are used in the sorting. (default: /spill/run,/spill/subrun)\n\
  [-o outfile]  output file name (default: 'out_list.txt')\n\
  infile        input file containing a list of HDF5 file paths (required)\n\n\
  This utility program re-order the file paths given in file 'infile' into a\n\
  sorted list, based on the increasing order of index datasets specified in\n\
  the argument of command-line option '-k'. When option '-k' is not used, the\n\
  default are datasets '/spill/run' and '/spill/subrun'. An example of its\n\
  usage is '-k /spill/run,/spill/subrun,/rec.hdr/cycle'. The index datasets\n\
  will be read and type-cast into internal buffers of type 'long long int',\n\
  before the sorting is applied. Sorting follows the order of datasets\n\
  appearing in the argument of option '-k'. Note the contents of outfile will\n\
  be overwritten if it exists. Requirements for the input HDF5 files:\n\
    1. must contain datasets '/spill/run' and '/spill/subrun' if command-line\n\
       option '-k' is not used\n\
    2. may contain multiple groups at root level\n\
    3. each group may contain multiple 2D datasets\n\
    4. 1st dimension of all datasets in the same group share the same size\n\
    5. datasets specified in argument '-k' must exist\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v|-d|-k paths|-o outfile] infile\n%s\n", progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int c, err_exit=0;
    char *infile=NULL, *outfile=NULL, *dset_list=NULL;
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
    while ((c = getopt(argc, argv, "hvdk:o:")) != -1)
        switch(c) {
            case 'k': dset_list=strdup(optarg);
                      break;
            case 'v': verbose = 1;
                      break;
            case 'd': debug = 1;
                      break;
            case 'o': outfile = strdup(optarg);
                      break;
            case 'h': usage(argv[0]);
                      return 0;
            default:  usage(argv[0]);
                      return 1;
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

    if (dset_list == NULL) { /* default datasets */
        it_op.add_level("/spill/run");
        it_op.add_level("/spill/subrun");
    }
    else {
        /* tokenize dest_list */
        char *token = strtok(dset_list, ",");
        while (token != NULL) {
            it_op.add_level(token);
            token = strtok(NULL, ",");
        }
        free(dset_list);
    }

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
        if (file_id < 0) HANDLE_ERROR("H5Fopen");

        if (debug) {
            it_op.reset();

            /* Iterate all objects to collect run and subrun IDs */
#if defined HAS_H5OVISIT3 && HAS_H5OVISIT3
            err = H5Ovisit3(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, get_IDs,
                            &it_op, H5O_INFO_ALL);
            if (err < 0) HANDLE_ERROR("H5Ovisit3");
#else
            err = H5Ovisit(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, get_IDs,
                           &it_op);
            if (err < 0) HANDLE_ERROR("H5Ovisit");
#endif
            if (it_op.err < 0) HANDLE_ERROR("H5Ovisit");

            if (verbose) {
                printf("File %zd:", i);
                for(auto eidx = 0u; eidx < it_op.names.size(); eidx++) {
                  printf(" %s ID = %lld", it_op.names[eidx].c_str(), it_op.event_index[eidx]);
                }
                printf("\n");
            }

            /* use sorted map in an increasing order of run and subrun */
            file_list[it_op.event_index] = in_list[i];
        }
        else {
            hsize_t one[2]={1,1}, offs[2]={0,0}, lens[2]={1,1};
            hid_t dset, mspace, fspace;

            /* reading one element only */
            mspace = H5Screate_simple(2, lens, NULL);
            if (mspace < 0) HANDLE_ERROR("H5Screate_simple");

            for (auto eidx = 0u; eidx < it_op.names.size(); eidx++) {
                int ndims;
                hsize_t dset_dims[2];

                /* open the dataset */
                dset = H5Dopen(file_id, it_op.names[eidx].c_str(), H5P_DEFAULT);
                if (dset < 0) HANDLE_ERROR(string("H5Dopen ")+it_op.names[eidx])

                /* check if this dataset is of zero size */
                fspace = H5Dget_space(dset);
                if (fspace < 0) HANDLE_ERROR("H5Dget_space");
                ndims = H5Sget_simple_extent_dims(fspace, dset_dims, NULL);
                if (ndims < 0)
                    HANDLE_ERROR(string("H5Sget_simple_extent_dims ")+it_op.names[eidx])
                if (dset_dims[0] == 0) {
                    cerr << "Error at line "<<__LINE__<<": file "<<in_list[i]<<
                         " zero-sized dataset '"<<it_op.names[eidx]<<"'"<<endl;
                    err_exit = -1;
                    goto fn_exit;
                }

                /* set the file space for reading the first element only */
                err = H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offs, NULL,
                                          one, lens);
                if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab");

                /* read data and type-cast to long long int */
                err = H5Dread(dset, H5T_NATIVE_LLONG, mspace, fspace,
                              H5P_DEFAULT, &it_op.event_index[eidx]);
                if (err < 0) HANDLE_ERROR("H5Dread");

                err = H5Dclose(dset);
                if (err < 0) HANDLE_ERROR("H5Dclose");

                err = H5Sclose(fspace);
                if (err < 0) HANDLE_ERROR("H5Sclose");
            }
            err = H5Sclose(mspace);
            if (err < 0) HANDLE_ERROR("H5Sclose");

            if (verbose) {
                printf("File %zd:", i);
                for(auto eidx = 0u; eidx < it_op.names.size(); eidx++) {
                    printf(" %s ID = %lld", it_op.names[eidx].c_str(), it_op.event_index[eidx]);
                }
                printf("\n");
            }

            /* check if key tuple (run, subrun, idx) has already existed */
            if (file_list.find(it_op.event_index) != file_list.end()) {
                cerr << "Error: key tuple (";
                for (size_t eidx = 0; eidx < it_op.names.size(); eidx++)
                    cerr << it_op.event_index[eidx] << " " ;
                cerr << ") already exists\n\n";;
                err_exit = -1;
                goto fn_exit;
            }

            /* use sorted map in an increasing order of run, subrun, and idx */
            file_list[it_op.event_index] = in_list[i];
        }
        check_h5_objects(in_list[i].c_str(), file_id);
        err = H5Fclose(file_id);
        if (err < 0) HANDLE_ERROR(string("H5Sclose")+in_list[i]);
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
