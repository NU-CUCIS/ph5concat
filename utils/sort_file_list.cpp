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

const std::vector<std::string> DEFAULT_LEVELS {"/spill/run", "/spill/subrun"};

struct compv {
  bool operator()(const std::vector<unsigned int>& lhs, const std::vector<unsigned int>& rhs)
  {
    assert(lhs.size() == rhs.size());
  
    for(auto i = 0u; i < lhs.size(); i++) {
      if(lhs[i] < rhs[i]) return true;
      if(lhs[i] > rhs[i]) return false;
    }
    return false;
  }
};

typedef map<std::vector<unsigned int>, string, compv> tuple_map;

/* parameters to be passed to the call back function */
struct op_data {
  std::vector<std::string> names;
  std::vector<bool> is_set;
  std::vector<unsigned int> event_index;
  herr_t   err;
  op_data(std::vector<std::string> _names) 
    : names(_names) 
  {    
    is_set      = std::vector<bool        >(names.size(), false);
    event_index = std::vector<unsigned int>(names.size(), 0    );
  }
  void add_level(std::string name)
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
  int index(std::string name) const
  {
    for(auto idx = 0; idx < names.size(); idx++) 
      if(strcmp(name.c_str(), names[idx].c_str()) == 0) return idx;
    return -1;
  }
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
    int level_idx;

    it_op->err = 0;

    if (info->type != H5O_TYPE_DATASET) return 0;

    /* remove group name */
    char *path = strdup(name);
    char *dset_name = basename(path);

    /* placeholder for level names */
    char *level_name;
    
    /* skip dataset that is not in list of index names
       Note: dset_name is compared to basename of the level index
             full path to ensure consistency within file.
	     Importantly, the current group does not contain
	     a matching dataset, this group will be skipped.
     */
    bool skip = true;
    for(auto eidx = 0u; eidx < it_op->names.size(); eidx++)
      /* basename requires non-const char* so make a copy */
      strncpy(level_name, it_op->names[eidx].c_str(), it_op->names[eidx].size());
      skip &= strcmp(dset_name, level_name) != 0;      
    if(skip) goto fn_exit;

    /* Open the dataset. Note that loc_id is not the dataset/group ID. */
    dset_id = H5Dopen(loc_id, name, H5P_DEFAULT);
    if (dset_id < 0) CALLBACK_ERROR("H5Dopen",name);

    /* Get dimension sizes */
    space_id = H5Dget_space(dset_id);
    if (space_id < 0) CALLBACK_ERROR("H5Dget_space",name);
    ndims = H5Sget_simple_extent_dims(space_id, dset_dims, NULL);
    if (ndims < 0) CALLBACK_ERROR("H5Sget_simple_extent_dims",name);
    err = H5Sclose(space_id);
    if (err < 0) CALLBACK_ERROR("H5Sclose",name);

    /* skip zero-size dataset */
    if (dset_dims[0] == 0) {
        err = H5Dclose(dset_id);
        if (err < 0) CALLBACK_ERROR("H5Dclose",name);
        goto fn_exit;
    }

    /* allocate read buffer */
    buf = new unsigned int [dset_dims[0]];

    /* read the entire dataset */
    err = H5Dread(dset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    if (err < 0) CALLBACK_ERROR("H5Dread",name);
    err = H5Dclose(dset_id);
    if (err < 0) CALLBACK_ERROR("H5Dclose",name);

    level_idx = it_op->index(dset_name);

    /* get the first value of level */
    if(! it_op->is_set[level_idx]) {
      it_op->is_set     [level_idx] = true;
      it_op->event_index[level_idx] = buf[0];
    }

    /* check level in all groups for consistency */
    for(i=0; i<dset_dims[0]; i++) {
      if(buf[i] != it_op->event_index[level_idx]) {
	printf("Error: inconsistent %s ID %u, expecting %u\n",
	       dset_name, it_op->event_index[level_idx]);
	err_exit = -1;
	goto fn_exit;
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
    struct op_data it_op(DEFAULT_LEVELS);
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
		printf(" %s ID = %d", it_op.names[eidx].c_str(), it_op.event_index[eidx]);
	      }
	      printf("\n");
	    }

            /* use sorted map in an increasing order of run and subrun */

            file_list[it_op.event_index] = in_list[i];
        }
        else {
            unsigned int run, subrun;
            hsize_t one[2]={1,1}, offs[2]={0,0}, lens[2]={1,1};
            hid_t dset_id, mem_space_id, file_space_id;
	    
            mem_space_id = H5Screate_simple(2, lens, NULL);
            if (mem_space_id < 0) HANDLE_ERROR("H5Screate_simple");

	    for(auto eidx = 0u; eidx < it_op.names.size(); eidx++) {
	      /* open the dataset */
	      dset_id = H5Dopen(file_id, it_op.names[eidx].c_str(), H5P_DEFAULT);
	      if (dset_id < 0) HANDLE_ERROR("H5Dopen");
	      
	      
	      file_space_id = H5Dget_space(dset_id);
	      if (file_space_id < 0) HANDLE_ERROR("H5Dget_space");

	      /* set the window */
	      err = H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offs, NULL,
					one, lens);
	      if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab");

	      /* read data */
	      err = H5Dread(dset_id, H5T_NATIVE_UINT, mem_space_id, file_space_id,
			    H5P_DEFAULT, &it_op.event_index[eidx]);
	      if (err < 0) HANDLE_ERROR("H5Dread");

	      err = H5Dclose(dset_id);
	      if (err < 0) HANDLE_ERROR("H5Dclose");
	      H5Sclose(file_space_id);
	    }
            H5Sclose(mem_space_id);

            if (verbose) {
	      printf("File %zd:", i);
	      for(auto eidx = 0u; eidx < it_op.names.size(); eidx++) {
		printf(" %s ID = %d", it_op.names[eidx].c_str(), it_op.event_index[eidx]);
	      }
	      printf("\n");
	    }

            /* use sorted map in an increasing order of run and subrun */
            file_list[it_op.event_index] = in_list[i];
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
