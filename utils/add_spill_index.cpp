/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator
 * Laboratory
 * See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strcmp(), strdup() */
#include <unistd.h> /* getopt() */
#include <assert.h> /* assert() */

#include <string>
#include <vector>

#include <hdf5.h>

using namespace std;

static int verbose;

#define RETURN_ERROR(func_name, obj_name) {				\
	printf("Error in %s line %d: calling %s for object %s\n",__FILE__,__LINE__,func_name,obj_name); \
	return -1;							\
    }

#define HANDLE_ERROR(msg) {				\
    printf("Error at line %d: %s\n",__LINE__, msg);	\
    err_exit = -1;					\
    goto fn_exit;					\
}

#define CALLBACK_ERROR(func_name, dset_name) {			\
    printf("Error in %s line %d: calling %s for group %s dataset %s\n",__FILE__,__LINE__,func_name, grp_name, dset_name); \
    return -1;								\
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
	printf("%s %4zd: type %s, name %s\n",
	       filename, ii, type_name.c_str(), obj_name);
    }
    free(objs);
    return howmany;
}

/*----< allocate_buffer() >--------------------------------------------------*/
/* support adding additional index values of integral type */
int 
allocate_buffer(void *&buffer, const size_t & len, const hid_t & h5type)
{
    htri_t equality = H5Tequal(h5type, H5T_INTEGER);
    
    if (equality > 0) { /* h5type is integral */
	buffer = malloc(len * H5Tget_size(h5type));
	return 0;
    }
    else if (equality == 0) { /* h5type is not integral */
	printf("Unsupported datatype\n");
	return -1;
    }
    else { /* equality check failed */
	RETURN_ERROR("H5Tequal", "h5type");
	/* returns -1 */
    }
    
    /*
    switch (H5Tget_class(h5type)) {
    case H5T_NATIVE_CHAR:
	buffer = (char*          ) malloc(len * sizeof(char          ));
	break;			 				     
    case H5T_NATIVE_SCHAR:	 				     
	buffer = (signed char*   ) malloc(len * sizeof(signed char   ));
	break;			 				     
    case H5T_NATIVE_UCHAR:	 				     
	buffer = (unsigned char* ) malloc(len * sizeof(unsigned char ));
	break;			 				     
    case H5T_NATIVE_SHORT:	 				     
	buffer = (short*         ) malloc(len * sizeof(short         ));
	break;
    case H5T_NATIVE_USHORT:
	buffer = (unsigned short*) malloc(len * sizeof(unsigned short));
	break;
    case H5T_NATIVE_INT:
	buffer = (int*           ) malloc(len * sizeof(int           ));
	break;
    case H5T_NATIVE_UINT:
        buffer = (unsigned int*  ) malloc(len * sizeof(unsigned int  ));
	break;
    case H5T_NATIVE_HSIZE:
	buffer = (hsize_t*       ) malloc(len * sizeof(hsize_t       ));
	break;
    case H5T_NATIVE_HSSIZE:
	buffer = (hssize_t*      ) malloc(len * sizeof(hssize_t      ));
	break;
    case H5T_NATIVE_HERR:
	buffer = (herr_t*        ) malloc(len * sizeof(herr_t        ));
	break;
    default:
	char dtype_name[1024];
	herr_t err;
	err = H5Tenum_nameof(h5type, dtype_name, 1024);
	if(err < 0) RETURN_ERROR("H5Tenum_nameof", "h5type");
	else {
	    printf("Error: %s not supported\n", dtype_name);
	    goto fn_exit;
	}
    }
    */
}



/*----< usage() >------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  [-h]          print this command usage message\n\
  [-v]          verbose mode (default: off)\n\
  -s source     full path to dataset who's value will be injected into the /spill group\n\
                multiple values allowed (required)\n\
  file_name     input/output HDF5 file name (required)\n\n\
  This utility program adds a new dataset in the /spill group of the input file.\n\
  The new dataset should be single-valued, is determined by the source dataset,\n\
  and is intended to be used as an additional index use for partition key generation.\n\
  Requirements for the HDF5 file:\n\
    1. must contain group /spill\n\
    2. contains multiple groups at root level\n\
    3. each group may contain multiple 2D datasets\n\
    4. all datasets in the same group share the 1st dimension size\n\
    5. each group must contain the source datasets\n\
    6. the second dimension size of the source datasets must be 1\n\
    7. data type of the source dataset(s) must be H5T_STD_U32LE\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v] -k base_name file_name\n%s\n", progname, USAGE);
}


int main(int argc, char **argv)
{
    char msg[1024], *infile=NULL, *grp_name;
    hid_t file_id=-1, fapl_id=-1;
    char * source=NULL;
    int c, err_exit;
    herr_t err;
    hsize_t one[2]={1,1}, offs[2]={0,0}, lens[2]={1,1};
    hid_t dset_id, mem_space_id, file_space_id;
    int *inbuf, *outbuf;
    bool dry_run;
    hid_t spill_id, space_id, run_id, dcpl_id, newsrc_id;
    hsize_t dset_dims[2], maxdims[2];
    int ndims;
    
    verbose = 0; /* default is quiet */

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvns:")) != -1) {
	switch(c) {
	    case 'h': usage(argv[0]);
		      err_exit = -1;
		      goto fn_exit;
	    case 'n': dry_run = true;
		      break;
	    case 's': source = strdup(optarg);
		      break;
	    case 'v': verbose = 1;
	              break;
	    default: break;
	}
    }
    if (argv[optind] == NULL) { /* input file name is mandatory */
        printf("input file name is missing\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    if (source == NULL) { /* source dataset is required */
	printf("source dataset is required\n");
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
    if (file_id < 0) {
        sprintf(msg, "Can't open input file %s\n", infile);
        HANDLE_ERROR(msg);
    }

    /* create memspace */
    mem_space_id = H5Screate_simple(2, lens, NULL);
    if (mem_space_id < 0) HANDLE_ERROR("H5Screate_simple");
    
    /* open the source dataset */
    dset_id = H5Dopen(file_id, source, H5P_DEFAULT);
    if (dset_id < 0) RETURN_ERROR("H5Dopen", source);

    /* allocate buffer for data to be read into */
    inbuf = new int[1];

    /* set the window */
    file_space_id = H5Dget_space(dset_id);
    if (file_space_id < 0) RETURN_ERROR("H5Dget_space", source);
	        
    err = H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, offs, NULL,
			      one, lens);
    if (err < 0) RETURN_ERROR("H5Sselect_hyperslab", source);

    /* read data */
    err = H5Dread(dset_id, H5T_NATIVE_INT, mem_space_id, file_space_id,
		  H5P_DEFAULT, inbuf);
    
    if (err < 0) RETURN_ERROR("H5Dread", source);
	
    err = H5Dclose(dset_id);
    if (err < 0) RETURN_ERROR("H5Dclose", source);
    H5Sclose(file_space_id);

    if(verbose) { /* print values we retrieved from the source datasets */
	printf("%s = %d\n", source, *inbuf);
    }
    
    /* open spill group */
    spill_id = H5Gopen(file_id, "/spill", H5P_DEFAULT);
    if (spill_id < 0) RETURN_ERROR("H5Gopen", "/spill");
			
    /* open run dataset in spill group to determine dataset size*/
    run_id = H5Dopen(spill_id, "run", H5P_DEFAULT);
    if (run_id < 0) RETURN_ERROR("H5Dopen", "run");

    /* Get dimension sizes of run dataset */
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

    if (verbose) printf("/spill/run dims (%d, 1)\n", dset_dims[0]);
    
    if (!dry_run) {
	char * newsrc_name = basename(source);

	/* allocate out buffer */
	outbuf = (int*) malloc(dset_dims[0] * sizeof(int));
	if (verbose) printf("Filling /spill/%s (%d, 1) with %d\n", newsrc_name, dset_dims[0], *inbuf);

	/* fill out buffer with value determined from source */
	for(auto ib = 0u; ib < dset_dims[0]; ib++) outbuf[ib] = inbuf[0];

	/* create new dataset */
	newsrc_id = H5Dcreate2(spill_id, newsrc_name, H5T_NATIVE_INT, space_id,
			       H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
	if (newsrc_id < 0) RETURN_ERROR("H5Dcreate2", newsrc_name);
	
	/* write to new dataset */
	if (dset_dims[0] > 0) {
	    err = H5Dwrite(newsrc_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
			   H5P_DEFAULT, outbuf);
	    if (err < 0) RETURN_ERROR("H5Dwrite", newsrc_name);
	}
	
	/* close dataset */
	err = H5Dclose(newsrc_id);
	if (err < 0) RETURN_ERROR("H5Dclose", newsrc_name);

    }

    
    /* close spill group */
    err = H5Gclose(spill_id);
    if (err < 0) RETURN_ERROR("H5Gclose", "/spill");

    /* release buffers */
    if (outbuf != NULL) free(outbuf);
    if (inbuf  != NULL) free(inbuf );
 fn_exit:
    if (file_id >= 0) {
        check_h5_objects(infile, file_id);
        err = H5Fclose(file_id);
        if (err < 0) printf("Error at line %d: H5Fclose\n",__LINE__);
    }
    if (fapl_id >= 0) {
        err = H5Pclose(fapl_id);
        if (err < 0) printf("Error at line %d: H5Pclose\n",__LINE__);
    }

}

