/*
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* strdup() */
#include <unistd.h> /* getopt() */
#include <assert.h> /* assert() */
#include <sys/types.h>
#include <sys/stat.h> /* stat() */
#include <errno.h>

#include <mpi.h>

#include <hdf5.h>
#include <hdf5_hl.h>

static int verbose;

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

#define HANDLE_ERROR(msg,name) { \
    printf("Error at line %d: func %s on %s\n",__LINE__,msg,name); \
    err_exit = -1; \
    goto fn_exit; \
}

/* parameters to be passed to the call back function */
typedef struct {
    hid_t      out_fd;  /* output file ID */
    herr_t     err;
} op_data;

/*----< metadata() >---------------------------------------------------------*/
/* call back function used in H5Ovisit() */
static
herr_t metadata(hid_t             loc_id,        /* object ID */
                const char       *name,          /* object name */
                const H5O_info_t *info,          /* object metadata */
                void             *operator_data) /* data passed from caller */
{
    int err_exit=0;
    herr_t err=0;
    op_data *it_op = (op_data*)operator_data;

    it_op->err = 0;

    if (info->type == H5O_TYPE_GROUP) {
        hid_t grp;

        /* Skip root group */
        if (!strcmp(name, ".")) return 0;

        // if (verbose) printf("GROUP %s\n",name);

        grp = H5Gcreate2(it_op->out_fd, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (grp < 0) HANDLE_ERROR("H5Gcreate2", name)

        err = H5Gclose(grp);
        if (err < 0) HANDLE_ERROR("H5Gclose", name)
    }
    else if (info->type == H5O_TYPE_DATASET) {
        /* create a dataset in output file with the same name */
        hid_t dset, out_dset, space_id, type_id, dcpl_id;

        // if (verbose) printf("Dataset %s\n",name);

        /* Open input dataset */
        dset = H5Dopen(loc_id, name, H5P_DEFAULT);
        if (dset < 0) HANDLE_ERROR("H5Dopen", name)

        /* Retrieve input dataset's data type */
        type_id = H5Dget_type(dset);
        if (type_id < 0) HANDLE_ERROR("H5Dget_type", name)

        /* Open input dataset's space */
        space_id = H5Dget_space(dset);
        if (space_id < 0) HANDLE_ERROR("H5Dget_space", name)

        /* Retrieve input dataset's creation property list */
        dcpl_id = H5Dget_create_plist(dset);
        if (dcpl_id < 0) HANDLE_ERROR("H5Dget_create_plist", name)

        /* Create a dataset in the output file. */
        out_dset = H5Dcreate2(it_op->out_fd, name, type_id, space_id,
                              H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        if (out_dset < 0) HANDLE_ERROR("H5Dcreate2 ", name)

        err = H5Dclose(out_dset);
        if (err < 0) HANDLE_ERROR("H5Dclose", name)

        err = H5Sclose(space_id);
        if (err < 0) HANDLE_ERROR("H5Sclose",name)

        err = H5Pclose(dcpl_id);
        if (err < 0) HANDLE_ERROR("H5Pclose dcpl_id", name)

        err = H5Tclose(type_id);
        if (err < 0) HANDLE_ERROR("H5Tclose", name)

        err = H5Dclose(dset);
        if (err < 0) HANDLE_ERROR("H5Dclose", name)
    }
    else {
        HANDLE_ERROR("unexpected data object type", name)
    }

fn_exit:
    it_op->err = err_exit;
    return (err_exit == 0) ? 0 : 1;
    /* return a positive value causes the visit iterator to
     * immediately return that positive value, indicating
     * short-circuit success.
     */
}

/*----< copy_data() >--------------------------------------------------------*/
/* call back function used in H5Ovisit() */
static
herr_t copy_data(hid_t             loc_id,        /* object ID */
                 const char       *name,          /* object name */
                 const H5O_info_t *info,          /* object metadata */
                 void             *operator_data) /* data passed from caller */
{
    int i, err_exit=0, ntimes;
    herr_t err=0;
    op_data *it_op = (op_data*)operator_data;

    it_op->err = 0;

    if (info->type == H5O_TYPE_GROUP) {
        return 0;
    }
    else if (info->type == H5O_TYPE_DATASET) {
        /* copy the dataset from input file to output file */
        char *buf;
        size_t type_size, buf_len;
        hid_t dset, out_dset, space_id, type_id, dcpl_id, dxpl_id;
        hsize_t dims[2], maxdims[2];
        hsize_t chunk_off[2], in_chunk_dims[2], chunk_nbytes;
        H5D_layout_t in_layout;

        // if (verbose) printf("Dataset %s\n",name);

        /* Open input dataset */
        dset = H5Dopen(loc_id, name, H5P_DEFAULT);
        if (dset < 0) HANDLE_ERROR("H5Dopen", name)

        dxpl_id = H5Pcreate(H5P_DATASET_XFER);
        if (dxpl_id < 0) HANDLE_ERROR("H5Pcreate", name)
        err = H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT);
        if (err < 0) HANDLE_ERROR("H5Pset_dxpl_mpio", name)

        /* Open output dataset */
        out_dset = H5Dopen(it_op->out_fd, name, H5P_DEFAULT);
        if (dset < 0) HANDLE_ERROR("H5Dopen", name)

        /* Retrieve input dataset's data type */
        type_id = H5Dget_type(dset);
        if (type_id < 0) HANDLE_ERROR("H5Dget_type", name)

        /* Retrieve data type size */
        type_size = H5Tget_size(type_id);
        if (type_size == 0) HANDLE_ERROR("H5Tget_size", name)

        /* Open input dataset's space */
        space_id = H5Dget_space(dset);
        if (space_id < 0) HANDLE_ERROR("H5Dget_space", name)

        /* Retrieve number of dimensions */
        int ndims = H5Sget_simple_extent_ndims(space_id);
        if (ndims < 0) HANDLE_ERROR("H5Sget_simple_extent_dims", name)

        /* retrieve dimension sizes */
        ndims = H5Sget_simple_extent_dims(space_id, dims, maxdims);
        if (ndims < 0) HANDLE_ERROR("H5Sget_simple_extent_dims", name)

        buf_len = type_size;
        for (i=0; i<ndims; i++) buf_len *= dims[i];
        buf = (char*) calloc(buf_len, 1);

        /* Retrieve input dataset's creation property list */
        dcpl_id = H5Dget_create_plist(dset);
        if (dcpl_id < 0) HANDLE_ERROR("H5Dget_create_plist", name)

        /* Retrieve input dataset's data layout */
        in_layout = H5Pget_layout(dcpl_id);
        if (in_layout < 0) HANDLE_ERROR("H5Pget_layout", name)

        /* Retrieve input dataset's data chunk dimensions */
        if (in_layout == H5D_CHUNKED) {
            unsigned int filter_mask;

            err = H5Pget_chunk(dcpl_id, 2, in_chunk_dims);
            if (err < 0) HANDLE_ERROR("H5Pget_chunk", name)

            assert(1 == H5Pget_nfilters(dcpl_id));

            /* read and write one chunk at a time */
            ntimes = dims[0] / in_chunk_dims[0];
            if (dims[0] % in_chunk_dims[0]) ntimes++;

            chunk_off[0] = 0;
            chunk_off[1] = 0;
            for (i=0; i<ntimes; i++) {
                err = H5Dget_chunk_storage_size(dset, chunk_off, &chunk_nbytes);
                if (err < 0) HANDLE_ERROR("H5Dget_chunk_storage_size", name)

#if 0
                H5DO APIs are high-level APIs, requiring -lhdf5_hl

                /* read a chunk in compressed form */
                err = H5DOread_chunk(dset, dxpl_id, chunk_off, &filter_mask, buf);
                if (err < 0) HANDLE_ERROR("H5Dread_chunk", name)

                /* write a chunk in compressed form */
                err = H5DOwrite_chunk(out_dset, dxpl_id, filter_mask, chunk_off,
                                     chunk_nbytes, buf);
                if (err < 0) HANDLE_ERROR("H5Dwrite_chunk", name)
#else
                /* read a chunk in compressed form */
                err = H5Dread_chunk(dset, dxpl_id, chunk_off, &filter_mask, buf);
                if (err < 0) HANDLE_ERROR("H5Dread_chunk", name)

                /* write a chunk in compressed form */
                err = H5Dwrite_chunk(out_dset, dxpl_id, filter_mask, chunk_off,
                                     chunk_nbytes, buf);
                if (err < 0) HANDLE_ERROR("H5Dwrite_chunk", name)
#endif

                /* move on to next chunk */
                chunk_off[0] += in_chunk_dims[0];
            }
        }
        else if (in_layout == H5D_CONTIGUOUS) {
            /* read the entire dataset */
            err = H5Dread(dset, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
            if (err < 0) HANDLE_ERROR("H5Dread ", name)

            /* write the entire dataset */
            err = H5Dwrite(out_dset, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
            if (err < 0) HANDLE_ERROR("H5Dwrite", name)
        }
        else assert(0);

        free(buf);

        err = H5Pclose(dxpl_id);
        if (err < 0) HANDLE_ERROR("H5Pclose", name)

        err = H5Pclose(dcpl_id);
        if (err < 0) HANDLE_ERROR("H5Pclose dcpl_id", name)

        err = H5Dclose(out_dset);
        if (err < 0) HANDLE_ERROR("H5Dclose", name)

        err = H5Sclose(space_id);
        if (err < 0) HANDLE_ERROR("H5Sclose",name)

        err = H5Tclose(type_id);
        if (err < 0) HANDLE_ERROR("H5Tclose", name)

        err = H5Dclose(dset);
        if (err < 0) HANDLE_ERROR("H5Dclose", name)
    }
    else {
        HANDLE_ERROR("unexpected data object type", name)
    }

fn_exit:
    it_op->err = err_exit;
    return (err_exit == 0) ? 0 : 1;
    /* return a positive value causes the visit iterator to
     * immediately return that positive value, indicating
     * short-circuit success.
     */
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
    if (howmany > 1) printf("open objects: %zd\n", howmany);

    for (ii=0; ii<howmany; ii++) {
         char *type_name="";
         ot = H5Iget_type(objs[ii]);
              if (ot == H5I_FILE)      continue; /* type_name = "H5I_FILE"; */
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
  [-h]         print this command usage message\n\
  [-v]         verbose mode (default: off)\n\
  -i file      name of input text file containing file names to be merged (required)\n\n\
  -o file      name of output HDF5 file (required)\n\n\
  This utility program merge multiple HDF5 files into one\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v] -i input_file -o output_file\n%s\n", progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    FILE *fd;
    int i, c, err_exit=0, nprocs, rank;
    int num_files, my_num_files, my_start_file;
    char *infname=NULL, *outfname=NULL, line[1024], **infile_names;
    herr_t err;
    hid_t *in_fd, fapl_id;
    op_data it_op;
    MPI_Info info=MPI_INFO_NULL;

    verbose = 0; /* default is quiet */

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvi:o:")) != -1)
        switch(c) {
            case 'v': verbose = 1;
                      break;
            case 'h': usage(argv[0]);
                      return 0;
            case 'i': infname = strdup(optarg);
                      break;
            case 'o': outfname = strdup(optarg);
                      break;
            default : usage(argv[0]);
                      return 1;
        }

    if (infname == NULL) { /* input file name is mandatory */
        printf("input file is missing\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    if (outfname == NULL) { /* output file name is mandatory */
        printf("output file is missing\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    if (verbose && rank == 0) {
        printf("input  file: %s\n", infname);
        printf("output file: %s\n", outfname);
    }
    memset(&it_op, 0, sizeof(op_data));

    fd = fopen(infname, "r");
    if (fd == NULL) {
        if (rank == 0)
            printf("Error: open fails on file %s (%s)\n",infname,strerror(errno));
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(1);
    }

    /* count number of input files */
    num_files = 0;
    while (fgets(line, 1024, fd)) {
        if (strlen(line) == 0)
            continue; /* skip empty lines */
        if (line[0] == '#')
            continue; /* skip comment line (start with #) */
        num_files++;
    }
    if (verbose && rank == 0) printf("Number of input files = %d\n",num_files);

    infile_names = (char**) malloc(num_files * sizeof(char*));
    infile_names[0] = (char*) malloc(num_files * 1024);
    for (i=1; i<num_files; i++)
        infile_names[i] = infile_names[i-1] + 1024;

    /* read input file names */
    rewind(fd);
    i = 0;
    while (fgets(line, 1024, fd)) {
        char *tail;
        if (strlen(line) == 0)
            continue; /* skip empty lines */
        if (line[0] == '#')
            continue; /* skip comment line (start with #) */
        /* remove blanks at tail. Note fgets stores newline to the buffer */
        tail = line + strlen(line) - 1;
        while (*tail == ' ' || *tail == '\t' || *tail == '\n') tail--;
        tail[1] = '\0';
        /* save file name to in_list */
        strcpy(infile_names[i], line);
        i++;
    }
    assert(i == num_files);
    fclose(fd);
    if (num_files < nprocs) {
        if (rank == 0)
            printf("Error: no. MPI processes (%d) must not be less than no. input files (%d)\n",
                   nprocs, num_files);
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(-1);
    }

    in_fd = (hid_t*) malloc(num_files * sizeof(hid_t));

    /* Create the output file using MPI-IO driver */
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id < 0) HANDLE_ERROR("H5Pcreate", outfname)
    err = H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info);
    if (err < 0) HANDLE_ERROR("H5Pset_fapl_mpio", outfname)
    /* set collective mode for metadata writes */
    err = H5Pset_coll_metadata_write(fapl_id, false);
    if (err < 0) HANDLE_ERROR("H5Pset_coll_metadata_write", outfname)

    /* collectively create output file */
    it_op.out_fd = H5Fcreate(outfname, H5F_ACC_EXCL, H5P_DEFAULT, fapl_id);
    if (it_op.out_fd < 0)
        HANDLE_ERROR("H5Fcreate in exclusive mode ", outfname)

    /* all processes open each input file and collectively create groups and
     * datasets in the shared outout file.
     */
    for (i=0; i<num_files; i++) {
        char *fname = infile_names[i];

        if (verbose) printf("rank %2d: open file %s\n",rank,fname);

        /* open input file in read-only mode */
        in_fd[i] = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (in_fd[i] < 0) HANDLE_ERROR("Can't open input file", fname)

        /* Iterate all objects in the input file i */
#if defined HAS_H5OVISIT3 && HAS_H5OVISIT3
        err = H5Ovisit3(in_fd[i], H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, metadata,
                        &it_op, H5O_INFO_ALL);
        if (err < 0) HANDLE_ERROR("H5Ovisit3", fname)
#else
        err = H5Ovisit(in_fd[i], H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, metadata,
                       &it_op);
        if (err < 0) HANDLE_ERROR("H5Ovisit", fname)
#endif
        if (it_op.err < 0) HANDLE_ERROR("H5Ovisit", fname)
    }

    check_h5_objects(outfname, it_op.out_fd);
    err = H5Fclose(it_op.out_fd);
    if (err < 0) HANDLE_ERROR("H5Fclose ",outfname)

    /* collectively open output file */
    it_op.out_fd = H5Fopen(outfname, H5F_ACC_RDWR, fapl_id);
    if (it_op.out_fd < 0)
        HANDLE_ERROR("H5Fopen in RDWR mode ", outfname)

    err = H5Pclose(fapl_id);
    if (err < 0) HANDLE_ERROR("H5Pclose", outfname)

    /* Disable metadata caching. Otherwise H5Dwrite_chunk or H5Fclose will hang
     */
    H5AC_cache_config_t config;
    config.version = H5AC__CURR_CACHE_CONFIG_VERSION;
    err = H5Fget_mdc_config(it_op.out_fd, &config);
    if (err < 0) HANDLE_ERROR("H5Fget_mdc_config", outfname)
    config.evictions_enabled = 0;
    config.incr_mode = H5C_incr__off;
    config.decr_mode = H5C_decr__off;
    config.flash_incr_mode = H5C_flash_incr__off;
    err = H5Fset_mdc_config(it_op.out_fd, &config);
    if (err < 0) HANDLE_ERROR("H5Fset_mdc_config", outfname)

    /* calculate how many files to be assigned to each process for copying */
    my_num_files = num_files / nprocs;
    my_start_file = my_num_files * rank;
    if (rank < num_files % nprocs) {
        my_start_file += rank;
        my_num_files++;
    }
    else
        my_start_file += num_files % nprocs;

    /* each process reads the groups and datasets of assigned input files and
     * copy them to the shared output file.
     */
    for (i=0; i<num_files; i++) {
        char *fname = infile_names[i];

        if (my_start_file <= i && i < my_start_file + my_num_files) {
            if (verbose) printf("rank %2d: copy file %s\n",rank,fname);

            /* Iterate all objects in the input file i */
#if defined HAS_H5OVISIT3 && HAS_H5OVISIT3
            err = H5Ovisit3(in_fd[i], H5_INDEX_CRT_ORDER, H5_ITER_NATIVE,
                            copy_data, &it_op, H5O_INFO_ALL);
            if (err < 0) HANDLE_ERROR("H5Ovisit3", fname)
#else
            err = H5Ovisit(in_fd[i], H5_INDEX_CRT_ORDER, H5_ITER_NATIVE,
                           copy_data, &it_op);
            if (err < 0) HANDLE_ERROR("H5Ovisit", fname)
#endif
            if (it_op.err < 0) HANDLE_ERROR("H5Ovisit", fname)
        }

        check_h5_objects(fname, in_fd[i]);
        err = H5Fclose(in_fd[i]);
        if (err < 0) HANDLE_ERROR("H5Fclose ", fname)
    }

    check_h5_objects(outfname, it_op.out_fd);
    err = H5Fclose(it_op.out_fd);
    if (err < 0) HANDLE_ERROR("H5Fclose ",outfname)

fn_exit:
    free(in_fd);
    free(infile_names[0]);
    free(infile_names);
    if (infname != NULL) free(infname);
    if (outfname != NULL) free(outfname);

    MPI_Finalize();
    return (err_exit != 0);
}
