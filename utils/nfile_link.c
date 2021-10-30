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
#include <errno.h>

#include <hdf5.h>

static int verbose;

#define HANDLE_ERROR(msg,name) { \
    printf("Error at line %d: func %s on %s\n",__LINE__,msg,name); \
    err_exit = -1; \
    goto fn_exit; \
}

/* parameters to be passed to the call back function */
typedef struct {
    char       in_fname[1024];
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
        /* Skip root group */
        if (!strcmp(name, ".")) return 0;

        if (verbose) printf("Create external link for GROUP %s\n",name);

        err = H5Lcreate_external(it_op->in_fname, name, it_op->out_fd, name,
                                 H5P_DEFAULT, H5P_DEFAULT);
        if (err < 0) HANDLE_ERROR("H5Lcreate_external", name)
    }
    /* ignore other object types */

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
  This utility program merge multiple HDF5 files as links into one master file\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v] -i input_file -o output_file\n%s\n", progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    FILE *fd;
    int i, c, err_exit=0, num_files;
    char *infname=NULL, *outfname=NULL, line[1024], **infile_names;
    herr_t err;
    op_data it_op;

    verbose = 0; /* default is quiet */

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
        exit(1);
    }
    if (outfname == NULL) { /* output file name is mandatory */
        printf("output file is missing\n");
        usage(argv[0]);
        exit(1);
    }
    if (verbose) {
        printf("input  file: %s\n", infname);
        printf("output file: %s\n", outfname);
    }
    memset(&it_op, 0, sizeof(op_data));

    fd = fopen(infname, "r");
    if (fd == NULL) {
        printf("Error: open fails on file %s (%s)\n",infname,strerror(errno));
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
    if (verbose) printf("Number of input files = %d\n",num_files);

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

    /* create output file */
    it_op.out_fd = H5Fcreate(outfname, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if (it_op.out_fd < 0)
        HANDLE_ERROR("H5Fcreate in exclusive mode ", outfname)

    for (i=0; i<num_files; i++) {
        char *fname = infile_names[i];

        if (verbose) printf("open file %s\n",fname);

        /* open input file in read-only mode */
        hid_t in_fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (in_fd < 0) HANDLE_ERROR("Can't open input file", fname)

        strcpy(it_op.in_fname, fname);

        /* Iterate all objects in the input file i */
#if defined HAS_H5OVISIT3 && HAS_H5OVISIT3
        err = H5Ovisit3(in_fd, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, metadata,
                        &it_op, H5O_INFO_ALL);
        if (err < 0) HANDLE_ERROR("H5Ovisit3", fname)
#else
        err = H5Ovisit(in_fd, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, metadata,
                       &it_op);
        if (err < 0) HANDLE_ERROR("H5Ovisit", fname)
#endif
        if (it_op.err < 0) HANDLE_ERROR("H5Ovisit", fname)

        check_h5_objects(fname, in_fd);
        err = H5Fclose(in_fd);
        if (err < 0) HANDLE_ERROR("H5Fclose ", fname)
    }

    check_h5_objects(outfname, it_op.out_fd);
    err = H5Fclose(it_op.out_fd);
    if (err < 0) HANDLE_ERROR("H5Fclose ",outfname)


fn_exit:

    free(infile_names[0]);
    free(infile_names);
    if (infname != NULL) free(infname);
    if (outfname != NULL) free(outfname);

    return (err_exit != 0);
}
