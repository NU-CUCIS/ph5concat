/*
 * Copyright (C) 2021, Northwestern University
 * See COPYRIGHT notice in top-level directory.
 */

/* This program is to be run in sequential.
 * It merges multiple HDF5 files into one, by calling H5Ocopy().
 * Because HDF5 1.12.x requires H5Ocopy() to be called collectively, this
 * program can only be run on one process. Note the object names must be unique
 * among all the input files.
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
    hid_t      out_fd;  /* output file ID */
    herr_t     err;
} op_data;

/*----< copy_obj() >---------------------------------------------------------*/
/* call back function used in H5Giterate() */
static
herr_t copy_obj(hid_t       group,         /* group ID */
                const char *name,          /* object name */
                void       *operator_data) /* data passed from caller */
{
    int err_exit=0;
    herr_t err=0;
    op_data *it_op = (op_data*)operator_data;

    it_op->err = 0;

    if (verbose) printf("Copying %s\n",name);

    err = H5Ocopy(group, name, it_op->out_fd, name, H5P_DEFAULT, H5P_DEFAULT);
    if (err < 0) HANDLE_ERROR("H5Ocopy", name)

fn_exit:
    it_op->err = err_exit;
    return (err_exit == 0) ? 0 : 1;
    /* return a positive value causes the visit iterator to
     * immediately return that positive value, indicating
     * short-circuit success.
     */
}

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  [-h]     print this command usage message\n\
  [-v]     verbose mode (default: off)\n\
  -i file  name of input text file containing file names to be merged (required)\n\
  -o file  name of output HDF5 file (required)\n\n\
  This utility program merge multiple HDF5 files into one\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v] -i input_file -o output_file\n%s\n", progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    FILE *fd;
    int i, c, err_exit=0, num_files;
    char *infname=NULL, *outfname=NULL, line[1024], **infile_names=NULL;
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
        err_exit = -1;
        goto fn_exit;
    }
    if (outfname == NULL) { /* output file name is mandatory */
        printf("output file is missing\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    printf("input  file: %s\n", infname);
    printf("output file: %s\n", outfname);
    memset(&it_op, 0, sizeof(op_data));

    fd = fopen(infname, "r");
    if (fd == NULL) {
        printf("Error: open fails on file %s (%s)\n",infname,strerror(errno));
        exit(1);
    }

    /* count number of input files */
    num_files = 0;
    while (fgets(line, 1024, fd)) {
        if (line[strlen(line)-1] == '\n')
            line[strlen(line)-1] = '\0';
        if (strlen(line) == 0)
            continue; /* skip empty lines */
        if (line[0] == '#')
            continue; /* skip comment line (start with #) */
        if (line[0] == '\n')
            continue; /* skip new line */
        num_files++;
    }
    printf("Number of input files = %d\n",num_files);

    infile_names = (char**) malloc(num_files * sizeof(char*));
    infile_names[0] = (char*) malloc(num_files * 1024);
    for (i=1; i<num_files; i++)
        infile_names[i] = infile_names[i-1] + 1024;

    /* read input file names */
    rewind(fd);
    i = 0;
    while (fgets(line, 1024, fd)) {
        char *tail;
        if (line[strlen(line)-1] == '\n')
            line[strlen(line)-1] = '\0';
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

    /* collectively create output file */
    it_op.out_fd = H5Fcreate(outfname, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
    if (it_op.out_fd < 0)
        HANDLE_ERROR("H5Fcreate in exclusive mode ", outfname)

    for (i=0; i<num_files; i++) {
        if (verbose) printf("copying file %s\n",infile_names[i]);

        /* open input file in read-only mode */
        hid_t in_fd = H5Fopen(infile_names[i], H5F_ACC_RDONLY, H5P_DEFAULT);
        if (in_fd < 0) HANDLE_ERROR("Can't open input file", infile_names[i])

        /* Iterate all objects in root group of input file i */
        err = H5Giterate(in_fd, ".", NULL, copy_obj, &it_op);
        if (err < 0) HANDLE_ERROR("H5Giterate", infile_names[i])

        if (it_op.err < 0) HANDLE_ERROR("H5Ovisit", infile_names[i])

        err = H5Fclose(in_fd);
        if (err < 0) HANDLE_ERROR("H5Fclose ", infile_names[i])
    }

    err = H5Fclose(it_op.out_fd);
    if (err < 0) HANDLE_ERROR("H5Fclose ",outfname)

fn_exit:
    if (infile_names != NULL) {
        free(infile_names[0]);
        free(infile_names);
    }
    if (infname != NULL) free(infname);
    if (outfname != NULL) free(outfname);

    return (err_exit != 0);
}
