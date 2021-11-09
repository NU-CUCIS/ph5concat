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
#include <libgen.h> /* dirname() */
#include <errno.h>

#include <hdf5.h>

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

typedef struct {
    char      name[128];      /* name of group */
    int       num_dsets;      /* number of datasets in this group */
    hsize_t   sum_row_size;   /* size sum of a row of all datasets */
    hsize_t   sum_dset_size;  /* size sum of all datasets */
    hsize_t   num_edges;      /* number of edges in this group */
    hsize_t   num_hits;       /* number of hits in this group */
} group_meta;

/* parameters to be passed to the call back function */
typedef struct {
    int          gid;         /* array index point to groups[] */
    hsize_t      num_groups;  /* number of groups */
    group_meta  *groups;      /* metadata of groups */
    int          num_links;   /* number of unique external links */
    char       **ext_links;   /* path names of external links */
    char        *evt_group;   /* event ID dataset's group name */
    char        *evt_dset;    /* event ID dataset name */
    hsize_t      num_events;  /* number of event IDs */
    hsize_t     *evt_size;    /* data size per event */
    int          hasEdges;    /* whether input file contains edges */
    int          hasHits;     /* whether input file contains hints */
    herr_t       err;
} op_data;

/*----< dataset_metadata() >-------------------------------------------------*/
static
herr_t dataset_metadata(hid_t       loc_id,
                        const char *name,
                        op_data    *it_op)
{
    int i, err_exit=0;
    char *type_name, *layout_name;
    size_t type_size;
    herr_t err=0;
    hid_t dset, space_id, type_id, r_dcpl;
    hsize_t nelems, dims[2], maxdims[2], chunk_dims[2];
    H5D_layout_t layout;
    H5T_class_t type_class;

    /* retrieve group name and dataset name */
    // char *name_copy = strdup(name);
    // char *group_name = strtok(name_copy, "/");
    const char *dataset_name = strrchr(name, '/');
    if (dataset_name == NULL) dataset_name = name;
    else dataset_name += 1;

    if (verbose) printf("check dataset %s\n",name);

    /* Open input dataset */
    dset = H5Dopen(loc_id, name, H5P_DEFAULT);
    if (dset < 0) HANDLE_ERROR("H5Dopen", name)

    /* Retrieve input dataset's data type */
    type_id = H5Dget_type(dset);
    if (type_id < 0) HANDLE_ERROR("H5Dget_type", name)

    /* Retrieve data type size */
    type_size = H5Tget_size(type_id);
    if (type_size == 0) HANDLE_ERROR("H5Tget_size", name)

    /* Retrieve data type class */
    type_class = H5Tget_class(type_id);
    if (type_class < 0) HANDLE_ERROR("H5Tget_class", name);

    if (type_class == H5T_STRING)
        type_name = "H5T_STRING";
    else if (H5Tequal(type_id, H5T_IEEE_F32LE))
        type_name = "H5T_IEEE_F32LE";
    else if (H5Tequal(type_id, H5T_STD_I32LE))
        type_name = "H5T_STD_I32LE";
    else if (H5Tequal(type_id, H5T_STD_I64LE))
        type_name = "H5T_STD_I64LE";
    else
        type_name = "unknown";

    /* Open input dataset's space */
    space_id = H5Dget_space(dset);
    if (space_id < 0) HANDLE_ERROR("H5Dget_space", name)

    /* Retrieve number of dimensions */
    int ndims = H5Sget_simple_extent_ndims(space_id);
    if (ndims < 0) HANDLE_ERROR("H5Sget_simple_extent_dims", name)

    /* retrieve dimension sizes */
    ndims = H5Sget_simple_extent_dims(space_id, dims, maxdims);
    if (ndims < 0) HANDLE_ERROR("H5Sget_simple_extent_dims", name)

    err = H5Sclose(space_id);
    if (err < 0) HANDLE_ERROR("H5Sclose",name)

    /* skip datasets event_id, event_id.seq, event_id.seq_cnt */
    if (it_op->evt_dset == NULL ||
        strncmp(it_op->evt_dset, dataset_name, strlen(it_op->evt_dset)))
        it_op->groups[it_op->gid].sum_row_size += dims[1] * type_size;

    nelems = 1;
    for (i=0; i<ndims; i++) nelems *= dims[i];

    /* accumulate all dataset sizes before compression */
    it_op->groups[it_op->gid].sum_dset_size += nelems * type_size;

    /* collect edge dimensions from dataset */
    if (strncmp(dataset_name,"edge_index", strlen("edge_index")) == 0) {
        /* string compare only the prefix name */
        if (ndims == 1)
            it_op->groups[it_op->gid].num_edges += dims[0];
        else if (ndims == 2)
            it_op->groups[it_op->gid].num_edges += dims[1];
        it_op->hasEdges = 1;
    }

    /* collect number of hits from each group */
    if (dataset_name[0] == 'y') {
        /* string compare only the prefix name */
        it_op->groups[it_op->gid].num_hits += dims[0];
        it_op->hasHits = 1;
    }

    /* Retrieve input dataset's creation property list */
    r_dcpl = H5Dget_create_plist(dset);
    if (r_dcpl < 0) HANDLE_ERROR("H5Dget_create_plist", name)

    /* Retrieve input dataset's data layout */
    layout = H5Pget_layout(r_dcpl);
    if (layout < 0) HANDLE_ERROR("H5Pget_layout", name)

    /* Retrieve input dataset's data chunk dimensions */
    if (layout == H5D_CHUNKED) {
        err = H5Pget_chunk(r_dcpl, ndims, chunk_dims);
        if (err < 0) HANDLE_ERROR("H5Pget_chunk", name)
        layout_name = "H5D_CHUNKED";
    }
    else if (layout == H5D_CONTIGUOUS)
        layout_name = "H5D_CONTIGUOUS";
    else if (layout == H5D_COMPACT)
        layout_name = "H5D_COMPACT";

    /* close input dataset creation property */
    err = H5Pclose(r_dcpl);
    if (err < 0) HANDLE_ERROR("H5Pclose r_dcpl", name)

    if (verbose) {
        printf("\tDATASET \"%s\"\n",dataset_name);
        printf("\t\tDATATYPE %s\n",type_name);
        if (ndims == 1) {
            printf("\t\tDATASPACE ( %lld ) / ", dims[0]);
            printf("\t\tDATASPACE ( %lld ) / ", dims[0]);
            if (maxdims[0] == H5S_UNLIMITED)
                printf("( H5S_UNLIMITED )\n");
            else
                printf("( %lld )\n", dims[0]);
            if (layout == H5D_CHUNKED)
                printf("\t\tSTORAGE_LAYOUT ( %lld )\n", chunk_dims[0]);
            else
                printf("\t\tSTORAGE_LAYOUT %s\n",layout_name);
        }
        else if (ndims == 2) {
            printf("\t\tDATASPACE ( %lld, %lld ) / ", dims[0], dims[1]);
            if (maxdims[0] == H5S_UNLIMITED)
                printf("( H5S_UNLIMITED,");
            else
                printf("( %lld", dims[0]);
            if (maxdims[1] == H5S_UNLIMITED)
                printf(" H5S_UNLIMITED)\n");
            else
                printf(" %lld)\n", dims[1]);
            if (layout == H5D_CHUNKED)
                printf("\t\tSTORAGE_LAYOUT ( %lld, %lld)\n", chunk_dims[0],
                       chunk_dims[1]);
            else
                printf("\t\tSTORAGE_LAYOUT %s\n",layout_name);
        }
    }

    err = H5Dclose(dset);
    if (err < 0) HANDLE_ERROR("H5Dclose", name)

fn_exit:
    return (err_exit == 0) ? 0 : -1;
}

/*----< ext_metadata() >-----------------------------------------------------*/
/* call back function used in H5Giterate() */
static
herr_t ext_metadata(hid_t       loc_id,        /* object ID */
                    const char *name,          /* object name */
                    void       *operator_data) /* data passed from caller */
{
    int err_exit=0;
    herr_t err=0;
    H5G_stat_t statbuf;

    op_data *it_op = (op_data*)operator_data;
    it_op->err = 0;

    if (verbose) printf("%s name = %s\n",__func__,name);

    /* retrieve info by following the external link */
    err = H5Gget_objinfo(loc_id, name, 1, &statbuf);
    if (err < 0) HANDLE_ERROR("H5Gget_objinfo", name)

    if (statbuf.type == H5G_DATASET) {
        // printf(" Object with name %s is a dataset \n", name);
        err = dataset_metadata(loc_id, name, it_op);
        if (err < 0) HANDLE_ERROR("dataset_metadata", name)
    }
    else if (statbuf.type == H5G_GROUP) {
        if (verbose) printf(" Object with name %s is a group \n", name);
    }
    else if (verbose)
        printf(" Unable to identify an object ");

fn_exit:
    it_op->err = err_exit;
    return (err_exit == 0) ? 0 : 1;
    /* return a positive value causes the visit iterator to
     * immediately return that positive value, indicating
     * short-circuit success.
     */
}

static
int search_link_path(char       **str_array,
                     int          num,
                     const char  *str)
{
    int i;
    for (i=0; i<num; i++)
        if (strcmp(str_array[i], str) == 0)
            return 1; /* found */

    return 0; /* not found */
}

/*----< link_metadata() >----------------------------------------------------*/
/* call back function used in H5Ovisit() */
static
herr_t link_metadata(hid_t             loc_id,        /* object ID */
                     const char       *name,          /* name of link */
                     const H5L_info_t *info,          /* object metadata */
                     void             *operator_data) /* data passed from caller */
{
    int err_exit=0;
    herr_t err=0;
    hid_t grp;
    H5G_info_t grp_info;

    op_data *it_op = (op_data*)operator_data;
    it_op->err = 0;

    if (info->type == H5L_TYPE_EXTERNAL) {
        /* each external link points to a group in an external file */
        grp = H5Oopen(loc_id, name, H5P_DEFAULT);
        if (grp < 0) HANDLE_ERROR("H5Oopen", name)

        unsigned int flags;
        char *linkbuf;
        const char *link_path;
        linkbuf = (char*) calloc(info->u.val_size+1, 1);
        err = H5Lget_val(loc_id, name, linkbuf, info->u.val_size, H5P_DEFAULT);
        if (err < 0) HANDLE_ERROR("H5Lget_val", name);
        err = H5Lunpack_elink_val(linkbuf, info->u.val_size, &flags, &link_path, NULL);
        if (err < 0) HANDLE_ERROR("H5Lunpack_elink_val", name);
        if (verbose) printf("group %s: external link: %s\n",name,link_path);

        /* search if link_path already exists */
        if (it_op->num_links == 0) {
            it_op->ext_links[0] = strdup(link_path);
            it_op->num_links++;
        }
        else if (search_link_path(it_op->ext_links, it_op->num_links, link_path) == 0)
            it_op->ext_links[it_op->num_links++] = strdup(link_path);
        free(linkbuf);

        /* inquire number of datasets in the external group */
        err = H5Gget_info(grp, &grp_info);
        if (err < 0) HANDLE_ERROR("H5Gget_info", name);

        if (verbose)
            printf("H5L_TYPE_EXTERNAL %s has %lld datasets\n",name,
                   grp_info.nlinks);

        it_op->gid++;
        strcpy(it_op->groups[it_op->gid].name, name);
        it_op->groups[it_op->gid].num_dsets = grp_info.nlinks;

        /* iterate all datasets in the group to collect their metadata */
        err = H5Giterate(loc_id, name, NULL, ext_metadata, it_op);
        if (err < 0) HANDLE_ERROR("H5Giterate", name)

        err = H5Oclose(grp);
        if (err < 0) HANDLE_ERROR("H5Oclose", name)
    }
    else if (info->type == H5L_TYPE_HARD) {
        /* not an external link */
        H5O_info_t obj_info;

        grp = H5Oopen(loc_id, name, H5P_DEFAULT);
        if (grp < 0) HANDLE_ERROR("H5Oopen", name)

        err = H5Oget_info(grp, &obj_info, H5O_INFO_BASIC);
        if (err < 0) HANDLE_ERROR("H5Oget_info", name)

        if (obj_info.type == H5O_TYPE_GROUP) {
            err = H5Gget_info(grp, &grp_info);
            if (err < 0) HANDLE_ERROR("H5Gget_info", name);

            if (verbose)
                printf("GROUP %s (contains %llu datasets)\n",
                       name, grp_info.nlinks);

            it_op->gid++;
            strcpy(it_op->groups[it_op->gid].name, name);
            it_op->groups[it_op->gid].num_dsets = grp_info.nlinks;
        }
        else if (obj_info.type == H5O_TYPE_DATASET) {
            if (verbose) printf("H5O_TYPE_DATASET %s\n",name);

            err = dataset_metadata(loc_id, name, it_op);
            if (err < 0) HANDLE_ERROR("dataset_metadata", name)
        }
        err = H5Oclose(grp);
        if (err < 0) HANDLE_ERROR("H5Oclose", name)
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
    if (howmany > 1) printf("open objects:\n");

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

/*----< cmpfunc() >----------------------------------------------------------*/
static
int cmpfunc (const void * a, const void * b) {

    if (*((hsize_t *)a) > *((hsize_t *)b)) return (1);
    if (*((hsize_t *)a) < *((hsize_t *)b)) return (-1);
    return (0);
}

/*----< histogram() >--------------------------------------------------------*/
static
void histogram(hsize_t  num,       /* number of elements in array values */
               hsize_t *values,    /* already sorted in an increasing order */
               int      bin_width,
               char    *out_fname)
{
    FILE *fp;
    hsize_t i;

    /* calculate number of bins */
    hsize_t num_bins = values[num-1] / bin_width;
    if (values[num-1] % bin_width) num_bins++;

    /* Initialize frequency array */
    hsize_t* frequency = (hsize_t*) calloc(num_bins, sizeof(hsize_t));

    /* Create histogram */
    for (i=0; i<num; i++) {
        int bin = values[i] / bin_width;
        frequency[bin]++;
    }

    /* Write histogram array to CSV file */
    fp = fopen(out_fname, "w");
    for (i = 0; i < num_bins; i++)
        fprintf(fp, "%llu, %llu\n", i, frequency[i]);

    fclose(fp);

    free(frequency);
}

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  [-h]         print this command usage message\n\
  [-v]         verbose mode (default: off)\n\
  [-e name]    event ID dataset name\n\
  [-b num]     bin width for creating event size histogram (default: 100000)\n\
  [-t num]     bin width for creating hit count histogram (default: 100)\n\
  [-d num]     bin width for creating edge count histogram (default: 100)\n\
  file         input file name (required)\n\n\
  This utility program prints the statistics of an input HDF5 file\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v|-e name|-b num|-t num|-d num] file\n%s\n", progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int i, j, k, c, err_exit=0;
    int event_bin_width=100000, hit_bin_width=100, edge_bin_width=100;
    char *fname=NULL, *evt_dset=NULL;
    herr_t err;
    hid_t fd=-1;
    hsize_t non_empty_evts;
    float max_evt_size, min_evt_size, median_evt_size;
    H5G_info_t grp_info;
    op_data it_op;
    struct stat file_stat;

    verbose = 0; /* default is quiet */

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hve:b:t:d:")) != -1)
        switch(c) {
            case 'v': verbose = 1;
                      break;
            case 'h': usage(argv[0]);
                      return 0;
            case 'e': evt_dset = strdup(optarg);
                      break;
            case 'b': event_bin_width = atoi(optarg);
                      break;
            case 't': hit_bin_width = atoi(optarg);
                      break;
            case 'd': edge_bin_width = atoi(optarg);
                      break;
            default : usage(argv[0]);
                      return 1;
        }

    if (argv[optind] == NULL) { /* input file name is mandatory */
        printf("input file is missing\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    fname = strdup(argv[optind]);
    if (verbose) printf("input file: %s\n", fname);
    memset(&it_op, 0, sizeof(op_data));

    /* open input file in read-only mode */
    fd = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fd < 0) HANDLE_ERROR("Can't open input file", fname)

    err = H5Gget_info_by_name(fd, "/", &grp_info, H5P_DEFAULT);
    if (err < 0) HANDLE_ERROR("H5Gget_info_by_name - root group",  fname)

    it_op.num_events = 0;
    it_op.num_groups = grp_info.nlinks;
    if (verbose) printf("Number of groups = %llu\n", it_op.num_groups);

    it_op.num_links = 0;
    it_op.ext_links = (char**) calloc(it_op.num_groups, sizeof(char*));

    it_op.groups = (group_meta*) calloc(it_op.num_groups, sizeof(group_meta));
    it_op.gid = -1;

    /* collect metadata of event ID dataset */
    it_op.evt_size    = NULL;
    it_op.evt_group   = NULL;
    it_op.evt_dset    = NULL;
    if (evt_dset != NULL) {
        char *ptr;
        /* extract group name */
        if (evt_dset[0] == '/')
            it_op.evt_group = strdup(evt_dset + 1);
        else
            it_op.evt_group = strdup(evt_dset);
        ptr = strchr(it_op.evt_group, '/');
        *ptr = '\0';

        /* extract dataset name */
        it_op.evt_dset = ptr + 1;

        /* check if the group of event ID dataset exists */
        htri_t grp_exist = H5Lexists(fd, it_op.evt_group, H5P_DEFAULT);
        if (grp_exist < 0) HANDLE_ERROR("H5Lexists", it_op.evt_group)
        if (grp_exist == 0) {
            printf("Warning: group '%s' does not exist.\n", it_op.evt_group);
            printf("Skip collecting metadata of event IDs.\n");
            free(evt_dset);
            evt_dset = NULL;
        }
        else { /* check if the event ID dataset exists */
            htri_t dset_exist = H5Lexists(fd, evt_dset, H5P_DEFAULT);
            if (dset_exist < 0) HANDLE_ERROR("H5Lexists", evt_dset)
            if (dset_exist == 0) {
                printf("Warning: dataset '%s' does not exist.\n", evt_dset);
                printf("Skip collecting metadata of event IDs.\n");
                free(evt_dset);
                evt_dset = NULL;
            }
        }
    }
    if (evt_dset != NULL) {
        hsize_t dims[2];

        /* open event ID dataset */
        hid_t dset = H5Dopen(fd, evt_dset, H5P_DEFAULT);
        if (dset < 0) HANDLE_ERROR("H5Dopen", evt_dset)

        /* retrieve dimension sizes */
        hid_t space_id = H5Dget_space(dset);
        if (space_id < 0) HANDLE_ERROR("H5Dget_space", evt_dset)
        int ndims = H5Sget_simple_extent_dims(space_id, dims, NULL);
        if (ndims < 0) HANDLE_ERROR("H5Sget_simple_extent_dims", evt_dset)

        /* number of event IDs is the size of dimension 0 */
        it_op.num_events = dims[0];
        if (verbose)
            printf("Number of partitioning keys = %llu\n", it_op.num_events);

        err = H5Sclose(space_id);
        if (err < 0) HANDLE_ERROR("H5Sclose",evt_dset)
        err = H5Dclose(dset);
        if (err < 0) HANDLE_ERROR("H5Dclose", evt_dset)

        /* allocate evt_size, amount of data per event across all groups */
        it_op.evt_size = (hsize_t*) calloc(it_op.num_events, sizeof(hsize_t));;
    }

    /* Iterate all objects and perform chunking adjustment */
    err = H5Lvisit(fd, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, link_metadata,
                   &it_op);
    if (err < 0) HANDLE_ERROR("H5Lvisit", fname)

    /* calculate event statistics */
    if (evt_dset != NULL) {
        for (k=0; k<it_op.num_groups; k++) {
            hsize_t dims[2];
            char dset_name[1024];

            /* skip event ID dataset's group */
            if (!strcmp(it_op.groups[k].name, it_op.evt_group))
                continue;

            /* construct dataset name for seq_cnt */
            sprintf(dset_name, "/%s/%s.seq_cnt", it_op.groups[k].name,
                    it_op.evt_dset);

            /* check if the dataset exists */
            htri_t dset_exist = H5Lexists(fd, dset_name, H5P_DEFAULT);
            if (dset_exist < 0) HANDLE_ERROR("H5Lexists", dset_name)
            if (dset_exist == 0) {
                printf("Error: dataset '%s' does not exist.\n", dset_name);
                HANDLE_ERROR("H5Lexists", dset_name)
            }

            /* open dataset */
            hid_t dset = H5Dopen(fd, dset_name, H5P_DEFAULT);
            if (dset < 0) HANDLE_ERROR("H5Dopen", dset_name)

            /* Open input dataset's space */
            hid_t space_id = H5Dget_space(dset);
            if (space_id < 0) HANDLE_ERROR("H5Dget_space", dset_name)

            /* retrieve dimension sizes */
            int ndims = H5Sget_simple_extent_dims(space_id, dims, NULL);
            if (ndims < 0) HANDLE_ERROR("H5Sget_simple_extent_dims", dset_name)

            err = H5Sclose(space_id);
            if (err < 0) HANDLE_ERROR("H5Sclose", dset_name)

            /* read the entire dataset .seq_cnt */
            assert(dims[1] == 2);
            int64_t *seq_cnt = (int64_t*) malloc(dims[0] * 2 * sizeof(int64_t));

            /* seq_cnt is a 2D array [number of sequences][count] */
            err = H5Dread(dset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, seq_cnt);
            if (err < 0) HANDLE_ERROR("H5Dread ", dset_name)

            /* accumulate the event data size */
            for (i=0,j=0; i<dims[0]; i++) {
                it_op.evt_size[seq_cnt[j]] += seq_cnt[j+1]
                                           * it_op.groups[k].sum_row_size;
                j += 2;
            }
            free(seq_cnt);

            err = H5Dclose(dset);
            if (err < 0) HANDLE_ERROR("H5Dclose", dset_name)
        }

        max_evt_size = 0;
        min_evt_size = LONG_MAX;
        non_empty_evts = 0;
        for (i=0; i<it_op.num_events; i++) {
            if (it_op.evt_size[i] == 0) continue;
            non_empty_evts++;
            max_evt_size = MAX(max_evt_size, it_op.evt_size[i]);
            min_evt_size = MIN(min_evt_size, it_op.evt_size[i]);
        }
        /* Sort event sizes array to create histogram */
        qsort(it_op.evt_size, it_op.num_events, sizeof(hsize_t), cmpfunc);
        median_evt_size = it_op.evt_size[(it_op.num_events+1)/2];
    }

fn_exit:
    if (fd >= 0) {
        check_h5_objects(fname, fd);
        err = H5Fclose(fd);
        if (err < 0) printf("Error at line %d: H5Fclose\n",__LINE__);
    }

    hsize_t total_num_dsets=0, *group_sizes, *edge_counts, *hit_counts;
    hsize_t max_edge_count, min_edge_count, median_edge_count,
            max_hit_count,  min_hit_count,  median_hit_count;
    float total_dset_size=0, max_grp_size, min_grp_size, median_grp_size; 
    group_sizes = (hsize_t*) malloc(it_op.num_groups * sizeof(hsize_t));
    edge_counts = (hsize_t*) malloc(it_op.num_groups * sizeof(hsize_t));
    hit_counts = (hsize_t*) malloc(it_op.num_groups * sizeof(hsize_t));
    max_grp_size = 0;
    min_grp_size = INT_MAX;
    max_edge_count = 0;
    min_edge_count = INT_MAX;
    max_hit_count = 0;
    min_hit_count = INT_MAX;
    for (i=0; i<it_op.num_groups; i++) {
        total_num_dsets += it_op.groups[i].num_dsets;
        total_dset_size += it_op.groups[i].sum_dset_size;
        max_grp_size = MAX(max_grp_size, it_op.groups[i].sum_dset_size);
        min_grp_size = MIN(min_grp_size, it_op.groups[i].sum_dset_size);
        group_sizes[i] = it_op.groups[i].sum_dset_size;
        if (it_op.hasEdges) {
            max_edge_count = MAX(max_edge_count, it_op.groups[i].num_edges);
            min_edge_count = MIN(min_edge_count, it_op.groups[i].num_edges);
            edge_counts[i] = it_op.groups[i].num_edges;
        }
        if (it_op.hasHits) {
            max_hit_count = MAX(max_hit_count, it_op.groups[i].num_hits);
            min_hit_count = MIN(min_hit_count, it_op.groups[i].num_hits);
            hit_counts[i] = it_op.groups[i].num_hits;
        }
    }
    /* sort and find median of group sizes */
    qsort(group_sizes, it_op.num_groups, sizeof(hsize_t), cmpfunc);
    median_grp_size = group_sizes[(it_op.num_groups+1)/2];
    histogram(it_op.num_groups, group_sizes, event_bin_width, "group_size_hist.csv");

    if (it_op.hasEdges) {
        /* sort and find median of edge counts */
        qsort(edge_counts, it_op.num_groups, sizeof(hsize_t), cmpfunc);
        median_edge_count = edge_counts[(it_op.num_groups+1)/2];
    }
   
    if (it_op.hasHits) {
        /* sort and find median of hint counts */
        qsort(hit_counts, it_op.num_groups, sizeof(hsize_t), cmpfunc);
        median_hit_count = hit_counts[(it_op.num_groups+1)/2];
    }

    /* obtain file size, including external links */
    off_t fsize = 0;
    if (stat(fname, &file_stat) < 0)
        fprintf(stderr,"Error: stat on file %s (%s)\n",fname,strerror(errno));
    else {
        fsize = file_stat.st_size;
        for (i=0; i<it_op.num_links; i++) {
            char path[1024];
            if (it_op.ext_links[i][0] == '/') /* absolute path name */
                strcpy(path, it_op.ext_links[i]);
            else { /* relative path to input file name */
                strcpy(path, fname);
                dirname(path); /* remove file name, keep only directory name */
                sprintf(path+strlen(path), "/%s", it_op.ext_links[i]);
            }
            if (stat(path, &file_stat) < 0) {
                printf("Error: stat on file %s (%s)\n", path, strerror(errno));
                break;
            }
            fsize += file_stat.st_size;
            free(it_op.ext_links[i]);
        }
    }

    if (err_exit == 0) {
        printf("input HDF5 file name              = %s\n",fname);
        printf("number of HDF5 external links     = %12d\n",it_op.num_links);
        printf("input file size, including links  = %12ld B = %8.1f MiB = %5.1f GiB\n",
               fsize, (float)fsize/1048576.0, (float)fsize/1073741824.0);
        printf("total dataset size                = %12.f B = %8.1f MiB = %5.1f GiB\n",
               total_dset_size, total_dset_size/1048576.0, total_dset_size/1073741824.0);
        printf("number of groups in the file      = %12llu\n", it_op.num_groups);
        printf("total number of datasets          = %12llu\n", total_num_dsets);
        printf("MAXIMUM single group data size    = %12.f B = %8.1f MiB = %5.1f GiB\n",
               max_grp_size, max_grp_size/1048576.0, max_grp_size/1073741824.0);
        printf("MINIMUM single group data size    = %12.f B = %8.1f MiB = %5.1f GiB\n",
               min_grp_size, min_grp_size/1048576.0, min_grp_size/1073741824.0);
        printf("MEDIAN  single group data size    = %12.f B = %8.1f MiB = %5.1f GiB\n",
               median_grp_size, median_grp_size/1048576.0, median_grp_size/1073741824.0);
        if (it_op.hasEdges) {
            hsize_t num_bins;
            char *name = "edge_count_hist.csv";
            printf("MAXIMUM edge count in a graph     = %12llu\n", max_edge_count);
            printf("MINIMUM edge count in a graph     = %12llu\n", min_edge_count);
            printf("MEDIAN  edge count in a graph     = %12llu\n", median_edge_count);
            histogram(it_op.num_groups, edge_counts, edge_bin_width, name);
            num_bins = edge_counts[it_op.num_groups-1] / edge_bin_width;
            if (edge_counts[it_op.num_groups-1] % edge_bin_width) num_bins++;
            printf("Edge count histogram file         = %s\n", name);
            printf("    histogram bin width           = %12d\n", edge_bin_width);
            printf("    number of histogram bins      = %12llu\n", num_bins);
        }
        if (it_op.hasHits) {
            hsize_t num_bins;
            char *name = "hit_count_hist.csv";
            printf("MAXIMUM hit  count in a graph     = %12llu\n", max_hit_count);
            printf("MINIMUM hit  count in a graph     = %12llu\n", min_hit_count);
            printf("MEDIAN  hit  count in a graph     = %12llu\n", median_hit_count);
            histogram(it_op.num_groups, hit_counts, hit_bin_width, name);
            num_bins = hit_counts[it_op.num_groups-1] / hit_bin_width;
            if (hit_counts[it_op.num_groups-1] % hit_bin_width) num_bins++;
            printf("Hit count histogram file          = %s\n", name);
            printf("    histogram bin width           = %12d\n", hit_bin_width);
            printf("    number of histogram bins      = %12llu\n", num_bins);
        }
        if (evt_dset != NULL) {
            hsize_t num_bins;
            char *name = "event_size_hist.csv";
            printf("event ID dataset                  = %s\n", evt_dset);
            printf("Total number of event IDs         = %12llu\n", it_op.num_events);
            printf("Total number of non-empty events  = %12llu\n", non_empty_evts);
            printf("MAXIMUM single event data size    = %12.f B = %8.1f KiB = %5.1f MiB\n",
                   max_evt_size, max_evt_size/1024.0, max_evt_size/1048576.0);
            printf("MINIMUM single event data size    = %12.f B = %8.1f KiB = %5.1f MiB\n",
                   min_evt_size, min_evt_size/1024.0, min_evt_size/1048576.0);
            printf("MEDIAN  single event data size    = %12.f B = %8.1f KiB = %5.1f MiB\n",
                   median_evt_size, median_evt_size/1024.0, median_evt_size/1048576.0);
            histogram(it_op.num_events, it_op.evt_size, event_bin_width, name);
            num_bins = it_op.evt_size[it_op.num_events-1] / event_bin_width;
            if (it_op.evt_size[it_op.num_events-1] % event_bin_width) num_bins++;
            printf("Event size histogram file         = %s\n", name);
            printf("    histogram bin width           = %12d\n", event_bin_width);
            printf("    number of histogram bins      = %12llu\n", num_bins);
        }
    }

    free(group_sizes);
    free(it_op.ext_links);
    if (it_op.groups != NULL) free(it_op.groups);
    if (evt_dset != NULL) free(evt_dset);
    if (it_op.evt_group != NULL) free(it_op.evt_group);
    if (it_op.evt_size != NULL) free(it_op.evt_size);
    if (fname != NULL) free(fname);

    return (err_exit != 0);
}
