/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
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

#include <hdf5.h>

static int verbose;

#if defined PROFILE && PROFILE
#define SET_TIMER(ts) { ts = wtime(); }
#define GET_TIMER(ts, te, t) { te = wtime(); t += te - ts; ts = te; }

static double dopen_time, dcreate_time, dread_time, dwrite_time, dclose_time;

#include <sys/time.h>
static
double wtime(void)
{
    double          now_time;
    struct timeval  etstart;
    struct timezone tzp;

    if (gettimeofday(&etstart, &tzp) == -1) {
        fprintf(stderr, "Error: calling gettimeofday() not successful.\n");
        return 0.0;
    }

    now_time = ((double)etstart.tv_sec) +              /* in seconds */
               ((double)etstart.tv_usec) / 1000000.0;  /* in microseconds */
    return now_time;
}
#else
#define SET_TIMER(ts)
#define GET_TIMER(ts, te, t)
#endif

#define HANDLE_ERROR(msg,name) { \
    printf("Error at line %d: func %s on %s\n",__LINE__,msg,name); \
    err_exit = -1; \
    goto fn_exit; \
}

/* parameters to be passed to the call back function */
struct op_data {
    hid_t        fd_out;          /* output file descriptor */
    int          zero_as_scalar;  /* define zero-sized variables as scalars */
    int          true_1d;         /* define true 1D datasets */
    int          raw_chunk_cache; /* to enable chunk caching for raw data */
    hsize_t      chunk_unit_1D;   /* chunk for true-1D datasets */
    hsize_t      chunk_unit_2D;   /* chunk for true-2D datasets */
    hsize_t      chunk_unit_MB;   /* chunk base size in MiB for 1D datasets */
    unsigned int gzip_level;      /* GZIP compression level */
    hsize_t      num_1D_dset;     /* number of 1D datasets */
    hsize_t      num_2D_dset;     /* number of 2D datasets */
    hsize_t      num_nonzero_1D;  /* number of non-zero 1D datasets */
    hsize_t      num_nonzero_2D;  /* number of non-zero 2D datasets */
    hsize_t      num_change_chunk;/* number of datasets changed chunking */
    size_t       io_buf_size;     /* I/O buffer size */
    void        *io_buf;          /* I/O buffer */
    double       timing[4];       /* breakdown timings */
    herr_t       err;
};

/*----< rechunk() >----------------------------------------------------------*/
/* call back function used in H5Ovisit() */
static
herr_t rechunk(hid_t             loc_id,        /* object ID */
               const char       *name,          /* object name */
               const H5O_info_t *info,          /* object metadata */
               void             *operator_data) /* data passed from caller */
{
    int err_exit=0;
    herr_t err=0;
    struct op_data *it_op = (struct op_data*)operator_data;
    hid_t fd_out = it_op->fd_out;
#if defined PROFILE && PROFILE
    double ts=0.0, te=0.0;
#endif

    it_op->err = 0;

    if (info->type == H5O_TYPE_GROUP) {
        if (verbose) printf("input group name=%s\n",name);

        /* Skip root group */
        if (!strcmp(name, ".")) return 0;

        /* create a group in output file with the same name */
        hid_t grp_id = H5Gcreate1(fd_out, name, 0);
        if (grp_id < 0) HANDLE_ERROR("H5Gcreate1", name)

        /* close group */
        err = H5Gclose(grp_id);
        if (err < 0) HANDLE_ERROR("H5Gclose", name)
    }
    else if (info->type == H5O_TYPE_DATASET) {
        /* create a dataset in output file with the same name */
        size_t type_size;
        hid_t in_dset, out_dset, space_id, type_id, r_dcpl, w_dcpl;
        hsize_t dims[2], in_chunk_dims[2], out_chunk_dims[2];
        H5D_layout_t in_layout;
        H5T_class_t type_class;

        if (verbose) printf("input daatset name=%s\n",name);

        SET_TIMER(ts)
        if (! it_op->raw_chunk_cache) {
            /* create access property for read */
            hid_t r_dapl = H5Pcreate(H5P_DATASET_ACCESS);
            if (r_dapl < 0) HANDLE_ERROR("H5Pcreate r_dapl", name)

            /* if a dataset is read by whole chunks and there is no need to
             * access the chunks more than once on the disk, the chunk cache is
             * not needed and can be set to 0 if there is a shortage of memory.
             */
            err = H5Pset_chunk_cache(r_dapl, 0, 0, 1.0);
            if (err < 0) HANDLE_ERROR("H5Pset_chunk_cache r_dapl", name)

            /* Open input dataset */
            in_dset = H5Dopen(loc_id, name, r_dapl);
            if (in_dset < 0) HANDLE_ERROR("H5Dopen", name)

            err = H5Pclose(r_dapl);
            if (err < 0) HANDLE_ERROR("H5Pclose r_dapl", name)
        }
        else {
            /* Open input dataset */
            in_dset = H5Dopen(loc_id, name, H5P_DEFAULT);
            if (in_dset < 0) HANDLE_ERROR("H5Dopen", name)
        }

        /* Retrieve input dataset's creation property list */
        r_dcpl = H5Dget_create_plist(in_dset);
        if (r_dcpl < 0) HANDLE_ERROR("H5Dget_create_plist", name)

        /* Retrieve input dataset's data layout */
        in_layout = H5Pget_layout(r_dcpl);
        if (in_layout < 0) HANDLE_ERROR("H5Pget_layout", name)

        /* Retrieve input dataset's data chunk dimensions */
        if (in_layout == H5D_CHUNKED) {
            err = H5Pget_chunk(r_dcpl, 2, in_chunk_dims);
            if (err < 0) HANDLE_ERROR("H5Pget_chunk", name)

            assert(1 == H5Pget_nfilters(r_dcpl));
        }

        /* close input dataset creation property */
        err = H5Pclose(r_dcpl);
        if (err < 0) HANDLE_ERROR("H5Pclose r_dcpl", name)

        /* Retrieve input dataset's data type */
        type_id = H5Dget_type(in_dset);
        if (type_id < 0) HANDLE_ERROR("H5Dget_type", name)

        /* Retrieve data type size */
        type_size = H5Tget_size(type_id);
        if (type_size == 0) HANDLE_ERROR("H5Tget_size", name)

        /* Retrieve data type class */
        type_class = H5Tget_class(type_id);
        if (type_class < 0) HANDLE_ERROR("H5Tget_class", name);

        /* Open input dataset's space */
        space_id = H5Dget_space(in_dset);
        if (space_id < 0) HANDLE_ERROR("H5Dget_space", name)

        /* Retrieve number of dimensions */
        int ndims = H5Sget_simple_extent_ndims(space_id);
        if (ndims < 0) HANDLE_ERROR("H5Sget_simple_extent_dims", name)
        assert(ndims == 2);

        /* retrieve dimension sizes */
        ndims = H5Sget_simple_extent_dims(space_id, dims, NULL);
        if (ndims < 0) HANDLE_ERROR("H5Sget_simple_extent_dims", name)

        assert(dims[0] >= 0);
        assert(dims[1] >= 0);

        if (dims[1] == 1)
            it_op->num_1D_dset++;
        else
            it_op->num_2D_dset++;
        GET_TIMER(ts, te, dopen_time)

        /* create dataset creation property for new dataset */
        w_dcpl = H5Pcreate(H5P_DATASET_CREATE);
        if (w_dcpl < 0) HANDLE_ERROR("H5Pcreate", name)

        /* Never fill in the dataset in advance except H5T_STRING, as
         * HDF5 requires fill mode enabled for H5T_STRING datasets.
         */
        if (type_class != H5T_STRING) {
            err = H5Pset_fill_time(w_dcpl, H5D_FILL_TIME_NEVER);
            if (err < 0) HANDLE_ERROR("H5Pset_fill_time", name)
        }

        /* For new dataset, determine if use chunk or compact layout */
        if (dims[0] == 0 || dims[1] == 0) { /* 0-size dataset */
            /* use COMPACT layout */
            err = H5Pset_layout(w_dcpl, H5D_COMPACT);
            if (err < 0) HANDLE_ERROR("H5Pset_layout", name)

            /* do not reuse data space */
            err = H5Sclose(space_id);
            if (err < 0) HANDLE_ERROR("H5Sclose", name)

            if (it_op->zero_as_scalar)
                /* Define zero-sized datasets as scalars. */
                space_id = H5Screate_simple(0, dims, NULL);
            else
                /* use fixed-size dimensions */
                space_id = H5Screate_simple(2, dims, NULL);

            if (space_id < 0) HANDLE_ERROR("H5Screate_simple", name)

            if (in_layout != H5D_COMPACT) it_op->num_change_chunk++;
        }
        else { /* not a 0-size dataset */
            if (it_op->true_1d && dims[1] == 1) { /* 1D dataset */
                /* define a true 1D dataset */
                err = H5Sclose(space_id);
                if (err < 0) HANDLE_ERROR("H5Sclose", name)
                space_id = H5Screate_simple(1, dims, NULL);
                if (space_id < 0) HANDLE_ERROR("H5Screate_simple", name)
            }

            /* use CHUNK layout */
            err = H5Pset_deflate(w_dcpl, it_op->gzip_level);
            if (err < 0) HANDLE_ERROR("H5Pset_deflate", name)

            if (dims[1] > 1) { /* 2D dataset */
                out_chunk_dims[0] = it_op->chunk_unit_2D;
                if (dims[0] < it_op->chunk_unit_2D)
                    out_chunk_dims[0] = dims[0];
            }
            else if (dims[0] > it_op->chunk_unit_1D) { /* large 1D dataset */
                if (it_op->chunk_unit_MB > 0)
                    out_chunk_dims[0] = it_op->chunk_unit_MB * 1048576 / type_size;
                else
                    out_chunk_dims[0] = it_op->chunk_unit_1D;
            }
            else /* small 1D dataset: one chunk */
                out_chunk_dims[0] = dims[0];

            out_chunk_dims[1] = dims[1];
            if (it_op->true_1d && dims[1] == 1) /* 1D dataset */
                /* chunk along 1D */
                err = H5Pset_chunk(w_dcpl, 1, out_chunk_dims);
            else
                err = H5Pset_chunk(w_dcpl, 2, out_chunk_dims);

            if (err < 0) HANDLE_ERROR("H5Pset_chunk", name)
        }

        if (! it_op->raw_chunk_cache) {
            /* write access property */
            hid_t w_dapl = H5Pcreate(H5P_DATASET_ACCESS);
            if (w_dapl < 0) HANDLE_ERROR("H5Pcreate w_dapl", name)

            /* As all chunks will be accessed only once, disable write chunk
             * caching in hope to skip some HDF5 internal memcpy due to
             * caching.
             */
            err = H5Pset_chunk_cache(w_dapl, 0, 0, 1.0);
            if (err < 0) HANDLE_ERROR("H5Pset_chunk_cache w_dapl", name)

            /* create output dataset */
            out_dset = H5Dcreate(fd_out, name, type_id, space_id, H5P_DEFAULT,
                                 w_dcpl, w_dapl);
            if (out_dset < 0) HANDLE_ERROR("H5Dcreate", name)

            err = H5Pclose(w_dapl);
            if (err < 0) HANDLE_ERROR("H5Pclose w_dapl", name)
        }
        else {
            /* create output dataset */
            out_dset = H5Dcreate(fd_out, name, type_id, space_id, H5P_DEFAULT,
                                 w_dcpl, H5P_DEFAULT);
            if (out_dset < 0) HANDLE_ERROR("H5Dcreate", name)
        }

        err = H5Pclose(w_dcpl);
        if (err < 0) HANDLE_ERROR("H5Pclose w_dcpl", name)

        err = H5Sclose(space_id);
        if (err < 0) HANDLE_ERROR("H5Sclose", name)
        GET_TIMER(ts, te, dcreate_time)

        if (dims[0] == 0 || dims[1] == 0) { /* return now for 0-size dataset */
            err = H5Tclose(type_id);
            if (err < 0) HANDLE_ERROR("H5Tclose", name)
            err = H5Dclose(in_dset);
            if (err < 0) HANDLE_ERROR("H5Dclose", name)
            err = H5Dclose(out_dset);
            if (err < 0) HANDLE_ERROR("H5Dclose", name)
            return 0;
        }

        if (dims[1] == 1)
            it_op->num_nonzero_1D++;
        else
            it_op->num_nonzero_2D++;

        /* when chunk setting is the same, do direct read and write chunks to
         * avoid data decompression and compression
         */
        if (in_chunk_dims[0] == out_chunk_dims[0] &&
            in_chunk_dims[1] == out_chunk_dims[1]) {
            hsize_t ii, ntimes, start[2], chunk_size;
            unsigned int filter_mask;

            /* read and write one chunk at a time */
            ntimes = dims[0] / in_chunk_dims[0];
            if (dims[0] % in_chunk_dims[0]) ntimes++;

            start[0] = 0;
            start[1] = 0;
            for (ii=0; ii<ntimes; ii++) {
                err = H5Dget_chunk_storage_size(in_dset, start, &chunk_size);
                if (err < 0) HANDLE_ERROR("H5Dget_chunk_storage_size", name)

                /* read a chunk in compressed form */
                err = H5Dread_chunk(in_dset, H5P_DEFAULT, start, &filter_mask,
                                    it_op->io_buf);
                if (err < 0) HANDLE_ERROR("H5Dread_chunk", name)
                GET_TIMER(ts, te, dread_time)

                /* write a chunk in compressed form */
                err = H5Dwrite_chunk(out_dset, H5P_DEFAULT, filter_mask, start,
                                     chunk_size, it_op->io_buf);
                if (err < 0) HANDLE_ERROR("H5Dwrite_chunk", name)

                /* move on to next chunk */
                start[0] += in_chunk_dims[0];
                GET_TIMER(ts, te, dwrite_time)
            }
        }
        else if (dims[0] * dims[1] * type_size <= it_op->io_buf_size) {
            /* Read the entire input dataset */
            err = H5Dread(in_dset, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          it_op->io_buf);
            if (err < 0) HANDLE_ERROR("H5Dread", name)
            GET_TIMER(ts, te, dread_time)

            /* Write the entire output dataset */
            err = H5Dwrite(out_dset, type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           it_op->io_buf);
            if (err < 0) HANDLE_ERROR("H5Dwrite", name)

            it_op->num_change_chunk++;

            GET_TIMER(ts, te, dwrite_time)
        }
        else {
            hid_t memspace_id, in_space_id, out_space_id;
            hsize_t ii, one[2], start[2], count[2], ntimes, remain;

            one[0] = one[1] = 1;
            start[0] = 0;
            start[1] = 0;
            count[0] = it_op->io_buf_size / (dims[1] * type_size);
            count[1] = dims[1];
            remain = dims[0];

            assert(count[0] > 0);

            ntimes = dims[0] / count[0];
            if (dims[0] % count[0]) ntimes++;

            /* Setup memory space */
            memspace_id = H5Screate_simple(2, count, NULL);
            if (memspace_id < 0) HANDLE_ERROR("H5Screate_simple", name)

            /* Setup hyperslab file space */
            in_space_id = H5Dget_space(in_dset);
            if (in_space_id < 0) HANDLE_ERROR("H5Dget_space", name)

            out_space_id = H5Dget_space(out_dset);
            if (out_space_id < 0) HANDLE_ERROR("H5Dget_space", name)

            for (ii=0; ii<ntimes; ii++) {
                err = H5Sselect_hyperslab(in_space_id, H5S_SELECT_SET, start,
                                          NULL, one, count);
                if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab", name)

                /* Read the input dataset in batches */
                err = H5Dread(in_dset, type_id, memspace_id, in_space_id,
                              H5P_DEFAULT, it_op->io_buf);
                if (err < 0) HANDLE_ERROR("H5Dread", name)
                GET_TIMER(ts, te, dread_time)

                err = H5Sselect_hyperslab(out_space_id, H5S_SELECT_SET, start,
                                          NULL, one, count);
                if (err < 0) HANDLE_ERROR("H5Sselect_hyperslab", name)

                /* Write the output dataset in batches */
                err = H5Dwrite(out_dset, type_id, memspace_id, out_space_id,
                               H5P_DEFAULT, it_op->io_buf);
                if (err < 0) HANDLE_ERROR("H5Dwrite", name)

                start[0] += count[0];
                remain   -= count[0];
                if (remain < count[0]) {
                    count[0] = remain;
                    /* update memspace */
                    err = H5Sset_extent_simple(memspace_id, 2, count, NULL);
                    if (err < 0) HANDLE_ERROR("H5Sset_extent_simple", name)
                }
                GET_TIMER(ts, te, dwrite_time)
            }
            err = H5Sclose(memspace_id);
            if (err < 0) HANDLE_ERROR("H5Sclose", name)
            err = H5Sclose(in_space_id);
            if (err < 0) HANDLE_ERROR("H5Sclose", name)
            err = H5Sclose(out_space_id);
            if (err < 0) HANDLE_ERROR("H5Sclose", name)

            it_op->num_change_chunk++;
        }

        err = H5Tclose(type_id);
        if (err < 0) HANDLE_ERROR("H5Tclose", name)
        err = H5Dclose(in_dset);
        if (err < 0) HANDLE_ERROR("H5Dclose", name)
        err = H5Dclose(out_dset);
        if (err < 0) HANDLE_ERROR("H5Dclose", name)
        GET_TIMER(ts, te, dclose_time)
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
    if (err < 0) HANDLE_ERROR("H5Pget_cache", __func__)

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
    if (err < 0) HANDLE_ERROR("H5Pset_cache", __func__)
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
  [-d]         disable in-memory I/O (default: enable)\n\
  [-r]         disable chunk caching for raw data (default: enable)\n\
  [-s]         re-define zero-sized datasets as scalars (default: no)\n\
  [-t]         define true 1D dataset (default: no)\n\
  [-c size]    chunk size along dimension 1 for true-1D datasets (default: 1048576)\n\
  [-C size]    chunk size along dimension 0 for true-2D datasets (default: 128)\n\
  [-M size]    chunk base size in MiB for true-1D datasets (overwrites -c option)\n\
  [-z level]   GZIP compression level (default: 6)\n\
  [-b size]    I/O buffer size in bytes (default: 1 GiB)\n\
  [-o outfile] output file name (default: out.h5)\n\
  infile       input file name (required and cannot be the same as output file)\n\n\
  This utility program copied an input HDF5 file to a new file using an\n\
  adjusted chunk settings in the new file.\n\n\
  Requirements of the input HDF5 file:\n\
    1. contains multiple groups only at root level\n\
    2. each group contains multiple 2D datasets\n\
    3. true-1D datasets are those whose 2nd dimension size is 1\n\
    4. true-2D datasets are those whose 2nd dimension size is larger than 1\n\n\
  Output HDF5 file:\n\
    1. zero-sized datasets will be stored in H5D_COMPACT layout\n\
    2. for non-zero sized datasets, chunking is only applied to 1st dimension\n\
    3. for true-1D datasets, new chunk dimensions use value from -c option\n\
    4. for true-2D datasets, new chunk dimensions use value from -C option\n\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v|-d|-r|-s|-t] [-c size] [-C size] [-M size] [-z level] [-b size] [-o outfile] infile\n%s\n",
           progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int c, err_exit=0, in_memory_io=1, raw_chunk_cache=1;
    char *infile=NULL, *outfile=NULL;
    herr_t err;
    hsize_t meta_block_size;
    hid_t fd_in=-1, fd_out=-1, fapl_id=-1;
    H5G_info_t grp_info;
    struct op_data it_op;
#if defined PROFILE && PROFILE
    double ts=0.0, te=0.0, timing[4]={0,0,0,0};
    dopen_time = dcreate_time = dread_time = dwrite_time = dclose_time = 0.0;
#endif

    verbose = 0; /* default is quiet */

    it_op.zero_as_scalar = 0;       /* define zero-sized datasets as scalars */
    it_op.true_1d        = 0;       /* define true 1D datasets */
    it_op.gzip_level     = 6;       /* default GZIP compression level */
    it_op.chunk_unit_1D  = 1048576; /* default chunk size for 1D datasets */
    it_op.chunk_unit_2D  = 128;     /* default chunk size for 2D datasets */
    it_op.chunk_unit_MB  = 0;       /* chunk base size in MiB for 1D datasets */
    it_op.io_buf_size    = 1073741824; /* default I/O buffer size */
    it_op.io_buf         = NULL;       /* I/O buffer */
    it_op.num_1D_dset    = 0;     /* number of 1D datasets */
    it_op.num_2D_dset    = 0;     /* number of 2D datasets */
    it_op.num_nonzero_1D = 0;     /* number of non-zero 1D datasets */
    it_op.num_nonzero_2D = 0;     /* number of non-zero 2D datasets */
    it_op.num_change_chunk = 0;   /* number of datasets with same chunking */

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvdrstc:z:o:c:C:M:b:")) != -1)
        switch(c) {
            case 'h': usage(argv[0]);
                      err_exit = -1;
                      goto fn_exit;
            case 'v': verbose = 1;
                      break;
            case 'd': in_memory_io = 0;
                      break;
            case 'r': raw_chunk_cache = 0;
                      break;
            case 's': it_op.zero_as_scalar = 1;
                      break;
            case 't': it_op.true_1d = 1;
                      break;
            case 'c': it_op.chunk_unit_1D = (hsize_t)(atoi(optarg));
                      break;
            case 'C': it_op.chunk_unit_2D = (hsize_t)(atoi(optarg));
                      break;
            case 'M': it_op.chunk_unit_MB = (hsize_t)(atoi(optarg));
                      break;
            case 'z': it_op.gzip_level = (unsigned int)(atoi(optarg));
                      break;
            case 'b': it_op.io_buf_size = (size_t)(atoll(optarg));
                      break;
            case 'o': outfile = strdup(optarg);
                      break;
            default: break;
        }

    if (argv[optind] == NULL) { /* input file name is mandatory */
        printf("input file is missing\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    infile = strdup(argv[optind]);

    if (verbose) printf("input file: %s\n", infile);

    if (outfile == NULL)  /* set default output file name */
        outfile = strdup("out.h5");

    if (verbose) printf("output file: %s\n", outfile);

#ifdef HAVE_ACCESS
    /* if access() is available, use it to check whether file already exists */
    if (access(outfile, F_OK) == 0) { /* output file already existed */
        fprintf(stderr, "Error: output file already existed %s\n", outfile);
        exit(1);
    }
#endif

    SET_TIMER(ts)

    /* open input file for read */
    if (in_memory_io) { /* enable in-memory I/O */
        /* create file access property */
        fapl_id = H5Pcreate(H5P_FILE_ACCESS);
        if (fapl_id < 0) HANDLE_ERROR("H5Pcreate", infile)

        err = H5Pset_fapl_core(fapl_id, 33554432, 1);
        if (err < 0) HANDLE_ERROR("H5Pset_fapl_core", infile)
    }
    else
        fapl_id = H5P_DEFAULT;

    /* open input file in read-only mode */
    fd_in = H5Fopen(infile, H5F_ACC_RDONLY, fapl_id);
    if (fd_in < 0) HANDLE_ERROR("Can't open input file", infile)

    if (in_memory_io) { /* in-memory I/O is disabled */
        err = H5Pclose(fapl_id);
        if (err < 0) HANDLE_ERROR("H5Pclose", infile)
        fapl_id = H5P_DEFAULT;
    }

#if 0
    if (verbose) {
        /* retrieve the property list used by the input file */
        fapl_id = H5Fget_access_plist(fd_in);
        if (fapl_id < 0) HANDLE_ERROR("H5Fget_access_plist", infile)

        err = H5Pget_meta_block_size(fapl_id, &meta_block_size);
        if (err < 0) HANDLE_ERROR("H5Pget_meta_block_size", infile)
        printf("metadata block size used by the input file is %lld\n",
               meta_block_size);
        err = H5Pclose(fapl_id);
        if (err < 0) HANDLE_ERROR("H5Pclose", infile)
        fapl_id = H5P_DEFAULT;
    }
#endif

    err = H5Gget_info_by_name(fd_in, "/", &grp_info, H5P_DEFAULT);
    if (err < 0) HANDLE_ERROR("H5Gget_info_by_name - root group",  infile)

    GET_TIMER(ts, te, timing[0])

    /* create file access property for creating output file */
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id < 0) HANDLE_ERROR("H5Pcreate", outfile)

    /* increase metadata block size to 1 MiB */
    err = H5Pset_meta_block_size(fapl_id, 1048576);
    if (err < 0) HANDLE_ERROR("H5Pset_meta_block_size", outfile)

/* This setting degrades the performance.
    err = incr_raw_data_chunk_cache(fapl_id);
    if (err < 0) HANDLE_ERROR("incr_raw_data_chunk_cache", outfile)
*/
    it_op.raw_chunk_cache = raw_chunk_cache;

    /* create output file */
    fd_out = H5Fcreate(outfile, H5F_ACC_EXCL, H5P_DEFAULT, fapl_id);
    if (fd_out < 0){
        printf("Error: Fail to create output file %s\n", outfile);
        printf("Error: output file may already exist\n");
        err_exit = -1;
        goto fn_exit;
    }
    it_op.fd_out = fd_out;

    if (verbose) {
        err = H5Pget_meta_block_size(fapl_id, &meta_block_size);
        if (err < 0) HANDLE_ERROR("H5Pget_meta_block_size", outfile)
        printf("metadata block size used by the output file is set to %lld\n",
               meta_block_size);
    }

    /* close the access property */
    err = H5Pclose(fapl_id);
    if (err < 0) HANDLE_ERROR("H5Pclose", outfile)
    fapl_id = H5P_DEFAULT;

    GET_TIMER(ts, te, timing[1])

    /* Allocate I/O buffer */
    it_op.io_buf = (void*) malloc(it_op.io_buf_size);

    /* Iterate all objects and perform chunking adjustment */
#if defined HAS_H5OVISIT3 && HAS_H5OVISIT3
    err = H5Ovisit3(fd_in, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, rechunk, &it_op, H5O_INFO_ALL);
    if (err < 0) HANDLE_ERROR("H5Ovisit3", infile)
#else
    err = H5Ovisit(fd_in, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, rechunk, &it_op);
    if (err < 0) HANDLE_ERROR("H5Ovisit", infile)
#endif
    if (it_op.err < 0) HANDLE_ERROR("H5Ovisit", infile)

    GET_TIMER(ts, te, timing[2])

fn_exit:
    if (fd_in >= 0) {
        check_h5_objects(infile, fd_in);
        err = H5Fclose(fd_in);
        if (err < 0) printf("Error at line %d: H5Fclose\n",__LINE__);
    }
    if (fd_out >= 0) {
        check_h5_objects(outfile, fd_out);
        err = H5Fclose(fd_out);
        if (err < 0) printf("Error at line %d: H5Fclose\n",__LINE__);
    }
    if (fapl_id >= 0 && fapl_id != H5P_DEFAULT) {
        err = H5Pclose(fapl_id);
        if (err < 0) printf("Error at line %d: H5Pclose\n",__LINE__);
    }

    if (it_op.io_buf != NULL) free(it_op.io_buf);
    if (outfile != NULL) free(outfile);
    if (infile  != NULL) free(infile);
    GET_TIMER(ts, te, timing[3])

#if defined PROFILE && PROFILE
    if (err_exit == 0) {
        printf("In-memory I/O                        = %s\n",
               (in_memory_io)?"enabled":"disabled");
        printf("Chunk caching for raw data           = %s\n",
               (raw_chunk_cache)?"enabled":"disabled");
        printf("I/O buffer size                      = %zd bytes\n", it_op.io_buf_size);
        printf("number of groups in the file         = %llu\n", grp_info.nlinks);
        printf("total number of 1D datasets          = %llu\n", it_op.num_1D_dset);
        printf("total number of 2D datasets          = %llu\n", it_op.num_2D_dset);
        printf("no. non-zero 1D datasets             = %llu\n", it_op.num_nonzero_1D);
        printf("no. non-zero 2D datasets             = %llu\n", it_op.num_nonzero_2D);
        printf("no. datasets chunking changed        = %llu\n", it_op.num_change_chunk);
        printf("0th  dim  chunk size for 2D datasets = %llu\n", it_op.chunk_unit_2D);
        printf("1st  dim  chunk size for 1D datasets = %llu\n", it_op.chunk_unit_1D);
        printf("MiB-based chunk size for 1D datasets = %llu MiB\n", it_op.chunk_unit_MB);
        printf("-------------------------------------------------------\n");
        printf("Input  file open   time              = %7.2f sec\n", timing[0]);
        printf("Output file create time              = %7.2f sec\n", timing[1]);
        printf("Re-chunk time                        = %7.2f sec\n", timing[2]);
        printf("  Re-chunk dataset open   time       = %7.2f sec\n", dopen_time);
        printf("  Re-chunk dataset create time       = %7.2f sec\n", dcreate_time);
        printf("  Re-chunk dataset read   time       = %7.2f sec\n", dread_time);
        printf("  Re-chunk dataset write  time       = %7.2f sec\n", dwrite_time);
        printf("  Re-chunk dataset close  time       = %7.2f sec\n", dclose_time);
        printf("Output file close time               = %7.2f sec\n", timing[3]);
        printf("-------------------------------------------------------\n");
        printf("Total time                           = %7.2f sec\n",
               timing[0] + timing[1] + timing[2] + timing[3]);
    }
#endif
    return (err_exit != 0);
}
