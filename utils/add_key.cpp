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
#include <libgen.h> /* dirname(), basename() */
#include <unistd.h> /* getopt() */
#include <assert.h> /* assert() */
#include <regex.h>

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <bits/stdc++.h>
using namespace std;

#include <hdf5.h>

static int verbose;

#ifndef MAX
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#endif

#if defined PROFILE && PROFILE
#define SET_TIMER(ts) { ts = wtime(); }
#define GET_TIMER(ts, te, t) { te = wtime(); t += te - ts; ts = te; }

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

#define RETURN_ERROR(func_name, obj_name) { \
    cerr<<"Error at line "<<__LINE__<<": calling "<<func_name<<" for object "<<obj_name<<endl; \
    return -1; \
}

#define HANDLE_ERROR(msg) { \
    cerr<<"Error at line "<<__LINE__<<": "<<msg<<endl; \
    err_exit = -1; \
    goto fn_exit; \
}

#define CALLBACK_ERROR(func_name, dset_name) { \
    cerr<<"Error at line "<<__LINE__<<": calling "<<func_name<<" for group "<<grp_name<<" dataset "<<dset_name<<endl; \
    return -1; \
}

/* lookup table based on n-tuple datasets */
typedef vector<long long> keyv;

struct hashv : public unary_function<keyv, size_t> {
    size_t operator()(const keyv& key) const {
        size_t h = 0;
        for(size_t ikey = 0; ikey < key.size(); ikey++)
            h ^= key[ikey];
        return h;
    }
};

typedef unordered_map<keyv, int64_t, hashv> table;

/* string-based basename */
string basename(string full_path)
{
    auto pos = full_path.find_last_of("/\\");
    if(pos != string::npos)
      return full_path.substr(pos+1);
    else
      return full_path;
}

static
table build_lookup_table(hsize_t    nelems,
                         hsize_t    ndsets,
                         long long *buf)
{
    hsize_t ii, jj;
    table ret;

    for (ii=0; ii<nelems; ii++) {
        keyv key(ndsets);
        for (jj=0; jj<ndsets; jj++)
            key[jj] = buf[nelems*jj + ii];
        ret[key] = ii;
    }

    return ret;
}

/* parameters to be passed to the call back function */
struct op_data {
    size_t   num_groups;
    char   **grp_names;     /* [num_groups] group names */
    herr_t   err;
    regex_t * regex;

    bool ignore(const char * group)
    {
        if(regex == NULL) return false;
        return regexec(regex, group, 0, NULL, 0) == 0;
    }
};

/*----< gather_grp_names() >-------------------------------------------------*/
/* call back function used in H5Ovisit() */
static
herr_t gather_grp_names(hid_t             loc_id,/* object ID */
                        const char       *name,  /* object name */
                        const H5O_info_t *info,  /* object metadata */
                        void             *op)    /* data passed from caller */
{
    struct op_data *it_op = (struct op_data*)op;
    it_op->err = 0;

    /* skip if this is a dataset */
    if (info->type == H5O_TYPE_DATASET) return 0;

    /* skip groups that match user's pattern */
    bool ignore = it_op->ignore(name);

    /* info->type == H5O_TYPE_GROUP */
    if (verbose)
        printf(ignore? "group name=%s (Ignoring)\n" : "group name=%s\n",name);

    /* skip root group '.' the first object */
    if (!strcmp(name, ".") || ignore)  return 0;

    it_op->grp_names[it_op->num_groups++] = strdup(name);
    return 0;
}

/*----< add_seq() >---------------------------------------------------------*/
static
herr_t add_seq(hid_t           fid,
               int             overwrite,
               vector<string>  index_levels,
               const char     *grp_name,
               const char     *part_key_base,
               table           lookup_table,
               int             create_seq,
               int             create_seq_cnt,
               size_t         *num_nonzero_groups,
               size_t         *max_len,    /* max seq size */
               size_t         *min_len,    /* min seq size */
               size_t         *avg_len,    /* avg seq size */
               double         *timing)
{
    const char *dset_name;
    char seq_name[1024], seq_cnt_name[1024];
    int ndims, isExist;
    size_t ii, jj, ndsets, chunk_size;
    int64_t *seq_buf;
    herr_t err;
    hid_t grp_id=-1, dset, memspace, seq_id, seq_cnt_id;
    hid_t dcpl_id, in_fspace=-1, out_fspace;
    hsize_t dset_dims[2], maxdims[2], chunk_dims[2];
    hsize_t start[2], count[2], ones[2]={1,1};
    htri_t src_exist;
    string dset_basename;
    map<int64_t, int64_t> mp;

#if defined PROFILE && PROFILE
    double ts=0.0, te=0.0;
#endif

    SET_TIMER(ts)

    /* open group */
    grp_id = H5Gopen(fid, grp_name, H5P_DEFAULT);
    if (grp_id < 0) RETURN_ERROR("H5Gopen", grp_name)

    /* first index dataset, for example 'run' in tuple (run, subrun, evt) */
    dset_name = index_levels[0].c_str();
    if (verbose)
        printf("%s: dataset %s/%s\n",__func__,grp_name,dset_name);

    /* check if the first index dataset exists */
    src_exist = H5Lexists(grp_id, dset_name, H5P_DEFAULT);
    if (src_exist < 0) RETURN_ERROR("H5Lexists", string(grp_name)+dset_name)
    if (src_exist == 0) {
        if (verbose)
            printf("Warn: dataset '%s/%s' not exist, skip this group.\n",
                   grp_name,dset_name);
        goto fn_exit;
    }

    /* open the first index dataset in this group */
    dset = H5Dopen(grp_id, dset_name, H5P_DEFAULT);
    if (dset < 0) CALLBACK_ERROR("H5Dopen", dset_name)

    /* Note size of dimension 0 of all datasets in a group should be the same */
    in_fspace = H5Dget_space(dset);
    if (in_fspace < 0) CALLBACK_ERROR("H5Dget_space", dset_name)
    ndims = H5Sget_simple_extent_dims(in_fspace, dset_dims, maxdims);
    if (ndims != 2) CALLBACK_ERROR("H5Sget_simple_extent_dims", dset_name)

    /* if the dataset is zero-sized, just create dataset */
    if (dset_dims[0] == 0) {
        err = H5Dclose(dset);
        if (err < 0) CALLBACK_ERROR("H5Dclose", dset_name)
        goto seq_create;
    }

    /* non-zero size group */
    (*num_nonzero_groups)++;

    /* collect statistics about the dataset sizes across groups */
    if (*num_nonzero_groups == 1) {
        *max_len = dset_dims[0];
        *min_len = dset_dims[0];
        *avg_len = dset_dims[0];
    }
    else {
        *max_len = MAX(*max_len, dset_dims[0]);
        *min_len = MIN(*min_len, dset_dims[0]);
        *avg_len += dset_dims[0];
    }

    /* number of index datasets, for example 3 in tuple (run, subrun, evt) */
    ndsets = index_levels.size();

seq_create:
    /* dataset create property */
    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (dcpl_id < 0) CALLBACK_ERROR("H5Pcreate", seq_name)

    /* Never fill in the chunk in advance. */
    err = H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);
    if (err < 0) CALLBACK_ERROR("H5Pset_fill_time", seq_name)

    /* setting H5D_ALLOC_TIME_LATE takes no effect for parallel I/O on
     * compressed datasets, which HDF5 always uses H5D_ALLOC_TIME_EARLY.
     */
    err = H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_LATE);
    if (err < 0) CALLBACK_ERROR("H5Pset_alloc_time", seq_name)

    /* set the chunk size to 256K x 1 elements */
    chunk_dims[0] = 256 * 1024;
    chunk_dims[1] = 1;
    err = H5Pset_chunk(dcpl_id, 2, chunk_dims);
    if (err < 0) CALLBACK_ERROR("H5Pset_chunk", seq_name)

    /* enable compression */
    err = H5Pset_deflate(dcpl_id, 6);
    if (err < 0) CALLBACK_ERROR("H5Pset_deflate", seq_name)

    if (create_seq) {
        /* seq dataset will be of the same dimension lengths */
        hsize_t maxdims[2] = {H5S_UNLIMITED, 1};
        hsize_t seq_dims[2] = {dset_dims[0], 1};
        out_fspace = H5Screate_simple(ndims, seq_dims, maxdims);
        if (err < 0) CALLBACK_ERROR("H5Screate_simple", seq_name)

        /* name of the new partition key dataset */
        sprintf(seq_name, "%s.seq", part_key_base);

        /* check if the dataset already exists */
        isExist = H5Lexists(grp_id, seq_name, H5P_DEFAULT);
        if (isExist >= 0) {
            if (overwrite) {
                if (verbose) printf("Delete dataset %s/%s\n",grp_name,seq_name);
                err = H5Ldelete(grp_id, seq_name, H5P_DEFAULT);
                if (err < 0) CALLBACK_ERROR("H5Ldelete", seq_name)
            }
            else {
                cerr<<"Error: dataset " << grp_name << "/" << seq_name << " already exists"<<endl;
                return -1;
            }
        }

        /* create the new partition key dataset */
        if (verbose) printf("Create dataset %s/%s\n",grp_name,seq_name);
        seq_id = H5Dcreate2(grp_id, seq_name, H5T_STD_I64LE, out_fspace,
                            H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        if (seq_id < 0) CALLBACK_ERROR("H5Dcreate2", seq_name)
    }

    /* skip reading/writeing if zero-size dataset */
    if (dset_dims[0] == 0) {
        if (create_seq_cnt) goto seq_cnt_create;
        goto fn_exit;
    }

#define READ_CHUNK (8 * 1024 * 1024)

    /* set the read and write buffer chunk size */
    chunk_size = MIN(READ_CHUNK, dset_dims[0]);
    if (verbose) printf("Read chunk size = %zd\n",chunk_size);

    /* allocate write buffer for partition key dataset, seq */
    seq_buf = (int64_t*) malloc(chunk_size * sizeof(int64_t));

    if (ndsets > 1) {
        /* index datasets are stored as separate datasets, e.g. run, subrun,
         * and evt.
         */
        long long *index_buf=NULL, *index_buf_ptr;
        size_t remain=dset_dims[0];
        hsize_t off=0;
        hid_t *dset_ids;

        /* array storing dataset IDs */
        dset_ids = (hid_t*) malloc(ndsets * sizeof(hid_t));
        dset_ids[0] = dset;

        /* open remaining index datasets, the 1st one has been opened earlier */
        for (ii=1; ii<ndsets; ii++) {
            /* check if dataset exists */
            src_exist = H5Lexists(grp_id, index_levels[ii].c_str(), H5P_DEFAULT);
            if (src_exist < 0) RETURN_ERROR("H5Lexists", string(grp_name)+index_levels[ii])
            if (src_exist == 0) {
                if (verbose)
                    printf("Warn: dataset '%s/%s' not exist, skip this group.\n",
                           grp_name,index_levels[ii].c_str());
                goto fn_exit;
            }

            /* open dataset */
            dset_ids[ii] = H5Dopen(grp_id, index_levels[ii].c_str(), H5P_DEFAULT);
            if (dset < 0) CALLBACK_ERROR("H5Dopen", index_levels[ii])
        }

        /* allocate read buffer */
        index_buf = (long long*) malloc(chunk_size * ndsets * sizeof(long long));

        /* loop over the index datasets until all data has been read */
        while (remain > 0) {
            /* set file space for the read chunk */
            start[0] = off;
            start[1] = 0;
            count[0] = MIN(chunk_size, remain);
            count[1] = 1;
            err = H5Sselect_hyperslab(in_fspace, H5S_SELECT_SET, start, NULL, ones, count);
            if (err < 0) CALLBACK_ERROR("H5Sselect_hyperslab", seq_name)

            /* set memory space, a contiguous buffer in memory */
            memspace = H5Screate_simple(1, count, NULL);

            /* read a chunk from each index datasets */
            index_buf_ptr = index_buf;
            for (ii=0; ii<ndsets; ii++) {
                /* read next chunk of the dataset */
                err = H5Dread(dset_ids[ii], H5T_NATIVE_LLONG, memspace, in_fspace, H5P_DEFAULT,
                              index_buf_ptr);
                if (err < 0) CALLBACK_ERROR("H5Dread", index_levels[ii])

                index_buf_ptr += count[0];
            }
            GET_TIMER(ts, te, timing[0])

            /* table look up the seq values */
            for (ii=0; ii<count[0]; ii++) {
                keyv key(ndsets);
                for (jj=0; jj<ndsets; jj++)
                    key[jj] = index_buf[count[0]*jj + ii];
                seq_buf[ii] = lookup_table[key];
                mp[seq_buf[ii]]++;
            }
            GET_TIMER(ts, te, timing[1])

            if (create_seq) {
                /* write a chunk to seq dataset */
                err = H5Sselect_hyperslab(out_fspace, H5S_SELECT_SET, start, NULL, ones, count);
                if (err < 0) CALLBACK_ERROR("H5Sselect_hyperslab", seq_name)

                err = H5Dwrite(seq_id, H5T_STD_I64LE, memspace, out_fspace, H5P_DEFAULT, seq_buf);
                if (err < 0) CALLBACK_ERROR("H5Dwrite", seq_name)
            }
            off    += count[0];
            remain -= count[0];
            GET_TIMER(ts, te, timing[2])
        }
        free(index_buf);

        /* close the index dataset */
        for (ii=0; ii<ndsets; ii++) {
            if ((err = H5Dclose(dset_ids[ii])) < 0)
                CALLBACK_ERROR("H5Dclose", index_levels[ii].c_str());
        }
        free(dset_ids);
    }
    else {
        /* index datasets are stored as a single dataset, for example
         * dataset /event_table/event_id which is a 2D array whose 2nd
         * dimension is 3, each row representing a tuple (run, subrun, evt)
         */
        long long **index_buf;
        size_t remain=dset_dims[0];
        hsize_t off=0, memlen;

        /* allocate read buffer */
        index_buf = (long long**) malloc(chunk_size * sizeof(long long*));
        index_buf[0] = (long long*) malloc(chunk_size * dset_dims[1] * sizeof(long long));
        for (ii=1; ii<chunk_size; ii++) index_buf[ii] = index_buf[ii-1] + dset_dims[1];

        if (verbose) {
            printf("Read index dataset /%s/%s\n",grp_name,dset_name);
            printf("     index dataset of shape %lu x %lu\n",dset_dims[0], dset_dims[1]);
        }

        /* loop over the index dataset until all data has been read */
        while (remain > 0) {
            /* set file space for the read chunk */
            start[0] = off;
            start[1] = 0;
            count[0] = MIN(chunk_size, remain);
            count[1] = dset_dims[1];
            err = H5Sselect_hyperslab(in_fspace, H5S_SELECT_SET, start, NULL, ones, count);

            /* set memory space, a contiguous buffer in memory */
            memlen = count[0] * count[1];
            memspace = H5Screate_simple(1, &memlen, NULL);

            /* read a chunk from each index datasets */
            err = H5Dread(dset, H5T_NATIVE_LLONG, memspace, in_fspace, H5P_DEFAULT, index_buf[0]);
            if (err < 0) RETURN_ERROR("H5Dread", dset_name)
            err = H5Sclose(memspace);
            if (err < 0) CALLBACK_ERROR("H5Sclose", dset_name)
            GET_TIMER(ts, te, timing[0])

            /* table look up the seq values */
            for (ii=0; ii<count[0]; ii++) {
                keyv key(dset_dims[1]);
                for (jj=0; jj<dset_dims[1]; jj++)
                    key[jj] = index_buf[ii][jj];
                seq_buf[ii] = lookup_table[key];
                mp[seq_buf[ii]]++;
            }
            GET_TIMER(ts, te, timing[1])

            if (create_seq) {
                /* write a chunk to seq dataset */
                count[1] = 1;
                err = H5Sselect_hyperslab(out_fspace, H5S_SELECT_SET, start, NULL, ones, count);
                memspace = H5Screate_simple(1, count, NULL);
                err = H5Dwrite(seq_id, H5T_STD_I64LE, memspace, out_fspace, H5P_DEFAULT, seq_buf);
                if (err < 0) CALLBACK_ERROR("H5Dwrite", seq_name)
                err = H5Sclose(memspace);
                if (err < 0) CALLBACK_ERROR("H5Sclose", seq_name)
            }
            off    += count[0];
            remain -= count[0];
            GET_TIMER(ts, te, timing[2])
        }

        free(index_buf[0]);
        free(index_buf);

        /* close the index dataset */
        if ((err = H5Dclose(dset)) < 0) RETURN_ERROR("H5Dclose", dset_name)
    }
    free(seq_buf);

    err = H5Sclose(in_fspace);
    if (err < 0) CALLBACK_ERROR("H5Sclose", dset_name)
    if (create_seq) {
        err = H5Sclose(out_fspace);
        if (err < 0) CALLBACK_ERROR("H5Sclose", seq_name)
        /* close dataset */
        err = H5Dclose(seq_id);
        if (err < 0) CALLBACK_ERROR("H5Dclose", seq_name)
    }

seq_cnt_create:
    if (create_seq_cnt) {
        /* create the new partition key sequence-count dataset */
        size_t nkeys = mp.size();
        int64_t *seq_cnt = (int64_t*) malloc(nkeys * 2 * sizeof(int64_t));
        // Traverse through map and print frequencies
        ii = 0;
        for (auto x : mp) {
            seq_cnt[ii++] = x.first;
            seq_cnt[ii++] = x.second;
        }

        /* name of new partition key dataset */
        sprintf(seq_cnt_name, "%s.seq_cnt", part_key_base);

        /* check if the dataset already exists */
        isExist = H5Lexists(grp_id, seq_cnt_name, H5P_DEFAULT);
        if (isExist >= 0) {
            if (overwrite) {
                if (verbose) printf("Delete dataset %s/%s\n",grp_name,seq_cnt_name);
                err = H5Ldelete(grp_id, seq_cnt_name, H5P_DEFAULT);
                if (err < 0) CALLBACK_ERROR("H5Ldelete", seq_cnt_name)
            }
            else {
                cerr<<"Error: dataset " << grp_name << "/" << seq_cnt_name << " already exists"<<endl;
                return -1;
            }
        }

        /* seq dataset will be of the same dimension lengths */
        hsize_t maxdims[2] = {H5S_UNLIMITED, 2};
        hsize_t seq_cnt_dims[2] = {nkeys, 2};
        out_fspace = H5Screate_simple(2, seq_cnt_dims, maxdims);
        if (err < 0) CALLBACK_ERROR("H5Screate_simple", seq_cnt_name)

        /* create the new partition key dataset */
        if (verbose) printf("Create dataset %s/%s\n",grp_name,seq_cnt_name);
        seq_cnt_id = H5Dcreate2(grp_id, seq_cnt_name, H5T_STD_I64LE, out_fspace,
                                H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        if (seq_cnt_id < 0) CALLBACK_ERROR("H5Dcreate2", seq_cnt_name)

        if (nkeys > 0) {
            err = H5Dwrite(seq_cnt_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, seq_cnt);
            if (err < 0) CALLBACK_ERROR("H5Dwrite", seq_cnt_name)
        }

        free(seq_cnt);

        err = H5Sclose(out_fspace);
        if (err < 0) CALLBACK_ERROR("H5Sclose", seq_name)
        /* close dataset */
        err = H5Dclose(seq_cnt_id);
        if (err < 0) CALLBACK_ERROR("H5Dclose", seq_cnt_name)
    }

fn_exit:
    /* close create property */
    err = H5Pclose(dcpl_id);
    if (err < 0) RETURN_ERROR("H5Pclose", dset_name)

    GET_TIMER(ts, te, timing[2])

    /* close group */
    if (grp_id >= 0) {
        err = H5Gclose(grp_id);
        if (err < 0) CALLBACK_ERROR("H5Gclose", grp_name)
    }

    GET_TIMER(ts, te, timing[3])

    return 0;
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
    if (err < 0) HANDLE_ERROR("H5Pget_cache")

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
    if (err < 0) HANDLE_ERROR("H5Pset_cache")
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

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  [-h]            print this command usage message\n\
  [-v]            verbose mode (default: off)\n\
  [-f]            overwrite key datasets if already exist\n\
  [-r pattern]    groups matching pattern are not injected with key dataset\n\
  [-k indx_names] dataset names separated by comma, to be used to generate\n\
                  partition keys (default: /spill/run,/spill/subrun,/spill/evt)\n\
  [-c]            create sequence-count datasets as partition keys, instead\n\
                  of sequence-only datasets. The suffix of the key datasets\n\
                  will be `seq_cnt`. (default: off)\n\
  [-a]            create both sequence and sequence-count datasets.\n\
                  (default: off)\n\
  file_name       input/output HDF5 file name (required)\n\n\
  This utility program adds a new dataset in each group of the input file.\n\
  The new dataset, referred as the partition key dataset and to be named as\n\
  'last_name.seq', where 'last_name' is the name of last dataset provided in\n\
  the argument 'indx_names' of command-line option '-k'. The partition key\n\
  dataset is to be used for data partitioning purpose in parallel read\n\
  operations of the HDF5 file. Its contents are generated based on those\n\
  index datasets whose names are provided in the argument of command-line\n\
  option '-k'. The default is '/spill/run,/spill/subrun,/spill/evt' if option\n\
  '-k' is not used. These datasets together provide a unique identifiers,\n\
  which will be used to generate an array of integers to be stored in the\n\
  file in an increasing order. When reading the HDF5 file with such key in\n\
  parallel, the data partitioning strategy can assign the dataset elements\n\
  with the same 3-tuple of (run, subrun, evt) to the same MPI process. Thus,\n\
  the partition key dataset, named 'evt.seq' in this example, created in each\n\
  group stores a list of unique IDs corresponding to the unique 3-tuples. The\n\
  values in dataset 'evt.seq' are consistent across all groups. If the index\n\
  datasets are stored as a 2D array, for example '/event_table/event_id'\n\
  whose 2nd dimension is of size 3 storing datasets run, subrun, and evt,\n\
  then the command-line option can be simply '-k /event_table/event_id'.\n\
  When option '-c' is used, the partition key datasets will be created as 2D\n\
  N x 2 arrays, where N is the number of unique partition key values and the\n\
  two elements in each row are the key value and the count of its repeat in\n\
  the index datasets. This option is expected to create the key datasets of\n\
  much smaller sizes.  Requirements for the input HDF5 file:\n\
    1. the group must be the same in option '-k'. If '-k' option is not used,\n\
       the default datasets '/spill/run,/spill/subrun,/spill/evt' must exist.\n\
    2. contains multiple groups at root level\n\
    3. each group may contain multiple 2D datasets\n\
    4. all datasets in the same group share the 1st dimension size\n\
    5. other groups may not contain datasets provided in option '-k'. For\n\
       those groups, adding the key partition datasets is skip.\n\
    6. second dimension size of datasets provided in option '-k' must be 1\n\
    7. datasets provided in option '-k' will be read and type-cast into\n\
       internal buffers of type 'long long int' before sorting is applied.\n\
       Users are warned for possible data type overflow, if there is any.\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    printf("Usage: %s [-h|-v|-f|-r pattern] -k indx_names file_name\n%s\n", progname, USAGE);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int c, err_exit=0, ndims, create_seq=1, create_seq_cnt=0, overwrite=0;
    char *infile=NULL, dset_name[1024], part_key_base[1024];
    char *token, *dset_list=NULL, *pattern=NULL, *grp_name=NULL, str_tmp[1024];
    size_t ii, jj, num_orig_groups=0, num_nonzero_groups=0, ndsets;
    size_t max_seq=0, min_seq=0, avg_seq=0;
    long long *buf=NULL, *buf_ptr;
    double timing[6], stime[6];
    herr_t err;
    hid_t fid=-1, fapl_id=-1, dset, space_id;
    hsize_t dset_dims[2], maxdims[2];
    H5G_info_t grp_info;
    table lookup_table;
    regex_t regex;
    vector<string> index_levels;

    struct op_data it_op;
    it_op.num_groups = 0;
    it_op.grp_names  = NULL;
    it_op.regex      = NULL;

#if defined PROFILE && PROFILE
    double ts=0.0, te=0.0;
#endif
    for (c=0; c<6; c++) timing[c] = stime[c] = 0.0;

    verbose = 0; /* default is quiet */

    /* command-line arguments */
    while ((c = getopt(argc, argv, "hvfacr:k:")) != -1)
        switch(c) {
            case 'v': verbose = 1;
                      break;
            case 'f': overwrite = 1;
                      break;
            case 'k': dset_list = strdup(optarg);
                      break;
            case 'r': pattern = strdup(optarg);
                      break;
            case 'c': create_seq_cnt = 1;
                      create_seq = 0;
                      break;
            case 'a': create_seq_cnt = 1;
                      create_seq = 1;
                      break;
            case 'h': usage(argv[0]);
                      return 0;
            default:  usage(argv[0]);
                      return 1;
        }

    if (argv[optind] == NULL) { /* input file name is mandatory */
        printf("input file is missing\n");
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }

    if (pattern && regcomp(&regex, pattern, 0)) { /* check valid regex */
        printf("Error: input pattern %s does not compile", pattern);
        usage(argv[0]);
        err_exit = -1;
        goto fn_exit;
    }
    else if(pattern) { /* valid regex */
        it_op.regex = &regex;
    }

    infile = strdup(argv[optind]);
    if (verbose) printf("input file: %s\n", infile);

    /* construct index datasets into individual strings */
    if (dset_list == NULL) /* default datasets */
        dset_list = strdup("/spill/run,/spill/subrun,/spill/evt");

    /* tokenize dest_list */
    token = strtok(dset_list, ",");
    while (token != NULL) {
        index_levels.push_back(token);
        token = strtok(NULL, ",");
    }

    /* number of index datasets */
    ndsets = index_levels.size();
    assert(ndsets > 0);
    if (verbose) cout << "number of index datasets = " << ndsets << endl;

    strcpy(str_tmp, index_levels[0].c_str());
    grp_name = strdup(dirname(str_tmp));
    /* strip the group name */
    strcpy(str_tmp, index_levels[0].c_str());
    index_levels[0].assign(basename(str_tmp));
    for (ii=1; ii<ndsets; ii++) {
        strcpy(str_tmp, index_levels[ii].c_str());
        if (strcmp(grp_name, dirname(str_tmp))) {
            printf("Error: group names are inconsistent\n");
            exit(1);
        }
        /* strip the group name */
        strcpy(str_tmp, index_levels[ii].c_str());
        index_levels[ii].assign(basename(str_tmp));
    }

    /* use the last dataset name, as base for key */
    strcpy(part_key_base, index_levels[ndsets-1].c_str());

    SET_TIMER(ts)

    /* create file access property for read and write */
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id < 0) HANDLE_ERROR("H5Pcreate")

    /* increase cache size for raw data */
    err = incr_raw_data_chunk_cache(fapl_id);
    if (err < 0) HANDLE_ERROR("incr_raw_data_chunk_cache")

    /* open input file */
    fid = H5Fopen(infile, H5F_ACC_RDWR, fapl_id);
    if (fid < 0) HANDLE_ERROR("H5Fopen")

    GET_TIMER(ts, te, timing[0])
#if defined PROFILE && PROFILE
    printf("Open input file = %7.2f sec\n",timing[0]);
#endif

    /* retrieve the number of groups */
    err = H5Gget_info_by_name(fid, "/", &grp_info, H5P_DEFAULT);
    if (err < 0) RETURN_ERROR("H5Gget_info_by_name", "root group")
    num_orig_groups = grp_info.nlinks;
    it_op.grp_names = (char**) malloc(num_orig_groups * sizeof(char*));

    /* Iterate all objects to collect names of all group */
#if defined HAS_H5OVISIT3 && HAS_H5OVISIT3
    err = H5Ovisit3(fid, H5_INDEX_NAME, H5_ITER_NATIVE, gather_grp_names,
                    &it_op, H5O_INFO_ALL);
    if (err < 0) HANDLE_ERROR("H5Ovisit3")
#else
    err = H5Ovisit(fid, H5_INDEX_NAME, H5_ITER_NATIVE, gather_grp_names,
                   &it_op);
    if (err < 0) HANDLE_ERROR("H5Ovisit")
#endif
    if (it_op.err < 0) HANDLE_ERROR("H5Ovisit")

    GET_TIMER(ts, te, timing[1])
#if defined PROFILE && PROFILE
    printf("Collect group names = %7.2f sec\n",timing[1]);
#endif

    if (ndsets > 1) {
        /* open the first index dataset in group /spill to collect the 1st
         * dimension size. Note all datasets in group /spill should have the
         * same 1st dimension size.
         */
        sprintf(dset_name, "%s/%s", grp_name,index_levels[0].c_str());
        dset = H5Dopen(fid, dset_name, H5P_DEFAULT);
        if (dset < 0) RETURN_ERROR("H5Dopen", dset_name)

        /* collect dimension sizes */
        space_id = H5Dget_space(dset);
        if (space_id < 0) RETURN_ERROR("H5Dget_space", dset_name)
        ndims = H5Sget_simple_extent_dims(space_id, dset_dims, maxdims);
        if (ndims < 0) RETURN_ERROR("H5Sget_simple_extent_dims", dset_name)
        err = H5Sclose(space_id);
        if (err < 0) RETURN_ERROR("H5Sclose", dset_name)

        /* allocate buffers for all index datasets */
        buf = (long long*) malloc(dset_dims[0] * ndsets * sizeof(long long));
        buf_ptr = buf;

        /* read the entire 1st index dataset */
        if (verbose) printf("Read index dataset %s\n", dset_name);
        err = H5Dread(dset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
        if (err < 0) RETURN_ERROR("H5Dread", dset_name)
        if ((err = H5Dclose(dset)) < 0) RETURN_ERROR("H5Dclose", dset_name)

        buf_ptr += dset_dims[0];

        /* read the remaining index datasets into buffer */
        for (ii=1; ii<ndsets; ii++) {
            sprintf(dset_name, "%s/%s", grp_name,index_levels[ii].c_str());
            if (verbose) printf("Read index dataset %s\n", dset_name);

            /* open index dataset */
            dset = H5Dopen(fid, dset_name, H5P_DEFAULT);
            if (dset < 0) RETURN_ERROR("H5Dopen", dset_name)

            /* read the entire dataset */
            err = H5Dread(dset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          buf_ptr);
            if (err < 0) RETURN_ERROR("H5Dread", dset_name)

            /* close dataset */
            if ((err = H5Dclose(dset))  < 0) RETURN_ERROR("H5Dclose", dset_name)

            buf_ptr += dset_dims[0];
        }
    }
    else { /* the index datasets are stored together in a 2D array */
        long long *tmp;

        /* open the dataset to collect the 1st dimension size. */
        sprintf(dset_name, "%s/%s", grp_name, index_levels[0].c_str());
        dset = H5Dopen(fid, dset_name, H5P_DEFAULT);
        if (dset < 0) RETURN_ERROR("H5Dopen", dset_name)

        /* collect dimension sizes */
        space_id = H5Dget_space(dset);
        if (space_id < 0) RETURN_ERROR("H5Dget_space", dset_name)
        ndims = H5Sget_simple_extent_dims(space_id, dset_dims, maxdims);
        if (ndims < 0) RETURN_ERROR("H5Sget_simple_extent_dims", dset_name)
        err = H5Sclose(space_id);
        if (err < 0) RETURN_ERROR("H5Sclose", dset_name)

        tmp = (long long*) malloc(dset_dims[0] * dset_dims[1] * sizeof(long long));

        /* read the entire index dataset */
        if (verbose) printf("Read index dataset [%lu][%lu] %s\n", dset_dims[0], dset_dims[1], dset_name);
        err = H5Dread(dset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
        if (err < 0) RETURN_ERROR("H5Dread", dset_name)
        if ((err = H5Dclose(dset)) < 0) RETURN_ERROR("H5Dclose", dset_name)

        buf = (long long*) malloc(dset_dims[0] * dset_dims[1] * sizeof(long long));
        buf_ptr = buf;

        /* copy into contiguous space */
        for (ii=0; ii<dset_dims[1]; ii++) {
            for (jj=0; jj<dset_dims[0]; jj++)
                buf_ptr[jj] = tmp[ii + jj * dset_dims[1]];

            buf_ptr += dset_dims[0];
        }
        free(tmp);

        ndsets = dset_dims[1];
    }

    GET_TIMER(ts, te, timing[2])
#if defined PROFILE && PROFILE
    printf("Read partition index datasets in group '%s' = %7.2f sec\n",grp_name,timing[2]);
#endif

    /* build a lookup table based on the index datasets */
    lookup_table = build_lookup_table(dset_dims[0], ndsets, buf);

    if (buf != NULL) free(buf);

    GET_TIMER(ts, te, timing[3])
#if defined PROFILE && PROFILE
    printf("Construct lookup table = %7.2f sec\n",timing[3]);
#endif

    /* Iterate all groups and create a new key dataset in each group */
    for (ii=0; ii<it_op.num_groups; ii++) {
        err = add_seq(fid, overwrite, index_levels, it_op.grp_names[ii],
                      part_key_base, lookup_table, create_seq, create_seq_cnt,
                      &num_nonzero_groups, &max_seq, &min_seq, &avg_seq, stime);
        if (err < 0) RETURN_ERROR("add_seq", it_op.grp_names[ii])
    }
    if (num_nonzero_groups > 0)
        avg_seq /= num_nonzero_groups;

    assert(num_orig_groups >= num_nonzero_groups);

    GET_TIMER(ts, te, timing[4])

fn_exit:
    if (it_op.grp_names != NULL) {
        for (ii=0; ii<num_orig_groups; ii++)
            free(it_op.grp_names[ii]);
        free(it_op.grp_names);
    }
    if (fid >= 0) {
        check_h5_objects(infile, fid);
        err = H5Fclose(fid);
        if (err < 0) printf("Error at line %d: H5Fclose\n",__LINE__);
    }
    if (fapl_id >= 0) {
        err = H5Pclose(fapl_id);
        if (err < 0) printf("Error at line %d: H5Pclose\n",__LINE__);
    }

    GET_TIMER(ts, te, timing[5])

#if defined PROFILE && PROFILE
    if (err_exit == 0) {
        printf("-------------------------------------------------------\n");
        printf("Input file name                   = %s\n", infile);
        printf("number of groups in the file      = %zd\n", num_orig_groups);
        printf("number of non-zero size groups    = %zd\n", num_nonzero_groups);
        printf("Partition index dataset group     = %s\n",grp_name);
        printf("Partition index dataset names     = %s",index_levels[0].c_str());
        for (ii=1; ii<index_levels.size(); ii++)
            printf(", %s", index_levels[ii].c_str());
        printf("\n");
        printf("Partition key dataset name        = %s.seq\n",part_key_base);
        printf("  max length among all groups     = %zd\n", max_seq);
        printf("  min length among all groups     = %zd\n", min_seq);
        printf("  avg length among all groups     = %zd\n", avg_seq);
        printf("-------------------------------------------------------\n");
        printf("Open input file                   = %7.2f sec\n", timing[0]);
        printf("Collect group names               = %7.2f sec\n", timing[1]);
        printf("Read 'spill' 3-tuples             = %7.2f sec\n", timing[2]);
        printf("Construct lookup table            = %7.2f sec\n", timing[3]);
        printf("Add partition key datasets        = %7.2f sec\n", timing[4]);
        printf("  Read partition index datasets   = %7.2f sec\n", stime[0]);
        printf("  Hash table lookup               = %7.2f sec\n", stime[1]);
        printf("  Write partition key datasets    = %7.2f sec\n", stime[2]);
        printf("  Close partition key datasets    = %7.2f sec\n", stime[3]);
        printf("File close                        = %7.2f sec\n", timing[5]);
        printf("-------------------------------------------------------\n");
        printf("End-to-end                        = %7.2f sec\n",
               timing[0] + timing[1] + timing[2] + timing[3] +
               timing[4] + timing[5]);
    }
#endif
    if (infile != NULL) free(infile);
    if (grp_name != NULL) free(grp_name);

    return (err_exit != 0);
}
