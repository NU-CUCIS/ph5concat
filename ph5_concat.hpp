/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
 * See COPYRIGHT notice in top-level directory.
 */

#ifndef _HPP_PH5_CONCAT
#define _HPP_PH5_CONCAT

#include <string>
#include <vector>
#include <unordered_map>
using namespace std;

#include <hdf5.h>
#include <mpi.h>

#define MAX_STR_LEN 64
#define PROFILE_HDF5_INTERNAL 0
#define MIN_DATASET_SIZE     128*1024*1024

#define RETURN_ERROR(func_name, obj_name) { \
    printf("Error in %s line %d: calling %s for object %s\n",__FILE__,__LINE__,func_name,obj_name); \
    return -1; \
}

#define HANDLE_ERROR(fname) { \
    cout<<"["<<__FILE__<<"]["<<__FUNCTION__<<"]["<<__LINE__<<"] "<<fname<<" failed."<<endl; \
    err_exit = -1; \
    goto fn_exit; \
}

#define HANDLE_DSET_ERROR(fname, dname) { \
    cout<<"["<<__FILE__<<"]["<<__FUNCTION__<<"()][line:"<<__LINE__<<"] "<<fname<<" for dataset \""<<dname<<"\" failed."<<endl; \
    err_exit = -1; \
    goto fn_exit; \
}

#if defined PROFILE && PROFILE
    #define SET_TIME(ts) { MPI_Barrier(comm); ts = MPI_Wtime(); }
    #define GET_TIME(ts, t) { t = MPI_Wtime() - ts; }
#else
    #define START_TIMER(ts)
    #define GET_TIME(ts, t)
#endif

/*
 * Metadata of each dataset.
 */
struct DSInfo_t {
    string        name;          /* name of dataset */
    hsize_t       local_dims[2]; /* dim sizes, aggregated from local files */
    hsize_t       global_dims[2];/* dim sizes, aggregated from all files */
    hsize_t       chunk_dims[2]; /* chunk dimension sizes */
    hid_t         type_id;       /* datatype ID */
    H5T_class_t   type_class;    /* datatype class */
    size_t        type_size;     /* size of data element */
    H5D_layout_t  layout;        /* COMPACT, CONTIGUOUS, CHUNK */
    hid_t         out_dset_id;   /* output HDF5 dataset ID */
    hsize_t       global_off;    /* start offset relative to global dataset */
    hsize_t       size_in_bytes; /* size of globally aggregated dataset from all files */
    hsize_t       chunk_size;    /* size in bytes of one chunk */
    size_t        num_writes;    /* number of rounds of collective writes */

    bool          is_key_base;   /* is this dataset partitioning key base ? */
    bool          is_key_seq;    /* is this dataset a partitioning key seq ? */
    hsize_t       seq_len;       /* number of unique IDs in key.seq dataset */
    long long    *seq_buf;       /* key.seq buffer */

    vector<hid_t>   in_dset_ids; /* dataset IDs of opened input files */
    vector<hsize_t> in_dim0;     /* 1st dim sizes in assigned input files */


    hsize_t cur_chunk_offset; /* used in concatenation phase: */
    hsize_t cur_offset;       /* used in concatenation phase: offsets for next round of collective write */
    vector<hsize_t> global_num_rows; /* TODO: [nprocs] boundaries of data partitioning */
};
typedef struct DSInfo_t DSInfo_t;

struct GrpInfo {
    string       name;        /* name of dataset */
    hid_t        id;          /* HDF5 group ID */
    size_t       num_dsets;   /* number of dataset objects in this group */
    size_t       shared_dim0; /* size of common 1st dimension (all datasets in
                               * the same group share the 1st dimension size)
                               */
    DSInfo_t    *dsets;       /* dataset objects in this group */
    DSInfo_t    *key_base;    /* point to dataset used to generate key */
    DSInfo_t    *seq_dset;    /* point to key.seq */
};
typedef struct GrpInfo GrpInfo;

/* lookup hash table */
typedef unordered_map<int, int64_t> table;

/*
 * Main class.
 */
class Concatenator {
public:
    Concatenator(int nprocs, int rank, MPI_Comm comm, MPI_Info info,
                 size_t num_input_files, string const& output,
                 bool posix_open, bool in_memory_io, bool chunk_caching,
                 size_t compress_threshold, bool one_process_create,
                 unsigned int zip_level, bool enforce_contiguous,
                 size_t buffer_size, int io_strategy, string const& part_key_base);
    ~Concatenator();
    int construct_metadata(vector<string> const &inputs);
    int file_create();

    /* File-based partitioning (all datasets are read independently, but
     * written collectively) */
    int concat_datasets(bool process_large_dsets);

    /* Dataset-based partitioning (small datasets are read independently, large
     * datasets are read collectively.  Writes are done collectively. */
    int concat_small_datasets(vector<string> const &inputs);
    int concat_large_datasets(vector<string> const &inputs);

    /* finalize and write partitioning key dataset to file */
    int write_partition_key_dataset();

    int close_output_file();
    int close_input_files();
    int collect_metadata(hid_t obj_id, char const *name, H5O_info_t const *obj_info);

    size_t inq_io_buffer_size() { return io_buffer_size; }
    size_t inq_num_groups()          { return num_groups; }
    size_t inq_original_num_groups() { return original_num_groups; }
    size_t inq_num_groups_have_key() { return num_groups_have_key; }
    size_t inq_num_datasets()                { return total_num_datasets; }
    size_t inq_original_total_num_datasets() { return original_total_num_datasets; }

    /* members */

    /* timers */
    double c_1d_2d;
    double o_1d;
    double o_2d;
    double r_1d;
    double r_2d;
    double w_1d;
    double w_2d;
    double o_f;
    double close_in_dsets;
    double close_out_dsets;
    int num_allreduce;   /* number of calls to MPI_Allreduce */
    int num_exscan;      /* number of calls to MPI_Exscan */

private:
    /* members */
    MPI_Comm comm;             /* MPI communicator */
    MPI_Info info;             /* MPI-IO hints */
    int      nprocs;           /* number of MPI processes */
    int      rank;             /* MPI rank of this process */
    size_t   num_input_files;
    size_t   io_buffer_size;   /* I/O buffer size */
    unsigned int zip;          /* GZIP compression level */
    bool     enforce_contiguous;    /* enforce contiguous storage layout for all datasets */
    bool     posix_open;    /* use POSIX I/O or MPI-IO to open input files */
    bool     in_memory_io;  /* enable HDF5 in-memory I/O to read input files */
    bool     chunk_caching; /* enable HDF5 caching for raw data chunks */
    size_t   compress_threshold; /* whether to compress small datasets */
    bool     one_process_create; /* enable one-process-create-then-all-open */
    bool     add_partition_key;  /* whether to create partition key */
    int      io_strategy;        /* 1 or 2 (parallel I/O strategy) */

    string output_file_name;
    string part_key_base; /* dataset used to create partition key */

    size_t   num_groups;          /* number of groups */
    size_t   original_num_groups; /* some groups may contain zero-sized data */
    size_t   num_groups_have_key; /* number of groups contain key base */
    GrpInfo *groups;              /* array of group objects */

    int      spill_grp_no;  /* spill's array index in group[] */
    table   *lookup_table;  /* [num_input_files] hash tables, one for each file,
                               served as lookup tables built based on the user
                               indicated partition key base dataset in group
                               /spill */

    size_t  total_num_datasets; /* total number of non-zero datasets */
    size_t  original_total_num_datasets;
    size_t  chunk_size_threshold;
    size_t  in_memory_cache_size;
    size_t  output_meta_cache_size;
    size_t  raw_chunk_cache_size;
    hsize_t max_local_size_in_bytes;
    unordered_map<string, hid_t> input_files;
    hid_t output_file_id;
    hid_t dxpl_id; /* HDF5 data transfer property list identifier */

    char *buffer;  /* I/O buffer */

    MPI_Request *async_reqs;
    MPI_Status *async_statuses;

    /* methods */

    /* accumulate dataset dimension sizes across all input files */
    void accumulate_dimensions();

    /* Compute number of rounds of collective writes for all the datasets based
     * on the given I/O buffer size.
     */
    void calculate_num_writes();

    /* calculate data chunk size for a dataset */
    void calculate_chunk_size(DSInfo_t &dset_info);

    /* create an HDF5 dataset object */
    int create_dataset(hid_t group_id, DSInfo_t &dset_info, bool toFill);

    /* Create partitioning key dataset and add to its group */
    int create_partition_key(GrpInfo &grp);

    /* Use partition base dataset /spill/key_base to generate the lookup table
     * which will be used to generate the partition key dataset in each group.
     */
    int generate_partition_key(GrpInfo &grp);

    /* Read function for all datasets in strategy 2.
     * concat_datasets() calls this function in file_partition.cpp. */
    int read_dataset2(DSInfo_t &dset_info, size_t file_no, hsize_t all_rows,
                      hsize_t round_off, hsize_t round_len, hsize_t mem_off);

    /* Write function for all datasets in strategy 2.
     * concat_datasets() calls this function in file_partition.cpp. */
    int write_dataset_2D(DSInfo_t &dset_info, hsize_t *offs, hsize_t *lens,
                         void *wbuf);

    int open_input_files(vector<string> files, bool collective_io);

    int numerology(DSInfo_t &dset_info, hsize_t *dset_size,
                   hsize_t *offsets, hsize_t *counts);

    /* Read function for 2D datasets in strategy 1.
     * concat_large_datasets() calls this function in dataset_partition.cpp. */
    int read_2d_dataset(hid_t dset_id, DSInfo_t &dset_info, hsize_t *counts,
                        hsize_t *offsets, hsize_t *dset_size);

    /* Write function for 2D datasets in strategy 1.
     * concat_large_datasets() calls this function in dataset_partition.cpp. */
    int write_2d_dataset(DSInfo_t &dset_info, hsize_t *counts,
                         hsize_t *offsets, hsize_t *dset_size);

    /* Read function for 1D datasets in strategy 1.
     * concat_small_datasets() calls this function in dataset_partition.cpp. */
    int read_dataset(int file_no, DSInfo_t &dset_info);

    /* Write function for 1D datasets in strategy 1.
     * concat_small_datasets() calls this function in dataset_partition.cpp. */
    int write_dataset(DSInfo_t &dset_info, hsize_t *offs, hsize_t *lens);

};
#endif
