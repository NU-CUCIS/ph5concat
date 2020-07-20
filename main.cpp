/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
 * See COPYRIGHT notice in top-level directory.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <cerrno>
#include <cstring> /* strerror() */
#include <string>
#include <vector>
#include <unistd.h>
#include <string.h> /* strdup() */
#include <libgen.h> /* dirname() */
#include <sys/time.h>
#include <time.h>
#include "ph5_concat.hpp"

#if defined PROFILE && PROFILE
    #define SET_TIMER(ts) { MPI_Barrier(MPI_COMM_WORLD); ts = MPI_Wtime(); }
    #define GET_TIMER(ts, t) { t = MPI_Wtime() - ts; }
    #define PRN_TIMER(t, msg) { \
        if (!opt.quiet && rank == 0) { \
            printf("%s takes %.4f seconds\n", msg, t); \
            fflush(stdout); \
        } \
    }
#else
    #define SET_TIMER(ts)
    #define GET_TIMER(ts, t)
    #define PRN_TIMER(t, msg) { \
        if (!opt.quiet && rank == 0) { \
            printf("%s ---- done\n", msg); \
            fflush(stdout); \
        } \
    }
#endif

class Options {
public:
    Options(int argc, char **argv);
    bool quiet;
    bool err_exit;
    bool one_process_create;
    bool posix_open;
    bool in_memory_io;
    bool chunk_caching;
    int  io_strategy;
    size_t compress_threshold;
    unsigned int zip_level;
    size_t buffer_size;
    vector<string> input_files;
    string input_dirname;
    string output_file;
    string part_key_base;
};

/*----< usage() >------------------------------------------------------------*/
static void
usage(char *progname)
{
#define USAGE   "\
  Concatenate multiple HDF5 files into an output file.\n\n\
  [-h]         print this command usage message\n\
  [-q]         enable quiet mode (default: disable)\n\
  [-d]         disable in-memory I/O (default: enable)\n\
  [-r]         disable chunk caching for raw data (default: enable)\n\
  [-s]         one process creates followed by all processes open file (default: off)\n\
  [-p]         use MPI-IO to open input files (default: POSIX)\n\
  [-t num]     use parallel I/O strategy 1 or 2 (default: 2)\n\
  [-m size]    disable compression for datasets of size smaller than 'size' MiB\n\
  [-k name]    name of dataset used to generate partitioning keys\n\
  [-z level]   GZIP compression level (default: 6)\n\
  [-b size]    I/O buffer size per process (default: 128 MiB)\n\
  [-o outfile] output file name (default: out.h5)\n\
  [-i infile]  input file containing HEP data files (default: list.txt)\n\n\
  *ph5concat version _PH5CONCAT_VERSION_ of _PH5CONCAT_RELEASE_DATE_\n"

    cout<<"Usage: "<<progname<<
    " [-h|-q|-d|-r|-s|-p] [-t num] [-m size] [-k name] [-z level] [-b size] [-o outfile] [-i infile]\n"
    << USAGE << endl;
}

Options::Options(int argc, char **argv) :
                 quiet(false),
                 err_exit(false),
                 one_process_create(false),
                 posix_open(true),
                 in_memory_io(true),
                 chunk_caching(true),
                 io_strategy(2),
                 compress_threshold(0),
                 zip_level(6),
                 buffer_size(128*1048576),
                 output_file("./out.h5"),
                 part_key_base("")
{
    int opt;
    char *in_filename = NULL;
    string line;
    ifstream fd;

    while ((opt = getopt(argc, argv, "hqspdrt:m:k:i:o:z:b:")) != -1) {
        switch (opt) {
            case 's':
                one_process_create = true;
                break;
            case 'p':
                posix_open = false;
                break;
            case 'd':
                in_memory_io = false;
                break;
            case 'r':
                chunk_caching = false;
                break;
            case 't':
                io_strategy = atoi(optarg);
                if (io_strategy != 1 && io_strategy != 2) {
                    printf("Error: supported I/O strategies are 1 or 2\n");
                    err_exit = true;
                }
                break;
            case 'm':
                compress_threshold = strtoul(optarg, NULL, 0);
                break;
            case 'i':
                in_filename = strdup(optarg);
                break;
            case 'o':
                output_file = string(optarg);
                break;
            case 'k':
                part_key_base = string(optarg);
                break;
            case 'z':
                zip_level = strtoul(optarg, NULL, 0);
                break;
            case 'b':
                buffer_size = strtoul(optarg, NULL, 0) * 1048576;
                break;
            case 'q':
                quiet = true;
                break;
            case 'h':
            default:
                usage(argv[0]);
                err_exit = true;
        }
    }
    if (in_filename == NULL)
        in_filename = strdup("./list.txt"); /* default input file name */

    /* open input file and catch error */
    try {
        fd.open(in_filename);
        if (!fd)
            throw ios_base::failure(strerror(errno));
    }
    catch (ifstream::failure& e) {
        cerr << "Error: opening file \""<<in_filename<<"\" ("
             << e.what() << ")" << endl;
        err_exit = true;
        free(in_filename);
        return;
    }

    /* read input file contents */
    while (getline(fd, line)) {
        if (line.length() == 0)
            continue; /* skip empty lines */
        if (line.at(0) == '#')
            continue; /* skip comment line (start with #) */
        input_files.push_back(line);
    }
    fd.close();

    char *pathcopy = strdup(input_files[0].c_str());
    input_dirname.assign(dirname(pathcopy));
    free(pathcopy);
    free(in_filename);
}

/*----< main() >-------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int err=0, nprocs, rank;
    vector<string> myinputs;
    size_t offset, length, remainder;
#if defined PROFILE && PROFILE
    double ts, step_time[8], max_time[10];
#endif

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* all processes read and parse the command-line options */
    Options opt(argc, argv);
    if (opt.err_exit) {
        MPI_Finalize();
        return 0;
    }

    if (!opt.quiet && rank == 0) {
        printf("Number of input HDF5 files: %zd\n",opt.input_files.size());
        printf("Input directory name: %s\n",opt.input_dirname.c_str());
        printf("Output file name: %s\n",opt.output_file.c_str());
        printf("Output datasets are compressed with level %u\n",opt.zip_level);
        fflush(stdout);
    }

    if (opt.input_files.size() < (size_t)nprocs) {
        cout<<"The number of input files should be larger than or equal to the number of processes."<<endl;
        MPI_Finalize();
        return 1;
    }

    /* Evenly assign the input files among all workers. */
    length = opt.input_files.size() / nprocs;
    remainder = opt.input_files.size() % nprocs;
    if (static_cast<unsigned int>(rank) < remainder)
        length++;

    offset = rank * (opt.input_files.size() / nprocs);
    offset += (static_cast<unsigned int>(rank) < remainder) ? rank : remainder;

    myinputs = vector<string>(opt.input_files.begin() + offset,
                              opt.input_files.begin() + offset + length);

#if defined DEBUG && DEBUG
    for (vector<string>::const_iterator it = myinputs.begin();
         it != myinputs.end(); it++)
        cout<<"R"<<rank<<" "<<it->c_str()<<endl;
    cout<<"R"<<rank<<" will work on "<<myinputs.size()<<" files."<<endl;
#endif

    Concatenator concat(nprocs, rank, MPI_COMM_WORLD, MPI_INFO_NULL,
                        myinputs.size(),  // number of assigned input files
                        opt.output_file,  // output file name
                        opt.posix_open,
                        opt.in_memory_io,
                        opt.chunk_caching,
                        opt.compress_threshold,
                        opt.one_process_create,
                        opt.zip_level,
                        opt.buffer_size,
                        opt.io_strategy,
                        opt.part_key_base);

    if (opt.zip_level > 0) {
        /* Check if gzip fileter is available */
        unsigned int filter_info;
        htri_t avail = H5Zfilter_avail(H5Z_FILTER_DEFLATE);
        if (avail < 0) {
            cout<<"Error failed when calling H5Zfilter_avail"<<endl;
            // TODO: in C++, we should catch exception.
            goto prog_exit;
        }

        err = H5Zget_filter_info(H5Z_FILTER_DEFLATE, &filter_info);
        if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
            !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED)) {
            cout<<"gzip filter not available for encoding and decoding!"<<endl;
            goto prog_exit;
        }
    }

    /* Collect all metadata from input files */
#if defined HAS_H5GET_ALLOC_STATS && HAS_H5GET_ALLOC_STATS
    size_t cur_bytes, md_malloc;
    H5get_alloc_stats(NULL, &cur_bytes, NULL, NULL, NULL, NULL, NULL);
#endif
    SET_TIMER(ts)
    err = concat.construct_metadata(myinputs);
    if(err != 0){
        cout<<"construct_metadata() failed."<<endl;
        goto prog_exit;
    }
    GET_TIMER(ts, step_time[0])
    PRN_TIMER(step_time[0], "Read metadata from input files")
#if defined HAS_H5GET_ALLOC_STATS && HAS_H5GET_ALLOC_STATS
    H5get_alloc_stats(NULL, &md_malloc, NULL, NULL, NULL, NULL, NULL);
    md_malloc -= cur_bytes;
#endif

    /* Create a new file and define all groups and datasets */
    SET_TIMER(ts)
    err = concat.file_create();
    if (err < 0) {
        cout<<"file_create() failed."<<endl;
        goto prog_exit;
    }
    GET_TIMER(ts, step_time[1])
    PRN_TIMER(step_time[1], "Create output file + datasets")

    if (opt.io_strategy == 1) {
        /* Concatenate 1D datasets first */
        SET_TIMER(ts)
        err = concat.concat_small_datasets(myinputs);
        if (err < 0) {
            cout<<"concat_small_datasets() failed."<<endl;
            goto prog_exit;
        }
        GET_TIMER(ts, step_time[2])
        PRN_TIMER(step_time[2], "Concatenating 1D datasets")

        if (opt.part_key_base.compare("") != 0) {
            SET_TIMER(ts)
            /* write the partition keys */
            err = concat.write_partition_key_dataset();
            if (err < 0) {
                cout<<"write_partition_key_dataset() failed."<<endl;
                goto prog_exit;
            }
            GET_TIMER(ts, step_time[3])
            PRN_TIMER(step_time[3], "Write partition key datasets")
        }

        err = concat.close_input_files();
        if (err < 0) {
            cout<<"close_input_files() failed."<<endl;
            goto prog_exit;
        }

        /* Concatenate 2D datasets */
        SET_TIMER(ts)
        err = concat.concat_large_datasets(opt.input_files);
        if (err < 0) {
            cout<<"concat_large_datasets() failed."<<endl;
            goto prog_exit;
        }
        GET_TIMER(ts, step_time[4])
        PRN_TIMER(step_time[4], "Concatenating 2D datasets")
    }
    else if (opt.io_strategy == 2) {
        /* Concatenate 1D datasets first */
        SET_TIMER(ts)
        err = concat.concat_datasets(false);
        if (err < 0) {
            cout<<"concat_datasets() failed."<<endl;
            goto prog_exit;
        }
        GET_TIMER(ts, step_time[2])
        PRN_TIMER(step_time[2], "Concatenating 1D datasets")

        if (opt.part_key_base.compare("") != 0) {
            SET_TIMER(ts)
            /* write the partition keys */
            err = concat.write_partition_key_dataset();
            if (err < 0) {
                cout<<"write_partition_key_dataset() failed."<<endl;
                goto prog_exit;
            }
            GET_TIMER(ts, step_time[3])
            PRN_TIMER(step_time[3], "Writ partition key datasets")
        }

        /* Concatenate 2D datasets */
        SET_TIMER(ts)
        err = concat.concat_datasets(true);
        if (err < 0) {
            cout<<"concat_datasets() failed."<<endl;
            goto prog_exit;
        }
        GET_TIMER(ts, step_time[4])
        PRN_TIMER(step_time[4], "Concatenating 2D datasets")
    }

    /* close all input files */
    SET_TIMER(ts)
    concat.close_input_files();
    GET_TIMER(ts, step_time[5])
    PRN_TIMER(step_time[5], "Close input files")

    /* close output file */
    SET_TIMER(ts)
    concat.close_output_file();
    GET_TIMER(ts, step_time[6])
    PRN_TIMER(step_time[6], "Close output files")

#if defined PROFILE && PROFILE
    /* calculate total time spent on each process */
    step_time[7] = 0.0;
    for (int i=0; i<7; i++) step_time[7] += step_time[i];

    /* find the max timings among all processes */
    MPI_Reduce(step_time, max_time, 8, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    for (int i=0; i<8; i++) step_time[i] = max_time[i];

    double local_time[10];
    local_time[0]  = concat.c_1d_2d;
    local_time[1]  = concat.o_1d;
    local_time[2]  = concat.r_1d;
    local_time[3]  = concat.w_1d;
    local_time[4]  = concat.o_2d;
    local_time[5]  = concat.r_2d;
    local_time[6]  = concat.w_2d;
    local_time[7]  = concat.o_f;
    local_time[8]  = concat.close_in_dsets;
    local_time[9]  = concat.close_out_dsets;
    MPI_Reduce(local_time, max_time, 10, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    concat.c_1d_2d         = max_time[0];
    concat.o_1d            = max_time[1];
    concat.r_1d            = max_time[2];
    concat.w_1d            = max_time[3];
    concat.o_2d            = max_time[4];
    concat.r_2d            = max_time[5];
    concat.w_2d            = max_time[6];
    concat.o_f             = max_time[7];
    concat.close_in_dsets  = max_time[8];
    concat.close_out_dsets = max_time[9];

    if (!opt.quiet && rank == 0) { /* only rank 0 reports timings */
        printf("-------------------------------------------------------------\n");
        printf("Input directory name:                    %s\n",opt.input_dirname.c_str());
        printf("Number of input HDF5 files:              %zd\n",opt.input_files.size());
        printf("Output HDF5 file name:                   %s\n",opt.output_file.c_str());
        printf("Parallel I/O strategy:                   %d\n", opt.io_strategy);
        printf("Use POSIX I/O to open file:              %s\n",opt.posix_open?"ON":"OFF");
        printf("POSIX In-memory I/O:                     %s\n",opt.in_memory_io?"ON":"OFF");
        printf("1-process-create-followed-by-all-open:   %s\n",opt.one_process_create?"ON":"OFF");
        printf("Chunk caching for raw data:              %s\n",opt.chunk_caching?"ON":"OFF");
        printf("GZIP level:                              %d\n",opt.zip_level);
        if (opt.compress_threshold > 0)
            printf("Disable compress for datasets of size < %4zd MiB\n",opt.compress_threshold);
        printf("Internal I/O buffer size:                %.1f MiB\n",(float)concat.inq_io_buffer_size()/1048576.0);
        if (opt.part_key_base.compare("") != 0) {
            printf("Dataset used to produce partition key:   %s\n",opt.part_key_base.c_str());
            printf("Name of partition key datasets:          %s.seq\n", opt.part_key_base.c_str());
        }
        printf("-------------------------------------------------------------\n");
        printf("Number of groups:                    %8zd\n", concat.inq_original_num_groups());
        printf("Number of non-zero-sized groups:     %8zd\n", concat.inq_num_groups());
        if (opt.part_key_base.compare("") != 0)
            printf("Number of groups have partition key: %8zd\n", concat.inq_num_groups_have_key());
        printf("Total number of datasets:            %8zd\n", concat.inq_original_total_num_datasets());
        printf("Total number of non-zero datasets:   %8zd\n", concat.inq_num_datasets());
        printf("-------------------------------------------------------------\n");
        printf("Number of MPI processes:             %8d\n", nprocs);
        printf("Number calls to MPI_Allreduce:       %8d\n", concat.num_allreduce);
        printf("Number calls to MPI_Exscan:          %8d\n", concat.num_exscan);
        printf("-------------------------------------------------------------\n");
        printf("H5Dcreate:                           %9.4f\n", concat.c_1d_2d);
        if (opt.io_strategy == 1) {
            printf("H5Dopen   for 1D datasets:           %9.4f\n", concat.o_1d);
            printf("H5Fopen   for 2D datasets:           %9.4f\n", concat.o_f);
            printf("H5Dopen   for 2D datasets:           %9.4f\n", concat.o_2d);
        }
        printf("H5Dread   for 1D datasets:           %9.4f\n", concat.r_1d);
        printf("H5Dwrite  for 1D datasets:           %9.4f\n", concat.w_1d);
        printf("H5Dread   for 2D datasets:           %9.4f\n", concat.r_2d);
        printf("H5Dwrite  for 2D datasets:           %9.4f\n", concat.w_2d);
        printf("H5Dclose  for  input datasets:       %9.4f\n", concat.close_in_dsets);
        printf("H5Dclose  for output datasets:       %9.4f\n", concat.close_out_dsets);
        printf("-------------------------------------------------------------\n");
        printf("Read metadata from input files:      %9.4f\n", step_time[0]);
        printf("Create output file + datasets:       %9.4f\n", step_time[1]);
        printf("Concatenate small datasets:          %9.4f\n", step_time[2]);
        if (opt.part_key_base.compare("") != 0)
            printf("Write to partition key datasets:     %9.4f\n", step_time[3]);
        printf("Concatenate large datasets:          %9.4f\n", step_time[4]);
        printf("Close  input files:                  %9.4f\n", step_time[5]);
        printf("Close output files:                  %9.4f\n", step_time[6]);
        printf("End-to-end:                          %9.4f\n", step_time[7]);
        printf("\n");
    }

#if defined HAS_H5GET_ALLOC_STATS && HAS_H5GET_ALLOC_STATS
    unsigned long long total_alloc_bytes, max_ull, min_ull, avg_ull;
    size_t curr_alloc_bytes, peak_alloc_bytes, max_block_size;
    size_t total_alloc_blocks_count, curr_alloc_blocks_count;
    size_t peak_alloc_blocks_count, zd[7], max_zd[7], min_zd[7], avg_zd[7];

    err = H5get_alloc_stats(&total_alloc_bytes, &curr_alloc_bytes,
			    &peak_alloc_bytes, &max_block_size,
			    &total_alloc_blocks_count,
                            &curr_alloc_blocks_count,
                            &peak_alloc_blocks_count);

    MPI_Reduce(&total_alloc_bytes, &min_ull, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_alloc_bytes, &max_ull, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_alloc_bytes, &avg_ull, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    min_ull /= 1048576;
    max_ull /= 1048576;
    avg_ull /= nprocs * 1048576;
    zd[0] = md_malloc;
    zd[1] = curr_alloc_bytes;
    zd[2] = peak_alloc_bytes;
    zd[3] = max_block_size;
    zd[4] = total_alloc_blocks_count;
    zd[5] = curr_alloc_blocks_count;
    zd[6] = peak_alloc_blocks_count;
    MPI_Reduce(zd, &min_zd, 7, MPI_UNSIGNED_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(zd, &max_zd, 7, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(zd, &avg_zd, 7, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    for (int i=0; i<4; i++) {
        min_zd[i] /= 1048576;
        max_zd[i] /= 1048576;
        avg_zd[i] /= nprocs * 1048576;
    }
    for (int i=4; i<7; i++) {
        min_zd[i] /= 1024;
        max_zd[i] /= 1024;
        avg_zd[i] /= nprocs * 1024;
    }

    if (!opt.quiet && rank == 0) { /* only rank 0 reports timings */
        printf("-------------------------------------------------------------\n");
        printf("Memory footprints (min, max, avg among all processes):\n");
        printf("construct metadata     (MiB) min=%8zd max=%8zd avg=%8zd\n",min_zd[0],max_zd[0],avg_zd[0]);
        printf("total_alloc_bytes      (MiB) min=%8llu max=%8llu avg=%8llu\n",min_ull,max_ull,avg_ull);
        // printf("curr_alloc_bytes       (MiB) min=%8zd max=%8zd avg=%8zd\n",min_zd[1],max_zd[1],avg_zd[1]);
        printf("peak_alloc_bytes       (MiB) min=%8zd max=%8zd avg=%8zd\n",min_zd[2],max_zd[2],avg_zd[2]);
        printf("max_block_size         (MiB) min=%8zd max=%8zd avg=%8zd\n",min_zd[3],max_zd[3],avg_zd[3]);
        printf("total_alloc_blocks_count (K) min=%8zd max=%8zd avg=%8zd\n",min_zd[4],max_zd[4],avg_zd[4]);
        // printf("curr_alloc_blocks_count  (K) min=%8zd max=%8zd avg=%8zd\n",min_zd[5],max_zd[5],avg_zd[5]);
        printf("peak_alloc_blocks_count  (K) min=%8zd max=%8zd avg=%8zd\n",min_zd[6],max_zd[6],avg_zd[6]);
        printf("-------------------------------------------------------------\n");
    }
#endif
#endif

prog_exit:
    MPI_Finalize();
    return (err != 0);
}
