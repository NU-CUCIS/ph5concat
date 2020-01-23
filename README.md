# Parallel Data Concatenation for High Energy Physics Data Analysis

This software package contains C++ programs for concatenating multiple HDF5
files into a single one by appending individual datasets one after another.

## Input HDF5 Files
* Each file contains multiple groups, each representing a "relational database
  table".
* Each group contains multiple datasets, each representing a column of the
  database table.
* Datasets in the same group are 2D arrays sharing the same size of 1st
  dimension (most significant). The size of 2nd dimension may be different.
* Some of the datasets are actually 1D arrays whose 2nd dimension if of size 1.
* Datasets can be of size zero, i.e. either dimension is of size 0.
* All the files have the same "schema", i.e. same structure of groups and
  datasets.
* A dataset in an input file may be of different 1st dimension size from the
  one in other files, while the 2nd dimension should be of the same size
  across files.

## Software Requirements
* A C++ compiler that support ISO C++0x standard or higher
* MPI C and C++ compilers
* An HDF5 library version  1.10.5 and later built with parallel I/O feature enabled

## Instructions to Build
1. Run command 'configure', for example
   ```
   ./configure --with-mpi=$HOME/MPICH/3.3 \
               --with-hdf5=$HOME/HDF5/1.10.5 \
               CFLAGS="-O2 -DNDEBUG" \
               CXXFLAGS="-O2 -DNDEBUG" \
               --enable-profiling
   ```
   * Option '--enable-profiling' enables timing measurement for internal
     functions and to report timing breakdowns on standard output.
2. Run command "make" to create the executable file named "ph5_concat"

## Command to Run
* Run command and command-line options are:
  ```
  mpiexec -n <np> ./ph5_concat [-h|-q|-d|-r|-s|-p|-x] [-t num] [-m size] [-k name] [-z level] [-b size] [-o outfile] [-i infile]

  [-h]         print this command usage message
  [-q]         enable quiet mode (default: disable)
  [-d]         disable in-memory I/O (default: enable)
  [-r]         disable chunk caching for raw data (default: enable)
  [-s]         one process creates followed by all processes open file (default: off)
  [-p]         use MPI-IO to open input files (default: POSIX)
  [-x]         disable calls to H5Dset_extent (default: enable)
  [-t num]     use parallel I/O strategy 1 or 2 (default: 2)
  [-m size]    disable compression for datasets of size smaller than 'size' MiB
  [-k name]    name of dataset used to generate partitioning keys
  [-z level]   GZIP compression level (default: 6)
  [-b size]    I/O buffer size per process (default: 128 MiB)
  [-o outfile] output file name (default: out.h5)
  [-i infile]  input file containing HEP data files (default: list.txt)
  ```
  + `<np>`: Number of MPI processes.
  + `partitioning keys`: when command-line option '-k' is used, two new
    datasets will be created to be used for the purpose of partitioning the
    datasets for parallel read operations. The two key datasets will have file
    name extensions '.key.seq' and '.key.cnt' appended to the dataset name
    provided in the option '-k'. When option '-k' is not used, the two
    partitioning key datasets will not be created.
  + I/O buffer size: when command-line option '-b' is used with value 0, this
    is equivalent to set the size to unlimited, i.e. the concatenator will
    allocate a buffer large enough to write each dataset in a single call to
    H5Dwrite. In this case, users may encounter out-of-memory errors.
  + I/O strategies: two I/O strategies (1 and 2) are currently supported. Both
    strategies share the same method for reading and writing the 1D datasets.
    At first, input files are assigned disjointly and evenly among all
    processes and each process reads 1D datasets only from the assigned files.
    Once 1D datasets are read, they are written to the output files using
    collective I/O. The difference between staretgies 1 and 2 are for reading
    and writing the 2D datasets. For strategy 1, each process reads only the
    assigned file (i.e. no shared-file reads) and then all processes
    collectively write to the output file. For strategy 2, all 2D datasets are
    read by all processes collectively (i.e. shared-file reads) and then all
    processes collectively write to the output file.

## An example output shown on screen
  ```
  % srun -n 128 ./ph5_concat -i ./nd_list_128.txt -o /scratch1/FS_1M_128/nd_out.h5 -b 512 -k evt -x

  Number of input HDF5 files: 128
  Input directory name: /global/cscratch1/sd/wkliao/FS_1M_8
  Output file name: /global/cscratch1/sd/wkliao/FS_1M_128/nd_out.h5
  Output datasets are compressed with level 6
  Read metadata from input files takes 1.2776 seconds
  Create output file + datasets takes 25.7466 seconds
  Concatenating 1D datasets takes 158.8101 seconds
  Writ partition key datasets takes 14.0372 seconds
  Concatenating 2D datasets takes 114.4464 seconds
  Close input files takes 0.0037 seconds
  Close output files takes 0.4797 seconds
  -------------------------------------------------------------
  Input directory name:                    /scratch/FS_1M_8
  Number of input HDF5 files:              128
  Output HDF5 file name:                   /scratch1/FS_1M_128/nd_out.h5
  Parallel I/O strategy:                   2
  Use POSIX I/O to open file:              ON
  POSIX In-memory I/O:                     ON
  1-process-create-followed-by-all-open:   OFF
  Chunk caching for raw data:              ON
  Enable calls to H5Dset_extent:           OFF
  GZIP level:                              6
  Internal I/O buffer size:                512.0 MiB
  Dataset used to produce partition key:   evt
  Name of partition key datasets:          evt.key.seq, evt.key.cnt
  -------------------------------------------------------------
  Number of groups:                         999
  Number of non-zero-sized groups:          108
  Number of groups have partition key:      108
  Total number of datasets:               17971
  Total number of non-zero datasets:       2795
  -------------------------------------------------------------
  Number of MPI processes:                  128
  Number calls to MPI_Allreduce:              3
  Number calls to MPI_Exscan:                 2
  -------------------------------------------------------------
  H5Dcreate:                             25.6772
  H5Dread   for 1D datasets:              1.8583
  H5Dwrite  for 1D datasets:            170.2729
  H5Dread   for 2D datasets:             19.9304
  H5Dwrite  for 2D datasets:             93.9265
  H5Dclose  for  input datasets:          0.0737
  H5Dclose  for output datasets:          0.0516
  -------------------------------------------------------------
  Read metadata from input files:         1.2782
  Create output file + datasets:         25.7466
  Small datasets total:                 158.8102
  Write to partition key datasets:       14.0372
  Large datasets total:                 114.4520
  Close  input files total:               0.0124
  Close output files total:               0.4799
  End-to-end:                           314.8095
  ```
## Questions/Comments:
* Sunwoo Lee <slz839@eecs.northwestern.edu>
* Wei-keng Liao <wkliao@eecs.northwestern.edu>

## Project funding supports:
This material is based upon work supported by the U.S. Department of Energy,
Office of Science, Office of Advanced Scientific Computing Research, Scientific
Discovery through Advanced Computing (SciDAC) program.

