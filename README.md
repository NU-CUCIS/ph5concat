# Parallel HDF5 Dataset Concatenation for High Energy Physics Data Analysis

This software package contains C++ programs for concatenating HDF5 datasets
across multiple files into a single file by appending individual datasets one
after another. In a typical neutrino particle collision experiment, the
detector collects data into files over a period of time. Each file is labeled
by IDs of 'run' and 'subrun'. Runs are divided into Subruns. A run can be, for
example, one day of data taking and a subrun can be about an hour of data
taking. Each file contains thousands of two-dimensional datasets, organized
into hundreds of group, containing data describing the properties of a given
particle type. To analyze the data, individual datasets are required to be
concatenated one after another across all files, preferably in an increasing
order of their run and subrun ID. As the data amount and number of files from a
given experiment can become very large, the performance scalability of the
parallel concatenater is important.

## Input HDF5 Files
* Each file contains data from a single subrun. All groups contain datasets
  named 'run', 'subrun', and 'evt'. The values in a 'run' dataset are the
  same, representing the ID of a run. Similarly for 'subrun', the values in
  a 'subrun' represent the ID of a subrun. However, the subrun IDs are unique
  among all input files.
* Each file contains multiple groups. The number of groups and group names are
  the same among all input files.
* Each group contains multiple datasets. The number of datasets in a group can
  be different from each other. The number of datasets, dataset names, and
  their memberships to the groups are the same among all input files.
* Datasets in the same group are 2D arrays sharing the same size of 1st
  (most significant) dimension. Their 2nd dimension size may be different.
* Datasets in different groups may be of different 1st dimension sizes.
* Some of the datasets are actually 1D arrays whose 2nd dimension is of size 1.
* Datasets can be of size zero, i.e. the 1st dimension being of size 0.
* All the files have the same "schema", i.e. same numbers of groups and
  datasets, and their names.
* The size of 1st dimension of a dataset in a group of an input file may be
  different from the one in the same group but a different files.
* The same datasets in the same group but in different input files share the
  size of the 2nd dimension.

## Compiler and Software Requirements
* A C++ compiler that support ISO C++0x standard or higher
* MPI C and C++ compilers
* An HDF5 library version 1.10.5 and later built with parallel I/O feature enabled

## Instructions to Build
0. If building from a git clone of this repository, then run command below
   first. Otherwise, if building from an official release, this step can be
   skipped.
   ```
   autoreconf -i
   ```
1. Run command 'configure', for example
   ```
   ./configure --with-mpi=$HOME/MPICH/3.3 \
               --with-hdf5=$HOME/HDF5/1.10.5 \
               CFLAGS="-O2 -DNDEBUG" \
               CXXFLAGS="-O2 -DNDEBUG" \
               LIBS="-ldl -lz" \
               --enable-profiling
   ```
   * Option '--enable-profiling' enables timing measurement for internal
     functions and to report timing breakdowns on standard output.
2. Run command "make" to create the executable file named "ph5_concat"

## Command to Run
* Command-line options are:
  ```
  mpiexec -n <np> ./ph5_concat [-h|-q|-d|-r|-s|-p] [-t num] [-m size] [-k base_name] [-z level] [-b size] [-o outfile] [-i infile]

  [-h]           print this command usage message
  [-q]           enable quiet mode (default: disable)
  [-d]           disable in-memory I/O (default: enable)
  [-r]           disable chunk caching for raw data (default: enable)
  [-s]           one process creates followed by all processes open file (default: off)
  [-p]           use MPI-IO to open input files (default: POSIX)
  [-t num]       use parallel I/O strategy 1 or 2 (default: 2)
  [-m size]      disable compression for datasets of size smaller than 'size' MiB
  [-k base_name] dataset name in group /spill to generate partitioning keys
  [-z level]     GZIP compression level (default: 6)
  [-b size]      I/O buffer size per process (default: 128 MiB)
  [-o outfile]   output file name (default: out.h5)
  [-i infile]    input file containing HEP data files (default: list.txt)
  ```
  + `<np>`: Number of MPI processes.
  + `partitioning keys`: when command-line option '-k' is used, a new dataset
    will be created in each group in the output file, which can be used for
    data partitioning in parallel read operations. The new dataset, referred
    as the partition key dataset, is named 'base_name.seq' where 'base_name' is
    the dataset name provided in the command-line option '-k'. Contents of the
    partition key dataset are generated based on the dataset 'base_name' in
    group '/spill'. This base dataset must contain a list of unique integer
    values, stored in an increasing order, not necessarily incremented by one.
    An example is the dataset '/spill/evt'. The data partitioning strategy for
    parallel reads is to assign the dataset elements with the same 3-tuple of
    'run', 'subrun', and the base dataset to the same MPI process. Thus the
    partition key dataset created in the output file stores a list of unique
    IDs corresponding to the unique 3-tuples. The unique IDs are consistent
    among datasets across all groups. When option '-k' is not used, the
    partition key dataset will not be created.
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

## Sample input and output files
* There are four sample input files provided in folder `examples`.
  + examples/sample_r11981_s06.h5
  + examples/sample_r11981_s07.h5
  + examples/sample_r11981_s08.h5
  + examples/sample_r11981_s09.h5
* Sample run commands
  ```
  mpiexec -n 2 ./ph5_concat -i examples/sample_list.txt -o sample_output.h5
  mpiexec -n 4 ./ph5_concat -i examples/sample_list.txt -o sample_output.h5 -k evt
  ```
  The output shown on screen is stored in `examples/sample_stdout.txt`.
* Sample output files
  + The output files from concatenating the 4 sample files are available in
    `examples/sample_output.h5` whose metadata dumped from command below is
    also available in `examples/sample_output.metadata`.
    ```
    h5dump -Hp sample_output.h5
    ```

## An example timing output from a run on Cori using 128 MPI processes.
  ```
  % srun -n 128 ./ph5_concat -i ./nd_list_128.txt -o /scratch1/FS_1M_128/nd_out.h5 -b 512 -k evt

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
  GZIP level:                              6
  Internal I/O buffer size:                512.0 MiB
  Dataset used to produce partition key:   evt
  Name of partition key datasets:          evt.seq
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
  Concatenate small datasets:           158.8102
  Write to partition key datasets:       14.0372
  Concatenate large datasets:           114.4520
  Close  input files:                     0.0124
  Close output files:                     0.4799
  End-to-end:                           314.8095
  ```

## Questions/Comments:
* Sunwoo Lee <slz839@eecs.northwestern.edu>
* Wei-keng Liao <wkliao@northwestern.edu>

## Project funding supports:
This material is based upon work supported by the U.S. Department of Energy,
Office of Science, Office of Advanced Scientific Computing Research, Scientific
Discovery through Advanced Computing ([SciDAC](https://www.scidac.gov)) program.
This work is a collaboration of [RAPIDS Institute](https://rapids.lbl.gov) and
[HEP Data Analytics on HPC](https://computing.fnal.gov/hep-on-hpc/).

