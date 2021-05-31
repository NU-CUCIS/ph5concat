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
order of their run and subrun IDs. As the data amount and number of files from
a given experiment can become very large, the performance scalability of such
parallel data concatenation is important.

## Input HDF5 Files
* Each file contains data from a single subrun. All groups contain at least 3
  datasets named 'run', 'subrun', and 'evt'. The values in a 'run' dataset are
  the same, representing the ID of a run. Similarly for 'subrun', the values in
  a 'subrun' represent the ID of a subrun. However, the subrun IDs are unique
  among all input files.
* Each file contains multiple groups. The number of groups and group names are
  the same among all input files.
* Each group contains multiple datasets. The number of datasets in a group can
  be different from each other. The number of datasets, dataset names, and
  their memberships to the groups are the same among all input files.
* All datasets are 2D arrays.
* Datasets in the same group share the same size of 1st (most significant)
  dimension. Their 2nd dimension sizes may be different.
* Datasets in different groups may be of different 1st dimension sizes.
* Some of the datasets are actually 1D arrays whose 2nd dimension is of size 1.
* Datasets can be of size zero where their 1st dimension is of size 0.
* All the files have the same "schema", i.e. same numbers of groups and
  datasets, and their names.
* The size of 1st dimension of a dataset in a group of an input file may be
  different from the one in the same group but a different file.
* The same datasets in the same group but in different input files share the
  size of the 2nd dimension.

## Output HDF5 File
* A single HDF5 output file will be created.
* The output file shares the same schema as the input files, i.e. the same
  numbers and names of groups datasets.
* The size of 1st dimension (most significant) of individual datasets are sum
  of the 1st dimension of the same dataset from all input files.
* HDF5 compression and data chunking settings can be customized by command-line
  options (see below.)

## Compiler and Software Requirements
* A C++ compiler that support ISO C++0x standard or higher
* MPI C and C++ compilers
* An HDF5 library version 1.10.5 and later built with parallel I/O feature
  enabled

## Instructions to Build
0. If building from a git clone of this repository, then the below command
   must run first. Otherwise, if building from an official release, this step
   can be skipped.
   ```
   git clone https://github.com/NU-CUCIS/ph5concat.git
   cd ph5concat
   autoreconf -i
   ```
1. Run command 'configure'. An example is given below.
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
  mpiexec -n <np> ./ph5_concat [-h|-q|-a|-d|-r|-s|-p] [-t num] [-m size] [-k base_name] [-z level] [-b size] [-o outfile] [-i infile]

  [-h]           print this command usage message
  [-q]           enable quiet mode (default: disable)
  [-a]           append concatenated data to an existing HDF5 file (default: no)
  [-d]           disable in-memory I/O (default: enable)
  [-r]           disable chunk caching for raw data (default: enable)
  [-s]           one process creates followed by all processes open file (default: off)
  [-p]           use MPI-IO to open input files (default: POSIX)
  [-t num]       use parallel I/O strategy 1 or 2 (default: 2)
  [-m size]      disable compression for datasets of size smaller than 'size' MiB
  [-k base_name] dataset name in group /spill to generate partitioning keys
  [-z level]     GZIP compression level (default: 6)
  [-c]           enforces the contiguous layout for all datasets (default: false)
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
    For 1D datasets, input files are first assigned disjointly and evenly among all
    processes. Each process reads each 1D dataset entirely from the assigned files
    and writes it to the output files using collective I/O. The difference between
    staretgies 1 and 2 are for the 2D datasets. In strategy 1, all processes open
    all input files collectively using MPI-IO and read all individual 2D datasets
    collectively (i.e. shared-file reads), followed by all processes collectively
    writing individual datasets to the output file. In this strategy, all reads
    and writes are collective for each dataset. In strategy 2, each process reads
    datasets only from the disjointly assigned file (i.e. no shared-file reads)
    and then all processes collectively write each of the datasets to the output
    file. In this strategy, reads are independent but writes are collective.
  + When using '-a' append mode, the output file must exist. The input files
    will be concatenated into the output file.
  + When using both options '-a' and '-k', the partitioning key datasets must
    have been created previously in the output file and their names must be the
    same as the name used in '-k' option.

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

## Development Team:
   * Sunwoo Lee <<slz839@eecs.northwestern.edu>>
   * Wei-keng Liao <<wkliao@northwestern.edu>>
   * Saba Sehrish <<ssehrish@fnal.gov>>
   * Marc Paterno <<paterno@fnal.gov>>
   * James Kowalkowski <<jbk@fnal.gov>>

## Questions/Comments:
   * Sunwoo Lee <<slz839@eecs.northwestern.edu>>
   * Wei-keng Liao <<wkliao@northwestern.edu>>

## Project funding supports:
This material is based upon work supported by the U.S. Department of Energy,
Office of Science, Office of Advanced Scientific Computing Research, Scientific
Discovery through Advanced Computing ([SciDAC](https://www.scidac.gov)) program.
This work is a collaboration of [RAPIDS Institute](https://rapids.lbl.gov) and
[HEP Data Analytics on HPC](https://computing.fnal.gov/hep-on-hpc/).

