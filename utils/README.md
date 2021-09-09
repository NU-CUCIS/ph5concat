# Utility Programs

* [rechunk](#rechunk) -- adjusts chunk settings of all datasets in a given
  HDF5 file.
* [add_spill_index](#add_spill_index) -- adds a new dataset in group '/spill'
  of the input HDF5 file.
* [add_key](#add_key) -- adds partitioning key datasets to all groups
* [sort_file_list](#sort_file_list) -- sorts file names based on the values
  of a list of given datasets.
* [check_seq_incr](#check_seq_incr) -- checks contents of partitioning
  key datasets in all groups for whether their values are organized in a
  monotonically nondecreasing order.

---
## rechunk

**rechunk** is a utility program running in sequential. Given an input HDF5
file, it copies all data groups and datasets in the input file to a new file
where the datasets in the new file use an adjust chunk setting. The new chunk
settings are described below in the command usage. The input HDF5 file should
come from [NOvA experiments](https://www.hep.ucl.ac.uk/nova/)

Command usage:
  ```
  % ./rechunk -h
  Usage: ./rechunk [-h|-v|-d|-r|-s|-t] [-c size] [-C size] [-M size] [-z level] [-b size] [-o outfile] infile

    [-h]         print this command usage message
    [-v]         verbose mode (default: off)
    [-d]         disable in-memory I/O (default: enable)
    [-r]         disable chunk caching for raw data (default: enable)
    [-s]         re-define zero-sized datasets as scalars (default: no)
    [-t]         define true 1D dataset (default: no)
    [-c size]    chunk size along 1st dimension for true-1D dataset (default: 1048576)
    [-C size]    chunk size along 1st dimension for true-2D dataset (default: 128)
    [-M size]    chunk base size in MiB for true-1D datasets (overwrites -c option)
    [-z level]   GZIP compression level (default: 6)
    [-b size]    I/O buffer size in bytes (default: 1 GiB)
    [-o outfile] output file name (default: out.h5)
    infile       input file name (required and cannot be the same as output file)

    This utility program copies an input HDF5 file to a new file where the
    datasets in the new file uses an adjusted chunk setting.

    Requirements of the input HDF5 file:
      1. contains multiple groups only at root level
      2. each group contains multiple 2D datasets
      3. true-1D datasets are those whose 2nd dimension size is 1
      4. true-2D datasets are those whose 2nd dimension size is larger than 1

    Output HDF5 file:
      1. zero-sized datasets will be stored in H5D_COMPACT layout
      2. for non-zero sized datasets, chunking is only applied to 1st dimension
      3. for true-1D datasets, new chunk dimensions use value from -c option
      4. for true-2D datasets, new chunk dimensions use value from -C option
    *ph5concat version 1.1.0 of March 1, 2020.
  ```
A shell script to run `rechunk` on multiple files in batch is given in
`batch.sh`. Example run and output:
  ```
  % ./rechunk -M 1 in.h5

  In-memory I/O                        = enabled
  Chunk caching for raw data           = enabled
  I/O buffer size                      = 1073741824 bytes
  number of groups in the file         = 999
  total number of 1D datasets          = 15965
  total number of 2D datasets          = 8
  no. non-zero 1D datasets             = 2573
  no. non-zero 2D datasets             = 6
  no. datasets chunking changed        = 21
  0th  dim  chunk size for 2D datasets = 128
  1st  dim  chunk size for 1D datasets = 1048576
  MiB-based chunk size for 1D datasets = 1 MiB
  -------------------------------------------------------
  Input  file open   time              =    0.05 sec
  Output file create time              =    0.04 sec
  Re-chunk time                        =   10.35 sec
    Re-chunk dataset open   time       =    1.35 sec
    Re-chunk dataset create time       =    1.81 sec
    Re-chunk dataset read   time       =    0.86 sec
    Re-chunk dataset write  time       =    4.82 sec
    Re-chunk dataset close  time       =    0.23 sec
  Output file close time               =    0.18 sec
  -------------------------------------------------------
  Total time                           =   10.62 sec
  ```

---
## add_spill_index

**add_spill_index** is a utility program to be run in sequential. Given an HDF5
file from NOvA experiments, it adds a new dataset in group '/spill'. The
argument of command-line option '-s' requires the full path of a dataset. The
path should be in the form of '/group/dset'. The new dataset to be created will
be '/spill/dset' whose contents will be single-valued, populated from the first
element of '/group/dset'. This new dataset is intended to be used as an
additional index dataset for generating the partition key datasets.

Command usage:
  ```
  % ./add_spill_index -h
  Usage: add_spill_index [-h|-v|-n] -s src_path file_name
    [-h]          print this command usage message
    [-v]          verbose mode (default: off)
    [-n]          dry-run mode (default: off)
    -s src_path   full path of dataset whose first element's value will be used
                  to populate the new dataset in group '/spill' (required)
    file_name     input/output HDF5 file name (required)

    This utility program adds a new dataset in group '/spill' of the input
    file. Argument 'src_path' should be in the form of '/group/dset'. The new
    dataset to be created will be '/spill/dset', whose contents will be
    single-valued, populated from the first element of '/group/dset'. This new
    dataset is intended to be used as an additional index dataset for
    generating partition key datasets. Requirements for the HDF5 file:
      1. must contain group '/spill'
      2. contains multiple groups at root level
      3. each group may contain multiple 2D datasets
      4. all datasets in the same group share the 1st dimension size
      5. dataset src_path must exist, except for the one in group '/spill'
      6. the second dimension size of the dataset 'dset' must be 1
    *ph5concat version 1.1.0 of March 1, 2020.
  ```
Example run:
  ```
  % ./add_spill_index -s /rec.energy.nue/cycle nd_data_165_files.h5
  ```
---
## add_key

**add_key** is a utility program to be run in sequential. Given an HDF5 file
from NOvA experiments, it adds a new dataset in each group of the input file.
The new dataset, referred as the partition key dataset and to be named as
'last_name.seq', where 'last_name' is the name of last dataset provided in the
argument 'indx_names' of command-line option '-k'. The partition key dataset is
to be used for data partitioning purpose in parallel read operations of the
HDF5 file. Its contents are generated based on those datasets in group
'/spill', whose names are provided in the argument of command-line option '-k'.
These datasets together provide a list of unique identifiers, which can be used
to generate an array integers stored in an increasing order. An example is '-k
run,subrun,evt'. The data partitioning strategy for parallel reads is to assign
the dataset elements with the same 3-tuple of (run, subrun, evt) to the same
MPI process. Thus, the partition key dataset, named 'evt.seq' in this example,
created in each group stores a list of unique IDs corresponding to the unique
3-tuples. The values in dataset 'evt.seq' are consistent across all groups.

At first, the index datasets specified in option '-k' in group '/spill' are
used to construct a hash table (a C++ `unordered_map` object), where the
N-tuple, where N is the number of index datasets, is a hash key. The hash
values are the indices of N-tuples. Because each of all the index datasets
contains a list of unique integers which have already been sorted in a
monotonically nondecreasing order, the hash values are unique for unique
N-tuples.  The constructed hash table is then used to create the partition key
datasets in all other groups. Note the contents of key datasets are consistent
cross all groups. In other words, two dataset elements have the same key
dataset value if their N-tuples are the same.

Note the integer numbers in the partition key datasets in all groups are sorted
in a monotonically nondecreasing order. Values can be repeated. However, in
group '/spill', the values are also sorted, but starting from 0 with increment
of 1, with no repeated value.

Command usage:
  ```
  % ./add_key -h
  Usage: ./add_key [-h|-v|-r pattern] -k indx_names file_name
    [-h]            print this command usage message
    [-v]            verbose mode (default: off)
    [-r pattern]    groups matching pattern are not injected with key dataset
    [-k indx_names] dataset names separated by comma, to be used to generate
                    partition keys (default: /spill/run,/spill/subrun,/spill/evt)
    [-c]            create sequence-count datasets as partition keys, instead
                    of sequence-only datasets. The suffix of the key datasets
                    will be `seq_cnt`. (default: off)
    [-a]            create both sequence and sequence-count datasets.
                    (default: off)
    file_name       input/output HDF5 file name (required)

    This utility program adds a new dataset in each group of the input file.
    The new dataset, referred as the partition key dataset and to be named as
    'last_name.seq', where 'last_name' is the name of last dataset provided in
    the argument 'indx_names' of command-line option '-k'. The partition key
    dataset is to be used for data partitioning purpose in parallel read
    operations of the HDF5 file. Its contents are generated based on those
    index datasets whose names are provided in the argument of command-line
    option '-k'. The default is '/spill/run,/spill/subrun,/spill/evt' if option
    '-k' is not used. These datasets together provide a unique identifiers,
    which will be used to generate an array of integers to be stored in the
    file in an increasing order. When reading the HDF5 file with such key in
    parallel, the data partitioning strategy can assign the dataset elements
    with the same 3-tuple of (run, subrun, evt) to the same MPI process. Thus,
    the partition key dataset, named 'evt.seq' in this example, created in each
    group stores a list of unique IDs corresponding to the unique 3-tuples. The
    values in dataset 'evt.seq' are consistent across all groups. If the index
    datasets are stored as a 2D array, for example '/event_table/event_id'
    whose 2nd dimension is of size 3 storing datasets run, subrun, and evt,
    then the command-line option can be simply '-k /event_table/event_id'.
    When option '-c' is used, the partition key datasets will be created as 2D
    N x 2 arrays, where N is the number of unique partition key values and the
    two elements in each row are the key value and the count of its repeat in
    the index datasets. This option is expected to create the key datasets of
    much smaller sizes.  Requirements for the input HDF5 file:
      1. the group must be the same in option '-k'. If '-k' option is not used,
         the default datasets '/spill/run,/spill/subrun,/spill/evt' must exist.
      2. contains multiple groups at root level
      3. each group may contain multiple 2D datasets
      4. all datasets in the same group share the 1st dimension size
      5. other groups may not contain datasets provided in option '-k'. For
         those groups, adding the key partition datasets is skip.
      6. second dimension size of datasets provided in option '-k' must be 1
      7. datasets provided in option '-k' will be read and type-cast into
         internal buffers of type 'long long int' before sorting is applied.
         Users are warned for possible data type overflow, if there is any.
    *ph5concat version 1.1.0 of March 1, 2020.
  ```
Example run and output:
  ```
  % ./add_key -k /spill/run,/spill/subrun,/spill/cycle,/spill/evt nd_data_165_files.h5
  Input file name                   = nd_data_165_files.h5
  number of groups in the file      = 999
  number of non-zero size groups    = 108
  Partition index dataset group     = /spill
  Partition index dataset names     = run, subrun, cycle, evt
  Partition key dataset name        = evt.seq
    max length among all groups     = 1695404032
    min length among all groups     = 193325
    avg length among all groups     = 30787342
  -------------------------------------------------------
  Open input file                   =    0.00 sec
  Collect group names               =    0.68 sec
  Read 'spill' 3-tuples             =    0.01 sec
  Construct lookup table            =    0.00 sec
  Add partition key dataset         =  721.86 sec
    Read partition index datasets   =  111.33 sec
    Hash table lookup               =  610.00 sec
    Write partition key datasets    =    0.53 sec
    Close partition key datasets    =    1.05 sec
  File close                        =    0.08 sec
  -------------------------------------------------------
  End-to-end                        =  722.81 sec
  ```
---
## sort_file_list

**sort_file_list** is a utility program to be run in sequential. Given a list
of NOvA file names, it sorts the file paths based on the values of first
element of a list of index datasets provided in the argument of command-line
option '-k'. When option '-k' is not used, the default datasets used for
sorting are '/spill/run' and '/spill/subrun'. The sorted file paths are output
in a text file. The sorting order follows the appearing order of datasets in
the argument of option '-k'. Note only the first element of the datasets are
used in the sorting. The output file is intended to be used as an input to the
parallel dataset concatenation program `ph5_concat`, so that the concatenated
data is organized in an increasing order of a given list of datasets.

Command usage:
  ```
  % ./sort_file_list -h
  Usage: ./sort_file_list [-h|-v|-d|-k paths|-o outfile] infile
    [-h]          print this command usage message
    [-v]          verbose mode (default: off)
    [-d]          debug mode (default: off)
    [-k paths]    full paths to datasets, separated by comma, to be used for
                  multi-index sorting. Only the first elements of the datasets
                  are used in the sorting. (default: /spill/run,/spill/subrun)
    [-o outfile]  output file name (default: 'out_list.txt')
    infile        input file containing a list of HDF5 file paths (required)

    This utility program re-order the file paths given in file 'infile' into a
    sorted list, based on the increasing order of index datasets specified in
    the argument of command-line option '-k'. When option '-k' is not used, the
    default are datasets '/spill/run' and '/spill/subrun'. An example of its
    usage is '-k /spill/run,/spill/subrun,/rec.hdr/cycle'. The index datasets
    will be read and type-cast into internal buffers of type 'long long int',
    before the sorting is applied. Sorting follows the order of datasets
    appearing in the argument of option '-k'. Note the contents of outfile will
    be overwritten if it exists. Requirements for the input HDF5 files:
      1. must contain datasets '/spill/run' and '/spill/subrun' if command-line
         option '-k' is not used
      2. may contain multiple groups at root level
      3. each group may contain multiple 2D datasets
      4. 1st dimension of all datasets in the same group share the same size
      5. datasets specified in argument '-k' must exist
    *ph5concat version 1.1.0 of March 1, 2020.
  ```
Example run and output:
  ```
  % cat sample_in_list.txt
  ND/neardet_r00011981_s06_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ND/neardet_r00011981_s07_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ND/neardet_r00011991_s01_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ND/neardet_r00011988_s01_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ND/neardet_r00011982_s04_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ND/neardet_r00011981_s24_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5

  % ./sort_file_list -o sample_out_list.txt sample_in_list.txt

  % cat sample_out_list.txt
  ND1/neardet_r00011981_s06_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ND1/neardet_r00011981_s07_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ND1/neardet_r00011981_s24_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ND1/neardet_r00011982_s04_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ND1/neardet_r00011988_s01_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ND1/neardet_r00011991_s01_t00_R19-02-23-miniprod5.i_v1_data.h5caf.h5
  ```
---
## check_seq_incr

**check_seq_incr** is a utility program to be run in sequential. Given a NOvA
HDF5 file that contains partition key datasets, it checks the contents of
partition key datasets in all groups. It first checks the partition dataset
specified at the command-line option '-d', a 1D array of integers, whose values
should start from 0 and increment by 1. The contents of partition key datasets
in all other groups are only required to be in a monotonically nondecreasing
order, i.e. may not start from 0, may contain repeated values, and may not
increment by 1. The name of partition dataset is a command-line parameter.

Command usage:
  ```
  % ./check_seq_incr -h
  Usage: ./check_seq_incr [-h|-v|-d name] infile
    [-h]       print this command usage message
    [-v]       verbose mode (default: off)
    [-d name]  partition key dataset name (default: 'evt.seq')
    infile     name of input HDF5 file (required)

    This utility program checks the contents of dataset 'name' in each group of
    the input file for whether the values are in a monotonically nondecreasing
    order. In particular, the partition dataset specified at command-line option
    '-d' is checked for whether the values start from 0 and increment by 1.
    The same datasets in all other groups are checked only for a monotonically
    nondecreasing order. Requirements for the input HDF5 file:
      1. contains multiple groups at root level
      2. each group may contain multiple 2D datasets
      3. the 1st  dimension of all datasets in the same group share same size
      4. each group must contain a partition key dataset named 'name'
    *ph5concat version 1.1.0 of March 1, 2020.
  ```

Example run:
  ```
  % ./check_seq_incr -d /event_table/evt.seq out.h5
  All datasets pass the check
  ```


