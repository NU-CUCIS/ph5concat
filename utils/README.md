# Utility Programs

## rechunk

* **rechunk** is a utility program running in sequential for adjusting the
  setting of data chunking for an input HDF5 file produced from the HEP NOvA
  experiments. The command usage is shown below.
  ```
  % ./rechunk -h
  Usage: ./rechunk [-h|-v|-d|-r|-s|-t] [-c size] [-C size] [-z level] [-b size] [-o outfile] infile

    [-h]         print this command usage message
    [-v]         verbose mode (default: off)
    [-d]         disable in-memory I/O (default: enable)
    [-r]         disable chunk caching for raw data (default: enable)
    [-s]         re-define zero-sized datasets as scalars (default: no)\n\
    [-t]         define true 1D dataset (default: no)s\n\
    [-c size]    chunk size along 1st dimension for true-1D dataset (default: 1048576)
    [-C size]    chunk size along 1st dimension for true-2D dataset (default: 128)
    [-z level]   GZIP compression level (default: 6)
    [-b size]    I/O buffer size in bytes (default: 1 GiB)
    [-o outfile] output file name (default: out.h5)
    infile       input file name (required and cannot be the same as output file)

    This utility program adjusts chunk setting of an input HDF5 file and save
    to an output file.

    Requirements of the input HDF5 file:
      1. contains multiple groups only at root level
      2. each group contains multiple 2D datasets
      3. true-1D datasets are those whose 2nd dimension size is 1
      4. true-2D datasets are those whose 2nd dimension size is larger than 1

    Output HDF5 file:
      1. zero-sized datasets will be stored in HDF5_COMPACT layout
      2. for non-zero sized datasets, chunking is only applied to 1st dimension
  ```
  A shell script to run `rechunk` on multiple files in batch is given in
  `batch.sh`. Example run and output:
  ```
  % ./rechunk in.h5

  In-memory I/O                  = enabled
  Chunk caching for raw data     = enabled
  I/O buffer size                = 1073741824 bytes
  number of groups in the file   = 999
  total number of 1D datasets    = 17935
  total number of 2D datasets    = 8
  no. non-zero 1D datasets       = 2779
  no. non-zero 2D datasets       = 6
  no. datasets chunking changed  = 43
  -------------------------------------------------------
  Input  file open   time        =    0.12 sec
  Output file create time        =    0.10 sec
  Re-chunk time                  =   20.44 sec
    Re-chunk dataset open   time =    1.53 sec
    Re-chunk dataset create time =    2.06 sec
    Re-chunk dataset read   time =    2.12 sec
    Re-chunk dataset write  time =   13.18 sec
    Re-chunk dataset close  time =    0.10 sec
  Output file close time         =    0.43 sec
  -------------------------------------------------------
  Total time                     =   21.09 sec
  ```

---
## add_key

* **add_key** is a utility program to be run in sequential. It adds each group
  two datasets to be used to calculate data partitioning for parallel programs
  that read the file. These two datasets are referred as partitioning key
  datasets. Their contents are calculated based on another dataset in the group
  chosen by the user, which is referred as the 'base' dataset. Together with
  'run' and 'subrun' datasets, the base dataset can uniquely identify data
  entries in a group that need to be processed as an unit. Thus the 3-tuple
  (run, subrun, base) works similarly to the primary key in relational
  database. When used for data partitioning, all data entries sharing the same
  3-tuple are assigned to the same MPI process.

  The two key datasets to be added have file names '.key.seq' and '.key.cnt'
  appended to the base name. For example, if the base dataset name is 'evt',
  then the two key datasets are 'evt.key.seq' and 'evt.key.cnt'. Dataset
  'key.seq' contains integer numbers corresponding to the unique 3-tuples in
  a group, starting from value 0 and ending with the number of unique 3-tuples
  minus one. Essentially, the values are a sequence of numbers incremented by
  one for every two consecutive entries whose 3-tuples are different. Dataset
  'key.cnt' stores the counters of unique 3-tuples. Its ith element contains
  the number of entries that share the same 3-tuple (also the ith unique
  3-tuple value). The length of 'key.cnt' can be shorter than 'key.seq', while
  the length of 'key.seq' is the same as all other datasets in the group.

The command usage is shown below.
  ```
  % ./add_key -h
  Usage: add_key [-h|-v] -k dataset_name file_name
    [-h]             print this command usage message
    [-v]             verbose mode (default: off)
    -k dataset_name  name of dataset used to generate partitioning keys (required)
    file_name        HDF5 file name (required)

    This utility program adds partitioning key datasets, key.seq and key.cnt,
    to an HDF5 file.

    Requirements for the HDF5 file:
      1. contains multiple groups only at root level
      2. each group may  contain multiple 2D datasets
    Requirements for the partitioning base daatset:
      1. the second dimension size of base dataset must be 1
      2. if base dataset is missing in a group, the key datasets will not be
         generated for that group
      3. currently supports base dataset of type either 4-byte unsigned int or
         2-byte short int
  ```
  Example run and output:
  ```
  % ./add_key -k evt nd_data_165_files.h5

  Dry-run mode                      = NO
  Input file name                   = nd_data_165_files.h5
  Partition key base dataset name   = evt
  number of groups in the file      = 999
  number of groups contain key base = 999
  number of non-zero key bases      = 108
  Partition key evt.key.seq :
    max length among all groups     = 1695404032
    min length among all groups     = 193325
    avg length among all groups     = 30787342
  Partition key evt.key.cnt :
    max length among all groups     = 410679
    min length among all groups     = 152244
    avg length among all groups     = 367746
  -------------------------------------------------------
  Open input file                   =    0.00 sec
  Collect partition key base names  =    0.68 sec
  Add partition key datasets        =  721.86 sec
    Create zero-size key datasets   =    0.19 sec
    Read run, subrun, base datasets =  111.33 sec
    Generate key.seq and key.cnt    =   30.19 sec
    Create key.seq                  =    0.06 sec
    Write  key.seq                  =  505.88 sec
    Close  key.seq                  =   55.88 sec
    Create key.cnt                  =    0.22 sec
    Write  key.cnt                  =    0.36 sec
    Close  key.cnt                  =   17.78 sec
  File close                        =    0.08 sec
  -------------------------------------------------------
  End-to-end                        =  722.81 sec
  ```
