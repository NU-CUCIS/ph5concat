# Utility Programs

* [rechunk](#rechunk)
* [add_key](#add_key)
* [sort_file_list](#sort_file_list)

---
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
  a new dataset named `base_name.seq` where 'base_name' is provided through
  command-line option '-k'. This dataset is referred as partitioning key
  dataset, to be used to calculate data partitioning for parallel programs that
  read the file. Its contents are calculated based on three existing datasets
  in the same group: `run`, `subrun`, and `base_name`. The 3-tuple (run[i],
  subrun[i], base_name[i]) forms a unique ID. Data elements with the same ID
  are assigned to the same MPI processes when the file is read in parallel. The
  3-tuple (run, subrun, base) works similarly to the primary key in relational
  database. Generation of `base_name.seq` is described as followed.

  Datasets `run`, `subrun`, and `base_name` in group `spill` are used to
  construct a C++ `unordered_map` as a hash table, where (run[i], subrun[i],
  base_name[i]) is a hash key. The hash values are the indices of 3-tuples.
  Because the contents of dataset `base_name` contain a list of unique
  integers, stored in an increasing order, the hash values are unique for
  unique 3-tuples. The constructed hash table is then used to create
  `base_name.seq` in all other groups.  Note the values in `base_name.seq` are
  consistent among datasets cross all groups. In other words, two dataset
  elements have the same `base_name.seq` value if their 3-tuples are the same.

  Note the integer numbers in `base_name.seq` in all groups are in a
  monotonically nondecreasing order. Values can be repeated.

The command usage is shown below.
  ```
  % ./add_key -h
  Usage: ./add_key [-h|-v] -k base_name file_name
    [-h]          print this command usage message
    [-v]          verbose mode (default: off)
    [-n]          dry run without creating key datasets (default: disabled)
    -k base_name  dataset name in group /spill to generate partitioning keys (required)
    file_name     input/output HDF5 file name (required)

    This utility program adds a new dataset in each group of the input file.
    The new dataset, referred as the partition key dataset and to be named as
    'base_name.seq', can be used for data partitioning in parallel read
    operations. Its contents are generated based on the dataset 'base_name' in
    group '/spill'. This base dataset must contain a list of unique integer
    values, stored in an increasing order, not necessarily incremented by one.
    An example is the dataset '/spill/evt'. The data partitioning strategy for
    parallel reads is to assign the dataset elements with the same 3-tuple of
    'run', 'subrun', and the base dataset to the same MPI process. Thus the
    partition key dataset created in the output file stores a list of unique
    IDs corresponding to the unique 3-tuples. The unique IDs are consistent
    among datasets across all groups. Requirements for the HDF5 file:
      1. must contain datasets /spill/run and /spill/subrun
      2. contains multiple groups at root level
      3. each group may contain multiple 2D datasets
      4. all datasets in the same group share the 1st dimension size
      5. each group must contain datasets run, subrun, and 'base_name'
      6. the second dimension size of the 3 daatsets must be 1
      7. data type of the 3 datasets msut be H5T_STD_U32LE
    *ph5concat version 1.1.0 of March 1, 2020.
  ```
  Example run and output:
  ```
  % ./add_key -k evt nd_data_165_files.h5

  Dry-run mode                      = NO
  Input file name                   = nd_data_165_files.h5
  Partition key base dataset name   = evt
  number of groups in the file      = 999
  number of non-zero size groups    = 108
  Partition key name                = evt.seq
    max length among all groups     = 1695404032
    min length among all groups     = 193325
    avg length among all groups     = 30787342
  -------------------------------------------------------
  Open input file                   =    0.00 sec
  Collect group names               =    0.68 sec
  Read 'spill' 3-tuples             =    0.01 sec
  Construct lookup table            =    0.00 sec
  Add partition key dataset         =  721.86 sec
    Read run, subrun, base datasets =  111.33 sec
    Hash table lookup               =  610.00 sec
    Write partition key seq         =    0.53 sec
    Close partition key seq         =    1.05 sec
  File close                        =    0.08 sec
  -------------------------------------------------------
  End-to-end                        =  722.81 sec
  ```
---
## sort_file_list

* **sort_file_list** is a utility program to be run in sequential. It reads
  a text input file containing a list of HDF5 file names and creates a new
  output text file containing the same set of HDF5 file names organized in an
  increasing order based on the run IDs and subrun ID stored in the HDF5 files.
  The output file can be used as input to the parallel dataset concatenation
  program `ph5_concat`, so that the data in the concatenated file follows the
  increasing  order of run and subrun IDs.

The command usage is shown below.
  ```
  % ./sort_file_list -h
  Usage: sort_file_list [-h|-v|-o outfile] infile
    [-h]          print this command usage message
    [-v]          verbose mode (default: off)
    [-d]          debug mode (default: off)
    [-o outfile]  output file name (default: 'out_list.txt')
    infile        input file name contains a list of HDF5 file names (required)

    This utility program re-order the files in infile into a sorted list based
    on the increasing order of 'run' and 'subrun' IDs. Requirements for the input
    HDF5 files:
      1. must contain datasets /spill/run and /spill/subrun
      2. contains multiple groups at root level
      3. each group may contain multiple 2D datasets
      4. all datasets in the same group share the 1st dimension size
      5. each group must contain datasets run and subrun'
      6. data type of datasets run and subrun must be H5T_STD_U32LE
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
