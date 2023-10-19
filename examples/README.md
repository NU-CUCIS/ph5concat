# Examples of running 'ph5concat'

This folder contains examples of running `ph5concat`, including a set of small
input files, run commands, concatenated output file, and metadata of input and
out files.

# Input files
* There are four input HDF5 files, each containing a small number of dataset
  extracted from neutrino simulation data for illustrative purpose.
  + sample_r11981_s06.gz
  + sample_r11981_s07.gz
  + sample_r11981_s08.gz
  + sample_r11981_s09.gz

* These input files are compressed. Run command 'make' will uncompress them
  into HDF5 files. They can also be uncompressed manually by commands below.
  ```console
  % gzip -dc sample_r11981_s06.gz > sample_r11981_s06.h5
  % gzip -dc sample_r11981_s07.gz > sample_r11981_s07.h5
  % gzip -dc sample_r11981_s08.gz > sample_r11981_s08.h5
  % gzip -dc sample_r11981_s09.gz > sample_r11981_s09.h5
  ```

* Metadata of input files
  * Metadata of individual files can be retrieved by running command 'h5ls -r'
    or 'h5dump -H'.
  * The metadata of files 'sample_r11981_s06.h5' and 'sample_r11981_s07.h5' is
    shown below.
    <p align="left">
    <img align="center" src="./s06_s07.tiff" width="1000">
    </p>
  * In these examples, each file contains 4 groups at the root level, namely
    'neutrino', 'rec.me.trkkalman', 'rec.training.cvnmaps', and 'spill'. The
    number of groups and their names must be identical among all input files to
    be concatenated.
  * In an input file, the number of datasets in a group can be different from
    another group.
  * Given a group, the number of datasets and their names it contains must be
    identical among all input files. However, the size of first dimension of
    datasets can be different.


* Run commands
  ```console
  % mpiexec -n 2 ../ph5_concat -i sample_list.txt -o sample_output.h5
  % mpiexec -n 4 ../ph5_concat -i sample_list.txt -o sample_output.h5 -k evt
  ```
  When completed, the output shown on screen is available in
  [sample_stdout.txt](./sample_stdout.txt).

* Concatenated output file
  + The concatenated output file is provided in `sample_output.h5.gz`.
  + The metadata retrieved from 'h5dump' command is also available in
    `sample_output.metadata`.
    ```console
    % gzip -dc sample_output.h5.gz > sample_output.h5
    % h5dump -Hp sample_output.h5
    ```
  * A short version of the metadata is shown below.
    <p align="left">
    <img align="center" src="./concated.tiff" width="400">
    </p>

## An example output from a run concatenating 128 files
Below is an example timing output from a larger run on Cori using 128 MPI
processes.
```console
  % srun -n 128 ./ph5_concat -i ./nd_list_128.txt -o /scratch1/FS_1M_128/nd_out.h5 -b 512 -k evt

  Number of input HDF5 files: 128
  Input directory name: /global/cscratch1/sd/wkliao/FS_1M_8
  Output file name: /global/cscratch1/sd/wkliao/FS_1M_128/nd_out.h5
  Output datasets are compressed with level 6
  Read metadata from input files takes 1.2776 seconds
  Create output file + datasets takes 25.7466 seconds
  Concatenating 1D datasets takes 158.8101 seconds
  Write partition key datasets takes 14.0372 seconds
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

