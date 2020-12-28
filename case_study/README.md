# Case Studies

---
## PandAna Parallel Read

**pandana_read.c** is a MPI-based  C program developed to study the performance
of parallel read from NOvA file produced by **ph5_concat**. It uses the same
data partitioning pattern as [PandAna](https://bitbucket.org/mpaterno/pandana),
a Python package that can be used to select NOvA events of interest.

### Data Partitioning Pattern
PandAna's read data partitioning pattern divides the event IDs exclusively into
contiguous ranges of event IDs evenly among all processes. The implementation
includes the followings.
1. Find the total number of unique event IDs. This is essentially the length of
   dataset **'/spill/evt.seq'** whose contents are 0, 1, 2, 3, ...., N-1, if
   its length is N.
2. Calculate the partitioning boundaries, so each process is responsible for a
   exclusive and contiguous range of event IDs. If N is not divisible by P, the
   number of processes, then the remainder IDs are assigned to the processes of
   lower ranks.
3. Note this calculation does not require reading the contents of
   **'/spill/evt.seq'**, but only inquires the length of **'/spill/evt.seq'**.
   ```
   my_count = N / nprocs;
   my_start = my_count * rank;
   if (rank < N % nprocs) {
       my_start += rank;
       my_count++;
   }
   else {
       my_start += N % nprocs;
   }
   my_end = my_start + my_count - 1;
   ```
4. Each process is responsible for the range from 'my_start' to 'my_end' inclusively.

### Parallel reads
Each group, **G**, in the concatenated file contains dataset 'evt.seq' whose values 
correspond to '/spill/evt.seq' and can be used to find the data element ranges of
all other datasets in the same group, to be read by all processes.
   * Note all datasets in the same group share the number of rows, i.e. the
     size of first dimension.
   * The contents of **'/G/evt.seq'** are monotonically non-decreasing. It is
     possible to have repeated event IDs in consecutive elements.

Parallel reads consist of the following steps.
1. Read **'/G/evt.seq'** and calculate array index ranges for all processes. 
   This can be done in 3 options.
   * Option 1. Root process reads the entire '/G/evt.seq' and then broadcasts 
     to the remaining processes. All processes use the contents of
     '/G/evt.seq' to calculate their responsible index ranges. 
   * Option 2. All processes collectively read the whole '/G/evt.seq'. All
     processes calculate their own responsible index ranges.
   * Option 3. Only root process reads '/G/evt.seq'. Root calculates
     responsible index ranges for all processes, and calls MPI_Scatter to 
     scatter the boundaries of ranges (start and end) to all other processes.
2. Calculate the responsible index ranges by checking the contents of 
   **'/G/evt.seq'** to find the starting and ending indices that point to range 
   of event IDs fall into its responsible range.
   * Two binary searches should be used, one to search for starting index and
     the other for ending index. This avoid sequentially checking the array
     contents.
3. All processes read the requested datasets in group G collectively, using
   the starting and ending indices (hyperslab), one dataset at a time.
   * Note read ranges are not overlapping among all processes.


### Run Command usage:
  ```
  % ./pandana_read -h
  Usage: ./pandana_read [-h|-v] [-p number] [-s number] [-m number] [-l file_name] [-i file_name]
    [-h]           print this command usage message
    [-v]           verbose mode (default: off)
    [-d]           debug mode (default: off)
    [-p number]    performance profiling method (0, 1, or 2)
                   0: report file open, close, read timings (default)
                   1: report number of chunks read per process
                   2: report read times for individual datasets
    [-s number]    read method for evt.seq (0, 1, or 2)
                   0: root process reads evt.seq and broadcasts (default)
                   1: all processes read the entire evt.seq collectively
                   2: root process reads evt.seq and scatters boundaries
                   3: A single MPI collective read all evt.seq and scatters boundaries
    [-m number]    read method for other datasets (0 or 1)
                   0: use H5Dread (default)
                   1: use MPI_file_read_all one dataset at a time
                   2: use MPI_file_read_all to read all datasets in one group at a time
                   1: use MPI_file_read_all
    [-l file_name] name of file containing dataset names to be read
    [-i file_name] name of input HDF5 file
    *ph5concat version 1.1.0 of March 1, 2020.
  ```

### Example Run:
A sample input file named 'dset.txt' is provided in this folder which includes
names of a list of datasets to be read from the concatenated HDF5 file.
Example run and output:
  ```
  % mpiexec -n 4 ./pandana_read -p 1 -s 2 -m 0 -l dset.txt -i nd_165_files_with_evtseq.h5
  Number of MPI processes = 4
  Input dataset name file 'dset.txt'
  Input concatenated HDF5 file 'nd_165_files_with_evtseq.h5'
  Number of datasets to read = 123
  Number of groups = 15
  Maximum number of datasets among groups = 13
  Read evt.seq method: root process reads evt.seq and scatters boundaries
  Read datasets method: H5Dread
  Number of unique evt IDs (size of /spill/evt.seq) = 410679
  ----------------------------------------------------
  MAX and MIN among all 4 processes
  MAX open_time=0.00 read_seq_time=0.40 read_dset_t=1.02 close_time=0.00
  MIN open_time=0.00 read_seq_time=0.40 read_dset_t=1.02 close_time=0.00
  ----------------------------------------------------
  Read amount MAX=4.77 MiB MIN=0.39 MiB (per dataset, per process)
  Amount of evt.seq datasets  259.68 MiB = 0.25 GiB (compressed  11.59 MiB = 0.01 GiB)
  Amount of  other  datasets  924.93 MiB = 0.90 GiB (compressed 181.23 MiB = 0.18 GiB)
  Sum amount of all datasets 1184.61 MiB = 1.16 GiB (compressed 192.82 MiB = 0.19 GiB)
  total number of chunks in all 123 datasets (exclude /spill/evt.seq): 1228
  Aggregate number of chunks read by all processes: 1552
          averaged among processes: 388.00
          averaged among processes among datasets: 3.15
  Out of 1228 chunks, number of chunks read by two or more processes: 320
  Out of 1228 chunks, most shared chunk is read by number of processes: 3
  ----------------------------------------------------


  rank   0: number of chunks read=490 (max=4 min=1 avg=3.98 among 108 datasets, exclude evt.seq)
  rank   1: number of chunks read=392 (max=6 min=1 avg=3.19 among 108 datasets, exclude evt.seq)
  rank   2: number of chunks read=332 (max=5 min=2 avg=2.70 among 108 datasets, exclude evt.seq)
  rank   3: number of chunks read=338 (max=6 min=1 avg=2.75 among 108 datasets, exclude evt.seq)
  ```
---
