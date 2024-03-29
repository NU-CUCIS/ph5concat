# ph5concat Release Notes

---
## Version 1.1.0 (Sep. 13, 2023)

* New features
  + Now support datasets of type H5T_STRING
  + New command-line option '-c' to use the HDF5 contiguous storage layout.

* Utility programs
  + New utility program check_seq_incr
  + add_key
    * option `-k` to take dataset names separated by comma, to be used to
      generate partition keys (default: /spill/run,/spill/subrun,/spill/evt)
    * add option `-r` to specify groups matching pattern are not injected with
      key dataset
    * add option `-c` to create sequence-count datasets as partition keys,
      instead of sequence-only datasets. The suffix of the key datasets will be
      `seq_cnt`. (default: off)
    * add option `-a` to create both sequence and sequence-count datasets.
      (default: off)
    * add option `-f` to overwrite key datasets if already exist.
  + add_spill_index
    * now keeps the data type same as the one provided in option '-s' in
      sort_file_list.
    * now can have any number of index datasets in option '-k'.
    * No longer requires index datasets to H5T_STD_U32LE.
  + New utility program nu_stat - to collect and print some statistics of data
    objects stored in an HDF5 file.
  + New utility program nfile_ocopy - to merge multiple HDF5 files into one, by
    calling H5Ocopy() for each group.


* Other updates
  + Chunk size setting now is in the unit of array elements, instead of MB.
  + Change the default chunk size to 256K array elements.
  + Add a case study based on PandAna read operations.

* Bugs fixed
  + Fix compile errors on Mac OS.
  + Fix when the number of processes is more than number of input files
  + Fix when a dataset is empty in one file but not in others


## Version 1.0.0 (Feb. 14, 2020)

* New features
  + Add appending mode. This allows to concatenate new input files and append
    to the end of an existing file. This feature enables users to concatenate
    files in multiple runs, which can be useful when there is a long list of
    input files to be concatenated. See 4e94ccc
  + Add support for datasets of string data type. See 5cdd0b7

* New command-line options
  + `-c` to store all datasets in contiguous storage layout. See 0343ecb

* Partitioning key datasets
  + ph5_concat now produces consistent partition key among all groups.
    See 4836809

* Utility `add_key`
  + Add option `-f` to overwrite key datasets if exist. See 0e687d5
  + Add option `-a` to create both `seq` and `seq_cnt` datasets. See c7d2990
  + Add option `-c` for creating sequence-count datasets. See 038d067
  + Add command-line option `-k` to allow users to provide the name of
    partition key base dataset. See 4836809

* New utility programs
  + `nfile_ocopy` merge files in parallel using C API `H5Ocopy`. See 57d8005
  + `nfile_merge` merge files in parallel by reading and writing compressed
    chunks directly. See 0da0034
  + `nfile_link` adds files as external links to a master file. See 0a520ac
  + `nu_stat` collects and display statistics of the groups and datasets in the
    file. See [README.md](utils/README.md), PR #4 and 708d5ff
  + `add_spill_index` adds a new index dataset in group '/spill'. Thanks Derek
    Doyle's contributions. See c8c2701
  + `check_seq_incr` checks the contents of partition key datasets in all
    groups. See [README.md](utils/README.md) for detailed information of usage
    and a0e51c9.
  + `sort_file_list` is a utility program to sort the file paths based on the
    values of first element of a list of index datasets provided in the
    argument of command-line option '-k'. See [README.md](utils/README.md)
    for detailed information of usage and 5da60ad.

* Case study
  + PandAna read of concatenated file. See 8005932

Bug fix:
  + Fix reading a subarray of datasets of variable-length string type. See
    dcfd754
  + Fix when a dataset is empty in one file but not in another. See 09e270c
  + Fix compression threshold calculation. See a9ea841

Other changes
  + Remove the HDF5 metadata memory footprint profiling, as H5get_alloc_stats
    has been removed from HDF5 since 1.14.0. See PR #5.
  + Remove the requirement of group `/spill`. See d4e0852
  + Change chunk setting for small datasets. Set the chunk size to 256K
    element, no matter if the dataset dim[0] is smaller than 256K element or
    not. See bb118bd
  + Improve NOvA MC indexing compatibility. Thanks Derek Doyle. See PR #2.
  + Extend the utility programs, add_keys and sort_file_list, to handle the
    case when the key dataset is missing from the input files. See PR #2.
  + Add memory footprint profiling. See PR #2.
  + Adding instrumentation for monitoring memory usage. Thanks to Derek Doyle's
    contributions. See 68b3f4f
  + Chunk size: change from 1-MiB based to 256K element based. See 8e85e29
    Chunk size setting is changed from 1-MiB based to 256K element based. This
    change is to allow all datasets in the same group to have the same chunk
    dimension size along the first dimension (most significant one). If the
    partition key dataset, e.g. evt.seq, is used to calculated the partitioning
    boundaries for all processes, then all datasets in the group can use the
    same partitioning boundaries. Note datasets in a group can be of different
    data types, which makes the 1-MiB based chunk size a poor choice for
    parallel reads (i.e. hard to implement a consistent chunk partitioning).
    See bd255d9 and 8e85e29
  + Change the chunk setting for non-empty datasets to be always 256K element
    based, no matter if their sizes are smaller than 256K elements. Previously,
    the chunk sizes are set to their current sizes if smaller than 256K
    elements.  See bb118bd

---

