# ph5concat Release Notes

---
## Version _PH5CONCAT_VERSION_ (_PH5CONCAT_RELEASE_DATE_)

* Chunk size setting is changed from 1-MiB based to 256K element based. This
  change is to allow all datasets in the same group to have the same chunk
  dimension size along the first dimension (most significant one). If the
  partition key dataset, e.g. evt.seq, is used to calculated the partitioning
  boundaries for all processes, then all datasets in the group can use the same
  partitioning boundaries. Note datasets in a group can be of different data
  types, which makes the 1-MiB based chunk size a poor choice for parallel
  reads (i.e. hard to implement a consistent chunk partitioning).
  See bd255d97f9e92e41b4bd11ca0bf444cc4b462dd9 and
  8e85e2910680a0d3a60eb8cfaadbe43d019557fa

* Extend the utility programs, add_keys and sort_file_list, to handle the case
  when the key dataset is missing from the input files. See pull request #2.

* Add memory footprint profiling. See pull request #2.

* Change the chunk setting for non-empty datasets to be always 256K element
  based, no matter if their sizes are smaller than 256K elements. Previously,
  the chunk sizes are set to their current sizes if smaller than 256K elements.
  See bb118bd9c36e1f4cffb1127b54d52e64a01deed5

* Add appending mode. This allows to concatenate new input files and append to
  the end of an existing file. This feature enables users to concatenate files
  in multiple runs, which can be useful when there is a long list of input
  files to be concatenated.

---


