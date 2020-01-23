#!/bin/sh
#
# Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

INPUT_DIR=ND_data
OUTPUT_DIR=ND_out

file_no=1
for entry in "$INPUT_DIR"/*
do
  if [ ${entry: -3} == ".h5" ]
  then
     infile_name=`basename $entry`
     echo "$infile_name"

     infile_size=`stat -c %s $entry`
     infile_size=$((infile_size / 1048575))
     echo "infile_size=$infile_size MiB"

     outfile_name="$OUTPUT_DIR/nd_out_$file_no.h5"
     rm -f $outfile_name

     ./rechunk -o $outfile_name $entry

     outfile_size=`stat -c %s $outfile_name`
     outfile_size=$((outfile_size / 1048575))
     echo "$outfile_name: size $outfile_size MiB"
     echo "---------------------------------------------------"

     file_no=$((file_no + 1))

  fi
done

