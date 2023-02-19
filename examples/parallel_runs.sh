#!/bin/bash -l
#
# Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

# echo "srcdir=$srcdir"
# echo "TESTMPIRUN=$TESTMPIRUN"
# echo "HDF5_LIB_PATH=$HDF5_LIB_PATH"

export LD_LIBRARY_PATH=${HDF5_LIB_PATH}:${LD_LIBRARY_PATH}

rm -f ./out.h5
CMD="${TESTMPIRUN} -n 2 ../ph5_concat -i ${srcdir}/sample_list.txt -o out.h5"
echo "CMD=$CMD"
$CMD
echo "========================================================================"
echo ""

rm -f ./out.h5
CMD="${TESTMPIRUN} -n 2 ../ph5_concat -i ${srcdir}/sample_list.txt -o out.h5 -k evt"
echo "CMD=$CMD"
$CMD
echo "========================================================================"
echo ""

rm -f ./out.h5
CMD="${TESTMPIRUN} -n 1 ../ph5_concat -i ${srcdir}/vl_string.txt -o out.h5 -b 1"
echo "CMD=$CMD"
$CMD
echo "========================================================================"

