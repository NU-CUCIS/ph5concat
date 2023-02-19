#
# Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

echo "srcdir=$srcdir"

rm -f ./out.h5
mpiexec -n 2 ../ph5_concat -i ${srcdir}/sample_list.txt -o out.h5

rm -f ./out.h5
mpiexec -n 2 ../ph5_concat -i ${srcdir}/sample_list.txt -o out.h5 -k evt

rm -f ./out.h5
mpiexec -n 1 ../ph5_concat -i ${srcdir}/vl_string.txt -o out.h5 -b 1

