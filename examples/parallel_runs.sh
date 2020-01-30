#
# Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
# See COPYRIGHT notice in top-level directory.
#

# Exit immediately if a command exits with a non-zero status.
set -e

rm -f ./out.h5
mpiexec -n 2 ../ph5_concat -i sample_list.txt

rm -f ./out.h5
mpiexec -n 2 ../ph5_concat -i sample_list.txt -k evt

