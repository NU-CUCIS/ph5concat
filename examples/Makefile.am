#
# Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

TESTS = parallel_runs.sh
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/parallel_runs.sh
SH_LOG_COMPILER =

INPUT_FILES = sample_r11981_s06.h5 \
              sample_r11981_s07.h5 \
              sample_r11981_s08.h5 \
              sample_r11981_s09.h5

CLEANFILES = core core.*

EXTRA_DIST = $(INPUT_FILES) \
             sample_list.txt \
             parallel_runs.sh \
             sample_output.h5 \
             sample_output.metadata \
             sample_stdout.txt

