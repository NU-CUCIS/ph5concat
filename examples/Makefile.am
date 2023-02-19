#
# Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

SUFFIXES = .h5 .gz

TESTS = parallel_runs.sh
TEST_EXTENSIONS = .sh
LOG_COMPILER = $(srcdir)/parallel_runs.sh
SH_LOG_COMPILER =

INPUT_FILES = sample_r11981_s06.gz \
              sample_r11981_s07.gz \
              sample_r11981_s08.gz \
              sample_r11981_s09.gz \
              vl_string.gz

BUILT_SOURCES = $(INPUT_FILES:.gz=.h5)

.gz.h5:
	gzip -dc $< > $@

CLEANFILES = core core.* out.h5

DISTCLEANFILES = $(INPUT_FILES:.gz=.h5)

EXTRA_DIST = $(INPUT_FILES) \
             sample_list.txt \
             parallel_runs.sh \
             sample_output.h5 \
             sample_output.metadata \
             sample_stdout.txt \
             vl_string.txt

