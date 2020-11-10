#
# Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

SUFFIXES = .cpp .c .h .o

AM_CPPFLAGS = -I.
if PH5CONCAT_DEBUG
   AM_CPPFLAGS += -DDEBUG=1
endif
if PH5CONCAT_PROFILING
   AM_CPPFLAGS += -DPROFILE=1
endif

noinst_PROGRAMS = pandana_read

pandana_read_SOURCES = pandana_read.c

EXTRA_DIST = README.md \
             dset.txt

dist-hook:
	$(SED_I) -e "s|_PH5CONCAT_VERSION_|@PH5CONCAT_VERSION@|g" $(distdir)/pandana_read.c
	$(SED_I) -e "s|_PH5CONCAT_RELEASE_DATE_|@PH5CONCAT_RELEASE_DATE@|g" $(distdir)/pandana_read.c

.PHONY:
