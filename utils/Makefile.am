#
# Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

SUFFIXES = .c .h .o

AM_DEFAULT_SOURCE_EXT = .c

AM_CPPFLAGS = -I.
if PH5CONCAT_DEBUG
   AM_CPPFLAGS += -DDEBUG=1
endif
if PH5CONCAT_PROFILING
   AM_CPPFLAGS += -DPROFILE=1
endif

noinst_PROGRAMS = rechunk \
                  add_key

EXTRA_DIST = README.md \
             batch.sh

dist-hook:
	$(SED_I) -e "s|_PH5CONCAT_VERSION_|@PH5CONCAT_VERSION@|g" $(distdir)/rechunk.c
	$(SED_I) -e "s|_PH5CONCAT_RELEASE_DATE_|@PH5CONCAT_RELEASE_DATE@|g" $(distdir)/rechunk.c
	$(SED_I) -e "s|_PH5CONCAT_VERSION_|@PH5CONCAT_VERSION@|g" $(distdir)/add_key.c
	$(SED_I) -e "s|_PH5CONCAT_RELEASE_DATE_|@PH5CONCAT_RELEASE_DATE@|g" $(distdir)/add_key.c

.PHONY: