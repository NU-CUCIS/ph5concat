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

noinst_PROGRAMS = rechunk \
                  add_key \
                  sort_file_list \
                  check_seq_incr \
                  add_spill_index \
                  nu_stat

rechunk_SOURCES = rechunk.c
add_key_SOURCES = add_key.cpp
sort_file_list_SOURCES = sort_file_list.cpp
check_seq_incr_SOURCES = check_seq_incr.cpp
add_spill_index_SOURCES = add_spill_index.cpp
nu_stat_SOURCES = nu_stat.c

EXTRA_DIST = README.md \
             batch.sh

dist-hook:
	$(SED_I) -e "s|_PH5CONCAT_VERSION_|@PH5CONCAT_VERSION@|g" $(distdir)/rechunk.c
	$(SED_I) -e "s|_PH5CONCAT_RELEASE_DATE_|@PH5CONCAT_RELEASE_DATE@|g" $(distdir)/rechunk.c
	$(SED_I) -e "s|_PH5CONCAT_VERSION_|@PH5CONCAT_VERSION@|g" $(distdir)/add_key.cpp
	$(SED_I) -e "s|_PH5CONCAT_RELEASE_DATE_|@PH5CONCAT_RELEASE_DATE@|g" $(distdir)/add_key.cpp
	$(SED_I) -e "s|_PH5CONCAT_VERSION_|@PH5CONCAT_VERSION@|g" $(distdir)/sort_file_list.cpp
	$(SED_I) -e "s|_PH5CONCAT_RELEASE_DATE_|@PH5CONCAT_RELEASE_DATE@|g" $(distdir)/sort_file_list.cpp
	$(SED_I) -e "s|_PH5CONCAT_VERSION_|@PH5CONCAT_VERSION@|g" $(distdir)/check_seq_incr.cpp
	$(SED_I) -e "s|_PH5CONCAT_RELEASE_DATE_|@PH5CONCAT_RELEASE_DATE@|g" $(distdir)/check_seq_incr.cpp
	$(SED_I) -e "s|_PH5CONCAT_RELEASE_DATE_|@PH5CONCAT_RELEASE_DATE@|g" $(distdir)/add_spill_index.cpp

.PHONY:
