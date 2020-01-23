#
# Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
# See COPYRIGHT notice in top-level directory.
#
# @configure_input@

SUFFIXES = .cpp .hpp .o

SUBDIRS = . utils

AM_CPPFLAGS = -I.
if PH5CONCAT_DEBUG
   AM_CPPFLAGS += -DDEBUG=1
endif
if PH5CONCAT_PROFILING
   AM_CPPFLAGS += -DPROFILE=1
endif

noinst_PROGRAMS = ph5_concat

ph5_concat_SOURCES = ph5_concat.cpp \
                     file_partition.cpp \
                     dataset_partition.cpp \
                     main.cpp

EXTRA_DIST = ph5_concat.hpp \
             COPYRIGHT \
             CREDITS.md \
             AUTHORS.md \
             README.md \
             RELEASE_NOTES.md \
             list.txt

dist-hook:
	$(SED_I) -e "s|_PH5CONCAT_VERSION_|@PH5CONCAT_VERSION@|g" $(distdir)/main.cpp
	$(SED_I) -e "s|_PH5CONCAT_RELEASE_DATE_|@PH5CONCAT_RELEASE_DATE@|g" $(distdir)/main.cpp
	$(SED_I) -e "1,10s|_PH5CONCAT_RELEASE_DATE_|@PH5CONCAT_RELEASE_DATE@|" $(distdir)/RELEASE_NOTES.md
	$(SED_I) -e "1,10s|_PH5CONCAT_VERSION_|$(PH5CONCAT_VERSION)|" $(distdir)/RELEASE_NOTES.md

.PHONY:
