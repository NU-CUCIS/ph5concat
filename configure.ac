dnl
dnl Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
dnl See COPYRIGHT notice in top-level directory.
dnl

dnl -*- Mode: shell-script-mode; -*-
dnl Process this file with GNU autoconf(1) to produce a configure script.
dnl

dnl autoconf v2.69 was released in 2012-04-24
AC_PREREQ([2.69])
AC_INIT([ph5concat],[1.1.0],[],[ph5concat],[])

AC_CONFIG_HEADERS([config.h])
AH_TOP([/*
 * Copyright (C) 2019, Northwestern University and Fermi National Accelerator Laboratory
 * See COPYRIGHT notice in top-level directory.
 */
])
AC_CONFIG_SRCDIR([ph5_concat.hpp])

AC_CONFIG_AUX_DIR([./scripts])

AM_INIT_AUTOMAKE([foreign])
AM_SILENT_RULES([yes])
AM_MAINTAINER_MODE([enable])

dnl parse the version numbers to 3 env variables
PH5CONCAT_VERSION_MAJOR=`echo ${PACKAGE_VERSION} | cut -d. -f1`
PH5CONCAT_VERSION_MINOR=`echo ${PACKAGE_VERSION} | cut -d. -f2`
PH5CONCAT_VERSION_SUB=`echo ${PACKAGE_VERSION} | cut -d. -f3`
PH5CONCAT_VERSION=${PACKAGE_VERSION}
AC_SUBST(PH5CONCAT_VERSION_MAJOR)
AC_SUBST(PH5CONCAT_VERSION_MINOR)
AC_SUBST(PH5CONCAT_VERSION_SUB)
AC_SUBST(PH5CONCAT_VERSION)

dnl Note that command 'date' is not portable across Unix platforms.
dnl But release date matters only to ph5concat developers who make the releases.
PH5CONCAT_RELEASE_DATE="`date '+%B %-d, %Y'`"
AC_SUBST(PH5CONCAT_RELEASE_DATE)
PH5CONCAT_RELEASE_DATE_FULL="`date '+%Y-%m-%d'`"
AC_SUBST(PH5CONCAT_RELEASE_DATE_FULL)
AC_SUBST(PACKAGE_VERSION)

AH_TOP([#ifndef _CONFIG_H
#define _CONFIG_H])
AH_BOTTOM([#endif])

AC_PROG_SED
dnl check sed command option -i and set SED_I (this requires RM defined)
UD_PROG_SED_I
AC_PROG_EGREP

MPI_INSTALL=
AC_ARG_WITH(mpi,
   [AS_HELP_STRING([--with-mpi=/path/to/implementation],
                   [The installation prefix path for MPI implementation.])
   ],[ dnl this clause is run when --with-mpi or --without-mpi is used
   if test "x${withval}" = xno ; then
      AC_MSG_ERROR([
      -----------------------------------------------------------------------
        ph5concat is built on top of MPI. Configure option --without-mpi or
        --with-mpi=no should not be used. Abort.
      -----------------------------------------------------------------------])
   elif test "x${withval}" = x ; then
      AC_MSG_ERROR(--with-mpi is set but the value is NULL)
   elif test "x${withval}" != xyes && test ! -d "${withval}" ; then
      # user may use --with-mpi without an argument, which results in withval
      # being "yes". This case is OK and we simply take no action, as ph5concat
      # requires MPI compilers and will check them.
      AC_MSG_ERROR(Directory '${withval}' specified in --with-mpi does not exist or is not a directory)
   fi
   MPI_INSTALL=${withval}
   ]
)

AC_ARG_VAR(MPICC,  [MPI C compiler, @<:@default: CC@:>@])
AC_ARG_VAR(MPICXX, [MPI C++ compiler, @<:@default: CXX@:>@])

ac_user_MPICC=$MPICC
ac_user_MPICXX=$MPICXX
if test "x$MPICC"  = x && test "x$CC"  != x ; then ac_user_MPICC=$CC   ; fi
if test "x$MPICXX" = x && test "x$CXX" != x ; then ac_user_MPICXX=$CXX ; fi

CANDIDATE_MPICC="${MPICC} mpicc mpicc_r"
CANDIDATE_MPICXX="${MPICXX} mpicxx mpic++ mpiCC mpcxx mpc++ mpicxx_r mpiCC_r mpcxx_r mpic++_r mpc++_r"
dnl add GNU MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpigcc mpgcc mpigcc_r mpgcc_r"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpig++ mpg++ mpig++_r mpg++_r"
dnl add IBM MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpcc_r mpcc mpixlc_r mpixlc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpCC_r mpCC mpixlcxx_r mpixlcxx mpixlC_r mpixlC"
dnl add IBM BGL MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC blrts_xlc mpxlc_r mpxlc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX blrts_xlC mpxlC_r mpxlC mpixlc++ mpxlcxx mpxlc++ mpxlCC mpixlc++_r mpxlcxx_r mpxlc++_r mpxlCC_r"
dnl add Fujitsu MPI compilers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpifccpx"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpiFCCpx"
dnl add Cray MPI compiler wrappers
CANDIDATE_MPICC="$CANDIDATE_MPICC cc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX CC"
dnl add Intel MPI compiler wrappers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpiicc icc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpiicpc mpiicxx mpiic++ mpiiCC icpc"
dnl add PGI MPI compiler wrappers
CANDIDATE_MPICC="$CANDIDATE_MPICC mpipgcc mppgcc"
CANDIDATE_MPICXX="$CANDIDATE_MPICXX mpipgCC mppgCC"

dnl find the full path of MPICC from CANDIDATE_MPICC and MPI_INSTALL
if test "x${ac_user_MPICC}" = x ; then
   dnl if MPICC or CC has not been set by users, then search from
   dnl CANDIDATE_MPICC, and find the full path of MPICC
   UD_MPI_PATH_PROGS([MPICC], [$CANDIDATE_MPICC])
else
   dnl check whether user specified MPICC is valid
   UD_MPI_PATH_PROG([MPICC], [$ac_user_MPICC])
fi

if test "x${MPICC}" = x ; then
   if test "x$ac_user_MPICC" = x ; then
      ERR_MSG="No MPI C compiler can be found"
   else
      ERR_MSG="Specified MPI C compiler \"$ac_user_MPICC\" cannot be found"
   fi
   if test "x$MPI_INSTALL" != x ; then
      ERR_MSG="$ERR_MSG under $MPI_INSTALL"
   fi
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     $ERR_MSG
     ph5concat requires a working MPI C compiler. Please specify the
     location of an MPI C compiler, either in the MPICC environment
     variable (not CC variable) or through --with-mpi configure flag.
     Abort.
   -----------------------------------------------------------------------])
fi
CC=${MPICC}
AC_PROG_CC

dnl find the full path of MPICXX from CANDIDATE_MPICXX and MPI_INSTALL
if test "x${ac_user_MPICXX}" = x ; then
   dnl if MPICXX or CXX has not been set by users, then search from
   dnl CANDIDATE_MPICXX, and find the full path of MPICXX
   UD_MPI_PATH_PROGS([MPICXX], [$CANDIDATE_MPICXX])
else
   dnl check whether user specified MPICXX is valid
   UD_MPI_PATH_PROG([MPICXX], [$ac_user_MPICXX])
fi
if test "x${MPICXX}" = x ; then
   if test "x$ac_user_MPICXX" = x ; then
      ERR_MSG="No MPI C++ compiler can be found"
   else
      ERR_MSG="Specified MPI C++ compiler \"$ac_user_MPICXX\" cannot be found"
   fi
   if test "x$MPI_INSTALL" != x ; then
      ERR_MSG="$ERR_MSG under $MPI_INSTALL"
   fi
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     $ERR_MSG
     ph5concat requires a working MPI C++ compiler. Please specify the
     location of an MPI C++ compiler, either in the MPICXX environment
     variable (not CXX variable) or through --with-mpi configure flag.
     Abort.
   -----------------------------------------------------------------------])
fi
CXX=${MPICXX}
AC_PROG_CXX

dnl Set output variable CPP to a command that runs the C preprocessor.
dnl Some C compilers require -E to be used as C preprocessor.
AC_PROG_CPP

dnl check if MPICXX works for basic MPI call: MPI_Comm_rank()
AC_LANG_PUSH(C++)
AC_CHECK_FUNC([MPI_Comm_rank], [],
   dnl maybe -lmpi is needed at link stage
   [AC_SEARCH_LIBS([MPI_Comm_rank], [mpi mpi++ mpich mpichcxx mpi_cxx], [],
                   [AC_MSG_ERROR([
   -----------------------------------------------------------------------
     Invalid MPI compiler specified or detected: "${MPICXX}"
     A working MPI C++ compiler is required. Please specify the location
     of one either in the MPICXX environment variable (not CXX variable)
     or through --with-mpi configure flag. Abort.
   -----------------------------------------------------------------------])
])])

AC_CHECK_FUNC([MPI_File_open], [],
   dnl maybe -lmpi++ is needed at link stage
   [AC_SEARCH_LIBS([MPI_File_open], [mpio], [],
                   [AC_MSG_ERROR([
   -----------------------------------------------------------------------
     The underneath MPI implementation does not support MPI-IO.
     ph5concat requires MPI-IO support to work properly. Abort.
   -----------------------------------------------------------------------])
])])
AC_LANG_POP(C++)

AC_CHECK_DECL([access], [], [], [[#include <unistd.h>]])
if test "x$ac_cv_have_decl_access" = xyes ; then
   AC_CHECK_FUNCS([access])
fi
AC_CHECK_DECL([unlink], [], [], [[#include <unistd.h>]])
if test "x$ac_cv_have_decl_unlink" = xyes ; then
   AC_CHECK_FUNCS([unlink])
fi

AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
                    [Enable ph5concat internal debug mode.
                     @<:@default: disabled@:>@])],
    [debug=${enableval}], [debug=no]
)
AM_CONDITIONAL(PH5CONCAT_DEBUG, [test "x$debug" = xyes])

if test "x${debug}" = xyes; then
   dnl add -g flag if not presented
   dnl remove all -O and -fast flags
   dnl add -O0 to all flags
   # check exit status of grep command is more portable than using -q
   str_found=`echo "${CXXFLAGS}" | ${EGREP} -- "-g"`
   if test "x$?" != x0 ; then
      CXXFLAGS="$CXXFLAGS -g"
   fi
   CXXFLAGS=`echo $CXXFLAGS | ${SED} 's/-O. *//g' | ${SED} 's/-fast *//g'`
   CXXFLAGS="$CXXFLAGS -O0"
   unset str_found
fi

AC_ARG_ENABLE([profiling],
   [AS_HELP_STRING([--enable-profiling],
                   [Enable internal time profiling. @<:@default: disabled@:>@])],
   [enable_profiling=${enableval}], [enable_profiling=no]
)
AM_CONDITIONAL(PH5CONCAT_PROFILING, [test "x$enable_profiling" = xyes])

AC_ARG_WITH([hdf5],
   [AS_HELP_STRING([--with-hdf5=/path/to/hdf5], [Specify HDF5 installation path(s):
    --with-hdf5=INC,LIB for include directory and library directory separated by a comma
    --with-hdf5=DIR for directory containing include/ and lib/ subdirectories
   ])], [
   case $withval in
     *,*)
        hdf5_inc="`echo $withval |cut -f1 -d,`"
        hdf5_lib="`echo $withval |cut -f2 -d, -s`"
        ;;
     *)
        if test -n "$withval"; then
          hdf5_inc="$withval/include"
          hdf5_lib="$withval/lib"
        fi
        ;;
   esac
   if test "x$hdf5_inc" != x ; then
      if test "x$CPPFLAGS" = x ; then
         CPPFLAGS="-I$hdf5_inc"
      elif ! echo "${CPPFLAGS}" | ${EGREP} -q -w -- "-I$hdf5_inc" ; then
         # append only if not already appear in CPPFLAGS
         CPPFLAGS="$CPPFLAGS -I$hdf5_inc"
      fi
   fi
   if test "x$hdf5_lib" != x ; then
      if test "x$LDFLAGS" = x ; then
         LDFLAGS="-L$hdf5_lib"
      elif ! echo "${LDFLAGS}" | ${EGREP} -q -w -- "-L$hdf5_lib" ; then
         # append only if not already appear in LDFLAGS
         LDFLAGS="$LDFLAGS -L$hdf5_lib"
      fi
   fi
])

# Only HDF5 1.10.x supports compression
AC_LANG_PUSH(C++)
AC_CHECK_HEADER([hdf5.h], [have_hdf5=yes], [have_hdf5=no])
if test "x$have_hdf5" = xno ; then
      AC_MSG_ERROR([
      -----------------------------------------------------------------------
      Missing HDF5-header files 'hdf5.h' required to build ph5concat. Use
      configure command-line option --with-hdf5=/path/to/implementation
      to specify the location of HDF5 installation. Abort.
      -----------------------------------------------------------------------])
fi
AC_MSG_CHECKING([whether HDF5 version is 1.10.4 or later])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <hdf5.h>
#if (H5_VERS_MAJOR*1000000 + H5_VERS_MINOR*1000 + H5_VERS_RELEASE < 1010004)
#error HDF5 version is older than 1.10.4
#endif
   ]])], [hdf5_ge_1_10_4=yes], [hdf5_ge_1_10_4=no])
AC_MSG_RESULT([$hdf5_ge_1_10_4])
if test x$hdf5_ge_1_10_4 = xno; then
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     ph5concat requires HDF5 1.10.4 and later. Abort.
   -----------------------------------------------------------------------])
fi
AC_MSG_CHECKING([whether HDF5 parallel I/O is enabled])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
#include <hdf5.h>
#if !defined(H5_HAVE_PARALLEL) || H5_HAVE_PARALLEL != 1
#error Parallel I/O is not enabled in HDF5
#endif
]])], [hdf5_parallel_io=yes], [hdf5_parallel_io=no])
AC_MSG_RESULT([$hdf5_parallel_io])
if test "x$hdf5_parallel_io" = xno ; then
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     ph5concat requires parallel I/O feature enabled in HDF5. Abort.
   -----------------------------------------------------------------------])
fi
AC_LANG_POP(C++)

if test "x$LIBS" = x ; then
   LIBS="-lhdf5"
else
   LIBS="$LIBS -lhdf5"
fi

saved_cross_compiling=$cross_compiling
cross_compiling="yes"
# Below two AC_CHECK_SIZEOF will try to link and run. Setting cross_compiling
# allows those checks to be done at compile time.
AC_CHECK_SIZEOF([size_t])
AC_CHECK_SIZEOF([hsize_t], [], [#include <hdf5.h>])
cross_compiling=$saved_cross_compiling
unset saved_cross_compiling
if test "x$ac_cv_sizeof_size_t" = x0 ; then
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     Cannot detect size of type 'size_t' Abort.
   -----------------------------------------------------------------------])
elif test "x$ac_cv_sizeof_size_t" = x8 ; then
   AC_DEFINE(MPI_SIZE_T, MPI_UNSIGNED_LONG_LONG, [Define MPI datatype corresponding to size_t])
else
   AC_DEFINE(MPI_SIZE_T, MPI_UNSIGNED_INT, [Define MPI datatype corresponding to size_t])
fi
if test "x$ac_cv_sizeof_hsize_t" = x0 ; then
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     Cannot detect size of type 'hsize_t' Abort.
   -----------------------------------------------------------------------])
elif test "x$ac_cv_sizeof_hsize_t" = x8 ; then
   AC_DEFINE(MPI_HSIZE_T, MPI_UNSIGNED_LONG_LONG, [Define MPI datatype corresponding to hsize_t])
else
   AC_DEFINE(MPI_HSIZE_T, MPI_UNSIGNED_INT, [Define MPI datatype corresponding to hsize_t])
fi

AC_CONFIG_FILES(Makefile utils/Makefile examples/Makefile)

AC_OUTPUT
echo "------------------------------------------------------------------------------"
echo \
"
   ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

   Features:  Internal debug mode         - ${debug}
              Internal profiling mode     - ${enable_profiling}"

echo "\

   Compilers: MPICC    = ${MPICC}
              MPICXX   = ${MPICXX}"
if test "x${CPPFLAGS}" != x ; then
   echo "\
              CPPFLAGS = ${CPPFLAGS}"
fi
echo "\
              CFLAGS   = ${CFLAGS}"
if test "${CXXFLAGS}" != x ; then
   echo "\
              CXXFLAGS = ${CXXFLAGS}"
fi
if test "x${LDFLAGS}" != x ; then
   echo "\
              LDFLAGS  = ${LDFLAGS}"
fi
if test "x${LIBS}" != x ; then
   echo "\
              LIBS     = ${LIBS}"
fi
echo "\

   Now run 'make' to build the executable.

------------------------------------------------------------------------------"

