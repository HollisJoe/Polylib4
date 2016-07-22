##############################################################################
#
# Polylib - Polygon Management Library
#
# Copyright (c) 2010-2011 VCAD System Research Program, RIKEN. 
# All rights reserved.
# 
# Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
##############################################################################
#
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.57])

AC_INIT([Polylib], [4.0.0], [keno@riken.jp],[Polylib])

echo "----------------------------------------"
echo "Running Polylib Configure Script"
echo "----------------------------------------"


AC_CONFIG_SRCDIR([src/Polylib.cxx])


#
# Specify directory for auxillary build tools (e.g., install-sh, config.sub, config.guess) and M4 files.
#

AC_CONFIG_AUX_DIR(config)


#
# To prevent the source codes from breaking, we suppress to do configure in the source directory.
#

if test -e configure.ac ; then
echo "Please run configure from a separate build directory. It is not allowd to run configure in the source directory."
echo "Pleaes try ..."
echo " "
echo "$ cd BUILD_DIR"
echo "$ ../configure [options]"
exit
fi


#
# Configure should create config.h from config.h.in
#

AC_CONFIG_HEADERS([config.h])


#
# Use automake
#

AM_INIT_AUTOMAKE()



# add to switch libtool
#AC_CONFIG_MACRO_DIR([m4]) // SO

#
# revision No.
#
PL_REVISION="20160228_1200"
AC_SUBST(PL_REVISION)


#
# C++
#
PL_CXX="$CXX"
AC_SUBST(PL_CXX)
PL_CLAGS=""

PL_CC="$CC"
AC_SUBST(PL_CC)

#
# FORTRAN
#

PL_FC=$F90
AC_SUBST(PL_FC)

#
# Additional LDFLAGS for tests
#  excepts for tests, added LDFLAGS automatically.
#
PL_LDFLAGS_ADD="$LDFLAGS"
AC_SUBST(PL_LDFLAGS_ADD)

#
# Additional Library for tests
#  excepts for tests, added LIBS automatically.
#
PL_LIBS_ADD="$LIBS"
AC_SUBST(PL_LIBS_ADD)


enable_mpi="none"

# if a wrapper compiler is used, "enable_mpi" is empty, but parallel.
# wrappers like mpicxx, mpicc, mpiCC, mpic++, mpiicc, mpiicpc, mpiFCCpx are taking into consideration, if not add in case.
# CC is for Hopper

case "$PL_CXX" in
  mpicxx)   enable_mpi="mpi" ;;
  mpicc)    enable_mpi="mpi" ;;
  mpiCC)    enable_mpi="mpi" ;;
  mpic++)   enable_mpi="mpi" ;;
  mpiicc)   enable_mpi="mpi" ;;
  mpiicpc)  enable_mpi="mpi" ;;
  mpiFCCpx) enable_mpi="mpi" ;;
  CC)       enable_mpi="mpi" ;;
esac

#
# MPI environment
#
AC_SUBST(MPI_DIR)
AC_SUBST(MPI_CFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_PL_OPT)


# mpich
AC_ARG_WITH(mpich, [AC_HELP_STRING([--with-mpich=dir],[Specify MPICH install directory])], , with_mpich=none)
if test x"$with_mpich" != x"none" ; then
  MPI_DIR=$with_mpich;
  MPI_CFLAGS=-I$MPI_DIR/include
  MPI_LDFLAGS=-L$MPI_DIR/lib
  MPI_LIBS=-lmpich
  enable_mpi="mpich"
fi


# openmpi
if test x"$with_mpich" = x"none" ; then
  AC_ARG_WITH( ompi, [AC_HELP_STRING([--with-ompi=dir],[Specify OpenMPI install directory])], , with_ompi=none)
  if test x"$with_ompi" != x"none" ; then
    MPI_DIR=$with_ompi;
    MPI_CFLAGS=-I$MPI_DIR/include
    MPI_LDFLAGS=-L$MPI_DIR/lib
    MPI_LIBS=-lmpi
    enable_mpi="ompi"
  fi
fi

if test x"$enable_mpi" = x"none" ; then
  echo "  Warning: Neither MPICH nor OpenMPI library is not specified."
  echo "           If compilation error occurs about MPI, specify --with-mpich or --with-ompi option."
fi

#
# test execution command
#
AC_SUBST(TEST_RUN_SERIAL)
TEST_RUN_SERIAL=test_run_serial.sh

# mpi execution command & mpi compile option
#       $MPI_EXEC $MPI_PREFIX_NUM_PROC $MPI_TEST_NUM_PROC
#       add -DMPI_PL for safe
AC_SUBST(TEST_RUN_MPI)

MPI_PL_OPT=
TEST_RUN_MPI=test_run_serial.sh

if test x"$enable_mpi" = x"mpich" ; then
  MPI_PL_OPT="-DMPI_PL"
  TEST_RUN_MPI=test_run_mpiexec.sh
fi
if test x"$enable_mpi" = x"ompi" ; then
  MPI_PL_OPT="-DMPI_PL"
  TEST_RUN_MPI=test_run_mpirun.sh
fi
if test x"$enable_mpi" = x"mpi" ; then
  MPI_PL_OPT="-DMPI_PL"
  if test x"$PL_CXX" = x"mpiFCCpx" ; then
    TEST_RUN_MPI=test_run_mpiexec.sh
  else
    echo "  Warning: TEST_RUN_MPI is not specified."
    echo "           Tests will be error."
  fi
fi

#
# Library for Misc.
#    for stl unordered_map
#
#CXXFLAGS_ADD="-std=c++11"
CXXFLAGS_ADD=

#
# default prefix setting
#
AC_PREFIX_DEFAULT([/usr/local/Polylib])
if test x"$prefix" = x"NONE" ; then
  PL_INST_DIR=/usr/local/Polylib
else
  PL_INST_DIR="$prefix"
fi


#
# TextParser environment
#
AC_ARG_WITH(parser, [AC_HELP_STRING([--with-parser=dir],[Specify TextParser install directory])], , with_parser=none)
AC_SUBST(TP_DIR)
AC_SUBST(TP_CFLAGS)
AC_SUBST(TP_LDFLAGS)
AC_SUBST(TP_MPI_LDFLAGS)

if test x"$with_parser" != x"none" ; then
  TP_DIR=$with_parser;
  TP_CFLAGS=`$TP_DIR/bin/tp-config --cflags`
  TP_LDFLAGS=`$TP_DIR/bin/tp-config --libs`
else
  echo "  Error: TextParser library is not specified."
  echo "         Specify --with-parser option."
  exit 1
fi


#
# Npatch library
#
AC_ARG_WITH(npatch, [AC_HELP_STRING([--with-npatch=dir],[Specify Nagata patch library install directory])], , with_npatch=none)
AC_SUBST(NPT_OPT)
AC_SUBST(NPT_DIR)
AC_SUBST(NPT_CFLAGS)
AC_SUBST(NPT_LDFLAGS)

if test x"$with_npatch" != x"none" ; then
  NPT_OPT=-DUSE_NPATCH_LIB
  NPT_DIR=$with_npatch;
  NPT_CFLAGS=`$NPT_DIR/bin/npt-config --cflags`
  NPT_LDFLAGS=`$NPT_DIR/bin/npt-config --libs`
else
#  echo "  Error: Nagata patch library is not specified."
#  echo "         Specify --with-npatch option."
#  exit 1
  NPT_OPT=
  NPT_DIR=
  NPT_CFLAGS=
  NPT_LDFLAGS=
fi


#
# example
#
AC_ARG_ENABLE(example, [AC_HELP_STRING([--enable-example=(no|yes)],[make example [no]])], , enable_example=no)

#
# test
#
AC_ARG_ENABLE(test, [AC_HELP_STRING([--enable-test=(no|yes)],[make check [no]])], , enable_test=no)




#
# Polylib special flags
#


if test x"$enable_mpi" = x"none" ; then
  PL_LIBS="-lPOLY"
else
  PL_LIBS="-lPOLYmpi"
fi

AC_SUBST(PL_LIBS)


# SERIALTARGET = true, if serial
AM_CONDITIONAL([SERIALTARGET], [test "$enable_mpi" = "none"])


PL_CFLAGS="-I$PL_INST_DIR/include"
AC_SUBST(PL_CFLAGS)

PL_LDFLAGS="-L$PL_INST_DIR/lib"
AC_SUBST(PL_LDFLAGS)


#
# Type of REAL
#
AC_ARG_WITH(real, [AC_HELP_STRING([--with-real=(float|double)],[Specify real type [float]])], , with_real=float)
AC_SUBST(REAL_OPT)

if test x"$with_real" = x"double" ; then
  REAL_OPT=-D_REAL_IS_DOUBLE_
else
  REAL_OPT=
fi


#
# C interface test
#
AC_ARG_ENABLE(test_c, [AC_HELP_STRING([--enable-test_c=(no|yes)],[Specify testing c interface [no]])], , enable_test_c=no)


#
# Fortran interface test
#
AC_ARG_ENABLE(test_fortran, [AC_HELP_STRING([--enable-test_fortran=(no|yes)],[Specify testing fortran inter face [no]])], , enable_test_fortran=no)


#
# MAKE_SUB_DIRS
#
AC_SUBST(MAKE_SUB_DIRS)
if test x"$with_npatch" != x"none" ; then
MAKE_SUB_DIRS="src tools doc"
else
MAKE_SUB_DIRS="src doc"
fi

if test x"$enable_example" = x"yes" ; then
  MAKE_SUB_DIRS="$MAKE_SUB_DIRS"" examples"
fi

if test x"$enable_test" = x"yes" ; then
  MAKE_SUB_DIRS="$MAKE_SUB_DIRS"" tests"
fi

#
# MAKE_TEST_SUB_DIRS
#
AC_SUBST(MAKE_TEST_SUB_DIRS)
if test x"$with_npatch" != x"none" ; then
MAKE_TEST_SUB_DIRS="file_io_stl file_io_npt search_polygon attribute move_polygon multi_bbox load_reduce_mem"
else
MAKE_TEST_SUB_DIRS="file_io_stl search_polygon attribute move_polygon multi_bbox load_reduce_mem"
fi

if test x"$enable_test_c" = x"yes" ; then
  MAKE_TEST_SUB_DIRS="$MAKE_TEST_SUB_DIRS"" c_interface"
fi

if test x"$enable_test_fortran" = x"yes" ; then
  MAKE_TEST_SUB_DIRS="$MAKE_TEST_SUB_DIRS"" f_interface"
fi


#
# C++ build compiler flags
#
AC_SUBST(PL_BUILD_CXXFLAGS)
PL_BUILD_CXXFLAGS="$CXXFLAGS"" $MPI_CFLAGS"" $MPI_PL_OPT"" $REAL_OPT"" $NPT_OPT"" $CXXFLAGS_ADD"

#
# C++ build linkage flags for examples or tests
#
AC_SUBST(PL_BUILD_CXX_LDFLAGS)
PL_BUILD_CXX_LDFLAGS="$MPI_LDFLAGS"" $MPI_LIBS"

#
# C build compiler flags
#
AC_SUBST(PL_BUILD_CFLAGS)
PL_BUILD_CFLAGS="$CFLAGS"" $MPI_CFLAGS"" $MPI_PL_OPT"" $REAL_OPT"" $NPT_OPT"

#
# C build linkage flags for examples or tests
#
AC_SUBST(PL_BUILD_C_LDFLAGS)
PL_BUILD_C_LDFLAGS="$PL_BUILD_LDFLAGS"" -lstdc++ -lm"

if test x"$enable_mpi" = x"ompi" ; then
  PL_BUILD_C_LDFLAGS="$PL_BUILD_C_LDFLAGS"" -lmpi_cxx"
fi

#
# Fortran build compiler flags for tests
#
AC_SUBST(PL_BUILD_FFLAGS)
PL_BUILD_FFLAGS="$F90FLAGS"" $MPI_PL_OPT"" $REAL_OPT"" $NPT_OPT"

#
# Fortran build linkage flags for tests
#
AC_SUBST(PL_BUILD_FORTRAN_LDFLAGS)
PL_BUILD_FORTRAN_LDFLAGS="$PL_BUILD_LDFLAGS"" -lstdc++ -lm"

if test x"$enable_mpi" = x"ompi" ; then
  PL_BUILD_FORTRAN_LDFLAGS="$PL_BUILD_FORTRAN_LDFLAGS"" -lmpi_cxx"
fi


#
# Checks for programs.
#
AC_PROG_CXX
AC_PROG_CC
#AC_PROG_F77
AC_PROG_FC
AC_PROG_CPP
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# add to switch libtool // SO
#AC_PROG_LIBTOOL
#AM_PROG_LIBTOOL
#AC_SUBST(LIBTOOL_DEPS)

#
# Checks for libraries.
#


#
# Checks for header files.
#
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h])


#
# Checks for typedefs, structures, and compiler characteristics.
#
AC_C_INLINE
AC_TYPE_SIZE_T


#
# Checks for library functions.
#
AC_CHECK_FUNCS([gettimeofday memset strchr strstr strcpy])



AC_CONFIG_FILES([Makefile \
                 doc/Makefile \
                 examples/Makefile \
                 examples/windmill/Makefile \
                 src/Makefile \
                 tests/Makefile \
                 tests/file_io_stl/Makefile \
                 tests/file_io_npt/Makefile \
                 tests/search_polygon/Makefile \
                 tests/attribute/Makefile \
                 tests/move_polygon/Makefile \
                 tests/multi_bbox/Makefile \
                 tests/load_reduce_mem/Makefile \
                 tests/c_interface/Makefile \
                 tests/f_interface/Makefile \
                 tools/Makefile \
                 tools/stl_to_npt/Makefile \
                 tools/npt_to_stl/Makefile \
                 tools/npt_to_stl4/Makefile \
                 include/polyVersion.h \
                 polylib-config ])
AC_OUTPUT

chmod +x ./polylib-config

echo "---------------------------------------------"
echo "Finished Running Polylib Configure Script"
echo "---------------------------------------------"
