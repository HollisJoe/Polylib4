#! /bin/sh

prefix=@prefix@
exec_prefix=@exec_prefix@
includedir=@includedir@
libdir=@libdir@

usage()
{
    cat <<EOF

Usage: polylib-config [OPTION]

Known values for OPTION are:

  --cxx       print C++ compiler command
  --cflags    print C/C++ pre-processor and compiler flags
  --libs      print library linking information for C++ program
  --help      display this help and exit
  --version   output version information

EOF

    exit $1
}

if test $# -eq 0; then
    usage 1
fi

cflags=false
libs=false

while test $# -gt 0; do
    case "$1" in
    -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
    *) optarg= ;;
    esac

    case "$1" in
    --version)
cat <<EOF

Polylib : Polygon Management Library  Version : @VERSION@ : @PL_REVISION@

Copyright (c) 2010-2011 VCAD System Research Program, RIKEN. 
All rights reserved.

Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
All rights reserved.

EOF
      exit 0
      ;;

    --help)
      usage 0
      ;;

    --cxx)
      echo @CMAKE_CXX_FLAGS@ @MPI_PL_OPT@ @REAL_OPT@ @NPT_OPT@ -I@PL@/include -I@PL@/include/common
      ;;

    --cflags)
      echo @CMAKE_CFLAGS@ @MPI_PL_OPT@ @REAL_OPT@ @NPT_OPT@ -I@PL@/include -I@PL@/include/common
      ;;

    --libs)
      echo -L@PL@/lib -l@PL_LIB@  @NPT_DIR_L@ @NPT_LIB_L@  -L@TP_DIR@/lib -l@TP_LIB@
      ;;

    *)
      usage
      exit 1
      ;;
    esac
    shift
done

exit 0