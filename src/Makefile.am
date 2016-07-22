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

AUTOMAKE_OPTIONS = subdir-objects

if SERIALTARGET
lib_LIBRARIES = libPOLY.a
#libPOLY_a_CXXFLAGS = @TP_CFLAGS@ @NPT_CFLAGS@ @REAL_OPT@ -I$(top_srcdir)/include
libPOLY_a_CXXFLAGS = @PL_BUILD_CXXFLAGS@ @TP_CFLAGS@ @NPT_CFLAGS@  -I../include -I$(top_srcdir)/include
libPOLY_a_SOURCES = \
     Polylib.cxx \
     groups/PolygonGroup.cxx \
     groups/VTree.cxx \
     polygons/Triangle.cxx \
     polygons/NptTriangle.cxx \
     file_io/FileIO_func.cxx \
     file_io/PolygonIO.cxx \
     c_lang/CPolylib.cxx \
     f_lang/FPolylib.cxx \
     util/time.cxx
else
lib_LIBRARIES = libPOLYmpi.a
#libPOLYmpi_a_CXXFLAGS = @MPI_CFLAGS@ @TP_CFLAGS@ @NPT_CFLAGS@ @MPI_PL_OPT@ @REAL_OPT@ -I$(top_srcdir)/include
libPOLYmpi_a_CXXFLAGS = @PL_BUILD_CXXFLAGS@ @TP_CFLAGS@ @NPT_CFLAGS@ -I../include -I$(top_srcdir)/include
libPOLYmpi_a_SOURCES = \
     Polylib.cxx \
     Polylib_MPI.cxx \
     groups/PolygonGroup.cxx \
     groups/PolygonGroup_MPI.cxx \
     groups/VTree.cxx \
     polygons/Triangle.cxx \
     polygons/NptTriangle.cxx \
     file_io/FileIO_func.cxx \
     file_io/PolygonIO.cxx \
     c_lang/CPolylib.cxx \
     f_lang/FPolylib.cxx \
     util/time.cxx
endif


# polyVersion.h.in -> BUILD_DIR/include/polyVersion.h -> install
POLY_includedir = $(includedir)
POLY_include_HEADERS = \
  ../include/polyVersion.h

# SRC_DIR/include/* -> install
nobase_include_HEADERS = \ 
  ../include/Polylib.h \
  ../include/Polylib_func.h \
  ../include/c_lang/CPolylib.h \
  ../include/common/BBox.h \
  ../include/common/PolylibCommon.h \
  ../include/common/PolylibDefine.h \
  ../include/common/PolylibStat.h \
  ../include/common/tt.h \
  ../include/common/Vec2.h \
  ../include/common/Vec3.h \
  ../include/f_lang/FPolylib.h \
  ../include/f_lang/FPolylib_User.inc \
  ../include/f_lang/FPolylib_User_f77.inc \
  ../include/f_lang/FPolylib_define.inc \
  ../include/f_lang/FPolylib_define_f77.inc \
  ../include/f_lang/FPolylib_precision.inc \
  ../include/f_lang/FPolylib_precision_def.inc \
  ../include/f_lang/FPolylib_prototype.inc \
  ../include/f_lang/FPolylib_prototype_f77_double.inc \
  ../include/f_lang/FPolylib_prototype_f77_float.inc \
  ../include/f_lang/FPolylib_struct.inc \
  ../include/file_io/FileIO_func.h \
  ../include/file_io/PolygonIO.h \
  ../include/groups/PolygonGroup.h \
  ../include/groups/VTree.h \
  ../include/polygons/NptTriangle.h \
  ../include/polygons/Triangle.h \
  ../include/util/time.h

EXTRA_DIST = Makefile_hand depend.inc
