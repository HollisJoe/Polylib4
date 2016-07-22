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

bin_PROGRAMS=stl_to_npt

if SERIALTARGET
stl_to_npt_SOURCES=main.cxx
else
stl_to_npt_SOURCES=main.cxx gather_polygons_to_npt.cxx
endif
stl_to_npt_CXXFLAGS= @PL_BUILD_CXXFLAGS@ @NPT_CFLAGS@  @TP_CFLAGS@ -I$(top_srcdir)/include
stl_to_npt_LDADD= -L../../src @PL_LIBS@ @PL_BUILD_CXX_LDFLAGS@ @NPT_LDFLAGS@ @TP_LDFLAGS@

EXTRA_DIST = Makefile_hand sample sample/sphere.stl
