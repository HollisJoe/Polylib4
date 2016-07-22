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

bin_PROGRAMS=npt_to_stl

npt_to_stl_SOURCES=main.cxx
npt_to_stl_CXXFLAGS= @PL_BUILD_CXXFLAGS@ @NPT_CFLAGS@  @TP_CFLAGS@ -I$(top_srcdir)/include
npt_to_stl_LDADD= -L../../src @PL_LIBS@ @PL_BUILD_CXX_LDFLAGS@ @NPT_LDFLAGS@ @TP_LDFLAGS@

EXTRA_DIST = Makefile_hand sample sample/sphere.npt
