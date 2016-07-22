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

noinst_PROGRAMS=exampleWindmill

DISTCLEANFILES=
CLEANFILES=

exampleWindmill_SOURCES=main.cxx
exampleWindmill_CXXFLAGS= @PL_BUILD_CXXFLAGS@ @TP_CFLAGS@ @NPT_CFLAGS@ -I./ -I$(top_srcdir)/include
exampleWindmill_LDADD= -L../../src @PL_LIBS@ @PL_BUILD_CXX_LDFLAGS@ @NPT_LDFLAGS@ @TP_LDFLAGS@

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)
