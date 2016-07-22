
########### install files ###############

install(FILES  ${PROJECT_SOURCE_DIR}/doc/refman.pdf ${PROJECT_SOURCE_DIR}/doc/Polylib_UG.pdf DESTINATION ${PL}/doc)



#original Makefile.am contents follow:

###############################################################################
##
## Polylib - Polygon Management Library
##
## Copyright (c) 2010-2011 VCAD System Research Program, RIKEN. 
## All rights reserved.
## 
## Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
## All rights reserved.
##
###############################################################################
#
#doc_includedir = $(prefix)/doc
#doc_include_HEADERS = \
#    $(top_builddir)/doc/refman.pdf \
#    $(top_builddir)/doc/Polylib_UG.pdf