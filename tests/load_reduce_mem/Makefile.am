
WSOURCES=$(top_srcdir)/tests/load_reduce_mem/main.cxx
WCFLAGS= @PL_BUILD_CXXFLAGS@ @TP_CFLAGS@ @NPT_CFLAGS@ -I$(top_srcdir)/include -I$(top_srcdir)/tests/include
WLDADD= -L../../src @PL_LIBS@ @PL_BUILD_CXX_LDFLAGS@ @NPT_LDFLAGS@ @TP_LDFLAGS@

check_SCRIPT=check
check:
	@PL_CXX@ -o test $(WSOURCES) $(WCFLAGS) $(WLDADD)
	ln -f -s $(top_srcdir)/tests/@TEST_RUN_MPI@ run.sh
	ln -f -s $(top_srcdir)/tests/load_reduce_mem/polylib_config.tp polylib_config.tp
	ln -f -s $(top_srcdir)/tests/load_reduce_mem/plane_in.stl  plane_in.stl
	ln -f -s $(top_srcdir)/tests/load_reduce_mem/sphere_in.stl sphere_in.stl
	sh ./run.sh > log.txt

DISTCLEANFILES=
CLEANFILES=test log.txt *.tpp *.stla

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)

