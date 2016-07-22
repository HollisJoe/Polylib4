
WSOURCES=$(top_srcdir)/tests/file_io_npt/main.cxx
WCXXFLAGS= @PL_BUILD_CXXFLAGS@ @TP_CFLAGS@ @NPT_CFLAGS@ -I$(top_srcdir)/include
WLDADD= -L../../src @PL_LIBS@ @PL_LIBS@ @PL_BUILD_CXX_LDFLAGS@ @NPT_LDFLAGS@ @TP_LDFLAGS@

check_SCRIPT=check
check:
	@PL_CXX@ -o test $(WSOURCES) $(WCXXFLAGS) $(WLDADD)
	ln -f -s $(top_srcdir)/tests/@TEST_RUN_MPI@ run.sh
	ln -f -s $(top_srcdir)/tests/file_io_npt/polylib_config.tp polylib_config.tp
	ln -f -s $(top_srcdir)/tests/file_io_npt/sphere_in_bin.npt sphere_in_bin.npt
	ln -f -s $(top_srcdir)/tests/file_io_npt/sphere_in_txt.npt sphere_in_txt.npt
	sh ./run.sh > log.txt

DISTCLEANFILES=
CLEANFILES=test log.txt *.tpp *.npta *.nptb

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)
