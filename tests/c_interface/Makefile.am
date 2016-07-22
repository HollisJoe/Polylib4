
WSOURCES=$(top_srcdir)/tests/c_interface/main.c
#WCFLAGS= $(CFLAGS) @PL_BUILD_CFLAGS@ @NPT_CFLAGS@ -I$(top_srcdir)/include
WCFLAGS= @PL_BUILD_CFLAGS@ @NPT_CFLAGS@ -I$(top_srcdir)/include
WLDADD= -L../../src @PL_LIBS@ @PL_BUILD_C_LDFLAGS@ @NPT_LDFLAGS@ @TP_LDFLAGS@

check_SCRIPT=check
check:
	@PL_CC@ -o test $(WSOURCES) $(WCFLAGS) $(WLDADD)
	ln -f -s $(top_srcdir)/tests/@TEST_RUN_MPI@ run.sh
	ln -f -s $(top_srcdir)/tests/c_interface/polylib_config.tp polylib_config.tp
	ln -f -s $(top_srcdir)/tests/c_interface/plane_in.stl  plane_in.stl
	ln -f -s $(top_srcdir)/tests/c_interface/sphere_in.stl sphere_in.stl
	sh ./run.sh > log.txt

DISTCLEANFILES=
CLEANFILES=test log.txt *.tpp *.stla

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)

