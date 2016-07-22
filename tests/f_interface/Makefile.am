
WSOURCES=$(top_srcdir)/tests/f_interface/main.F90
WFFLAGS= $(F90FLAGS) @PL_BUILD_FFLAGS@ -I$(top_srcdir)/include/f_lang
WLDADD= -L../../src @PL_LIBS@ @PL_BUILD_FORTRAN_LDFLAGS@ @NPT_LDFLAGS@ @TP_LDFLAGS@

check_SCRIPT=check
check:
	@PL_FC@ -o test $(WSOURCES) $(WFFLAGS) $(WLDADD)
	ln -f -s $(top_srcdir)/tests/@TEST_RUN_MPI@ run.sh
	ln -f -s $(top_srcdir)/tests/f_interface/polylib_config.tp polylib_config.tp
	ln -f -s $(top_srcdir)/tests/f_interface/plane_in.stl  plane_in.stl
	ln -f -s $(top_srcdir)/tests/f_interface/sphere_in.stl sphere_in.stl
	sh ./run.sh > log.txt

DISTCLEANFILES=
CLEANFILES=test log.txt *.tpp *.stla

dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)
