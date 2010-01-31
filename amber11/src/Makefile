#   top-level $(MAKE)file for the AmberTools suite of programs
AMBERTOOLS=AmberTools, version 1.3
CONFIG_FILE=config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI

ifeq ($(AMBERHOME),)
export AMBERHOME=$(PWD:/src=)
endif

ifeq ($(TARGET_PLATFORM),)
ifeq ($(findstring darwin,$(OSTYPE)),darwin)
export TARGET_PLATFORM=MAC
else
export TARGET_PLATFORM=MSW
endif
endif

include $(CONFIG_FILE)

install: serial

serial: configured_serial
	@echo "Starting installation of ${AMBERTOOLS} at `date`".
	mkdir -p $(LIBDIR)
#	utility routines and libraries:
#	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
#	(cd byacc && $(MAKE) install )
	(cd carpack && $(MAKE) install );
	(cd clapack && $(MAKE) $(LAPACK) )
	(cd cblas && $(MAKE) $(BLAS) )
	(cd f2c &&  $(MAKE) $(F2C)  )
	(cd c9x-complex && $(MAKE) $(C9XCOMPLEX) )
#	(cd etc && $(MAKE) install )
#	(cd protonate && $(MAKE) install )
#	(cd chamber && $(MAKE) install )
#	(cd pbsa && $(MAKE) install )

#	antechamber:
	mkdir -p $(BINDIR)
	(cd antechamber && $(MAKE) install )
	(cd sqm && $(MAKE) install )

#   ptraj:
#	(cd ptraj && $(MAKE) install )

#   miscellaneous:
#	(cd reduce && $(MAKE) install )

#   leap and gleap:
#	(cd leap && $(MAKE) install )
#	(cd gleap && $(MAKE) install )

#	nab:
#	(cd cifparse && $(MAKE) install )
#	(cd nab && $(MAKE) install )
#	(cd nss && $(MAKE) install )

#   mdgx:
#	(cd mdgx1 && $(MAKE) $(MDGX) )

	@echo ""
	@echo "Installation of ${AMBERTOOLS} is complete at `date`."
	@echo ""

nabonly:
#	utility routines and libraries:
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd byacc && $(MAKE) install )
	(cd carpack && $(MAKE) install );
	(cd clapack && $(MAKE) $(LAPACK) )
	(cd cblas && $(MAKE) $(BLAS) )
	(cd f2c &&  $(MAKE) $(F2C)  )
	(cd c9x-complex && $(MAKE) $(C9XCOMPLEX) )

#	nab:
	(cd cifparse && $(MAKE) install )
	(cd nab && $(MAKE) install )
	(cd nss && $(MAKE) install )

	@echo
	@echo "Completed installation of NAB"
	@echo

parallel: configured_parallel
	@echo "Starting installation of ${AMBERTOOLS} MPI at `date`".
#	(cd sff && $(MAKE) clean )
#	(cd nab && $(MAKE) install )
	(cd ptraj && $(MAKE) parallel )

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
		echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE) !" ;\
		exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial !" ;\
		echo "  Rerun $(CONFIG_COMMAND) without the -mpi flag." ;\
		exit 2 ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is of type serial, not parallel !" ;\
			echo "  Rerun $(CONFIG_COMMAND) with the -mpi flag." ;\
			exit 2 ;\
		} ;\
	)

clean:
#	(cd ucpp-1.3 && $(MAKE) clean )
#	(cd byacc && $(MAKE) clean )
#	(cd cifparse && $(MAKE) clean )
#	(cd nab && $(MAKE) clean )
#	(cd sff && $(MAKE) clean )
#	(cd nss && $(MAKE) clean )
#	(cd leap && $(MAKE) clean )
#	(cd gleap && $(MAKE) clean )
	(cd carpack && $(MAKE) clean )
	(cd cblas && $(MAKE) clean )
	(cd clapack && $(MAKE) clean )
	(cd blas && $(MAKE) clean )
	(cd lapack && $(MAKE) clean )
	(cd f2c && $(MAKE) clean )
#	(cd c9x-complex && $(MAKE) clean )
#	(cd etc && $(MAKE) clean )
#	(cd protonate && $(MAKE) clean )
#	(cd chamber && $(MAKE) clean )
#	(cd pbsa && $(MAKE) clean )
	(cd antechamber && $(MAKE) clean )
	(cd sqm && $(MAKE) clean )
	(cd lib && $(MAKE) clean )
#	(cd ptraj && $(MAKE) clean )
#	(cd reduce && $(MAKE) clean )
#	(cd mdgx1 && $(MAKE) clean )
#	(cd paramfit && $(MAKE) clean )

realclean: clean uninstall_libraries
	-(cd netcdf/src && $(MAKE) clean)
	-(cd netcdf/lib && rm -f *.a *.la)
	-(cd netcdf/include && rm -f *.mod netcdf.*)
	-(cd netcdf/bin && rm -f *)
	-(cd pnetcdf/src && $(MAKE) clean)
	-(cd pnetcdf/lib && rm -f *.a *.la)
	-(cd pnetcdf/include && rm -f *.mod netcdf.*)
	-(cd pnetcdf/bin && rm -f *)
	-(cd fftw-3.2.2 && $(MAKE) clean)

# This target realcleans, removes the configuration file and cleans the tests.
distclean: realclean
	-cd ../test && $(MAKE) distclean
	-/bin/rm -f config.h

# This target distcleans and removes all executables.
superclean: uninstall distclean

uninstall: uninstall_libraries
	(cd ucpp-1.3 && $(MAKE) uninstall )
	(cd byacc && $(MAKE) uninstall )
	(cd cifparse && $(MAKE) uninstall )
	(cd nab && $(MAKE) uninstall )
	(cd sff && $(MAKE) uninstall )
	(cd nss && $(MAKE) uninstall )
	(cd leap && $(MAKE) uninstall )
	(cd gleap && $(MAKE) uninstall )
	(cd carpack && $(MAKE) uninstall )
	(cd cblas && $(MAKE) uninstall )
	(cd clapack && $(MAKE) uninstall )
	(cd blas && $(MAKE) uninstall )
	(cd lapack && $(MAKE) uninstall )
	(cd f2c && $(MAKE) uninstall )
	(cd c9x-complex && $(MAKE) uninstall )
	(cd antechamber && $(MAKE) uninstall )
	(cd pbsa && $(MAKE) uninstall )
	(cd lib && $(MAKE) -f Makefile_at uninstall )
	(cd resp && $(MAKE) uninstall )
	(cd ptraj && $(MAKE) uninstall )
	(cd reduce && $(MAKE) uninstall )
#	(cd mdgx1 && $(MAKE) uninstall )

uninstall_libraries:
	(cd ../lib && rm -f libnab.a libsff.a cifparse.a libpdb.a carpack.a \
              cblas.a clapack.a blas.a lapack.a f2c.a libsym.a libfftw3.a)
