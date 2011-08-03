#
#************************************************************************
#                             AMBER                                    **
#                                                                      **
#                Copyright (c) 1986, 1991, 1995                        **
#             Regents of the University of California                  **
#                       All Rights Reserved.                           ** 
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#
#
SHELL=/bin/sh
VPATH=
PREC=DPREC
include ../config.h

OBJ=  \
sys.o wallclock.o random.o veclib.o amopen.o mexit.o nxtsec.o rgroup.o \
rfree.o matinv.o

#  default rules for Fortran and C compilation:

.f.o:   $<
	$(FPP) $< > _$<
	$(FC) -c $(FFLAGS) -o $@ _$<

.c.o:
	$(CC) -c $(CFLAGS) -o $@ $<

#  special rules for free-format routines:

random.o:   random.f
	$(FPP) random.f > _random.f
	$(FC) -c $(FFLAGS) $(FREEFORMAT_FLAG) -o $@ _random.f

veclib.o:   veclib.f
	$(FPP) veclib.f > _veclib.f
	$(FC) -c $(FFLAGS) $(FREEFORMAT_FLAG) -o $@ _veclib.f

#
all:
		@echo "(lib: make specific .o files, driven by other Makefiles)"

new2oldparm$(SFX):  new2oldparm.o nxtsec.o
		$(FC) -o new2oldparm$(SFX) new2oldparm.o nxtsec.o $(LDFLAGS) $(FLDFLAGS)


sys.a:   sys.o wallclock.o random.o veclib.o
		$(AR) sys.a sys.o wallclock.o random.o veclib.o
		$(RANLIB) sys.a

sys2.a:            $(OBJ)
		$(AR) sys2.a $(OBJ)
		$(RANLIB) sys2.a

install: new2oldparm$(SFX)
		mv new2oldparm$(SFX) $(BINDIR)
		/bin/rm -f _*.f
		
clean:
		rm -f malloc_test new2oldparm$(SFX) *.o _*.f *.a

uninstall:
	-rm -f $(BINDIR)/new2oldparm$(SFX)

