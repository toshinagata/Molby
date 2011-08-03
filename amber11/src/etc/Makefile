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
SHELL=/bin/sh
NSRC= 	nucgen.f ngfil.f
NOBJ=	nucgen.o ngfil.o
ASRC=   ambpdb.f forFIRSTtype.f forFIRSTbond.f forFIRSThbond.f \
                 forFIRSTteth.f forFIRSToutp.f forFIRSTsugar.f
AOBJ=   ambpdb.o forFIRSTtype.o forFIRSTbond.o forFIRSThbond.o \
                 forFIRSTteth.o forFIRSToutp.o forFIRSTsugar.o
include ../config.h

PROGS= ambpdb$(SFX) nucgen$(SFX) elsize$(SFX) molsurf$(SFX) resp$(SFX)
PERLPROGS= calcpka.pl cpinutil.pl CPin.pm lmodprmtop

all:		$(PROGS) 

nukit$(SFX):	nukit.o
	$(FC) -o nukit nukit.o $(LDFLAGS) $(FLDFLAGS)

nucgen$(SFX):	$(NOBJ) syslib
	$(FC) -o nucgen$(SFX) $(NOBJ) ../lib/sys2.a $(LDFLAGS) $(FLDFLAGS)

ambpdb$(SFX): 	$(AOBJ) syslib 
	$(FC) -o ambpdb$(SFX) $(AOBJ) ../lib/sys2.a $(LDFLAGS) $(FLDFLAGS)

elsize$(SFX): elsize.o
	$(CC) $(CFLAGS) -o elsize$(SFX) elsize.o $(LM)

molsurf$(SFX):  molsurf.o 
	$(CC) $(CFLAGS) -o molsurf$(SFX) molsurf.o $(LM)
molsurf.o:  molsurf.h

resp$(SFX):	lapack.o resp.o
	$(FC) $(LDFLAGS)  lapack.o resp.o -o resp$(SFX) $(FLDFLAGS)

syslib::
	cd ../lib && $(MAKE) sys2.a

install:	$(PROGS)
	/bin/mv $(PROGS) $(BINDIR)
	/bin/cp $(PERLPROGS) $(BINDIR)

clean:
	/bin/rm -f *.o _*.f $(PROGS)

uninstall:
	-for file in $(PROGS) $(PERLPROGS) ; do \
	rm -f $(BINDIR)/$${file} ; \
	done

#  default rules for Fortran and C compilation:

.f.o:   $<
	$(FPP) $< > _$<
	$(FC) -c $(FFLAGS) -o $@ _$<

.c.o:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<
