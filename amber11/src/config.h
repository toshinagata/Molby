#  Amber configuration file
#  Created by ./configure, and modified by hand
#  Toshi Nagata, 2010/01/25
###############################################################################

# (1)  Location of the installation

ifeq ($(TARGET_PLATFORM),MAC)
ISYSROOT=-isysroot /Developer/SDKs/MacOSX10.4u.sdk -mmacosx-version-min=10.4 -arch ppc -arch i386
endif
BINDIR=$(AMBERHOME)/bin
LIBDIR=$(AMBERHOME)/lib
INCDIR=$(AMBERHOME)/include
DATDIR=$(AMBERHOME)/dat
#NABHOME=$(AMBERHOME)/dat

###############################################################################


#  (2) If you want to search additional libraries by default, add them
#      to the FLIBS variable here.  (External libraries can also be linked into
#      NAB programs simply by including them on the command line; libraries
#      included in FLIBS are always searched.)

ifeq ($(TARGET_PLATFORM),MAC)
FLIBS = $(LIBDIR)/libsym.a $(LIBDIR)/carpack.a $(LIBDIR)/f2c.a -framework Accelerate
FLIBSF= $(LIBDIR)/carpack.a $(LIBDIR)/f2c.a $(FLIBS_EXTRA)
else
FLIBS= $(LIBDIR)/libsym.a $(LIBDIR)/carpack.a $(LIBDIR)/clapack.a $(LIBDIR)/cblas.a $(LIBDIR)/f2c.a $(LIBDIR)/libmc.a
FLIBSF= $(LIBDIR)/carpack.a $(LIBDIR)/lapack.a $(LIBDIR)/blas.a $(LIBDIR)/f2c.a $(LIBDIR)/libmc.a
endif
#FLIBS_PTRAJ= $(LIBDIR)/carpack.a $(LIBDIR)/f2c.a $(FLIBS_EXTRA)
#FLIBS_FFTW2=
###############################################################################

#  (3)  Modify any of the following if you need to change, e.g. to use gcc
#        rather than cc, etc.

SHELL=/bin/sh

#  Set the C compiler, etc. 

#          For GNU:  CC-->gcc; LEX-->flex; YACC-->bison -y -t;
#          Note: If your lexer is "really" flex, you need to set
#          LEX=flex below.  For example, on many linux distributions,
#          /usr/bin/lex is really just a pointer to /usr/bin/flex,
#          so LEX=flex is necessary.  In general, gcc seems to need
#          flex.

ifeq ($(TARGET_PLATFORM),MAC)
CC=gcc-4.0 $(ISYSROOT)
CXX=g++-4.0 $(ISYSROOT)
CPLUSPLUS=g++-4.0 $(ISYSROOT)
CFLAGS= -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ $(AMBERBUILDFLAGS)
OCFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ $(AMBERBUILDFLAGS)
#NABFLAGS=
LDFLAGS= $(AMBERBUILDFLAGS) -framework Accelerate
FLDFLAGS= -nodefaultlibs -lgfortran-static -lgcc -lc -lm -lSystem -lSystemStubs -lgfortranbegin
else
CC=gcc
CXX=g++
CPLUSPLUS=g++
CFLAGS= -DUSE_AMBER_C9XCOMPLEX -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE  $(AMBERBUILDFLAGS) -DWINDOWS=1
OCFLAGS=-O3 -DUSE_AMBER_C9XCOMPLEX -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE $(AMBERBUILDFLAGS) -DWINDOWS=1
NABFLAGS=
LDFLAGS=-Wl,--stack=0x01000000 $(AMBERBUILDFLAGS)
FLDFLAGS=
endif

#LEX=   flex
#YACC=  $(BINDIR)/yacc
AR=    ar rv
M4=    m4
RANLIB=ranlib
MAKE=make

#  Set the C-preprocessor.  Code for a small preprocessor is in
#    uccp-1.3;  it gets installed as $(BINDIR)/ucpp;
#    this can generally be used (maybe not on 64-bit machines like altix).

CPP=    $(BINDIR)/ucpp -l

#  These variables control whether we will use compiled versions of BLAS
#  and LAPACK (which are generally slower), or whether those libraries are
#  already available (presumably in an optimized form).

ifeq ($(TARGET_PLATFORM),MAC)
LAPACK=skip
BLAS=skip
else
LAPACK=install
BLAS=install
endif
F2C=install

#  These variables determine whether builtin versions of certain components
#  can be used, or whether we need to compile our own versions.

UCPP=install
ifeq ($(TARGET_PLATFORM),MAC)
C9XCOMPLEX=skip
else
C9XCOMPLEX=libmc.a
endif

#  For Windows/cygwin, set SFX to ".exe"; for Unix/Linux leave it empty:

ifeq ($(TARGET_PLATFORM),MAC)
SFX=
else
SFX=.exe
endif

#  Information about Fortran compilation:

ifeq ($(TARGET_PLATFORM),MAC)
FC=gfortran $(ISYSROOT)
else
FC=gfortran
endif
FFLAGS= -O0 $(LOCALFLAGS) $(AMBERBUILDFLAGS)
FOPTFLAGS= -O3 $(LOCALFLAGS) $(AMBERBUILDFLAGS)
FREEFORMAT_FLAG= -ffree-form
LM=-lm
FPP=cpp -traditional -P  -DNO_SANDER_DIVCON -DBINTRAJ 
FPPFLAGS=-P  -DNO_SANDER_DIVCON -DBINTRAJ 

BUILD_SLEAP=install_sleap
XHOME= /usr/X11R6
XLIBS= -L/usr/X11R6/lib
MAKE_XLEAP=install_xleap

LIBDIVCON=
INCDIVCON=
MODULEDIR=-I
TESTSANDERDIVCON=skipsanderDIVCON

NETCDF=netcdf.mod
NETCDFLIB=../netcdf/lib/libnetcdf.a
PNETCDF=no
PNETCDFLIB=

HASFC=yes
MDGX=no

COMPILER=gnu
MKL=
MKL_PROCESSOR=
