SHELL=/bin/sh
SRCDIR=..
VPATH=
include ../config.h

.c.o:
	$(CC) -c $(OCFLAGS) -o $@ $<

install: carpack.a

ARPACKLIB= dgetv0.o dlaqrb.o dnaitr.o dnapps.o \
dnaup2.o dnaupd.o dnconv.o dneigh.o \
dneupd.o dngets.o dsaitr.o dsapps.o dsaup2.o \
dsaupd.o dsconv.o dseigt.o dsesrt.o dseupd.o \
dsgets.o dsortc.o dsortr.o dstatn.o dstats.o dstqrb.o dsarpack.o \
wallclock.o dmout.o ivout.o dvout.o

SARPACKLIB= svout.o sgetv0.o slaqrb.o snaitr.o snapps.o \
snaup2.o snaupd.o snconv.o sneigh.o smout.o \
sneupd.o sngets.o ssaitr.o ssapps.o ssaup2.o \
ssaupd.o ssconv.o sseigt.o ssesrt.o sseupd.o \
ssgets.o ssortc.o ssortr.o sstatn.o sstats.o sstqrb.o

carpack.a: $(ARPACKLIB)
	$(AR) carpack.a $(ARPACKLIB)
	$(RANLIB) carpack.a
	mv carpack.a $(LIBDIR)

clean:
	-rm -f *.o carpack.a

uninstall:
	-rm -f $(LIBDIR)/carpack.a

$(ARPACKLIB):  f2c.h
