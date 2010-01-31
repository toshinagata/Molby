include ../config.h
SHELL=/bin/sh
SRCDIR=..
VPATH=

.c.o:
	$(CC) -c $(OCFLAGS) -o $@ $<

OBJS= dasum.o daxpy.o dcopy.o ddot.o \
         dgbmv.o dgemm.o dgemv.o dger.o dnrm2.o drot.o \
         drotg.o dsbmv.o dscal.o dspmv.o dspr.o dspr2.o \
         dswap.o dsymm.o dsymv.o dsyr.o dsyr2.o dsyr2k.o \
         dsyrk.o dtbmv.o dtbsv.o dtpmv.o dtpsv.o dtrmm.o \
         dtrmv.o dtrsm.o dtrsv.o idamax.o lsame.o xerbla.o

install: $(OBJS)
	$(AR) cblas.a $(OBJS)
	$(RANLIB) cblas.a
	mv cblas.a $(LIBDIR)

skip:
	@echo "skipping compilation of BLAS"

xerbla: xerbla.o
	mv xerbla.o $(LIBDIR)

clean:
	-rm -f *.o cblas.a

uninstall:
	-rm -f $(LIBDIR)/cblas.a

$(OBJS): f2c.h
