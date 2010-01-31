include ../config.h
SHELL=/bin/sh
SRCDIR=..
VPATH=

.c.o:
	$(CC) -c $(OCFLAGS) -o $@ $<


OBJS= dsyev.o dsyevd.o dormtr.o dormql.o dormqr.o dlansy.o dsytrd.o dorgtr.o dsteqr.o \
	dsterf.o dlascl.o dlamch.o ilaenv.o dlassq.o dsytd2.o dorm2l.o \
	dlatrd.o dorgql.o dorgqr.o dlasrt.o dlaset.o dlanst.o dlaev2.o \
	dlasr.o dlapy2.o dlartg.o dlae2.o dorg2r.o dlarft.o dlarfb.o \
	dorg2l.o dlarfg.o dlarf.o dspev.o dlansp.o dsptrd.o dopgtr.o \
	dposv.o dppsv.o dspsv.o dgeev.o dsptrs.o dsptrf.o dpptrs.o dpptrf.o \
	dgebak.o dtrevc.o dlacpy.o dhseqr.o dorghr.o dgehrd.o dgebal.o \
	dlange.o dlabad.o dlahqr.o dlarfx.o dlanhs.o dlaln2.o dgehd2.o \
	dlahrd.o dlanv2.o dladiv.o dgerq2.o dorm2r.o dlarnv.o dgeqr2.o \
	dlaruv.o dspevd.o dstedc.o dopmtr.o dlaed0.o dlaed1.o dlaed7.o \
	dlaeda.o dlaed8.o dlaed9.o dlamrg.o dlaed2.o dlaed3.o dsygv.o \
	dsygst.o dsygs2.o dsyevr.o dstebz.o dstein.o dstemr.o dspevx.o \
	dlaed4.o dlaed5.o dlaed6.o dpotrf.o dpotrs.o dpotf2.o ieeeck.o \
	dsyevd.o dormtr.o dormqr.o dormql.o dorm2l.o dlarrc.o dlarre.o \
	dlarrj.o dlarrr.o dlarrv.o dlaebz.o dlagtf.o dlagts.o dlarrd.o \
	dlarrf.o dlarrk.o dlar1v.o dlasq2.o dlarra.o dlarrb.o disnan.o \
	dlaneg.o dlasq3.o dlaisnan.o dlasq4.o dlasq5.o dlasq6.o dsysv.o \
	dgelsd.o dgelsy.o dgeqrf.o dgesvd.o dtzrzf.o dorgbr.o dorglq.o dormbr.o \
	dormlq.o dormrz.o dbdsqr.o dgebrd.o dlaic1.o dgelqf.o dgeqp3.o dlalsd.o \
	dsytrf.o dsytrs.o dorgl2.o dorml2.o dormr3.o dgebd2.o dlabrd.o dgelq2.o \
	dlalsa.o dlaqp2.o dsytf2.o dlasda.o dlasdq.o dlasq1.o dlaqps.o dlasv2.o \
	dlarzb.o dlarzt.o dlasyf.o dlatrz.o dlas2.o dlals0.o dlasd6.o dlarfp.o \
	dlasdt.o dlarz.o dlasd7.o dlasd8.o dlasd4.o dlasd5.o dsyevx.o \
	dgetrf.o dgetrs.o dgetf2.o dlaswp.o iparmq.o iladlr.o iladlc.o dlahr2.o \
	dlaqr0.o dlaqr1.o dlaqr2.o dlaqr3.o dlaqr4.o dlaqr5.o dormhr.o dtrexc.o \
	dlaexc.o dlasy2.o

install: $(OBJS)
	$(AR) clapack.a $(OBJS)
	$(RANLIB) clapack.a
	mv clapack.a $(LIBDIR)

skip:
	@echo "skipping compilation of LAPACK"

clean:
	-/bin/rm -f *.o clapack.a

dlaruv.o:   dlaruv.c
	$(CC) -c -O0 -o dlaruv.o dlaruv.c

uninstall:
	rm -f $(LIBDIR)/clapack.a

$(OBJS):  f2c.h
