include ../config.h
SHELL=/bin/sh
SRCDIR=..
VPATH=

.c.o:
	$(CC) -c $(OCFLAGS) -o $@ $<

OBJS= abort_.o  d_atan.o  d_cosh.o  d_int.o   d_nint.o  d_sinh.o  \
    d_abs.o   d_atn2.o  d_dim.o   d_lg10.o  d_prod.o  d_sqrt.o  exit_.o \
    d_acos.o  d_exp.o   d_log.o   d_sign.o  d_tan.o   \
    d_asin.o  d_cos.o   d_mod.o   d_sin.o   d_tanh.o  \
    wsfe.o s_cmp.o s_copy.o s_paus.o s_rnge.o s_stop.o wsle.o \
    util.o wref.o wrtfmt.o sfe.o sue.o uio.o \
    open.o rdfmt.o rewind.o rsfe.o rsli.o ilnw.o inquire.o lread.o fmtlib.o \
    ftell_.o iio.o err.o close.o dfe.o dolio.o due.o endfile.o backspace.o \
    fmt.o sig_die.o pow_dd.o pow_di.o pow_ii.o lwrite.o i_len.o s_cat.o  \
	F77_aloc.o wsne.o typesize.o i_dnnt.o i_nint.o

install: $(OBJS)
	$(AR) f2c.a $(OBJS)
	$(RANLIB) f2c.a
	mv f2c.a $(LIBDIR)

skip:
	@echo "skipping compilation of F2C"

clean:
	-rm -f *.o f2c.a

uninstall:
	-rm -f $(LIBDIR)/f2c.a
