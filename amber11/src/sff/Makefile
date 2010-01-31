include ../config.h

.c.o:
	$(CC) -c -Dflex $(OCFLAGS) -o $@ $<

OBJS = binpos.o conjgrad.o lmodC.o memutil.o nblist.o newton.o nmode.o \
	prm.o rand2.o sasad.o sff.o time.o traceback.o xminC.o

install: libsff.a
	mv libsff.a $(LIBDIR)

libsff.a: $(OBJS)
	-rm -f libsff.a
	$(AR) libsff.a $(OBJS)	
	$(RANLIB) libsff.a

uninstall:
	-rm -f $(LIBDIR)/libsff.a

clean:
	-/bin/rm -f $(OBJS) libsff.a lex.mm_options.c

lex.mm_options.c: mm_options.l
	$(LEX) -t mm_options.l | sed -e 's/{stdout}/{DONTchangeSTDOUThere}/' -e 's/yy/mmo/g' -e 's/stdout/nabout/' -e 's/DONTchangeSTDOUThere/stdout/' > lex.mm_options.c

#  dependencies (hand-created for now):

binpos.o: memutil.h

conjgrad.o: memutil.h

lmodC.o: sff.h

memutil.o: sff.h

nblist.o: sff.h

newton.o: sff.h memutil.h

nmode.o: sff.h memutil.h

prm.o:   sff.h errormsg.h

rand2.c:  sff.h

sasad.o: sff.h

sff.o: sff.h memutil.h \
       debug.h errormsg.h traceback.h eff_box.c dssp.c eff.c intersect.c \
       gbsa.c egbr6.c sff2.c mask.c rattle.c lex.mm_options.c hcp.c

time.o:  sff.h

traceback.o:  errormsg.h

xminC.o:  sff.h
