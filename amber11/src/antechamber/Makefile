# makefile for antechamber

include ../config.h

PROGS = \
	am1bcc$(SFX) \
	acdoctor$(SFX) \
	antechamber$(SFX) \
	atomtype$(SFX) \
	bondtype$(SFX) \
	charmmgen$(SFX) \
	database$(SFX) \
	espgen$(SFX) \
	parmcal$(SFX) \
	parmchk$(SFX) \
	prepgen$(SFX) \
	residuegen$(SFX) \
	respgen$(SFX) \
	translate$(SFX) \
	top2mol2$(SFX)

TESTS = test-penalties test-int-part test-combination

#############################################################################

all: $(PROGS)

tests: $(TESTS)

install: all
	mv -f $(PROGS) $(BINDIR)
        
clean:
	-rm -f *.o $(PROGS) $(TESTS) bondtype.orig$(SFX)

uninstall:
	-for file in $(PROGS) ; do \
	rm -f $(BINDIR)/$${file} ; \
	done

#############################################################################

#_____________________________________________________________________________
# Dependencies, for ./makedepend:

acdoctor.o: common.h
acdoctor.o: define.h
acdoctor.o: atom.h
acdoctor.o: utility.c
acdoctor.o: common.c
acdoctor.o: ring.c
acdoctor.o: rotate.c
acdoctor.o: ac.c
acdoctor.o: charmm.c
acdoctor.o: mol2.c
acdoctor.o: mopcrt.c
acdoctor.o: divcrt.c
acdoctor.o: mopint.c
acdoctor.o: mopout.c
acdoctor.o: divout.c
acdoctor.o: gcrt.c
acdoctor.o: gzmat.c
acdoctor.o: gout.c
acdoctor.o: pdb.c
acdoctor.o: csd.c
acdoctor.o: mdl.c
acdoctor.o: alc.c
acdoctor.o: hin.c
acdoctor.o: prep.c
acdoctor.o: rst.c
acdoctor.o: jzmat.c
acdoctor.o: jcrt.c
acdoctor.o: jout.c
am1bcc.o: common.h
am1bcc.o: define.h
am1bcc.o: atom.h
am1bcc.o: utility.c
am1bcc.o: common.c
am1bcc.o: rotate.c
am1bcc.o: ac.c
am1bcc.o: pdb.c
antechamber.o: 
antechamber.o: common.h
antechamber.o: define.h
antechamber.o: atom.h
antechamber.o: utility.c
antechamber.o: common.c
antechamber.o: ring.c
antechamber.o: rotate.c
antechamber.o: ac.c
antechamber.o: charmm.c
antechamber.o: mol2.c
antechamber.o: mopcrt.c
antechamber.o: divcrt.c
antechamber.o: sqmcrt.c
antechamber.o: mopint.c
antechamber.o: mopout.c
antechamber.o: divout.c
antechamber.o: gcrt.c
antechamber.o: gzmat.c
antechamber.o: gout.c
antechamber.o: pdb.c
antechamber.o: csd.c
antechamber.o: mdl.c
antechamber.o: alc.c
antechamber.o: hin.c
antechamber.o: prep.c
antechamber.o: rst.c
antechamber.o: jzmat.c
antechamber.o: jcrt.c
antechamber.o: jout.c
antechamber.o: charge.c
atomtype.o: common.h
atomtype.o: define.h
atomtype.o: atom.h
atomtype.o: utility.c
atomtype.o: common.c
atomtype.o: ring.c
atomtype.o: rotate.c
atomtype.o: ac.c
atomtype.o: mol2.c
bondtype.o: 
bondtype.o: ncsu-penalties.h
bondtype.o: common.h
bondtype.o: define.h
bondtype.o: atom.h
bondtype.o: utility.c
bondtype.o: common.c
bondtype.o: ring.c
bondtype.o: ac.c
bondtype.o: mol2.c
charmmgen.o: common.h
charmmgen.o: define.h
charmmgen.o: atom.h
charmmgen.o: utility.c
charmmgen.o: common.c
charmmgen.o: rotate.c
charmmgen.o: ac.c
charmmgen.o: mol2.c
charmmgen.o: prep.c
crdgrow.o: common.h
crdgrow.o: define.h
crdgrow.o: atom.h
crdgrow.o: utility.c
crdgrow.o: common.c
crdgrow.o: rotate.c
crdgrow.o: pdb.c
crdgrow.o: ac.c
crdgrow.o: prep.c
database.o: 
database.o: utility.c
delphigen.o: common.h
delphigen.o: define.h
delphigen.o: atom.h
delphigen.o: common.c
delphigen.o: ac.c
delphigen.o: pdb.c
espgen.o: common.h
espgen.o: define.h
espgen.o: atom.h
match.o: 
ncsu-penalties.o: 
ncsu-penalties.o: ncsu-penalties.h
parmcal.o: common.h
parmcal.o: define.h
parmcal.o: utility.c
parmchk.o: common.h
parmchk.o: define.h
parmchk.o: atom.h
parmchk.o: utility.c
parmchk.o: common.c
parmchk.o: rotate.c
parmchk.o: ac.c
parmchk.o: mol2.c
parmchk.o: prep.c
parmjoin.o: common.h
parmjoin.o: define.h
parmjoin.o: utility.c
prepgen.o: common.h
prepgen.o: define.h
prepgen.o: atom.h
prepgen.o: utility.c
prepgen.o: common.c
prepgen.o: rotate.c
prepgen.o: ring.c
prepgen.o: ac.c
residuegen.o: common.h
residuegen.o: define.h
residuegen.o: atom.h
residuegen.o: utility.c
residuegen.o: common.c
residuegen.o: ac.c
respgen.o: common.h
respgen.o: define.h
respgen.o: atom.h
respgen.o: utility.c
respgen.o: common.c
respgen.o: ac.c
rotate.o: common.h
smiles.o: common.h
smiles.o: define.h
smiles.o: atom.h
smiles.o: globals.h
test-penalties.o: 
test-penalties.o: ncsu-penalties.h
top2mol2.o: 
translate.o: common.h
translate.o: define.h
translate.o: atom.h
translate.o: utility.c
translate.o: common.c
translate.o: rotate.c
translate.o: ac.c
translate.o: pdb.c
translate.o: mol2.c
translate.o: prep.c
translate.o: lsfit.c

.c.o:
	$(CC) -c $(CFLAGS)  -o $@ $<

#_____________________________________________________________________________
# Programs

am1bcc$(SFX): am1bcc.o
	$(CC) $(LDFLAGS) am1bcc.o $(LM) -o am1bcc$(SFX)  

acdoctor$(SFX): acdoctor.o
	$(CC) $(LDFLAGS) acdoctor.o $(LM) -o acdoctor$(SFX)  

antechamber$(SFX): antechamber.o
	$(CC) $(LDFLAGS) antechamber.o $(LM) -o antechamber$(SFX)  

atomtype$(SFX): atomtype.o
	$(CC) $(LDFLAGS) atomtype.o $(LM) -o atomtype$(SFX)  

bondtype$(SFX): bondtype.o ncsu-penalties.o
	$(CC) $(LDFLAGS) ncsu-penalties.o bondtype.o $(LM) -o bondtype$(SFX)  

bondtype.orig.o: bondtype.cpp
	$(CXX) -c $(CFLAGS) -o bondtype.orig.o bondtype.cpp

bondtype.orig$(SFX): bondtype.orig.o
	$(CXX) $(LDFLAGS) bondtype.orig.o $(LM) -o bondtype.orig$(SFX)  

charmmgen$(SFX): charmmgen.o
	$(CC) $(LDFLAGS) charmmgen.o $(LM) -o charmmgen$(SFX)

crdgrow$(SFX): crdgrow.o
	$(CC) $(LDFLAGS) crdgrow.o $(LM) -o crdgrow$(SFX)

database$(SFX): database.o
	$(CC) $(LDFLAGS) database.o $(LM) -o database$(SFX)

delphigen$(SFX): delphigen.o
	$(CC) $(LDFLAGS) delphigen.o -o delphigen$(SFX)

espgen$(SFX): espgen.o
	$(CC) $(LDFLAGS) espgen.o -o espgen$(SFX)

parmcal$(SFX): parmcal.o 
	$(CC) $(LDFLAGS) parmcal.o $(LM) -o parmcal$(SFX)

parmchk$(SFX): parmchk.o 
	$(CC) $(LDFLAGS) parmchk.o $(LM) -o parmchk$(SFX)

parmjoin$(SFX): parmjoin.o 
	$(CC) $(LDFLAGS) parmjoin.o -o parmjoin$(SFX)

prepgen$(SFX): prepgen.o 
	$(CC) $(LDFLAGS) prepgen.o $(LM) -o prepgen$(SFX)

residuegen$(SFX): residuegen.o 
	$(CC) $(LDFLAGS) residuegen.o $(LM) -o residuegen$(SFX)

respgen$(SFX): respgen.o 
	$(CC) $(LDFLAGS) respgen.o $(LM) -o respgen$(SFX)

translate$(SFX): translate.o 
	$(CC) $(LDFLAGS) translate.o $(LM) -o translate$(SFX)

top2mol2$(SFX): top2mol2.o
	$(CC) $(LDFLAGS) top2mol2.o $(LM) -o top2mol2$(SFX)

test-penalties: ncsu-penalties.o test-penalties.o
	$(CC) -o test-penalties $(LDFLAGS) ncsu-penalties.o test-penalties.o

test-int-part: ncsu-penalties.h ncsu-penalties.c
	$(CC) -o test-int-part ncsu-penalties.c $(CFLAGS) $(LDFLAGS) \
   -DENABLE_INT_PART_MAIN=1

test-combination: ncsu-penalties.h ncsu-penalties.c
	$(CC) -o test-combination ncsu-penalties.c $(CFLAGS) $(LDFLAGS) \
   -DENABLE_COMBINATION_MAIN=1

# DO NOT DELETE
