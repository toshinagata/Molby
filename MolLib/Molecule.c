/*
 *  Molecule.c
 *
 *  Created by Toshi Nagata on 06/03/11.
 *  Copyright 2006 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#include "MolLib.h"
#include <string.h>
#include <stddef.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

#include "Missing.h"
#include "Dcd.h"
#include "MD/MDCore.h"
#include "MD/MDPressure.h"

static Molecule *sMoleculeRoot = NULL;
static int sMoleculeUntitledCount = 0;

Int gSizeOfAtomRecord = sizeof(Atom);

/*  These are the pasteboard data type. Since the internal representation of the
	pasteboard data includes binary data that may be dependent on the software version,
    the revision number is appended to these strings on startup (See MyApp::OnInit())  */
char *gMoleculePasteboardType = "Molecule";
char *gParameterPasteboardType = "Parameter";

#pragma mark ====== Utility function ======

int
strlen_limit(const char *s, int limit)
{
	int len;
	for (len = 0; *s != 0 && (limit < 0 || len < limit); s++, len++);
	return len;
}

#pragma mark ======  Atom handling  ======

static Atom *
s_AtomDuplicate(Atom *dst, const Atom *src, Int copy_frame)
{
	if (dst == NULL) {
		dst = (Atom *)malloc(gSizeOfAtomRecord);
		if (dst == NULL)
			return NULL;
	}
	memmove(dst, src, gSizeOfAtomRecord);
	if (src->aniso != NULL) {
		dst->aniso = (Aniso *)malloc(sizeof(Aniso));
		if (dst->aniso != NULL)
			memmove(dst->aniso, src->aniso, sizeof(Aniso));
	}
	if (src->frames != NULL && copy_frame) {
		dst->frames = (Vector *)malloc(sizeof(Vector) * src->nframes);
		if (dst->frames != NULL) {
			memmove(dst->frames, src->frames, sizeof(Vector) * src->nframes);
			dst->nframes = src->nframes;
		} else {
			dst->nframes = 0;
		}
	}
	if (src->connect.count > ATOM_CONNECT_LIMIT) {
		dst->connect.u.ptr = NULL;
		dst->connect.count = 0;
		NewArray(&(dst->connect.u.ptr), &(dst->connect.count), sizeof(Int), src->connect.count);
		memmove(dst->connect.u.ptr, src->connect.u.ptr, sizeof(Int) * src->connect.count);
	}
	if (src->anchor != NULL) {
		dst->anchor = (PiAnchor *)malloc(sizeof(PiAnchor));
		if (dst->anchor != NULL)
			memmove(dst->anchor, src->anchor, sizeof(PiAnchor));
		if (dst->anchor->connect.count > ATOM_CONNECT_LIMIT) {
			dst->anchor->connect.u.ptr = NULL;
			dst->anchor->connect.count = 0;
			NewArray(&(dst->anchor->connect.u.ptr), &(dst->anchor->connect.count), sizeof(Int), src->anchor->connect.count);
			memmove(dst->anchor->connect.u.ptr, src->anchor->connect.u.ptr, sizeof(Int) * src->anchor->connect.count);
		}
		if (dst->anchor->ncoeffs > 0) {
			NewArray(&(dst->anchor->coeffs), &(dst->anchor->ncoeffs), sizeof(Double), src->anchor->ncoeffs);
			memmove(dst->anchor->coeffs, src->anchor->coeffs, sizeof(Double) * src->anchor->ncoeffs);
		}
	}
	return dst;
}

Atom *
AtomDuplicate(Atom *dst, const Atom *src)
{
	return s_AtomDuplicate(dst, src, 1);
}

Atom *
AtomDuplicateNoFrame(Atom *dst, const Atom *src)
{
	return s_AtomDuplicate(dst, src, 0);
}

void
AtomClean(Atom *ap)
{
	if (ap->aniso != NULL) {
		free(ap->aniso);
		ap->aniso = NULL;
	}
	if (ap->frames != NULL) {
		free(ap->frames);
		ap->frames = NULL;
		ap->nframes = 0;
	}
	if (ap->connect.count > ATOM_CONNECT_LIMIT) {
		ap->connect.count = 0;
		free(ap->connect.u.ptr);
		ap->connect.u.ptr = NULL;
	}
}

void
CubeRelease(Cube *cp)
{
	if (cp != NULL) {
		if (cp->dp != NULL)
			free(cp->dp);
		free(cp);
	}
}

void
BasisSetRelease(BasisSet *bset)
{
	int i;
	if (bset == NULL)
		return;
	if (bset->shells != NULL)
		free(bset->shells);
	if (bset->priminfos != NULL)
		free(bset->priminfos);
	if (bset->mo != NULL)
		free(bset->mo);
	if (bset->cns != NULL)
		free(bset->cns);
	if (bset->moenergies != NULL)
		free(bset->moenergies);
	if (bset->scfdensities != NULL)
		free(bset->scfdensities);
/*	if (bset->pos != NULL)
		free(bset->pos); */
	if (bset->nuccharges != NULL)
		free(bset->nuccharges);
	if (bset->cubes != NULL) {
		for (i = 0; i < bset->ncubes; i++) {
			CubeRelease(bset->cubes[i]);
		}
		free(bset->cubes);
	}
	free(bset);
}

Int *
AtomConnectData(AtomConnect *ac)
{
	if (ac == NULL)
		return NULL;
	return ATOM_CONNECT_PTR(ac);
}

void
AtomConnectResize(AtomConnect *ac, Int nconnects)
{
	Int *p;
	if (ac == NULL)
		return;
	if (nconnects <= ATOM_CONNECT_LIMIT) {
		if (ac->count > ATOM_CONNECT_LIMIT) {
			p = ac->u.ptr;
			memmove(ac->u.data, p, sizeof(Int) * nconnects);
			free(p);
		}
	} else {
		if (ac->count <= ATOM_CONNECT_LIMIT) {
			p = NULL;
			ac->count = 0;
			NewArray(&p, &(ac->count), sizeof(Int), nconnects);
			memmove(p, ac->u.data, sizeof(Int) * ac->count);
			ac->u.ptr = p;
		} else if (ac->count < nconnects) {
			/*  Reallocate  */
			AssignArray(&(ac->u.ptr), &(ac->count), sizeof(Int), nconnects - 1, NULL);
		}
	}
	ac->count = nconnects;
}

void
AtomConnectInsertEntry(AtomConnect *ac, Int idx, Int connect)
{
	Int n, *p;
	if (ac == NULL)
		return;
	if (idx > ac->count)
		idx = ac->count;
	else if (idx < 0) {
		/*  Insert after the last component that is smaller than connect
		    (i.e. keep them sorted)  */
		p = ATOM_CONNECT_PTR(ac);
		for (idx = 0; idx < ac->count; idx++) {
			if (p[idx] >= connect)
				break;
		}
	}
	AtomConnectResize(ac, ac->count + 1);
	n = ac->count - idx - 1;  /*  Number of entries to be moved towards the end  */
	p = ATOM_CONNECT_PTR(ac);
	if (n > 0) {
		memmove(p + idx + 1, p + idx, sizeof(Int) * n);
	}
	p[idx] = connect;
}

void
AtomConnectDeleteEntry(AtomConnect *ac, Int idx)
{
	Int n, *p;
	if (ac == NULL)
		return;
	if (idx < 0 || idx >= ac->count)
		return;
	n = ac->count - idx - 1;  /*  Number of entries to be moved towards the top  */
	p = ATOM_CONNECT_PTR(ac);
	if (n > 0) {
		memmove(p + idx, p + idx + 1, sizeof(Int) * n);
	}
	AtomConnectResize(ac, ac->count - 1);
}

int
AtomConnectHasEntry(AtomConnect *ac, Int ent)
{
	Int n, *p;
	if (ac == NULL)
		return 0;
	p = ATOM_CONNECT_PTR(ac);
	for (n = 0; n < ac->count; n++) {
		if (ent == p[n])
			return 1;
	}
	return 0;
}

#pragma mark ====== Accessor types ======

MolEnumerable *
MolEnumerableNew(Molecule *mol, int kind)
{
	MolEnumerable *mseq = (MolEnumerable *)calloc(sizeof(MolEnumerable), 1);
	if (mseq != NULL) {
		mseq->mol = MoleculeRetain(mol);
		mseq->kind = kind;
	}
	return mseq;
}

void
MolEnumerableRelease(MolEnumerable *mseq)
{
	if (mseq != NULL) {
		MoleculeRelease(mseq->mol);
		free(mseq);
	}
}

AtomRef *
AtomRefNew(Molecule *mol, int idx)
{
	AtomRef *aref = (AtomRef *)calloc(sizeof(AtomRef), 1);
	if (aref != NULL) {
		aref->mol = MoleculeRetain(mol);
		aref->idx = idx;
	}
	return aref;
}

void
AtomRefRelease(AtomRef *aref)
{
	if (aref != NULL) {
		MoleculeRelease(aref->mol);
		free(aref);
	}
}

#pragma mark ====== Creation of molecules ======

Molecule *
MoleculeNew(void)
{
	char name[40];
	Molecule *mp = (Molecule *)calloc(sizeof(Molecule), 1);
	if (mp == NULL)
		Panic("Cannot allocate new molecule record");
	snprintf(name, sizeof name, "Untitled %d", sMoleculeUntitledCount++);
	ObjectInit((Object *)mp, (Object **)&sMoleculeRoot, name);
	mp->mview = MainView_new();
	mp->mview->mol = mp;
	return mp;
}

Molecule *
MoleculeNewWithName(const char *name)
{
	Molecule *mp = MoleculeNew();
	MoleculeSetName(mp, name);
	return mp;
}

Molecule *
MoleculeInitWithAtoms(Molecule *mp, const Atom *atoms, int natoms)
{
	int i;
	if (mp == NULL)
		mp = MoleculeNew();
	if (natoms == 0)
		return mp;
	if (NewArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, natoms) == NULL)
		Panic("Cannot allocate memory for atoms");
	for (i = 0; i < natoms; i++)
		AtomDuplicate(mp->atoms + i, atoms + i);
	mp->nframes = -1;  /*  Should be recalculated later  */
	return mp;
}

Molecule *
MoleculeInitWithMolecule(Molecule *mp2, Molecule *mp)
{
	int i, n;
	MoleculeFlushFrames(mp);
	MoleculeInitWithAtoms(mp2, mp->atoms, mp->natoms);
	if (mp->nbonds > 0) {
		if (NewArray(&mp2->bonds, &mp2->nbonds, sizeof(Int)*2, mp->nbonds) == NULL)
			goto error;
		memmove(mp2->bonds, mp->bonds, sizeof(Int) * 2 * mp->nbonds);
	}
	if (mp->nangles > 0) {
		if (NewArray(&mp2->angles, &mp2->nangles, sizeof(Int)*3, mp->nangles) == NULL)
			goto error;
		memmove(mp2->angles, mp->angles, sizeof(Int) * 3 * mp->nangles);
	}
	if (mp->ndihedrals > 0) {
		if (NewArray(&mp2->dihedrals, &mp2->ndihedrals, sizeof(Int)*4, mp->ndihedrals) == NULL)
			goto error;
		memmove(mp2->dihedrals, mp->dihedrals, sizeof(Int) * 4 * mp->ndihedrals);
	}
	if (mp->nimpropers > 0) {
		if (NewArray(&mp2->impropers, &mp2->nimpropers, sizeof(Int)*4, mp->nimpropers) == NULL)
			goto error;
		memmove(mp2->impropers, mp->impropers, sizeof(Int) * 4 * mp->nimpropers);
	}
	if (mp->nresidues > 0) {
		if (NewArray(&mp2->residues, &mp2->nresidues, sizeof(mp->residues[0]), mp->nresidues) == NULL)
			goto error;
		memmove(mp2->residues, mp->residues, sizeof(mp->residues[0]) * mp->nresidues);
	}
	if (mp->cell != NULL) {
		mp2->cell = (XtalCell *)calloc(sizeof(XtalCell), 1);
		memmove(mp2->cell, mp->cell, sizeof(XtalCell));
	}
	if (mp->nsyms > 0) {
		NewArray(&(mp2->syms), &(mp2->nsyms), sizeof(Transform), mp->nsyms);
		memmove(mp2->syms, mp->syms, sizeof(Transform) * mp2->nsyms);
	}

	/*	mp2->useFlexibleCell = mp->useFlexibleCell; */
	if (mp->nframe_cells > 0) {
		if (NewArray(&mp2->frame_cells, &mp2->nframe_cells, sizeof(Vector) * 4, mp->nframe_cells) == NULL)
			goto error;
		memmove(mp2->frame_cells, mp->frame_cells, sizeof(Vector) * 4 * mp->nframe_cells);
	}
	
	if (mp->nmolprops > 0) {
		if (NewArray(&mp2->molprops, &mp2->nmolprops, sizeof(MolProp), mp->nmolprops) == NULL)
			goto error;
		n = MoleculeGetNumberOfFrames(mp);
		for (i = 0; i < mp2->nmolprops; i++) {
			mp2->molprops[i].propname = strdup(mp->molprops[i].propname);
			mp2->molprops[i].propvals = (Double *)malloc(sizeof(Double) * n);
			memcpy(mp2->molprops[i].propvals, mp->molprops[i].propvals, sizeof(Double) * n);
		}
	}
	
	/* FIXME: should bset (basis set info) and elpot be duplicated or not?  */

	if (mp->par != NULL)
		mp2->par = ParameterDuplicate(mp->par);
	if (mp->arena != NULL) {
		md_arena_new(mp2);
		md_arena_init_from_arena(mp2->arena, mp->arena);
	}
	
	return mp2;
  error:
	Panic("Cannot allocate memory for duplicate molecule");
	return NULL;  /*  Not reached  */
}

/*  Assign a unique name to this parameter record  */
void
MoleculeSetName(Molecule *mp, const char *name)
{
	ObjectSetName((Object *)mp, name, (Object *)sMoleculeRoot);
}

const char *
MoleculeGetName(Molecule *mp)
{
	return ObjectGetName((Object *)mp);
}

Molecule *
MoleculeWithName(const char *name)
{
	return (Molecule *)ObjectWithName(name, (Object *)sMoleculeRoot);
}

void
MoleculeSetPath(Molecule *mol, const char *fname)
{
	char *buf, *cwd;
	if (mol == NULL || fname == NULL)
		return;
	if (fname[0] == '/' || (isalpha(fname[0]) && fname[1] == ':')) {
		/*  Full path  */
		buf = strdup(fname);
	} else {
		cwd = getcwd(NULL, 0);
		asprintf(&buf, "%s/%s", cwd, fname);
		free(cwd);
	}
	if (mol->path != NULL) {
		if (strcmp(mol->path, buf) == 0) {
			/*  No change  */
			free(buf);
			return;
		}
		free((void *)(mol->path));
	}
	mol->path = buf;
	if (mol->arena != NULL) {
		md_close_output_files(mol->arena);
	}
}

const char *
MoleculeGetPath(Molecule *mol)
{
	if (mol == NULL)
		return NULL;
	return mol->path;
}

Molecule *
MoleculeRetain(Molecule *mp)
{
	ObjectIncrRefCount((Object *)mp);
	MoleculeRetainExternalObj(mp);
	return mp;
}

void
MoleculeClear(Molecule *mp)
{
	int i;
	if (mp == NULL)
		return;
	if (mp->arena != NULL) {
		md_arena_set_molecule(mp->arena, NULL);
		mp->arena = NULL;
	}
	if (mp->par != NULL) {
		ParameterRelease(mp->par);
		mp->par = NULL;
	}
	if (mp->atoms != NULL) {
		for (i = 0; i < mp->natoms; i++)
			AtomClean(mp->atoms + i);
		free(mp->atoms);
		mp->atoms = NULL;
		mp->natoms = 0;
	}
	if (mp->bonds != NULL) {
		free(mp->bonds);
		mp->bonds = NULL;
		mp->nbonds = 0;
	}
	if (mp->angles != NULL) {
		free(mp->angles);
		mp->angles = NULL;
		mp->nangles = 0;
	}
	if (mp->dihedrals != NULL) {
		free(mp->dihedrals);
		mp->dihedrals = NULL;
		mp->ndihedrals = 0;
	}
	if (mp->impropers != NULL) {
		free(mp->impropers);
		mp->impropers = NULL;
		mp->nimpropers = 0;
	}
	if (mp->residues != NULL) {
		free(mp->residues);
		mp->residues = NULL;
		mp->nresidues = 0;
	}
	if (mp->cell != NULL) {
		free(mp->cell);
		mp->cell = NULL;
	}
	if (mp->syms != NULL) {
		free(mp->syms);
		mp->syms = NULL;
		mp->nsyms = 0;
	}
	if (mp->selection != NULL) {
		IntGroupRelease(mp->selection);
		mp->selection = NULL;
	}
	if (mp->frame_cells != NULL) {
		free(mp->frame_cells);
		mp->frame_cells = NULL;
		mp->nframe_cells = 0;
	}
	if (mp->bset != NULL) {
		BasisSetRelease(mp->bset);
		mp->bset = NULL;
	}
	if (mp->mcube != NULL) {
		MoleculeDeallocateMCube(mp->mcube);
		mp->mcube = NULL;
	}
	if (mp->molprops != NULL) {
		for (i = 0; i < mp->nmolprops; i++) {
			free(mp->molprops[i].propname);
			free(mp->molprops[i].propvals);
		}
		free(mp->molprops);
		mp->molprops = NULL;
		mp->nmolprops = 0;
	}
	if (mp->par != NULL) {
		ParameterRelease(mp->par);
		mp->par = NULL;
	}
	if (mp->elpots != NULL) {
		free(mp->elpots);
		mp->elpots = NULL;
		mp->nelpots = 0;
	}
	if (mp->path != NULL) {
		free((void *)mp->path);
		mp->path = NULL;
	}
}

void
MoleculeRelease(Molecule *mp)
{
	if (mp == NULL)
		return;
	MoleculeReleaseExternalObj(mp);
	if (ObjectDecrRefCount((Object *)mp) == 0) {
		MoleculeClear(mp);
		mp->mview->mol = NULL;
		MainView_release(mp->mview);
		ObjectDealloc((Object *)mp, (Object **)&sMoleculeRoot);
	}
}

void
MoleculeExchange(Molecule *mp1, Molecule *mp2)
{
	Molecule mp_temp;
	struct MainView *mview1, *mview2;
	struct MDArena *arena1, *arena2;
	/*  mview and arena must be kept as they are  */
	mview1 = mp1->mview;
	mview2 = mp2->mview;
	arena1 = mp1->arena;
	arena2 = mp2->arena;
	/*  'natoms' is the first member to be copied  */
	int ofs = offsetof(Molecule, natoms);
	memmove((char *)(&mp_temp) + ofs, (char *)mp1 + ofs, sizeof(Molecule) - ofs);
	memmove((char *)mp1 + ofs, (char *)mp2 + ofs, sizeof(Molecule) - ofs);
	memmove((char *)mp2 + ofs, (char *)(&mp_temp) + ofs, sizeof(Molecule) - ofs);
	mp1->arena = arena1;
	mp2->arena = arena2;
	mp1->mview = mview1;
	mp2->mview = mview2;
/*	if (mp1->arena != NULL && mp1->arena->mol == mp2)
		mp1->arena->mol = mp1;
	if (mp1->arena != NULL && mp1->arena->xmol == mp2)
		mp1->arena->xmol = mp1;
	if (mp2->arena != NULL && mp2->arena->mol == mp1)
		mp2->arena->mol = mp2;
	if (mp2->arena != NULL && mp2->arena->xmol == mp1)
		mp2->arena->xmol = mp2; */
}

#pragma mark ====== Mutex ======

void
MoleculeLock(Molecule *mol)
{
	if (mol == NULL || mol->mutex == NULL)
		return;
	MoleculeCallback_lockMutex(mol->mutex);
}

void
MoleculeUnlock(Molecule *mol)
{
	if (mol == NULL || mol->mutex == NULL)
		return;
	MoleculeCallback_unlockMutex(mol->mutex);
}

#pragma mark ====== Modify count ======

void
MoleculeIncrementModifyCount(Molecule *mp)
{
	if (mp != NULL) {
		if (++(mp->modifyCount) == 1)
			MoleculeCallback_notifyModification(mp, 0);
	/*	fprintf(stderr, "MoleculeIncrementModifyCount: %d\n", mp->modifyCount); */
	}
}

void
MoleculeClearModifyCount(Molecule *mp)
{
	if (mp != NULL) {
		mp->modifyCount = 0;
	/*	fprintf(stderr, "MoleculeClearModifyCount: %d\n", mp->modifyCount); */
	}
}

#pragma mark ====== File handling functions ======

static const char *
guessMoleculeType(const char *fname)
{
	char buf[1024], *p;
	FILE *fp;
	const char *retval = NULL;
	fp = fopen(fname, "rb");
	if (fp != NULL) {
		memset(buf, 0, sizeof buf);
		if (fread(buf, 1, sizeof buf - 1, fp) > 0) {
			if (strncmp(buf, "PSF", 3) == 0)
				retval = "psf";
			else if (((p = strstr(buf, "ATOM")) != NULL && (p == buf || p[-1] == '\n' || p[-1] == '\r'))
			|| ((p = strstr(buf, "HETATM")) != NULL && (p == buf || p[-1] == '\n' || p[-1] == '\r')))
				retval = "pdb";
			else
				retval = "???";  /*  unknown  */
		}
		fclose(fp);
	}
	return retval;
}

static int
guessElement(Atom *ap)
{
	int atomicNumber = -1;
	if (ap->atomicNumber > 0)
		atomicNumber = ap->atomicNumber;
	else {
		atomicNumber = GuessAtomicNumber(ap->element, ap->weight);
		if (atomicNumber <= 0 && ap->aname[0] != 0)
			atomicNumber = GuessAtomicNumber(ap->aname, ap->weight);
	}
	if (atomicNumber >= 0) {
		ap->atomicNumber = atomicNumber;
		if (ap->weight <= 0)
			ap->weight = WeightForAtomicNumber(atomicNumber);
		if (ap->element[0] == 0)
			ElementToString(atomicNumber, ap->element);
	}
	return atomicNumber;
}

static int
sReadLineWithInterrupt(char *buf, int size, FILE *stream, int *lineNumber)
{
	static int lastLineNumber = 0;
	if (lineNumber != NULL) {
		if (*lineNumber == 0)
			lastLineNumber = 0;
		else if (*lineNumber >= lastLineNumber + 1000) {
			if (MyAppCallback_checkInterrupt() != 0)
				return -1;  /*  User interrupt  */
			lastLineNumber = *lineNumber;
		}
	}
	return ReadLine(buf, size, stream, lineNumber);
}

static int
s_append_asprintf(char **buf, const char *fmt, ...)
{
	int len;
	char *s;
	va_list va;
	va_start(va, fmt);
	vasprintf(&s, fmt, va);
	len = (*buf == NULL ? 0 : strlen(*buf));
	if (s == NULL)
		return len;
	len += strlen(s);
	if (*buf == NULL) {
		*buf = malloc(len + 1);
		**buf = 0;
	} else {
		*buf = realloc(*buf, len + 1);
	}
	strcat(*buf, s);
	free(s);
	return len;
}

int
MoleculeLoadFile(Molecule *mp, const char *fname, const char *ftype, char **errbuf)
{
	int retval;
	if (ftype == NULL || *ftype == 0) {
		const char *cp;
		cp = strrchr(fname, '.');
		if (cp != NULL)
			ftype = cp + 1;
		else {
			cp = guessMoleculeType(fname);
			if (strcmp(cp, "???") != 0)
				ftype = cp;
		}
	}
	if (strcasecmp(ftype, "psf") == 0) {
		retval = MoleculeLoadPsfFile(mp, fname, errbuf);
	} else if (strcasecmp(ftype, "pdb") == 0) {
		retval = MoleculeReadCoordinatesFromPdbFile(mp, fname, errbuf);
	} else if (strcasecmp(ftype, "tep") == 0) {
		retval = MoleculeLoadTepFile(mp, fname, errbuf);
	} else if (strcasecmp(ftype, "res") == 0 || strcasecmp(ftype, "ins") == 0) {
		retval = MoleculeLoadShelxFile(mp, fname, errbuf);
	} else if (strcasecmp(ftype, "fchk") == 0 || strcasecmp(ftype, "fch") == 0) {
		retval = MoleculeLoadGaussianFchkFile(mp, fname, errbuf);
	} else {
		s_append_asprintf(errbuf, "Unknown format %s", ftype);
		return 1;
	}
/*	if (retval != 0) {
		retval = MoleculeLoadPsfFile(mp, fname, errbuf, errbufsize);
	} */
	if (retval == 0)
		MoleculeSetPath(mp, fname);
	return retval;
}

int
MoleculeLoadMbsfFile(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	char buf[1024];
	int i, j, k, err, fn, nframes, nwarnings;
	int lineNumber;
	int ibuf[12];
	Int iibuf[4];
	double dbuf[12];
	int mview_ibuf[18];
	double mview_dbuf[10];
	char cbuf[12][8];
	const char **pp;
	char *bufp, *valp, *comp;
	Int *ip;
	Double *dp;
	Vector v;
	Atom *ap;
	const int kUndefined = -10000000;
	err = 0;
	*errbuf = NULL;
	nwarnings = 0;
	if (mp->natoms != 0 || mp->par != NULL || mp->arena != NULL) {
		s_append_asprintf(errbuf, "The molecule must be empty");
		return 1;
	}
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot open file");
		return 1;
	}
	for (i = 0; i < 10; i++)
		mview_dbuf[i] = kUndefined;
	for (i = 0; i < 18; i++)
		mview_ibuf[i] = kUndefined;
	/*	flockfile(fp); */
	lineNumber = 0;
	fn = 0;
	nframes = 0;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		if (strncmp(buf, "!:", 2) != 0)
			continue;   /*  Skip until section header is found  */
		bufp = buf;
		strsep(&bufp, " \t\n");
		if (strcmp(buf, "!:atoms") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx seg_name res_seq res_name name type charge weight element atomic_number occupancy temp_factor int_charge */
				if (sscanf(buf, "%d %6s %d %6s %6s %6s %lf %lf %6s %d %lf %lf %d", &ibuf[0], cbuf[0], &ibuf[1], cbuf[1], cbuf[2], cbuf[3], &dbuf[0], &dbuf[1], cbuf[4], &ibuf[2], &dbuf[2], &dbuf[3], &ibuf[3]) < 13) {
					s_append_asprintf(errbuf, "line %d: coordinates cannot be read for atom %d", lineNumber, mp->natoms + 1);
					goto err_exit;
				}
				ap = AssignArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, mp->natoms, NULL);
				strncpy(ap->segName, cbuf[0], 4);
				ap->resSeq = ibuf[1];
				strncpy(ap->resName, cbuf[1], 4);
				strncpy(ap->aname, cbuf[2], 4);
				ap->type = AtomTypeEncodeToUInt(cbuf[3]);
				ap->charge = dbuf[0];
				ap->weight = dbuf[1];
				strncpy(ap->element, cbuf[4], 2);
				ap->atomicNumber = ibuf[2];
				ap->occupancy = dbuf[2];
				ap->tempFactor = dbuf[3];
				ap->intCharge = ibuf[3];
			}
			continue;
		} else if (strcmp(buf, "!:atoms_symop") == 0) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx symop symbase */
				if (sscanf(buf, "%d %d %d", &ibuf[0], &ibuf[1], &ibuf[2]) < 3) {
					s_append_asprintf(errbuf, "line %d: symmetry operations cannot be read for atom %d", lineNumber, i + 1);
					goto err_exit;
				}
				if (i >= mp->natoms) {
					s_append_asprintf(errbuf, "line %d: too many atomic symmetry info\n", lineNumber);
					goto err_exit;
				}
				ap = ATOM_AT_INDEX(mp->atoms, i);
				ap->symop.sym = ibuf[1] / 1000000;
				ap->symop.dx = (ibuf[1] % 1000000) / 10000;
				ap->symop.dy = (ibuf[1] % 10000) / 100;
				ap->symop.dz = ibuf[1] % 100;
				ap->symbase = ibuf[2];
				i++;
			}
			continue;
		} else if (strcmp(buf, "!:atoms_fix") == 0) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx fix_force fix_pos */
				if (sscanf(buf, "%d %lf %lf %lf %lf", &ibuf[0], &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3]) < 5) {
					s_append_asprintf(errbuf, "line %d: fix atom info cannot be read for atom %d", lineNumber, i + 1);
					goto err_exit;
				}
				if (i >= mp->natoms) {
					s_append_asprintf(errbuf, "line %d: too many fix atom info\n", lineNumber);
					goto err_exit;
				}
				ap = ATOM_AT_INDEX(mp->atoms, i);
				ap->fix_force = dbuf[0];
				ap->fix_pos.x = dbuf[1];
				ap->fix_pos.y = dbuf[2];
				ap->fix_pos.z = dbuf[3];
				i++;
			}
			continue;
		} else if (strcmp(buf, "!:uff_types") == 0) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx uff_type */
				if (sscanf(buf, "%d %6s", &ibuf[0], cbuf[0]) < 2) {
					s_append_asprintf(errbuf, "line %d: uff type info cannot be read for atom %d", lineNumber, i + 1);
					goto err_exit;
				}
				if (i >= mp->natoms) {
					s_append_asprintf(errbuf, "line %d: too many uff type info\n", lineNumber);
					goto err_exit;
				}
				ap = ATOM_AT_INDEX(mp->atoms, i);
				strncpy(ap->uff_type, cbuf[0], 5);
				ap->uff_type[5] = 0;
				i++;
			}
		} else if (strcmp(buf, "!:mm_exclude") == 0) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx mm_exclude periodic_exclude */
				if (sscanf(buf, "%d %d %d", &ibuf[0], &ibuf[1], &ibuf[2]) < 3) {
					s_append_asprintf(errbuf, "line %d: mm_exclude flags cannot be read for atom %d", lineNumber, i + 1);
					goto err_exit;
				}
				if (i >= mp->natoms) {
					s_append_asprintf(errbuf, "line %d: too many mm_exclude flags\n", lineNumber);
					goto err_exit;
				}
				ap = ATOM_AT_INDEX(mp->atoms, i);
				ap->mm_exclude = (ibuf[1] != 0);
				ap->periodic_exclude = (ibuf[2] != 0);
				i++;
			}
			continue;
		} else if (strcmp(buf, "!:pi_anchor") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx count */
				if ((j = sscanf(buf, "%d %d", &ibuf[0], &ibuf[1])) < 2) {
					s_append_asprintf(errbuf, "line %d: bad format for pi_anchor", lineNumber);
					goto err_exit;
				}
				i = ibuf[0];
				ap = ATOM_AT_INDEX(mp->atoms, i);
				if (ap->anchor != NULL) {
					s_append_asprintf(errbuf, "line %d: warning: duplicate pi_anchor entry", lineNumber);
					AtomConnectResize(&ap->anchor->connect, 0);
					free(ap->anchor->coeffs);
					free(ap->anchor);
				}
				ap->anchor = (PiAnchor *)calloc(sizeof(PiAnchor), 1);
				if (ibuf[1] < 2 || ibuf[1] >= mp->natoms) {
					s_append_asprintf(errbuf, "line %d: bad number of components for pi_anchor", lineNumber);
					goto err_exit;
				}
				AtomConnectResize(&ap->anchor->connect, ibuf[1]);
				ip = AtomConnectData(&ap->anchor->connect);
				NewArray(&ap->anchor->coeffs, &ap->anchor->ncoeffs, sizeof(Double), ibuf[1]);
				j = ibuf[1];
				for (i = 0; i < j; i++) {
					if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
						s_append_asprintf(errbuf, "line %d: unexpected end of file while reading pi_anchors", lineNumber);
						goto err_exit;
					}
					if (sscanf(buf, "%d %lf", &ibuf[0], &dbuf[0]) < 2) {
						s_append_asprintf(errbuf, "line %d: bad format for pi_anchor", lineNumber);
						goto err_exit;
					}
					if (ibuf[0] < 0 || ibuf[0] >= mp->natoms) {
						s_append_asprintf(errbuf, "line %d: atom index out of range", lineNumber);
						goto err_exit;
					}
					if (dbuf[0] <= 0.0) {
						s_append_asprintf(errbuf, "line %d: the pi anchor weights should be positive", lineNumber);
						goto err_exit;
					}
					ip[i] = ibuf[0];
					ap->anchor->coeffs[i] = dbuf[0];
				}
			}
			continue;
		} else if (strcmp(buf, "!:positions") == 0) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx x y z */
				if ((j = sscanf(buf, "%d %lf %lf %lf %lf %lf %lf", &ibuf[0], &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5])) < 4) {
					s_append_asprintf(errbuf, "line %d: atom position cannot be read for atom %d frame %d", lineNumber, i + 1, nframes);
					goto err_exit;
				}
				if (j > 4 && nframes != 0) {
					s_append_asprintf(errbuf, "line %d: atom position sigma can only be given for frame 0", lineNumber);
					goto err_exit;
				}
				if (j > 4 && j != 7) {
					s_append_asprintf(errbuf, "line %d: atom position sigma cannot be read for atom %d frame %d", lineNumber, i + 1, nframes);
					goto err_exit;
				}
				if (i >= mp->natoms) {
					s_append_asprintf(errbuf, "line %d: too many atom position records\n", lineNumber);
					goto err_exit;
				}
				v.x = dbuf[0];
				v.y = dbuf[1];
				v.z = dbuf[2];
				ap = ATOM_AT_INDEX(mp->atoms, i);
				if (nframes > 0) {
					AssignArray(&ap->frames, &ap->nframes, sizeof(Vector), nframes, &v);
					if (nframes == 1)
						ap->frames[0] = ap->r;
				}
				ap->r = v;
				if (j == 7) {
					ap->sigma.x = dbuf[3];
					ap->sigma.y = dbuf[4];
					ap->sigma.z = dbuf[5];
				}
				i++;
			}
			nframes++;
			if (nframes >= 2) {
				mp->nframes = nframes;
				mp->cframe = nframes - 1;
			} else {
				mp->nframes = mp->cframe = 0;
			}
			continue;
		} else if (strcmp(buf, "!:bonds") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* from1 to1 from2 to2 from3 to3 from4 to4 */ 
				i = sscanf(buf, "%d %d %d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5], &ibuf[6], &ibuf[7]);
				if (i < 2 || i % 2 != 0) {
					s_append_asprintf(errbuf, "line %d: bad bond format", lineNumber);
					goto err_exit;
				}
				for (j = 0; j < i; j += 2) {
					iibuf[0] = ibuf[j];
					iibuf[1] = ibuf[j + 1];
					if (iibuf[0] < 0 || iibuf[0] >= mp->natoms || iibuf[1] < 0 || iibuf[1] >= mp->natoms || iibuf[0] == iibuf[1]) {
						s_append_asprintf(errbuf, "line %d: warning: bad bond specification (%d-%d) - skipped\n", lineNumber, iibuf[0], iibuf[1]);
						nwarnings++;
					} else if (AtomConnectHasEntry(&(ATOM_AT_INDEX(mp->atoms, iibuf[0])->connect), iibuf[1])) {
						s_append_asprintf(errbuf, "line %d: warning: bond %d-%d is already present - skipped\n", lineNumber, iibuf[0], iibuf[1]);
						nwarnings++;
					} else {
						AssignArray(&mp->bonds, &mp->nbonds, sizeof(Int) * 2, mp->nbonds, iibuf);
						AtomConnectInsertEntry(&(ATOM_AT_INDEX(mp->atoms, iibuf[0])->connect), -1, iibuf[1]);
						AtomConnectInsertEntry(&(ATOM_AT_INDEX(mp->atoms, iibuf[1])->connect), -1, iibuf[0]);
					}
				}
			}
			continue;
		} else if (strcmp(buf, "!:bond_orders") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* b1 b2 b3 b4 */
				i = sscanf(buf, "%lf %lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3]);
				if (i == 0) {
					s_append_asprintf(errbuf, "line %d: bad bond order format", lineNumber);
					goto err_exit;
				}
				for (j = 0; j < i; j++) {
					AssignArray(&mp->bondOrders, &mp->nbondOrders, sizeof(Double), mp->nbondOrders, &dbuf[j]);
				}
			}
			if (mp->nbondOrders > mp->nbonds) {
				s_append_asprintf(errbuf, "line %d: warning: the number of bond order info (%d) exceeds number of bonds (%d) - ignoring excess info\n", lineNumber, mp->nbondOrders, mp->nbonds);
				nwarnings++;
				mp->nbondOrders = mp->nbonds;
			} else if (mp->nbondOrders < mp->nbonds) {
				s_append_asprintf(errbuf, "line %d: warning: the number of bond order info (%d) is less than number of bonds (%d)\n", lineNumber, mp->nbondOrders, mp->nbonds);
				nwarnings++;
				j = mp->nbondOrders;
				AssignArray(&mp->bondOrders, &mp->nbondOrders, sizeof(Double), mp->nbonds - 1, NULL);
				for (i = j; i < mp->nbonds; i++)
					mp->bondOrders[i] = 0.0;
			}
			continue;
			
		} else if (strcmp(buf, "!:angles") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a1 b1 c1 a2 b2 c2 a3 b3 c3 */ 
				i = sscanf(buf, "%d %d %d %d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5], &ibuf[6], &ibuf[7], &ibuf[8]);
				if (i == 0 || i % 3 != 0) {
					s_append_asprintf(errbuf, "line %d: bad angle format", lineNumber);
					goto err_exit;
				}
				for (j = 0; j < i; j += 3) {
					iibuf[0] = ibuf[j];
					iibuf[1] = ibuf[j + 1];
					iibuf[2] = ibuf[j + 2];
					if (iibuf[0] < 0 || iibuf[0] >= mp->natoms || iibuf[1] < 0 || iibuf[1] >= mp->natoms || iibuf[2] < 0 || iibuf[2] >= mp->natoms || iibuf[0] == iibuf[1] || iibuf[1] == iibuf[2]) {
						s_append_asprintf(errbuf, "line %d: warning: bad angle specification (%d-%d-%d) - skipped\n", lineNumber, iibuf[0], iibuf[1], iibuf[2]);
						nwarnings++;
					} else if (MoleculeAreAtomsConnected(mp, iibuf[1], iibuf[0]) == 0 || MoleculeAreAtomsConnected(mp, iibuf[1], iibuf[2]) == 0) {
						s_append_asprintf(errbuf, "line %d: warning: angle with non-bonded atoms (%d-%d-%d) - skipped\n", lineNumber, iibuf[0], iibuf[1], iibuf[2]);
						nwarnings++;						
					} else if (MoleculeLookupAngle(mp, iibuf[0], iibuf[1], iibuf[2]) >= 0) {
						s_append_asprintf(errbuf, "line %d: warning: angle %d-%d-%d is already present - skipped\n", lineNumber, iibuf[0], iibuf[1], iibuf[2]);
						nwarnings++;
					} else {
						AssignArray(&mp->angles, &mp->nangles, sizeof(Int) * 3, mp->nangles, iibuf);
					}
				}
			}
			continue;
		} else if (strcmp(buf, "!:dihedrals") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a1 b1 c1 d1 a2 b2 c2 d2 */ 
				i = sscanf(buf, "%d %d %d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5], &ibuf[6], &ibuf[7]);
				if (i == 0 || i % 4 != 0) {
					s_append_asprintf(errbuf, "line %d: bad dihedral format", lineNumber);
					goto err_exit;
				}
				for (j = 0; j < i; j += 4) {
					iibuf[0] = ibuf[j];
					iibuf[1] = ibuf[j + 1];
					iibuf[2] = ibuf[j + 2];
					iibuf[3] = ibuf[j + 3];
					if (iibuf[0] < 0 || iibuf[0] >= mp->natoms || iibuf[1] < 0 || iibuf[1] >= mp->natoms || iibuf[2] < 0 || iibuf[2] >= mp->natoms || iibuf[3] < 0 || iibuf[3] >= mp->natoms || iibuf[0] == iibuf[1] || iibuf[1] == iibuf[2] || iibuf[2] == iibuf[3] || iibuf[0] == iibuf[2] || iibuf[1] == iibuf[3] || iibuf[0] == iibuf[3]) {
						s_append_asprintf(errbuf, "line %d: warning: bad dihedral specification (%d-%d-%d-%d) - skipped\n", lineNumber, iibuf[0], iibuf[1], iibuf[2], iibuf[3]);
						nwarnings++;
					} else if (MoleculeAreAtomsConnected(mp, iibuf[1], iibuf[0]) == 0 || MoleculeAreAtomsConnected(mp, iibuf[1], iibuf[2]) == 0 || MoleculeAreAtomsConnected(mp, iibuf[2], iibuf[3]) == 0) {
						s_append_asprintf(errbuf, "line %d: warning: dihedral with non-bonded atoms (%d-%d-%d-%d) - skipped\n", lineNumber, iibuf[0], iibuf[1], iibuf[2], iibuf[3]);
						nwarnings++;						
					} else if (MoleculeLookupDihedral(mp, iibuf[0], iibuf[1], iibuf[2], iibuf[3]) >= 0) {
						s_append_asprintf(errbuf, "line %d: warning: dihedral %d-%d-%d-%d is already present - skipped\n", lineNumber, iibuf[0], iibuf[1], iibuf[2], iibuf[3]);
						nwarnings++;
					} else {
						AssignArray(&mp->dihedrals, &mp->ndihedrals, sizeof(Int) * 4, mp->ndihedrals, iibuf);
					}
				}
			}
			continue;
		} else if (strcmp(buf, "!:impropers") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a1 b1 c1 d1 a2 b2 c2 d2 */ 
				i = sscanf(buf, "%d %d %d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5], &ibuf[6], &ibuf[7]);
				if (i == 0 || i % 4 != 0) {
					s_append_asprintf(errbuf, "line %d: bad improper format", lineNumber);
					goto err_exit;
				}
				for (j = 0; j < i; j += 4) {
					iibuf[0] = ibuf[j];
					iibuf[1] = ibuf[j + 1];
					iibuf[2] = ibuf[j + 2];
					iibuf[3] = ibuf[j + 3];
					if (iibuf[0] < 0 || iibuf[0] >= mp->natoms || iibuf[1] < 0 || iibuf[1] >= mp->natoms || iibuf[2] < 0 || iibuf[2] >= mp->natoms || iibuf[3] < 0 || iibuf[3] >= mp->natoms || iibuf[0] == iibuf[1] || iibuf[1] == iibuf[2] || iibuf[2] == iibuf[3] || iibuf[0] == iibuf[2] || iibuf[1] == iibuf[3] || iibuf[0] == iibuf[3]) {
						s_append_asprintf(errbuf, "line %d: warning: bad improper specification (%d-%d-%d-%d) - skipped\n", lineNumber, iibuf[0], iibuf[1], iibuf[2], iibuf[3]);
						nwarnings++;
					} else if (MoleculeAreAtomsConnected(mp, iibuf[2], iibuf[0]) == 0 || MoleculeAreAtomsConnected(mp, iibuf[2], iibuf[1]) == 0 || MoleculeAreAtomsConnected(mp, iibuf[2], iibuf[3]) == 0) {
						s_append_asprintf(errbuf, "line %d: warning: improper with non-bonded atoms (%d-%d-%d-%d) - skipped\n", lineNumber, iibuf[0], iibuf[1], iibuf[2], iibuf[3]);
						nwarnings++;						
					} else if (MoleculeLookupImproper(mp, iibuf[0], iibuf[1], iibuf[2], iibuf[3]) >= 0) {
						s_append_asprintf(errbuf, "line %d: warning: improper %d-%d-%d-%d is already present - skipped\n", lineNumber, iibuf[0], iibuf[1], iibuf[2], iibuf[3]);
						nwarnings++;
					} else {
						AssignArray(&mp->impropers, &mp->nimpropers, sizeof(Int) * 4, mp->nimpropers, iibuf);
					}
				}
			}
			continue;
		} else if (strcmp(buf, "!:xtalcell") == 0 && mp->cell == NULL) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a b c alpha beta gamma [sigmaflag] */ 
				if ((j = sscanf(buf, "%lf %lf %lf %lf %lf %lf %d", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5], &ibuf[0])) < 6) {
					s_append_asprintf(errbuf, "line %d: bad xtalcell format", lineNumber);
					goto err_exit;
				}
				MoleculeSetCell(mp, dbuf[0], dbuf[1], dbuf[2], dbuf[3], dbuf[4], dbuf[5], 0);
				if (j == 7 && ibuf[0] != 0) {
					if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
						s_append_asprintf(errbuf, "line %d: sigma for xtalcell are missing", lineNumber);
						goto err_exit;
					}
					if (sscanf(buf, "%lf %lf %lf %lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5]) < 6) {
						s_append_asprintf(errbuf,"line %d: bad xtalcell sigma format", lineNumber);
						goto err_exit;
					}
					if (mp->cell != NULL) {
						mp->cell->has_sigma = 1;
						for (i = 0; i < 6; i++) {
							mp->cell->cellsigma[i] = dbuf[i];
						}
					} else {
						s_append_asprintf(errbuf, "line %d: cell sigma are given while cell is not given", lineNumber);
					}
				}
			}
			continue;
		} else if (strcmp(buf, "!:symmetry_operations") == 0) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				Transform tr;
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a11 a12 a13; a21 a22 a23; a31 a32 a33; t1 t2 t3 */
				if (sscanf(buf, "%lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2]) < 3) {
					s_append_asprintf(errbuf, "line %d: bad symmetry_operation format", lineNumber);
					goto err_exit;
				}
				if (i < 3) {
					tr[i] = dbuf[0];
					tr[i + 3] = dbuf[1];
					tr[i + 6] = dbuf[2];
				} else {
					tr[9] = dbuf[0];
					tr[10] = dbuf[1];
					tr[11] = dbuf[2];
				}
				i++;
				if (i == 4) {
					AssignArray(&mp->syms, &mp->nsyms, sizeof(Transform), mp->nsyms, tr);
					i = 0;
				}
			}
			continue;
		} else if (strcmp(buf, "!:anisotropic_thermal_parameters") == 0) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* b11 b22 b33 b12 b13 b23 [has_sigma] */
				if ((j = sscanf(buf, "%lf %lf %lf %lf %lf %lf %d", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5], &ibuf[0])) < 6) {
					s_append_asprintf(errbuf, "line %d: anisotropic thermal parameters cannot be read for atom %d", lineNumber, i + 1);
					goto err_exit;
				}
				if (i >= mp->natoms) {
					s_append_asprintf(errbuf, "line %d: too many anisotropic thermal parameters\n", lineNumber);
					goto err_exit;
				}
				if (dbuf[0] == 0.0 && dbuf[1] == 0.0 && dbuf[2] == 0.0 && dbuf[3] == 0.0 && dbuf[4] == 0.0 && dbuf[5] == 0.0) {
					/*  Skip it  */
				} else {
					MoleculeSetAniso(mp, i, 0, dbuf[0], dbuf[1], dbuf[2], dbuf[3], dbuf[4], dbuf[5], NULL);
				}
				if (j == 7 && ibuf[0] != 0) {
					if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
						s_append_asprintf(errbuf, "line %d: anisotropic thermal parameters sigma missing", lineNumber);
						goto err_exit;
					}
					if (sscanf(buf, "%lf %lf %lf %lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5]) < 6) {
						s_append_asprintf(errbuf, "line %d: anisotropic thermal parameters sigma cannot be read for atom %d", lineNumber, i + 1);
						goto err_exit;
					}
					ap = ATOM_AT_INDEX(mp->atoms, i);
					if (ap->aniso == NULL) {
						s_append_asprintf(errbuf, "line %d: anisotropic thermal parameters sigma are given while the parameters are not given", lineNumber);
						goto err_exit;
					}
					ap->aniso->has_bsig = 1;
					for (j = 0; j < 6; j++)
						ap->aniso->bsig[j] = dbuf[j];
				}
				i++;
			}
			continue;
		} else if (strcmp(buf, "!:periodic_box") == 0) {
			Vector vs[5];
			Byte has_sigma = 0;
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* ax ay az; bx by bz; cx cy cz; ox oy oz; fx fy fz [sigma; sa sb sc s_alpha s_beta s_gamma] */
				if (i < 4) {
					if (sscanf(buf, "%lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2]) < 3) {
						s_append_asprintf(errbuf, "line %d: bad periodic_box format", lineNumber);
						goto err_exit;
					}
					vs[i].x = dbuf[0];
					vs[i].y = dbuf[1];
					vs[i].z = dbuf[2];
					i++;
					continue;
				}
				if ((j = sscanf(buf, "%d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3])) < 3) {
					s_append_asprintf(errbuf, "line %d: bad periodic_box format", lineNumber);
					goto err_exit;
				}
				if (j == 4 && ibuf[3] != 0)
					has_sigma = 1;
				cbuf[0][0] = ibuf[0];
				cbuf[0][1] = ibuf[1];
				cbuf[0][2] = ibuf[2];
				MoleculeSetPeriodicBox(mp, vs, vs + 1, vs + 2, vs + 3, cbuf[0], 0);
				if (has_sigma) {
					if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
						s_append_asprintf(errbuf, "line %d: sigma for cell parameters are missing", lineNumber);
						goto err_exit;
					}
					if (sscanf(buf, "%lf %lf %lf %lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5]) < 6) {
						s_append_asprintf(errbuf, "line %d: bad periodic_box sigma format", lineNumber);
						goto err_exit;
					}
					if (mp->cell != NULL) {
						mp->cell->has_sigma = 1;
						for (i = 0; i < 6; i++) {
							mp->cell->cellsigma[i] = dbuf[i];
						}
					} else {
						s_append_asprintf(errbuf, "line %d: cell sigma are given while cell is not given", lineNumber);
					}
				}
				break;
			}
			continue;
		} else if (strcmp(buf, "!:frame_periodic_boxes") == 0) {
			Vector vs[5];
			i = 0;
		/*	mp->useFlexibleCell = 1;  *//*  The presence of this block causes asserting this flag  */
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				if (sscanf(buf, "%lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2]) < 3) {
					s_append_asprintf(errbuf, "line %d: bad frame_periodic_box format", lineNumber);
					goto err_exit;
				}
				vs[i].x = dbuf[0];
				vs[i].y = dbuf[1];
				vs[i].z = dbuf[2];
				i++;
				if (i == 4) {
					AssignArray(&mp->frame_cells, &mp->nframe_cells, sizeof(Vector) * 4, mp->nframe_cells, vs);
					i = 0;
				}
			}
			if (mp->cframe < mp->nframe_cells) {
				/*  mp->cframe should already have been set when positions are read  */
				Vector *vp = &mp->frame_cells[mp->cframe * 4];
				static char defaultFlags[] = {1, 1, 1};
				char *flags = (mp->cell != NULL ? mp->cell->flags : defaultFlags);
				MoleculeSetPeriodicBox(mp, vp, vp + 1, vp + 2, vp + 3, flags, 0);
			}
			continue;
		} else if (strcmp(buf, "!:md_parameters") == 0) {
			MDArena *arena;
			if (mp->arena == NULL)
				mp->arena = md_arena_new(NULL);
			arena = mp->arena;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				bufp = buf;
				comp = strsep(&bufp, " \t");
				if (bufp != NULL) {
					while (*bufp == ' ' || *bufp == '\t')
						bufp++;
					valp = strsep(&bufp, "\n");
				} else valp = NULL;
				if (strcmp(comp, "alchem_flags") == 0) {
					j = (valp == NULL ? 0 : atoi(valp));
					if (j > 0) {
						valp = (char *)malloc(j);
						i = 0;
						while ((k = fgetc(fp)) >= 0) {
							ungetc(k, fp);
							if (k < '0' || k > '9') {
								s_append_asprintf(errbuf, "line %d: too few flags in alchem_flags block", lineNumber + 1);
								free(valp);
								goto err_exit;
							}
							ReadLine(buf, sizeof buf, fp, &lineNumber);
							bufp = buf;
							while (*bufp != 0) {
								if (*bufp >= '0' && *bufp <= '2') {
									if (i >= j) {
										s_append_asprintf(errbuf, "line %d: too many flags in alchem_flags block", lineNumber);
										free(valp);
										goto err_exit;
									}
									valp[i++] = *bufp - '0';
								} else if (*bufp != ' ' && *bufp != '\t' && *bufp != '\n') {
									s_append_asprintf(errbuf, "line %d: strange character (0x%02x) in alchem_flags block", lineNumber, (int)*bufp);
									free(valp);
									goto err_exit;
								}
								bufp++;
							}
							if (i == j)
								break;
						}
						md_set_alchemical_flags(arena, j, valp);
						free(valp);
					}
					continue;
				}
				/*  In the following, the redundant "!= NULL" is to suppress suprious warning  */
				if ((strcmp(comp, "log_file") == 0 && (pp = &arena->log_result_name) != NULL)
					|| (strcmp(comp, "coord_file") == 0 && (pp = &arena->coord_result_name) != NULL)
					|| (strcmp(comp, "vel_file") == 0 && (pp = &arena->vel_result_name) != NULL)
					|| (strcmp(comp, "force_file") == 0 && (pp = &arena->force_result_name) != NULL)
					|| (strcmp(comp, "debug_file") == 0 && (pp = &arena->debug_result_name) != NULL)) {
					if (*valp == 0 || strstr(valp, "(null)") == valp)
						*pp = NULL;
					else {
						valp = strdup(valp);
						if (valp != NULL) {
							char *valp1 = strchr(valp, '\n');
							if (valp1 != NULL)
								*valp1 = 0;
						}
						*pp = valp;
					}
				} else if ((strcmp(comp, "debug_output_level") == 0 && (ip = &arena->debug_output_level) != NULL)
						   || (strcmp(comp, "coord_output_freq") == 0 && (ip = &arena->coord_output_freq) != NULL)
						   || (strcmp(comp, "energy_output_freq") == 0 && (ip = &arena->energy_output_freq) != NULL)
						   || (strcmp(comp, "coord_frame") == 0 && (ip = &arena->coord_result_frame) != NULL)
						   || (strcmp(comp, "andersen_freq") == 0 && (ip = &arena->andersen_thermo_freq) != NULL)
						   || (strcmp(comp, "random_seed") == 0 && (ip = &arena->random_seed) != NULL)
						   || (strcmp(comp, "use_xplor_shift") == 0 && (ip = &arena->use_xplor_shift) != NULL)
						   || (strcmp(comp, "relocate_center") == 0 && (ip = &arena->relocate_center) != NULL)
						   || (strcmp(comp, "surface_potential_freq") == 0 && (ip = &arena->surface_potential_freq) != NULL)
						   || (strcmp(comp, "use_graphite") == 0 && (ip = &arena->use_graphite) != NULL)) {
					*ip = (valp == NULL ? 0 : atoi(valp));
				} else if ((strcmp(comp, "timestep") == 0 && (dp = &arena->timestep) != NULL)
						   || (strcmp(comp, "cutoff") == 0 && (dp = &arena->cutoff) != NULL)
						   || (strcmp(comp, "electro_cutoff") == 0 && (dp = &arena->electro_cutoff) != NULL)
						   || (strcmp(comp, "pairlist_distance") == 0 && (dp = &arena->pairlist_distance) != NULL)
						   || (strcmp(comp, "switch_distance") == 0 && (dp = &arena->switch_distance) != NULL)
						   || (strcmp(comp, "temperature") == 0 && (dp = &arena->temperature) != NULL)
						   || (strcmp(comp, "andersen_coupling") == 0 && (dp = &arena->andersen_thermo_coupling) != NULL)
						   || (strcmp(comp, "dielectric") == 0 && (dp = &arena->dielectric) != NULL)
						   || (strcmp(comp, "gradient_convergence") == 0 && (dp = &arena->gradient_convergence) != NULL)
						   || (strcmp(comp, "coordinate_convergence") == 0 && (dp = &arena->coordinate_convergence) != NULL)
						   || (strcmp(comp, "scale14_vdw") == 0 && (dp = &arena->scale14_vdw) != NULL)
						   || (strcmp(comp, "scale14_elect") == 0 && (dp = &arena->scale14_elect) != NULL)
						   || (strcmp(comp, "surface_probe_radius") == 0 && (dp = &arena->probe_radius) != NULL)
						   || (strcmp(comp, "surface_tension") == 0 && (dp = &arena->surface_tension) != NULL)
						   || (strcmp(comp, "alchemical_lambda") == 0 && (dp = &arena->alchem_lambda) != NULL)
						   || (strcmp(comp, "alchemical_delta_lambda") == 0 && (dp = &arena->alchem_dlambda) != NULL)) {
					*dp = (valp == NULL ? 0.0 : strtod(valp, NULL));
				}
			}
			continue;
		} else if (strcmp(buf, "!:pressure_control_parameters") == 0) {
			MDPressureArena *pressure;
			if (mp->arena == NULL)
				mp->arena = md_arena_new(mp);
			if (mp->arena->pressure == NULL)
				mp->arena->pressure = pressure_new();
			pressure = mp->arena->pressure;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				bufp = buf;
				comp = strsep(&bufp, " \t");
				if (bufp != NULL) {
					while (*bufp == ' ' || *bufp == '\t')
						bufp++;
					valp = strsep(&bufp, "\n");
				} else valp = NULL;
				if (strcmp(comp, "pressure") == 0) {
					if (sscanf(valp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5], &dbuf[6], &dbuf[7], &dbuf[8]) < 9) {
						s_append_asprintf(errbuf, "line %d: bad format", lineNumber);
						goto err_exit;
					}
					for (i = 0; i < 9; i++)
						pressure->apply[i] = dbuf[i];
				} else if (strcmp(comp, "pressure_cell_flexibility") == 0) {
					if (sscanf(valp, "%lf %lf %lf %lf %lf %lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5], &dbuf[6], &dbuf[7]) < 8) {
						s_append_asprintf(errbuf, "line %d: bad format", lineNumber);
						goto err_exit;
					}
					for (i = 0; i < 8; i++)
						pressure->cell_flexibility[i] = dbuf[i];
				} else if ((strcmp(comp, "pressure_freq") == 0 && (ip = &pressure->freq) != NULL)) {
					*ip = (valp == NULL ? 0 : atoi(valp));
				} else if ((strcmp(comp, "pressure_coupling") == 0 && (dp = &pressure->coupling) != NULL)
						   || (strcmp(comp, "pressure_fluctuate_cell_origin") == 0 && (dp = &pressure->fluctuate_cell_origin) != NULL)
						   || (strcmp(comp, "pressure_fluctuate_cell_orientation") == 0 && (dp = &pressure->fluctuate_cell_orientation) != NULL)) {
					*dp = (valp == NULL ? 0.0 : strtod(valp, NULL));
				}
			}
			continue;
		} else if (strcmp(buf, "!:velocity") == 0) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx vx vy vz */
				if (sscanf(buf, "%d %lf %lf %lf", &ibuf[0], &dbuf[0], &dbuf[1], &dbuf[2]) < 4) {
					s_append_asprintf(errbuf, "line %d: atom velocity cannot be read for atom %d", lineNumber, i + 1);
					goto err_exit;
				}
				if (i >= mp->natoms) {
					s_append_asprintf(errbuf, "line %d: too many atom velocity records\n", lineNumber);
					goto err_exit;
				}
				ap = ATOM_AT_INDEX(mp->atoms, i);
				ap->v.x = dbuf[0];
				ap->v.y = dbuf[1];
				ap->v.z = dbuf[2];
				i++;
			}
			continue;
		} else if (strcmp(buf, "!:force") == 0) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx fx fy fz */
				if (sscanf(buf, "%d %lf %lf %lf", &ibuf[0], &dbuf[0], &dbuf[1], &dbuf[2]) < 4) {
					s_append_asprintf(errbuf, "line %d: atom force cannot be read for atom %d", lineNumber, i + 1);
					goto err_exit;
				}
				if (i >= mp->natoms) {
					s_append_asprintf(errbuf, "line %d: too many atom force records\n", lineNumber);
					goto err_exit;
				}
				ap = ATOM_AT_INDEX(mp->atoms, i);
				ap->f.x = dbuf[0];
				ap->f.y = dbuf[1];
				ap->f.z = dbuf[2];
				i++;
			}
			continue;
		} else if (strcmp(buf, "!:parameter") == 0 || strcmp(buf, "!:parameters") == 0) {
			Parameter *par = mp->par;
			if (par == NULL) {
				mp->par = ParameterNew();
				par = mp->par;
			}
			bufp = NULL;
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				j = ParameterReadFromString(par, buf, &bufp, fname, lineNumber, 0);
				if (j < 0) {
					s_append_asprintf(errbuf, "%s", bufp);
					free(bufp);
					goto err_exit;
				}
				i += j;
			}
			if (bufp != NULL) {
				s_append_asprintf(errbuf, "%s", bufp);
				free(bufp);
			}
			continue;
		} else if (strcmp(buf, "!:trackball") == 0) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				if (mp->mview == NULL || mp->mview->track == NULL)
					continue;  /*  Skip (this should not happen though)  */
				/* scale; trx try trz; theta_deg x y z */
				if ((i == 0 && sscanf(buf, "%lf", &dbuf[0]) < 1)
					|| (i == 1 && sscanf(buf, "%lf %lf %lf",
										 &dbuf[1], &dbuf[2], &dbuf[3]) < 3)
					|| (i == 2 && sscanf(buf, "%lf %lf %lf %lf",
										 &dbuf[4], &dbuf[5], &dbuf[6], &dbuf[7]) < 4)) {
					s_append_asprintf(errbuf, "line %d: bad trackball format", lineNumber);
					goto err_exit;
				}
				if (i == 0)
					TrackballSetScale(mp->mview->track, dbuf[0]);
				else if (i == 1)
					TrackballSetTranslate(mp->mview->track, dbuf + 1);
				else if (i == 2)
					TrackballSetRotate(mp->mview->track, dbuf + 4);
				i++;
			}
			continue;
		} else if (strcmp(buf, "!:view") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				if (mp->mview == NULL)
					continue;  /*  Skip (this should not happen, though)  */
				bufp = buf;
				comp = strsep(&bufp, " \t");
				if (bufp != NULL) {
					while (*bufp == ' ' || *bufp == '\t')
						bufp++;
					valp = strsep(&bufp, "\n");
				} else valp = NULL;
				if (strcmp(comp, "show_unit_cell") == 0)
					mp->mview->showUnitCell = atoi(valp);
				else if (strcmp(comp, "show_periodic_box") == 0)
					mp->mview->showPeriodicBox = atoi(valp);
				else if (strcmp(comp, "show_expanded_atoms") == 0)
					mp->mview->showExpandedAtoms = atoi(valp);
				else if (strcmp(comp, "show_ellipsoids") == 0)
					mp->mview->showEllipsoids = atoi(valp);
				else if (strcmp(comp, "show_hydrogens") == 0)
					mp->mview->showHydrogens = atoi(valp);
				else if (strcmp(comp, "show_dummy_atoms") == 0)
					mp->mview->showDummyAtoms = atoi(valp);
				else if (strcmp(comp, "show_rotation_center") == 0)
					mp->mview->showRotationCenter = atoi(valp);
				else if (strcmp(comp, "show_graphite_flag") == 0)
					mp->mview->showGraphiteFlag = atoi(valp);
				else if (strcmp(comp, "show_periodic_image_flag") == 0)
					mp->mview->showPeriodicImageFlag = atoi(valp);
				else if (strcmp(comp, "show_graphite") == 0)
					mp->mview->showGraphite = atoi(valp);
				else if (strcmp(comp, "show_expanded_atoms") == 0)
					mp->mview->showExpandedAtoms = atoi(valp);
				else if (strcmp(comp, "atom_resolution") == 0 && (i = atoi(valp)) >= 6)
					mp->mview->atomResolution = i;
				else if (strcmp(comp, "bond_resolution") == 0 && (i = atoi(valp)) >= 4)
					mp->mview->bondResolution = i;
				else if (strcmp(comp, "atom_radius") == 0)
					mp->mview->atomRadius = strtod(valp, NULL);
				else if (strcmp(comp, "bond_radius") == 0)
					mp->mview->bondRadius = strtod(valp, NULL);
				else if (strcmp(comp, "show_periodic_image") == 0) {
					sscanf(valp, "%d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5]);
					for (i = 0; i < 6; i++)
						mp->mview->showPeriodicImage[i] = ibuf[i];
				}
			}
			continue;
		} else if (strcmp(buf, "!:property") == 0) {
			char dec[1024];
			i = 0;
			bufp = buf + 13;
			while (*bufp != 0 && *bufp != '\n' && bufp < (buf + sizeof buf - 3)) {
				if (*bufp == '%') {
					dec[i] = bufp[1];
					dec[i + 1] = bufp[2];
					dec[i + 2] = 0;
					dec[i++] = strtol(dec, NULL, 16);
					bufp += 3;
				} else {
					dec[i++] = *bufp++;
				}
				if (i >= 1000)
					break;
			}
			if (i == 0)
				continue;
			dec[i] = 0;
			i = MoleculeCreateProperty(mp, dec);
			if (i < 0) {
				s_append_asprintf(errbuf, "line %d: warning: duplicate molecular property %s - ignored\n", lineNumber, dec);
				nwarnings++;
				continue;
			}
			j = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				if (j >= nframes) {
					s_append_asprintf(errbuf, "line %d: warning: too many molecular property %s - ignored\n", lineNumber, dec);
					nwarnings++;
					break;
				}
				dbuf[0] = strtod(buf, NULL);
				mp->molprops[i].propvals[j] = dbuf[0];
				j++;
			}
			continue;
		} else if (strcmp(buf, "!:gaussian_primitives") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* sym nprims a_idx */
				if (sscanf(buf, "%6s %d %d", cbuf[0], &ibuf[0], &ibuf[1]) < 3) {
					s_append_asprintf(errbuf, "line %d: the gaussian primitive info cannot be read", lineNumber);
					goto err_exit;
				}
				if (strcasecmp(cbuf[0], "S") == 0) {
					ibuf[2] = 0;
				} else if (strcasecmp(cbuf[0], "P") == 0) {
					ibuf[2] = 1;
				} else if (strcasecmp(cbuf[0], "SP") == 0) {
					ibuf[2] = -1;
				} else if (strcasecmp(cbuf[0], "D") == 0) {
					ibuf[2] = 2;
				} else if (strcasecmp(cbuf[0], "D5") == 0) {
					ibuf[2] = -2;
				} else if (strcasecmp(cbuf[0], "F") == 0) {
					ibuf[2] = 3;
				} else if (strcasecmp(cbuf[0], "F7") == 0) {
					ibuf[2] = -3;
				} else if (strcasecmp(cbuf[0], "G") == 0) {
					ibuf[2] = 4;
				} else if (strcasecmp(cbuf[0], "G9") == 0) {
					ibuf[2] = -4;
				} else {
					s_append_asprintf(errbuf, "line %d: the gaussian primitive type %s is unknown", lineNumber, cbuf[0]);
					goto err_exit;
				}
				if (ibuf[0] <= 0) {
					s_append_asprintf(errbuf, "line %d: the number of primitive (%d) must be positive", lineNumber, ibuf[0]);
					goto err_exit;
				}
				if (ibuf[1] < 0 || ibuf[1] >= mp->natoms) {
					s_append_asprintf(errbuf, "line %d: the atom index (%d) is out of range", lineNumber, ibuf[1]);
					goto err_exit;
				}
				MoleculeAddGaussianOrbitalShell(mp, ibuf[1], ibuf[2], ibuf[0]);
				i = ibuf[0];
				while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
					if (buf[0] == '!')
						continue;
					if (buf[0] == '\n')
						break;
					if (sscanf(buf, "%lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2]) < 3) {
						s_append_asprintf(errbuf, "line %d: cannot read gaussian primitive coefficients", lineNumber);
						goto err_exit;
					}
					MoleculeAddGaussianPrimitiveCoefficients(mp, dbuf[0], dbuf[1], dbuf[2]);
					if (--i == 0)
						break;
				}
				if (buf[0] == '\n')
					break;
			}
			continue;
		} else if (strcmp(buf, "!:mo_info") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				if (sscanf(buf, "%6s %d %d", cbuf[0], &ibuf[0], &ibuf[1]) < 3) {
					s_append_asprintf(errbuf, "line %d: the MO info cannot be correctly read", lineNumber);
					goto err_exit;
				}
				if (strcasecmp(cbuf[0], "RHF") == 0) {
					ibuf[2] = 1;
				} else if (strcasecmp(cbuf[0], "ROHF") == 0) {
					ibuf[2] = 2;
				} else if (strcasecmp(cbuf[0], "UHF") == 0) {
					ibuf[2] = 0;
				} else {
					s_append_asprintf(errbuf, "line %d: unknown HF type: %s", lineNumber, cbuf[0]);
					goto err_exit;
				}
				if (ibuf[0] < 0 || ibuf[1] < 0) {
					s_append_asprintf(errbuf, "line %d: incorrect number of electrons", lineNumber);
					goto err_exit;
				}
				MoleculeSetMOInfo(mp, ibuf[2], ibuf[0], ibuf[1]);
			}
			continue;
		} else if (strcmp(buf, "!:mo_coefficients") == 0) {
			if (mp->bset == NULL || mp->bset->nshells == 0) {
				s_append_asprintf(errbuf, "line %d: the :gaussian_primitive section must come before :mo_coefficients", lineNumber);
				goto err_exit;
			}
			/*  Count the number of components  */
			dp = (Double *)malloc(sizeof(Double) * mp->bset->ncomps);
			i = 1;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				if (sscanf(buf, "MO %d %lf", &ibuf[0], &dbuf[6]) < 2) {
					s_append_asprintf(errbuf, "line %d: cannot read the MO index or energy", lineNumber);
					goto err_exit;
				}
				if (ibuf[0] != i) {
					s_append_asprintf(errbuf, "line %d: the MO index (%d) must be in ascending order", lineNumber, ibuf[0]);
					goto err_exit;
				}
				i = 0;
				while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
					j = sscanf(buf, "%lf %lf %lf %lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5]);
					if (j == 0) {
						s_append_asprintf(errbuf, "line %d: cannot read the MO coefficients", lineNumber);
						goto err_exit;
					}
					for (k = 0; k < j; k++, i++) {
						if (i >= mp->bset->ncomps) {
							s_append_asprintf(errbuf, "line %d: too many MO coefficients", lineNumber);
							goto err_exit;
						}
						dp[i] = dbuf[k];
					}
					if (i >= mp->bset->ncomps)
						break;
				}
				i = MoleculeSetMOCoefficients(mp, ibuf[0], dbuf[6], mp->bset->ncomps, dp);
				if (i != 0) {
					s_append_asprintf(errbuf, "line %d: cannot set MO coefficients", lineNumber);
					goto err_exit;
				}
				i = ibuf[0] + 1;  /*  For next entry  */
			}
			continue;
		} else if (strcmp(buf, "!:graphics") == 0) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				MainViewGraphic *gp = NULL;
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				if (mp->mview == NULL)
					continue;  /*  Skip  */
			redo:
				if (strcmp(buf, "line\n") == 0) {
					ibuf[0] = kMainViewGraphicLine;
				} else if (strcmp(buf, "poly\n") == 0) {
					ibuf[0] = kMainViewGraphicPoly;
				} else if (strcmp(buf, "cylinder\n") == 0) {
					ibuf[0] = kMainViewGraphicCylinder;
				} else if (strcmp(buf, "cone\n") == 0) {
					ibuf[0] = kMainViewGraphicCone;
				} else if (strcmp(buf, "ellipsoid\n") == 0) {
					ibuf[0] = kMainViewGraphicEllipsoid;
				} else {
					continue;  /*  Skip  */
				}
				gp = (MainViewGraphic *)calloc(sizeof(MainViewGraphic), 1);
				gp->kind = ibuf[0];
				i = 0;
				while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
					if (buf[0] == '!')
						continue;
					if (buf[0] == '\n')
						break;
					if (i == 0) {
						if (sscanf(buf, "%d %d", &ibuf[0], &ibuf[1]) < 2) {
							s_append_asprintf(errbuf, "line %d: the closed/visible flags cannot be read for graphic object", lineNumber);
							goto err_exit;
						}
						gp->closed = ibuf[0];
						gp->visible = ibuf[1];
					} else if (i == 1) {
						if (sscanf(buf, "%lf %lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3]) < 4) {
							s_append_asprintf(errbuf, "line %d: the color cannot be read for graphic object", lineNumber);
							goto err_exit;
						}
						for (j = 0; j < 4; j++)
							gp->rgba[j] = dbuf[j];
					} else if (i == 2) {
						j = atoi(buf);
						if (j < 0) {
							s_append_asprintf(errbuf, "line %d: the number of control points must be non-negative", lineNumber);
							goto err_exit;
						}
						if (j > 0)
							NewArray(&gp->points, &gp->npoints, sizeof(GLfloat) * 3, j);
					} else if (i >= 3 && i < gp->npoints + 3) {
						if (sscanf(buf, "%lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2]) < 3) {
							s_append_asprintf(errbuf, "line %d: the control point cannot be read for graphic object", lineNumber);
							goto err_exit;
						}
						j = (i - 3) * 3;
						gp->points[j++] = dbuf[0];
						gp->points[j++] = dbuf[1];
						gp->points[j] = dbuf[2];
					} else if (i == gp->npoints + 3) {
						j = atoi(buf);
						if (j < 0) {
							s_append_asprintf(errbuf, "line %d: the number of normals must be non-negative", lineNumber);
							goto err_exit;
						}
						if (j > 0)
							NewArray(&gp->normals, &gp->nnormals, sizeof(GLfloat) * 3, j);
					} else if (i >= gp->npoints + 4 && i < gp->npoints + gp->nnormals + 4) {
						if (sscanf(buf, "%lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2]) < 3) {
							s_append_asprintf(errbuf, "line %d: the normal vector cannot be read for graphic object", lineNumber);
							goto err_exit;
						}
						j = (i - gp->npoints - 4) * 3;
						gp->normals[j++] = dbuf[0];
						gp->normals[j++] = dbuf[1];
						gp->normals[j] = dbuf[2];
					} else break;
					i++;
				}
				MainView_insertGraphic(mp->mview, -1, gp);
				free(gp);
				if (buf[0] == '\n' || buf[0] == 0)
					break;
				goto redo;
			}
			continue;
		} else if (strncmp(buf, "!:@", 3) == 0) {
			/*  Plug-in implemented in the ruby world  */
			Int stringLen;
			char *stringBuf, *returnString;
			i = strlen(buf);
			NewArray(&stringBuf, &stringLen, sizeof(char), i + 1);
			strcpy(stringBuf, buf);
			k = lineNumber;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				/*  The comment lines are _not_ skipped  */
				if (buf[0] == '\n')
					break;
				j = strlen(buf);
				AssignArray(&stringBuf, &stringLen, sizeof(char), i + j, NULL);
				strncpy(stringBuf + i, buf, j);
				i += j;
			}
			if (MolActionCreateAndPerform(mp, SCRIPT_ACTION("si;s"),
										  "proc { |i| loadmbsf_plugin(i) rescue \"line #{i}: #{$i.to_s}\" }",
										  stringBuf, k, &returnString) != 0) {
				s_append_asprintf(errbuf, "line %d: cannot invoke Ruby plugin", lineNumber);
				goto err_exit;
			} else if (returnString[0] != 0) {
				s_append_asprintf(errbuf, "%s", returnString);
				goto err_exit;
			}
			free(stringBuf);
			continue;
		}
		/*  Unknown sections are silently ignored  */
	}

	MoleculeCleanUpResidueTable(mp);
	if (mp->arena != NULL)
		md_arena_set_molecule(mp->arena, mp);

	fclose(fp);

/*	if (mp->mview != NULL) {
		if (mview_ibuf[0] != kUndefined)
			mp->mview->showUnitCell = mview_ibuf[0];
		if (mview_ibuf[1] != kUndefined)
			mp->mview->showPeriodicBox = mview_ibuf[1];
		if (mview_ibuf[2] != kUndefined)
			mp->mview->showExpandedAtoms = mview_ibuf[2];
		if (mview_ibuf[3] != kUndefined)
			mp->mview->showEllipsoids = mview_ibuf[3];
		if (mview_ibuf[4] != kUndefined)
			mp->mview->showHydrogens = mview_ibuf[4];
		if (mview_ibuf[5] != kUndefined)
			mp->mview->showDummyAtoms = mview_ibuf[5];
		if (mview_ibuf[6] != kUndefined)
			mp->mview->showRotationCenter = mview_ibuf[6];
		if (mview_ibuf[7] != kUndefined)
			mp->mview->showGraphiteFlag = mview_ibuf[7];
		if (mview_ibuf[8] != kUndefined)
			mp->mview->showPeriodicImageFlag = mview_ibuf[8];
		if (mview_ibuf[9] != kUndefined)
			mp->mview->showGraphite = mview_ibuf[9];
		if (mview_ibuf[10] != kUndefined && mview_ibuf[10] >= 6)
			mp->mview->atomResolution = mview_ibuf[10];
		if (mview_ibuf[11] != kUndefined && mview_ibuf[11] >= 4)
			mp->mview->bondResolution = mview_ibuf[11];
		for (i = 0; i < 6; i++) {
			if (mview_ibuf[12 + i] != kUndefined)
				mp->mview->showPeriodicImage[i] = mview_ibuf[12 + i];
		}
		if (mview_dbuf[8] != kUndefined)
			mp->mview->atomRadius = mview_dbuf[8];
		if (mview_dbuf[9] != kUndefined)
			mp->mview->bondRadius = mview_dbuf[9];		
		if (mp->mview->track != NULL) {
			if (mview_dbuf[0] != kUndefined)
				TrackballSetScale(mp->mview->track, mview_dbuf[0]);
			if (mview_dbuf[1] != kUndefined)
				TrackballSetTranslate(mp->mview->track, mview_dbuf + 1);
			if (mview_dbuf[4] != kUndefined)
				TrackballSetRotate(mp->mview->track, mview_dbuf + 4);
		}
	}
*/

	return 0;

err_exit:
	fclose(fp);
	/*  The content of mp may be broken, so make it empty  */
	MoleculeClear(mp);
	return -1;	
}

int
MoleculeLoadPsfFile(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	char buf[1024];
	char *p;
	int section = -1;
	int i, j, err, fn;
	int lineNumber;
	Int ibuf[12];
	Vector *frames = NULL;
	Atom *ap;
	err = 0;
	*errbuf = NULL;
	if (mp == NULL)
		mp = MoleculeNew();
	else MoleculeClear(mp);
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot open file");
		return 1;
	}
/*	flockfile(fp); */
	lineNumber = 0;
	fn = 0;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		if (strncmp(buf, "PSF", 3) == 0) {
			section = 0;
			continue;
		} else {
			for (p = buf; *p != 0 && isspace(*p); p++) {}
			if (*p == 0) {
				section++;
				continue;
			}
		}
		if (strstr(buf, "!COORD") != NULL) {
			/*  Extended psf file with coordinates  */
			if (fn > 0) {
				/*  Allocate a temporary storage for frames  */
				size_t size = sizeof(Vector) * mp->natoms * fn;
				if (frames == NULL)
					frames = (Vector *)malloc(size);
				else
					frames = (Vector *)realloc(frames, size);
				if (frames == NULL)
					goto panic;
			}
			/*  Read coordinates  */
			for (i = 0; i < mp->natoms; i++) {
				double dval[3];
				Vector r;
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					err = 1;
					s_append_asprintf(errbuf, "line %d: premature end of file while reading coordinates (frame %d)", lineNumber, fn);
					goto exit;
				}
				if (sscanf(buf, "%lg %lg %lg", dval, dval + 1, dval + 2) != 3) {
					err = 1;
					s_append_asprintf(errbuf, "line %d: coordinates cannot be read for atom %d", lineNumber, i + 1);
					goto exit;
				}
				r.x = dval[0];
				r.y = dval[1];
				r.z = dval[2];
				if (fn == 0)
					ATOM_AT_INDEX(mp->atoms, i)->r = r;
				else
					frames[mp->natoms * (fn - 1) + i] = r;
			}
			fn++;
			continue;
		}
		
		if (section == 2) {
			/*  Atoms  */
			Int natoms;
			ReadFormat(buf, "I8", &natoms);
			if (natoms == 0)
				continue;
			if (NewArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, natoms) == NULL)
				goto panic;
			mp->nresidues = 0;
			for (i = 0; i < natoms; i++) {
				struct {
					char segName[5], resName[4], atomName[5], atomType[3], element[3];
					Int serial;
				} w;
				memset(&w, 0, sizeof(w));
				ap = ATOM_AT_INDEX(mp->atoms, i);
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					err = 1;
					s_append_asprintf(errbuf, "line %d: premature end of file while reading atoms", lineNumber);
					goto exit;
				}
				ReadFormat(buf, "I8 x1 S4 I5 x1 S3 x2 S4 x1 S4 F16 F10",
					&w.serial, w.segName, &ap->resSeq, w.resName, w.atomName, 
					w.atomType, &ap->charge, &ap->weight);
				strncpy(ap->segName, w.segName, 4);
				strncpy(ap->resName, w.resName, 3);
				strncpy(ap->aname, w.atomName, 4);
				ap->type = AtomTypeEncodeToUInt(w.atomType);
				/*  $element = ($name =~ /([A-Za-z]{1,2})/); # in Perl  */
				ap->atomicNumber = GuessAtomicNumber(w.atomName, ap->weight);
				ElementToString(ap->atomicNumber, w.element);
				strncpy(ap->element, w.element, 2);
			/*	w.element[0] = 0;
				for (p = w.atomName; *p != 0; p++) {
					if (isalpha(*p) && *p != '_') {
						w.element[0] = toupper(*p);
						if (isalpha(p[1]) && p[1] != '_') {
							w.element[1] = toupper(p[1]);
							w.element[2] = 0;
						} else {
							w.element[1] = 0;
						}
						break;
					}
				}
				strncpy(ap->element, w.element, 2);
				ap->atomicNumber = ElementToInt(w.element); */
				if (w.resName[0] == 0)
					strncpy(ap->resName, "XXX", 3);
				if (ap->resSeq > mp->nresidues)
					mp->nresidues = ap->resSeq;
			}
			if (mp->residues != NULL)
				free(mp->residues);
			if (NewArray(&mp->residues, &mp->nresidues, sizeof(char (*)[4]), mp->nresidues + 1) == 0)
				goto panic;
			for (i = 0; i < mp->natoms; i++) {
				j = mp->atoms[i].resSeq;
				if (mp->residues[j][0] == 0)
					strncpy(mp->residues[j], mp->atoms[i].resName, 4);
			}
			continue;
		} else if (section == 3) {
			/*  Bonds  */
			Int nbonds;
			Int *bp;
			ReadFormat(buf, "I8", &nbonds);
			if (nbonds == 0)
				continue;
			if (NewArray(&mp->bonds, &mp->nbonds, sizeof(Int) * 2, nbonds) == NULL)
				goto panic;
			bp = mp->bonds;
			for (i = 0; i < nbonds; i += 4) {
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					s_append_asprintf(errbuf, "line %d: premature end of file while reading bonds", lineNumber);
					err = 1;
					goto exit;
				}
				ReadFormat(buf, "I8I8I8I8I8I8I8I8", ibuf, ibuf + 1, ibuf + 2, ibuf + 3,
					ibuf + 4, ibuf + 5, ibuf + 6, ibuf + 7);
				for (j = 0; j < 4 && i + j < nbonds; j++) {
					Int b1, b2;
					Atom *ap;
					b1 = ibuf[j * 2] - 1;    /* Internal atom number is 0-based */
					b2 = ibuf[j * 2 + 1] - 1;
					if (b1 < 0 || b1 >= mp->natoms || b2 < 0 || b2 >= mp->natoms) {
						s_append_asprintf(errbuf, "line %d: The bond %d-%d includes non-existent atom", lineNumber, b1+1, b2+1);
						err = 1;
						goto exit;
					}
					*bp++ = b1;
					*bp++ = b2;
					ap = ATOM_AT_INDEX(mp->atoms, b1);
					AtomConnectInsertEntry(&ap->connect, -1, b2);
					ap = ATOM_AT_INDEX(mp->atoms, b2);
					AtomConnectInsertEntry(&ap->connect, -1, b1);
				}
			}
			continue;
		} else if (section == 4) {
			/*  Angles  */
			Int nangles;
			Int *gp;
			ReadFormat(buf, "I8", &nangles);
			if (nangles == 0)
				continue;
			if (NewArray(&mp->angles, &mp->nangles, sizeof(Int) * 3, nangles) == NULL)
				goto panic;
			gp = mp->angles;
			for (i = 0; i < nangles; i += 3) {
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					s_append_asprintf(errbuf, "line %d: premature end of file while reading angles", lineNumber);
					err = 1;
					goto exit;
				}
				ReadFormat(buf, "I8I8I8I8I8I8I8I8I8", ibuf, ibuf + 1, ibuf + 2, ibuf + 3,
					ibuf + 4, ibuf + 5, ibuf + 6, ibuf + 7, ibuf + 8);
				for (j = 0; j < 3 && i + j < nangles; j++) {
					Int a1, a2, a3;
					a1 = ibuf[j * 3] - 1;   /* Internal atom number is 0-based */
					a2 = ibuf[j * 3 + 1] - 1;
					a3 = ibuf[j * 3 + 2] - 1;
					if (a1 < 0 || a1 >= mp->natoms || a2 < 0 || a2 >= mp->natoms || a3 < 0 || a3 >= mp->natoms) {
						s_append_asprintf(errbuf, "line %d: The angle %d-%d-%d includes non-existent atom", lineNumber, a1+1, a2+1, a3+1);
						err = 1;
						goto exit;
					}
					*gp++ = a1;
					*gp++ = a2;
					*gp++ = a3;
				}
			}
			continue;
		} else if (section == 5 || section == 6) {
			/*  Dihedrals and Impropers  */
			Int ndihedrals;
			Int *dp;
			ReadFormat(buf, "I8", &ndihedrals);
			if (ndihedrals == 0)
				continue;
			if (section == 5) {
				if (NewArray(&mp->dihedrals, &mp->ndihedrals, sizeof(Int) * 4, ndihedrals) == NULL)
					goto panic;
				dp = mp->dihedrals;
			} else {
				if (NewArray(&mp->impropers, &mp->nimpropers, sizeof(Int) * 4, ndihedrals) == NULL)
					goto panic;
				dp = mp->impropers;
			}
			for (i = 0; i < ndihedrals; i += 2) {
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					fclose(fp);
					s_append_asprintf(errbuf, "line %d: premature end of file while reading %s", lineNumber, (section == 5 ? "dihedral" : "improper"));
					err = 1;
					goto exit;
				}
				ReadFormat(buf, "I8I8I8I8I8I8I8I8", ibuf, ibuf + 1, ibuf + 2, ibuf + 3, ibuf + 4, ibuf + 5, ibuf + 6, ibuf + 7);
				for (j = 0; j < 2 && i + j < ndihedrals; j++) {
					Int d1, d2, d3, d4;
					d1 = ibuf[j * 4] - 1;   /*  Internal atom number is 0-based  */
					d2 = ibuf[j * 4 + 1] - 1;
					d3 = ibuf[j * 4 + 2] - 1;
					d4 = ibuf[j * 4 + 3] - 1;
					if (d1 < 0 || d1 >= mp->natoms || d2 < 0 || d2 >= mp->natoms || d3 < 0 || d3 >= mp->natoms || d4 < 0 || d4 >= mp->natoms) {
						s_append_asprintf(errbuf, "line %d: The %s %d-%d-%d-%d angle includes non-existent atom", lineNumber, (section == 5 ? "dihedral" : "improper"), d1+1, d2+1, d3+1, d4+1);
						err = 1;
						goto exit;
					}
					*dp++ = d1;
					*dp++ = d2;
					*dp++ = d3;
					*dp++ = d4;
				}
			}
			continue;
		}
	}
	
	/*  Create frames for each atom if necessary  */
	if (fn > 1) {
		for (i = 0; i < mp->natoms; i++) {
			ap = ATOM_AT_INDEX(mp->atoms, i);
			NewArray(&ap->frames, &ap->nframes, sizeof(Vector), fn);
			if (ap->frames == NULL)
				goto panic;
			for (j = 0; j < fn; j++)
				ap->frames[j] = frames[mp->natoms * j + i];
		}
		free(frames);
		frames = NULL;
	}

  exit:
/*	funlockfile(fp); */
	fclose(fp);
	mp->nframes = -1;  /*  Should be recalculated later  */
	if (err)
		return 1;
	else if (section == -1)
		return -1;
	return 0;
  panic:
	Panic("low memory while reading structure file %s", fname);
	return 1; /* not reached */
}

/* ("-x", "y", "-z+0.5") -> (-1,0,0,0,0,1,0,0,0,0,-1,0.5)  */
static int
sMoleculeSymopStringsToTransform(char **symops, Transform tr)
{
	int i;
	char *symop;
	memset(tr, 0, sizeof(Transform));
	for (i = 0; i < 3; i++) {
		symop = symops[i];
		if (symop == NULL)
			return 1;
		while (*symop != 0) {
			int sn = 1;
			while (isspace(*symop))
				symop++;
			if (*symop == 0 || *symop == '\r' || *symop == 'n')
				break;
			if (*symop == '-') {
				sn = -1;
				symop++;
			} else if (*symop == '+') {
				sn = 1;
				symop++;
			}
			while (isspace(*symop))
				symop++;
			if (*symop == '.' || isdigit(*symop)) {
				/*  Numerical offset  */
				double d = strtod(symop, &symop);
				if (*symop == '/') {
					double dd = strtod(symop + 1, &symop);
					if (dd > 0)
						d /= dd;
					else
						return 1;  /*  Bad format  */
				}
				tr[9 + i] = d * sn;
			} else if (*symop == 'x' || *symop == 'X') {
				tr[i] = sn;
				symop++;
			} else if (*symop == 'y' || *symop == 'Y') {
				tr[i + 3] = sn;
				symop++;
			} else if (*symop == 'z' || *symop == 'Z') {
				tr[i + 6] = sn;
				symop++;
			} else return 1;  /*  Bad format  */
		} /* end while (*symop != 0) */
	}
	return 0;
}

static void
sMoleculeGenerateSymopWithTransform(Molecule *mp, Transform gtr, int num)
{
	int i, j;
	Transform tr;
	if (num <= 0)
		num = mp->nsyms;
	for (i = 0; i < num; i++) {
		memmove(tr, mp->syms[i], sizeof(Transform));
		TransformMul(tr, gtr, tr);
		for (j = 9; j < 12; j++) {
			if (tr[j] >= 1.0)
				tr[j] -= 1.0;
			else if (tr[j] <= 0.0)
				tr[j] += 1.0;
		}
		AssignArray(&mp->syms, &mp->nsyms, sizeof(Transform), mp->nsyms, tr);
	}
}

static char *
sChomp(char *buf)
{
	char *p = buf + strlen(buf) - 1;
	if (p >= buf && (*p == '\n' || *p == '\r')) {
		*p = 0;
		if (--p >= buf && (*p == '\n' || *p == '\r')) {
			*p = 0;
		}
	}
	return buf;
}

int
MoleculeLoadTepFile(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	char buf[1024];
	int section = -1;
	int lineNumber;
	int cellType;
	Int ibuf[12];
	Double fbuf[12];
	Int *bonds, nbonds;
	*errbuf = NULL;
	if (mp == NULL)
		mp = MoleculeNew();
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot open file");
		return 1;
	}
	lineNumber = 0;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		if (section == -1) {
			/*  Title  */
			section = 0;
			continue;
		}
		if (section == 0) {
			/*  XtalCell  */
			ReadFormat(buf, "I1F8F9F9F9F9F9", ibuf, fbuf, fbuf+1, fbuf+2, fbuf+3, fbuf+4, fbuf+5);
			cellType = ibuf[0];
			MoleculeSetCell(mp, fbuf[0], fbuf[1], fbuf[2], fbuf[3], fbuf[4], fbuf[5], 0);
			section = 1;
			continue;
		}
		if (section == 1) {
			/*  Symmetry  */
			Transform tr;
			if (cellType == 0) {
				ReadFormat(buf, "I1F14F3F3F3F15F3F3F3F15F3F3F3", ibuf, fbuf, fbuf+1, fbuf+2, fbuf+3, fbuf+4, fbuf+5, fbuf+6, fbuf+7, fbuf+8, fbuf+9, fbuf+10, fbuf+11);
				tr[0] = fbuf[1];
				tr[3] = fbuf[2];
				tr[6] = fbuf[3];
				tr[1] = fbuf[5];
				tr[4] = fbuf[6];
				tr[7] = fbuf[7];
				tr[2] = fbuf[9];
				tr[5] = fbuf[10];
				tr[8] = fbuf[11];
				tr[9] = fbuf[0];
				tr[10] = fbuf[4];
				tr[11] = fbuf[8];
			} else {
				char *symops[3], *brks;
				sChomp(buf);
				memset(tr, 0, sizeof(Transform));
				ReadFormat(buf, "I1", ibuf);
				symops[0] = strtok_r(buf + 1, ", ", &brks);
				symops[1] = strtok_r(NULL, ", ", &brks);
				symops[2] = strtok_r(NULL, ", ", &brks);
				if (sMoleculeSymopStringsToTransform(symops, tr)) {
					s_append_asprintf(errbuf, "line %d: bad symmetry specification", lineNumber);
					return 1;
				}
			}
			if (AssignArray(&mp->syms, &mp->nsyms, sizeof(Transform), mp->nsyms, tr) == 0)
				goto panic;
			if (ibuf[0] != 0)
				section = 2;
			continue;
		}
		if (section == 2) {	 /*  Atoms  */
			char name[8];
			Atom *ap;
			int atomType;
			int atomIndex = mp->natoms;
			ap = AssignArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, atomIndex, NULL);
			memset(ap, 0, gSizeOfAtomRecord);
			ReadFormat(buf, "S6x3x9x9F9F9F9F9", name, fbuf, fbuf+1, fbuf+2, fbuf+3);
			strncpy(ap->aname, name, 4);
			ap->r.x = fbuf[0];
			ap->r.y = fbuf[1];
			ap->r.z = fbuf[2];
			MoleculeXtalToCartesian(mp, &(ap->r), &(ap->r));
		/*	ap->atomicNumber = AtomNameToElement(ap->name);
			ElementToString(ap->atomicNumber, ap->element); */
		/*	sAtomSetElement(ap, -1, ap->name); */
			guessElement(ap);
			atomType = fbuf[3];
			if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
				s_append_asprintf(errbuf, "unexpected end of file");
				return 1;
			}
			ReadFormat(buf, "I1F8F9F9F9F9F9F9", ibuf, fbuf, fbuf+1, fbuf+2, fbuf+3, fbuf+4, fbuf+5, fbuf+6);
			atomType = fbuf[6];
			if ((atomType >= 0 && atomType <= 5) || (atomType >= 8 && atomType <= 10)) { 
				/*  Anisotropic thermal parameters  */
				MoleculeSetAniso(mp, atomIndex, atomType, fbuf[0], fbuf[1], fbuf[2], fbuf[3], fbuf[5], fbuf[4], NULL);
			}
			if (ibuf[0] != 0)
				section = 3;
			continue;
		}
	}
	fclose(fp);
	MoleculeGuessBonds(mp, 0.0, &nbonds, &bonds);
	if (nbonds > 0) {
		MoleculeAddBonds(mp, nbonds, bonds, NULL, 1);
		free(bonds);
	}
	mp->nframes = -1;  /*  Should be recalculated later  */
	return 0;
  panic:
	Panic("low memory while reading structure file %s", fname);
	return -1; /* not reached */
}

int
MoleculeLoadShelxFile(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	char buf[1024];
	char *p1, *p2;
	int n;
	int lineNumber;
	int latticeType;
	int currentResSeq = 0;
	char currentResName[6];
	Transform tr;
	int ibuf[12];
	float fbuf[12];
	Double dbuf[12];
	Int nsfacs = 0;
	Int nbonds, *bonds;
	char (*sfacs)[4] = NULL;

	*errbuf = NULL;
	if (mp == NULL)
		mp = MoleculeNew();
	currentResName[0] = 0;
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot open file");
		return 1;
	}
	lineNumber = 0;
	tr[0] = tr[4] = tr[8] = 1;
	tr[1] = tr[2] = tr[3] = tr[5] = tr[6] = tr[7] = tr[9] = tr[10] = tr[11] = 0;
	if (AssignArray(&mp->syms, &mp->nsyms, sizeof(Transform), 0, tr) == 0)
		goto panic;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		if (strncmp(buf, "CELL", 4) == 0) {
			/*  XtalCell  */
			sscanf(buf + 4, " %f %f %f %f %f %f %f", fbuf, fbuf+1, fbuf+2, fbuf+3, fbuf+4, fbuf+5, fbuf+6);
			MoleculeSetCell(mp, fbuf[1], fbuf[2], fbuf[3], fbuf[4], fbuf[5], fbuf[6], 0);
			continue;
		} else if (strncmp(buf, "SFAC", 4) == 0) {
			sChomp(buf);
			for (p1 = strtok_r(buf + 4, " ", &p2); p1 != NULL; p1 = strtok_r(NULL, " ", &p2)) {
				char *pp = (char *)AssignArray(&sfacs, &nsfacs, 4, nsfacs, NULL);
				if (pp == NULL)
					goto panic;
				strncpy(pp, p1, 3);
				pp[3] = 0;
			}
			continue;
		} else if (strncmp(buf, "LATT", 4) == 0) {
			sscanf(buf + 4, " %d", &latticeType);
			continue;
		} else if (strncmp(buf, "SYMM", 4) == 0) {
			char *symops[3], *brks;
			memset(tr, 0, sizeof(Transform));
		//	ReadFormat(buf + 4, "I1", ibuf);
			sChomp(buf);
			symops[0] = strtok_r(buf + 4, ",", &brks);
			symops[1] = strtok_r(NULL, ",", &brks);
			symops[2] = strtok_r(NULL, ",", &brks);
			if (sMoleculeSymopStringsToTransform(symops, tr)) {
				s_append_asprintf(errbuf, "line %d: bad symmetry specification", lineNumber);
				return 1;
			}
			if (AssignArray(&mp->syms, &mp->nsyms, sizeof(Transform), mp->nsyms, tr) == 0)
				goto panic;
			continue;
		} else if (strncmp(buf, "RESI", 4) == 0) {
			for (p1 = buf + 4; isspace(*p1); p1++);
			if (isalpha(*p1)) {
				for (p2 = p1 + 1; isalnum(*p2); p2++);
				*p2 = 0;
				strncpy(currentResName, p1, 4);
				currentResName[4] = 0;
				p1 = p2 + 1;
			} else currentResName[0] = 0;
			sscanf(buf + 4, " %d", &currentResSeq);
			continue;
		} else {
			/* Atom name: [A-Za-z]{1,2}[0-9]*  */
			for (p1 = buf; p1 < buf + 2 && (isalpha(*p1) && *p1 != '_'); p1++);
			if (p1 > buf) {
				while (isdigit(*p1))
					p1++;
			}
			if (p1 > buf && p1 <= buf + 4 && isspace(*p1)) {
				/*  Atom  */
				Atom *ap;
				char cont[4];
				int atomIndex = mp->natoms;
				ap = AssignArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, atomIndex, NULL);
				memset(ap, 0, gSizeOfAtomRecord);
				strncpy(ap->aname, buf, 4);
				n = sscanf(p1, " %d %f %f %f %f %f %f %2s", ibuf, fbuf, fbuf+1, fbuf+2, fbuf+3, fbuf+4, fbuf+5, cont);
				if (n == 8 && strcmp(cont, "=") == 0) {
					if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
						s_append_asprintf(errbuf, "line %d: unexpected end of file within the atom cards", lineNumber);
						return 1;
					}
					sscanf(buf, " %f %f %f %f", fbuf+6, fbuf+7, fbuf+8, fbuf+9);
					n = 10;   /*  Aniso  */
				} else n = 5; /*  Iso  */
				ap->r.x = fbuf[0];
				ap->r.y = fbuf[1];
				ap->r.z = fbuf[2];
				MoleculeXtalToCartesian(mp, &(ap->r), &(ap->r));
				ap->occupancy = fbuf[3];
				if (ap->aname[0] != 'Q' && ibuf[0] >= 1 && ibuf[0] <= nsfacs) {
					strncpy(ap->element, sfacs[ibuf[0] - 1], 2);
					ap->element[2] = 0;
				/*	sAtomSetElement(ap, -1, sfacs[ibuf[0] - 1]); */
				/*	strncpy(ap->element, sfacs[ibuf[0] - 1], 4);
					ap->atomicNumber = ElementToInt(ap->element); */
			/*	} else {
					sAtomSetElement(ap, -1, ap->name); */
				/*	ap->atomicNumber = AtomNameToElement(ap->name);
					ElementToString(ap->atomicNumber, ap->element); */
				}
				guessElement(ap);
				if (n == 10 || fbuf[4] >= 0.0) {
					int i, c, j;
					/*  Read in the standard deviations  */
					ReadLine(buf, sizeof buf, fp, &lineNumber);
					for (i = 0; i < 9; i++) {
						j = 3 + i * 8;
						c = buf[j + 8];
						buf[j + 8] = 0;
						dbuf[i] = strtod(buf + j, NULL);
						buf[j + 8] = c;
					}
					ap->sigma.x = dbuf[0];
					ap->sigma.y = dbuf[1];
					ap->sigma.z = dbuf[2];
				}
				if (n == 5)
					ap->tempFactor = fbuf[4] * 78.9568352087147; /* 8*pi*pi */
				else
					MoleculeSetAniso(mp, atomIndex, 8, fbuf[4], fbuf[5], fbuf[6], fbuf[9], fbuf[7], fbuf[8], dbuf);
				ap->resSeq = currentResSeq;
				strncpy(ap->resName, currentResName, 4);
			}
			continue;
		}
	}
	fclose(fp);

	/*  Add symmetry operations according to the lattice type  */
	switch (latticeType < 0 ? -latticeType : latticeType) {
		static Transform tr_i = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0.5, 0.5, 0.5};
		static Transform tr_c = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0.5, 0.5, 0};
		static Transform tr_a = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0.5, 0.5};
		static Transform tr_b = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0.5, 0, 0.5};
		static Transform tr_r1 = {0, 1, 0, -1, -1, 0, 0, 0, 1, 0, 0, 0};
		static Transform tr_r2 = {-1, -1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0};
		case 1:  /* P */
			break;
		case 2:  /* I */
			sMoleculeGenerateSymopWithTransform(mp, tr_i, 0);
			break;
		case 3:  /* Rhombohedral obverse on hexagonal axes  */
			n = mp->nsyms;
			sMoleculeGenerateSymopWithTransform(mp, tr_r1, n);
			sMoleculeGenerateSymopWithTransform(mp, tr_r2, n);
			break;
		case 4:  /* F */
			n = mp->nsyms;
			sMoleculeGenerateSymopWithTransform(mp, tr_a, n);
			sMoleculeGenerateSymopWithTransform(mp, tr_b, n);
			sMoleculeGenerateSymopWithTransform(mp, tr_c, n);
			break;
		case 5:  /* A */
			sMoleculeGenerateSymopWithTransform(mp, tr_a, 0);
			break;
		case 6:  /* B */
			sMoleculeGenerateSymopWithTransform(mp, tr_b, 0);
			break;
		case 7:  /* C */
			sMoleculeGenerateSymopWithTransform(mp, tr_c, 0);
			break;
	}
		
	if (latticeType > 0) {
		static Transform tr_inv = {-1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0};
		sMoleculeGenerateSymopWithTransform(mp, tr_inv, 0);
	}
	
	MoleculeGuessBonds(mp, 0.0, &nbonds, &bonds);
	if (nbonds > 0) {
		MoleculeAddBonds(mp, nbonds, bonds, NULL, 1);
		free(bonds);
	}
	mp->nframes = -1;  /*  Should be recalculated later  */
	return 0;
  panic:
	Panic("low memory while reading structure file %s", fname);
	return -1; /* not reached */
}

/*  Add one gaussian orbital shell information (not undoable)  */
int
MoleculeAddGaussianOrbitalShell(Molecule *mol, Int a_idx, Int sym, Int nprims)
{
	BasisSet *bset;
	ShellInfo *shellp;
	if (mol == NULL)
		return -1;  /*  Molecule is empty  */
	bset = mol->bset;
	if (bset == NULL) {
		bset = mol->bset = (BasisSet *)calloc(sizeof(BasisSet), 1);
		if (bset == NULL)
			return -2;  /*  Low memory  */
	}
	shellp = AssignArray(&bset->shells, &bset->nshells, sizeof(ShellInfo), bset->nshells, NULL);
	if (shellp == NULL)
		return -2;  /*  Low memory  */
	switch (sym) {
		case 0:  shellp->sym = kGTOType_S;  shellp->ncomp = 1; break;
		case 1:  shellp->sym = kGTOType_P;  shellp->ncomp = 3; break;
		case -1: shellp->sym = kGTOType_SP; shellp->ncomp = 4; break;
		case 2:  shellp->sym = kGTOType_D;  shellp->ncomp = 6; break;
		case -2: shellp->sym = kGTOType_D5; shellp->ncomp = 5; break;
		case 3:  shellp->sym = kGTOType_F;  shellp->ncomp = 10; break;
		case -3: shellp->sym = kGTOType_F7; shellp->ncomp = 7; break;
		case 4:  shellp->sym = kGTOType_G;  shellp->ncomp = 15; break;
		case -4: shellp->sym = kGTOType_G9; shellp->ncomp = 9; break;
		default:
			return -3;  /* Unsupported shell type  */
	}
	shellp->nprim = nprims;
	shellp->a_idx = a_idx;
	if (bset->shells < shellp) {
		shellp->m_idx = shellp[-1].m_idx + shellp[-1].ncomp;
		shellp->p_idx = shellp[-1].p_idx + shellp[-1].nprim;
	} else {
		shellp->m_idx = 0;
		shellp->p_idx = 0;
	}
	/*  Update the number of components (if not yet determined)  */
	if (bset->ncomps < shellp->m_idx + shellp->ncomp)
		bset->ncomps = shellp->m_idx + shellp->ncomp;
	return 0;
}

/*  Add a set of gaussian primitive coefficients (not undoable)  */
int
MoleculeAddGaussianPrimitiveCoefficients(Molecule *mol, Double exponent, Double contraction, Double contraction_sp)
{
	BasisSet *bset;
	PrimInfo *primp;
	if (mol == NULL)
		return -1;  /*  Molecule is empty  */
	bset = mol->bset;
	if (bset == NULL) {
		bset = mol->bset = (BasisSet *)calloc(sizeof(BasisSet), 1);
		if (bset == NULL)
			return -2;  /*  Low memory  */
	}
	primp = AssignArray(&bset->priminfos, &bset->npriminfos, sizeof(PrimInfo), bset->npriminfos, NULL);
	if (primp == NULL)
		return -2;  /*  Low memory  */
	primp->A = exponent;
	primp->C = contraction;
	primp->Csp = contraction_sp;
	return 0;
}

/*  Get the shell information from the component index  */
/*  The outLabel must have space for at least 23 non-Null characters  */
int
MoleculeGetGaussianComponentInfo(Molecule *mol, Int comp_idx, Int *outAtomIdx, char *outLabel, Int *outShellIdx)
{
	BasisSet *bset;
	ShellInfo *shellp;
	int si;
	if (mol == NULL || (bset = mol->bset) == NULL)
		return -1;  /*  No basis set info  */
	if (comp_idx < 0 || comp_idx >= bset->ncomps)
		return -2;  /*  Component index out of range  */
	for (si = 0, shellp = bset->shells; si < bset->nshells; si++, shellp++) {
		if (comp_idx >= shellp->ncomp) {
			comp_idx -= shellp->ncomp;
			continue;
		} else {
			static const char *type_p = "xyz";
			static const char *type_d = "xxyyzzxyxzyz";
			static const char *type_d5[] = {"xy","yz","zz", "xz", "xx-yy"};
			static const char *type_f = "xxxyyyzzzxxyxxzxyyyyzxzzyzzxyz";
			static const char *type_f7[] = {"x3-3xy2", "x2z-y2z", "x(5z2-r2)", "z(5z2-3r2)", "y(5z2-r2)", "xyz", "3x2y-y3"};
			static const char *type_g[] = {"x4", "y4", "z4", "x3y", "x3z", "xy3", "y3z", "xz3", "yz3", "x2y2", "x2z2", "y2z2", "x2yz", "x2yz", "xyz2"};
			static const char *type_g9[] = {"x4+y4-6x2y2", "xz(x2-3y2)", "(x2-y2)(7z2-r2)", "xz(7z2-3r2)", "35z4-30z2r2+3r4", "yz(7z2-3r2)", "xy(7z2-r2)", "yz(3x2-y2)", "xy(x2-y2)"};
			*outAtomIdx = shellp->a_idx;
			*outShellIdx = si;
			switch (shellp->sym) {
				case kGTOType_S:
					strcpy(outLabel, "S");
					break;
				case kGTOType_P:
					outLabel[0] = 'P';
					outLabel[1] = type_p[comp_idx];
					outLabel[2] = 0;
					break;
				case kGTOType_SP:
					if (comp_idx == 0)
						strcpy(outLabel, "S");
					else {
						outLabel[0] = 'P';
						outLabel[1] = type_p[comp_idx - 1];
						outLabel[2] = 0;
					}
					break;
				case kGTOType_D:
					outLabel[0] = 'D';
					strncpy(outLabel + 1, type_d + comp_idx * 2, 2);
					outLabel[3] = 0;
					break;
				case kGTOType_D5:
					outLabel[0] = 'D';
					strcpy(outLabel + 1, type_d5[comp_idx]);
					break;
				case kGTOType_F:
					outLabel[0] = 'F';
					strncpy(outLabel + 1, type_f + comp_idx * 3, 3);
					outLabel[4] = 0;
					break;
				case kGTOType_F7:
					outLabel[0] = 'F';
					strcpy(outLabel + 1, type_f7[comp_idx]);
					break;
				case kGTOType_G:
					outLabel[0] = 'G';
					strcpy(outLabel + 1, type_g[comp_idx]);
					break;
				case kGTOType_G9:
					outLabel[0] = 'G';
					strcpy(outLabel + 1, type_g9[comp_idx]);
					break;
				default:
					return -3;  /*  Unsupported orbital type (internal error) */
			}
			return 0;
		}
	}
	return -4;  /*  comp_idx out of range? (internal error)  */
}

/*  Set MO coefficients for idx-th MO (1-based)  */
int
MoleculeSetMOCoefficients(Molecule *mol, Int idx, Double energy, Int ncomps, Double *coeffs)
{
	BasisSet *bset;
	int i, n;
	if (mol == NULL)
		return -1;  /*  Molecule is empty  */
	bset = mol->bset;
	if (bset == NULL) {
		bset = mol->bset = (BasisSet *)calloc(sizeof(BasisSet), 1);
		if (bset == NULL)
			return -2;  /*  Low memory  */
	}
	if (bset->nmos == 0) {
		if (bset->nshells > 0) {
			/*  Shell info is already set: calculate the number of MOs from there  */
			for (i = n = 0; i < bset->nshells; i++)
				n += bset->shells[i].ncomp;
			bset->ncomps = n;
		} else if (ncomps > 0) {
			bset->ncomps = ncomps;
		}
		if (bset->rflag == 0)
			bset->nmos = bset->ncomps * 2;
		else
			bset->nmos = bset->ncomps;
		if (bset->nmos <= 0)
			return -3;  /*  Bad or inconsistent number of MOs  */
		bset->mo = (Double *)calloc(sizeof(Double), (bset->nmos + 1) * bset->ncomps);
		bset->moenergies = (Double *)calloc(sizeof(Double), bset->nmos + 1);
		if (bset->mo == NULL || bset->moenergies == NULL) {
			if (bset->mo != NULL)
				free(bset->mo);
			if (bset->moenergies != NULL)
				free(bset->moenergies);
			bset->mo = NULL;
			bset->moenergies = NULL;
			bset->nmos = 0;
			return -2;  /*  Low memory  */
		}
	}
	if (idx < 0)
		idx = -idx + bset->ncomps;
	if (idx < 0 || idx > bset->nmos)
		return -4;  /*  Bad MO index  */
	if (idx == 0)
		idx = bset->nmos;  /*  Arbitrary vector  */
	else
		idx--;
	if (energy != -1000000)
		bset->moenergies[idx] = energy;
	if (ncomps < bset->ncomps)
		return -5;  /*  Insufficient number of data provided  */
	memmove(bset->mo + (idx * bset->ncomps), coeffs, sizeof(Double) * bset->ncomps);
	if (bset->cns != NULL) {
		/*  Clear the cached values  */
		free(bset->cns);
		bset->cns = NULL;
		bset->ncns = 0;
	}
	return 0;
}

/*  Get MO coefficients for idx-th MO (1-based)  */
/*  Caution: *ncoeffs and *coeffs should be valid _before_ calling this function, i.e.  */
/*  *ncoeffs = 0 && *coeffs = NULL or *coeffs is a valid memory pointer and *ncoeffs  */
/*  properly designates the memory size as an array of Doubles.  */
int
MoleculeGetMOCoefficients(Molecule *mol, Int idx, Double *energy, Int *ncoeffs, Double **coeffs)
{
	BasisSet *bset;
	if (mol == NULL)
		return -1;  /*  Molecule is empty  */
	bset = mol->bset;
	if (bset == NULL || bset->ncomps <= 0)
		return -2;  /*  No basis set info  */
	if (idx < 0)
		idx = -idx + bset->ncomps;
	if (idx < 0 || idx > bset->nmos)
		return -3;  /*  MO index out of range  */
	if (idx == 0)
		idx = bset->nmos;  /*  Arbitrary vector  */
	else
		idx--;
	if (energy != NULL)
		*energy = bset->moenergies[idx];
	if (ncoeffs != NULL && coeffs != NULL) {
		if (*ncoeffs < bset->ncomps || *coeffs == NULL) {
			if (*coeffs != NULL)
				free(*coeffs);  /*  Caution: possible cause of SIGBUS if *coeff is not initialized properly */
			*coeffs = (Double *)calloc(sizeof(Double), bset->ncomps);
			*ncoeffs = bset->ncomps;
		}
		memmove(*coeffs, bset->mo + (idx * bset->ncomps), sizeof(Double) * bset->ncomps);
	}
	return 0;
}

/*  Set Basic MO Info. rflag: 0, UHF; 1, RHF; 2, ROHF; -1, clear
    ne_alpha: number of alpha electrons, ne_beta: number of beta electrons   */
int
MoleculeSetMOInfo(Molecule *mol, Int rflag, Int ne_alpha, Int ne_beta)
{
	BasisSet *bset;
	if (mol == NULL || mol->natoms == 0)
		return -1;  /*  Molecule is empty  */
	if (rflag < 0) {
		if (mol->bset != NULL) {
			BasisSetRelease(mol->bset);
			mol->bset = NULL;
		}
		return 0;
	}
	bset = mol->bset;
	if (bset == NULL) {
		bset = mol->bset = (BasisSet *)calloc(sizeof(BasisSet), 1);
		if (bset == NULL)
			return -2;  /*  Low memory  */
	}
	bset->natoms_bs = mol->natoms;
	bset->ne_alpha = ne_alpha;
	bset->ne_beta = ne_beta;
	bset->rflag = rflag;
	return 0;
}

static void
sSeparateTokens(char *inString, char **outPtr, int size)
{
	char *p;
	int i;
	for (i = 0; i < size; i++) {
		p = strtok((i == 0 ? inString : NULL), " \r\n");
		if (p == NULL)
			break;
		outPtr[i] = p;
	}
	while (i < size) {
		outPtr[i++] = NULL;
	}
}

static int
sReadNumberArray(void *basep, Int *countp, Int size, Int num, FILE *fp, int *lnp)
{
	char buf[256];
	Int i, n;
	*((void **)basep) = NULL;
	*countp = 0;
	if (AssignArray(basep, countp, size, num - 1, NULL) == NULL)
		return 4;  /*  Out of memory  */
	n = 0;
	while (ReadLine(buf, sizeof buf, fp, lnp) > 0) {
		char *tokens[16], *p;
		sSeparateTokens(buf, tokens, 16);
		for (i = 0; i < 16; i++) {
			if (tokens[i] == NULL)
				break;
			if (size == sizeof(Int)) {
				(*((Int **)basep))[n] = strtol(tokens[i], &p, 0);
			} else if (size == sizeof(Double)) {
				(*((Double **)basep))[n] = strtod(tokens[i], &p);
			} else return -1;  /*  Internal error  */
			if (tokens[i] == p || *p != 0)
				return 1;  /*  Non-digit character  */
			if (++n == num) {
				if (i < 15 && tokens[i + 1] != NULL)
					return 2;  /*  Too many data  */
				return 0;  /*  All data are successfully read  */
			}
		}
	}
	return 3;  /*  Unexpected EOF  */			
}

static int
sSetupGaussianCoefficients(BasisSet *bset)
{
	ShellInfo *sp;
	PrimInfo *pp;
	int i, j, k;
	Double *dp, d;
	
	/*  Cache the contraction coefficients for efficient calculation  */
	/*  Sum up the number of components for all primitives  */
	for (i = k = 0, sp = bset->shells; i < bset->nshells; i++, sp++) {
		sp->cn_idx = k;
		k += sp->nprim * sp->ncomp;
	}
	/*  Allocate memory for the cached values  */
	if (AssignArray(&bset->cns, &bset->ncns, sizeof(Double), k - 1, NULL) == NULL)
		return 1;
	/*  Iterate over all primitives  */
	dp = bset->cns;
	for (i = 0, sp = bset->shells; i < bset->nshells; i++, sp++) {
		for (j = 0, pp = bset->priminfos + sp->p_idx; j < sp->nprim; j++, pp++) {
			switch (sp->sym) {
				case kGTOType_S:
					// (8 alpha^3/pi^3)^0.25 exp(-alpha r^2)
					*dp++ = pp->C * pow(pp->A, 0.75) * 0.71270547;
					break;
				case kGTOType_P:
					// (128 alpha^5/pi^3)^0.25 [x|y|z]exp(-alpha r^2)
					d = pp->C * pow(pp->A, 1.25) * 1.425410941;
					*dp++ = d;
					*dp++ = d;
					*dp++ = d;
					break;
				case kGTOType_SP:
					*dp++ = pp->C * pow(pp->A, 0.75) * 0.71270547;
					d = pp->Csp * pow(pp->A, 1.25) * 1.425410941;
					*dp++ = d;
					*dp++ = d;
					*dp++ = d;
					break;
				case kGTOType_D:
					//  xx|yy|zz: (2048 alpha^7/9pi^3)^0.25 [xx|yy|zz]exp(-alpha r^2)
					//  xy|yz|zx: (2048 alpha^7/pi^3)^0.25 [xy|xz|yz]exp(-alpha r^2)
					d = pp->C * pow(pp->A, 1.75);
					dp[0] = dp[1] = dp[2] = d * 1.645922781;
					dp[3] = dp[4] = dp[5] = d * 2.850821881;
					dp += 6;
					break;
				case kGTOType_D5:
					//  3zz-rr:   (128 alpha^7/9pi^3)^0.25 (3zz-rr)exp(-alpha r^2)
					//  xy|yz|zx: (2048 alpha^7/pi^3)^0.25 [xy|xz|yz]exp(-alpha r^2)
					//  xx-yy:    (128 alpha^7/pi^3)^0.25 (xx-yy)exp(-alpha r^2)
					d = pp->C * pow(pp->A, 1.75);
					dp[0] = d * 0.822961390;
					dp[1] = dp[2] = dp[4] = d * 2.850821881;
					dp[3] = d * 1.425410941;
					dp += 5;
					break;
				/*  TODO: Support F/F7 and G/G9 type orbitals  */
			}
		}
	}
	return 0;
}

int
MoleculeLoadGaussianFchkFile(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	char buf[1024];
	int lineNumber;
	int natoms, nbasis, i, j, k, n, mxbond, retval, ncomps, nprims, nelec;
	BasisSet *bset;
	ShellInfo *sp;
	PrimInfo *pp;
	Int nary;
	Int *iary;
	Double *dary;
	Atom *ap;
/*	Vector *vp; */
	Double w;

	*errbuf = NULL;
	if (mp == NULL)
		mp = MoleculeNew();
	bset = (BasisSet *)calloc(sizeof(BasisSet), 1);
	if (bset == NULL)
		goto panic;
	mp->bset = bset;
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot open file");
		return 1;
	}
	lineNumber = 0;
	natoms = nbasis = -1;
	mxbond = 0;
	ncomps = 0;
	nelec = 0;
	nprims = 0;
	nary = 0;
	iary = NULL;
	dary = NULL;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		char *tokens[16];
		char *p = buf + 41;
		if (lineNumber == 2) {
			/*  job info line  */
			if (buf[10] == 'U')
				bset->rflag = 0;  /*  UHF  */
			else if (buf[11] == 'O')
				bset->rflag = 2;  /*  ROHF  */
			else bset->rflag = 1; /*  RHF  */
			continue;
		}
		while (p > buf && *p == ' ')
			p--;
		p[1] = 0;
		sSeparateTokens(buf + 42, tokens, 16);
		if (strcmp(buf, "Number of atoms") == 0) {
			if (tokens[1] == NULL || (natoms = atoi(tokens[1])) <= 0) {
				s_append_asprintf(errbuf, "Line %d: strange number of atoms: %s", lineNumber, tokens[1]);
				retval = 2;
				goto cleanup;
			}
			bset->natoms_bs = natoms;
			/*  Allocate atom records (all are empty for now)  */
			AssignArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, natoms - 1, NULL);
			/*  Also allocate atom position array for MO calculations  */
		/*	AssignArray(&bset->pos, &bset->natoms, sizeof(Vector), natoms - 1, NULL); */
			/*  Also allocate nuclear charge array  */
			bset->nuccharges = (Double *)calloc(sizeof(Double), natoms);
		} else if (strcmp(buf, "Number of electrons") == 0) {
			if (tokens[1] == NULL || (i = atoi(tokens[1])) < 0) {
				s_append_asprintf(errbuf, "Line %d: strange number of electrons: %s", lineNumber, tokens[1]);
				retval = 2;
				goto cleanup;
			}
			nelec = i;
		} else if (strcmp(buf, "Number of alpha electrons") == 0) {
			if (tokens[1] == NULL || (i = atoi(tokens[1])) < 0) {
				s_append_asprintf(errbuf, "Line %d: strange number of alpha electrons: %s", lineNumber, tokens[1]);
				retval = 2;
				goto cleanup;
			}
			bset->ne_alpha = i;
		} else if (strcmp(buf, "Number of beta electrons") == 0) {
			if (tokens[1] == NULL || (i = atoi(tokens[1])) < 0) {
				s_append_asprintf(errbuf, "Line %d: strange number of beta electrons: %s", lineNumber, tokens[1]);
				retval = 2;
				goto cleanup;
			}
			bset->ne_beta = i;
			if (bset->ne_alpha + bset->ne_beta != nelec) {
				s_append_asprintf(errbuf, "Line %d: sum of alpha (%d) and beta (%d) electrons does not match the number of electrons (%d)", lineNumber, (int)bset->ne_alpha, (int)bset->ne_beta, (int)nelec);
				retval = 2;
				goto cleanup;
			}
		} else if (strcmp(buf, "Number of basis functions") == 0) {
			if (tokens[1] == NULL || (nbasis = atoi(tokens[1])) <= 0) {
				s_append_asprintf(errbuf, "Line %d: strange number of basis functions: %s", lineNumber, tokens[1]);
				retval = 2;
				goto cleanup;
			}
		} else if (strcmp(buf, "Atomic numbers") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != natoms) {
				s_append_asprintf(errbuf, "Line %d: wrong or inconsistent number of atoms: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&iary, &nary, sizeof(Int), natoms, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read atomic numbers", lineNumber);
				retval = 2;
				goto cleanup;
			}
			for (i = 0, ap = mp->atoms; i < natoms; i++, ap = ATOM_NEXT(ap)) {
				ap->atomicNumber = iary[i];
				bset->nuccharges[i] = iary[i];
				ElementToString(ap->atomicNumber, ap->element);
				memmove(ap->aname, ap->element, 4);
				if ((w = WeightForAtomicNumber(ap->atomicNumber)) > 0.0)
					ap->weight = w;
			}
			free(iary);
			iary = NULL;
		} else if (strcmp(buf, "Nuclear charges") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != natoms) {
				s_append_asprintf(errbuf, "Line %d: wrong or inconsistent number of atoms: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), natoms, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read nuclear charges", lineNumber);
				retval = 2;
				goto cleanup;
			}
			for (i = 0, ap = mp->atoms; i < natoms; i++, ap = ATOM_NEXT(ap)) {
				bset->nuccharges[i] = dary[i];
			}
			free(iary);
			iary = NULL;
		} else if (strcmp(buf, "Current cartesian coordinates") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != natoms * 3) {
				s_append_asprintf(errbuf, "Line %d: wrong or inconsistent number of cartesian coordinates: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), natoms * 3, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read cartesian coordinates", lineNumber);
				retval = 2;
				goto cleanup;
			}
			for (i = 0, ap = mp->atoms; i < natoms; i++, ap = ATOM_NEXT(ap)) {
				ap->r.x = dary[i * 3] * kBohr2Angstrom;
				ap->r.y = dary[i * 3 + 1] * kBohr2Angstrom;
				ap->r.z = dary[i * 3 + 2] * kBohr2Angstrom;
			}
			free(dary);
			dary = NULL;
		} else if (strcmp(buf, "MxBond") == 0) {
			if (tokens[1] == NULL || (mxbond = atoi(tokens[1])) <= 0) {
				s_append_asprintf(errbuf, "Line %d: strange number of bonds per atom: %s", lineNumber, tokens[1]);
				retval = 2;
				goto cleanup;
			}
		} else if (strcmp(buf, "IBond") == 0) {
			Int *bonds;
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != natoms * mxbond) {
				s_append_asprintf(errbuf, "Line %d: wrong or inconsistent number of bonds: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&iary, &nary, sizeof(Int), natoms * mxbond, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read bond information", lineNumber);
				retval = 2;
				goto cleanup;
			}
			bonds = (Int *)malloc(sizeof(Int) * (mxbond * 2 + 1));
			for (i = 0; i < natoms; i++) {
				for (j = k = 0; j < mxbond; j++) {
					n = iary[i * mxbond + j] - 1;
					if (n > i) {
						/*  Connect atom i and atom n  */
						bonds[k++] = i;
						bonds[k++] = n;
					}
				}
				if (k > 0) {
					bonds[k] = kInvalidIndex;
					MoleculeAddBonds(mp, k / 2, bonds, NULL, 1);
				}
			}
			free(iary);
			free(bonds);
			iary = NULL;
		} else if (strcmp(buf, "Shell types") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0) {
				s_append_asprintf(errbuf, "Line %d: wrong number of shell types: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&iary, &nary, sizeof(Int), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read shell types", lineNumber);
				retval = 2;
				goto cleanup;
			}
			/*  Allocate ShellInfo table and store shell type information  */
			AssignArray(&bset->shells, &bset->nshells, sizeof(ShellInfo), nary - 1, NULL);
			for (i = n = 0, sp = bset->shells; i < nary; i++, sp++) {
				switch (iary[i]) {
					case 0:  sp->sym = kGTOType_S;  sp->ncomp = 1; break;
					case 1:  sp->sym = kGTOType_P;  sp->ncomp = 3; break;
					case -1: sp->sym = kGTOType_SP; sp->ncomp = 4; break;
					case 2:  sp->sym = kGTOType_D;  sp->ncomp = 6; break;
					case -2: sp->sym = kGTOType_D5; sp->ncomp = 5; break;
					case 3:  sp->sym = kGTOType_F;  sp->ncomp = 10; break;
					case -3: sp->sym = kGTOType_F7; sp->ncomp = 7; break;
					case 4:  sp->sym = kGTOType_G;  sp->ncomp = 15; break;
					case -4: sp->sym = kGTOType_G9; sp->ncomp = 9; break;
					default:
						s_append_asprintf(errbuf, "Line %d: unsupported shell type %d", lineNumber, iary[i]);
						retval = 2;
						goto cleanup;
				}
				sp->m_idx = n;
				n += sp->ncomp;
			}
			bset->ncomps = ncomps = n;
			free(iary);
			iary = NULL;
		} else if (strcmp(buf, "Number of primitives per shell") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != bset->nshells) {
				s_append_asprintf(errbuf, "Line %d: wrong size of the primitive table: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&iary, &nary, sizeof(Int), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read primitive table", lineNumber);
				retval = 2;
				goto cleanup;
			}
			for (i = n = 0, sp = bset->shells; i < bset->nshells; i++, sp++) {
				sp->nprim = iary[i];
				sp->p_idx = n;
				n += sp->nprim;
			}
			nprims = n;
			free(iary);
			iary = NULL;
		} else if (strcmp(buf, "Shell to atom map") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != bset->nshells) {
				s_append_asprintf(errbuf, "Line %d: wrong size of the shell-to-atom map: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&iary, &nary, sizeof(Int), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read shell-to-atom table", lineNumber);
				retval = 2;
				goto cleanup;
			}
			for (i = 0, sp = bset->shells; i < bset->nshells; i++, sp++) {
				sp->a_idx = iary[i] - 1;
			}
			free(iary);
			iary = NULL;
		} else if (strcmp(buf, "Primitive exponents") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != nprims) {
				s_append_asprintf(errbuf, "Line %d: wrong number of primitive exponents: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read primitive exponents", lineNumber);
				retval = 2;
				goto cleanup;
			}
			/*  Allocate PrimInfo table  */
			AssignArray(&bset->priminfos, &bset->npriminfos, sizeof(PrimInfo), nprims - 1, NULL);
			for (i = 0, pp = bset->priminfos; i < nprims; i++, pp++) {
				pp->A = dary[i];
			}
			free(dary);
			dary = NULL;
		} else if (strcmp(buf, "Contraction coefficients") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != bset->npriminfos) {
				s_append_asprintf(errbuf, "Line %d: wrong number of contraction coefficients: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read contraction coefficients", lineNumber);
				retval = 2;
				goto cleanup;
			}
			for (i = 0, pp = bset->priminfos; i < bset->npriminfos; i++, pp++) {
				pp->C = dary[i];
			}
			free(dary);
			dary = NULL;
		} else if (strcmp(buf, "P(S=P) Contraction coefficients") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != bset->npriminfos) {
				s_append_asprintf(errbuf, "Line %d: wrong number of P(S=P) contraction coefficients: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read P(S=P) contraction coefficients", lineNumber);
				retval = 2;
				goto cleanup;
			}
			for (i = 0, pp = bset->priminfos; i < bset->npriminfos; i++, pp++) {
				pp->Csp = dary[i];
			}
			free(dary);
			dary = NULL;
		} else if (strcmp(buf, "Alpha Orbital Energies") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != ncomps) {
				s_append_asprintf(errbuf, "Line %d: wrong or inconsistent number of alpha orbitals: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&bset->moenergies, &bset->nmos, sizeof(Double), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read alpha orbital energies", lineNumber);
				retval = 2;
				goto cleanup;
			}
		} else if (strcmp(buf, "Alpha MO coefficients") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != ncomps * ncomps) {
				s_append_asprintf(errbuf, "Line %d: wrong or inconsistent number of alpha MO coefficients: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&bset->mo, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read MO coefficients", lineNumber);
				retval = 2;
				goto cleanup;
			}
		} else if (strcmp(buf, "Beta Orbital Energies") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != ncomps) {
				s_append_asprintf(errbuf, "Line %d: wrong or inconsistent number of beta orbitals: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read beta orbital energies", lineNumber);
				retval = 2;
				goto cleanup;
			}
			bset->moenergies = (Double *)realloc(bset->moenergies, sizeof(Double) * 2 * ncomps);
			bset->nmos = ncomps * 2;
			bset->mo = (Double *)realloc(bset->mo, sizeof(Double) * 2 * ncomps * ncomps);
			memmove(bset->moenergies + ncomps, dary, sizeof(Double) * ncomps);
			memset(bset->mo + ncomps * ncomps, 0, sizeof(Double) * ncomps * ncomps);
			free(dary);
			dary = NULL;
		} else if (strcmp(buf, "Beta MO coefficients") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != ncomps * ncomps) {
				s_append_asprintf(errbuf, "Line %d: wrong or inconsistent number of beta MO coefficients: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read alpha MO coefficients", lineNumber);
				retval = 2;
				goto cleanup;
			}
			bset->mo = (Double *)realloc(bset->mo, sizeof(Double) * 2 * ncomps * ncomps);  /*  Should be unnecessary, just in case  */
			memmove(bset->mo + ncomps * ncomps, dary, sizeof(Double) * ncomps * ncomps);
			free(dary);
			dary = NULL;
		} else if (strcmp(buf, "Total SCF Density") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != ncomps * (ncomps + 1) / 2) {
				s_append_asprintf(errbuf, "Line %d: wrong or inconsistent number of SCF densities: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&bset->scfdensities, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				s_append_asprintf(errbuf, "Line %d: cannot read SCF densities", lineNumber);
				retval = 2;
				goto cleanup;
			}
		}
	}
	if (mp->natoms == 0) {
		s_append_asprintf(errbuf, "Atom information is missing");
		retval = 2;
		goto cleanup;
	}
	if (bset->shells == NULL || bset->priminfos == NULL) {
		s_append_asprintf(errbuf, "Gaussian primitive information is missing");
		retval = 2;
		goto cleanup;
	}
	if (bset->mo == NULL) {
		s_append_asprintf(errbuf, "MO coefficients were not found");
		retval = 2;
		goto cleanup;
	}
	if (sSetupGaussianCoefficients(bset) != 0) {
		s_append_asprintf(errbuf, "Internal error during setup MO calculation");
		retval = 2;
		goto cleanup;
	}
	mp->nframes = -1;
	retval = 0;
cleanup:
	fclose(fp);
	if (iary != NULL)
		free(iary);
	if (dary != NULL)
		free(dary);
	if (retval != 0) {
		if (mp->bset != NULL) {
			BasisSetRelease(mp->bset);
			mp->bset = NULL;
		}
	}
	return retval;
panic:
	Panic("low memory while reading fchk file %s", fname);
	return -1; /* not reached */	
}

int
MoleculeLoadGamessDatFile(Molecule *mol, const char *fname, char **errbuf)
{
	FILE *fp;
	int newmol = 0;
	char buf[1024];
	int lineNumber, i, j, k, len, natoms = 0;
	int nframes = 0;
	int n1;
	int retval = 0;
	int ival[8];
	double dval[8];
	char sval[16];
	Vector *vbuf = NULL;
	IntGroup *ig;
	int optimizing = 0, status = 0;
	
	*errbuf = NULL;
	if (mol == NULL) {
		mol = MoleculeNew();
	}
	if (mol->natoms == 0)
		newmol = 1;

	fp = fopen(fname, "rb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot open file");
		return 1;
	}
	
	/*  ESP is cleared (not undoable!)  */
	if (mol->elpots != NULL) {
		free(mol->elpots);
		mol->elpots = NULL;
		mol->nelpots = 0;
	}
	
	lineNumber = 0;
	while ((status = sReadLineWithInterrupt(buf, sizeof buf, fp, &lineNumber)) > 0) {
	redo:
		n1 = 0;
		if (strncmp(buf, " $DATA", 6) == 0) {
			/*  Initial geometry  */
			if (!newmol) {
				vbuf = (Vector *)calloc(sizeof(Vector), mol->natoms);
			}
			i = 0;
			ReadLine(buf, sizeof buf, fp, &lineNumber);  /*  Title  */
			ReadLine(buf, sizeof buf, fp, &lineNumber);  /*  Symmetry  */
			while ((status = sReadLineWithInterrupt(buf, sizeof buf, fp, &lineNumber)) > 0) {
				if (strncmp(buf, " $END", 5) == 0)
					break;
				if (sscanf(buf, "%12s %lf %lf %lf %lf", sval, &dval[0], &dval[1], &dval[2], &dval[3]) < 5) {
					s_append_asprintf(errbuf, "Line %d: bad format in $DATA section", lineNumber);
					retval = 2;
					goto exit_loop;
				}
				if (newmol) {
					Atom a;
					memset(&a, 0, sizeof(a));
					strncpy(a.aname, sval, 4);
					a.r.x = dval[1];
					a.r.y = dval[2];
					a.r.z = dval[3];
					a.atomicNumber = (Int)dval[0];
					strncpy(a.element, ElementToString(a.atomicNumber, sval), 3);
					a.type = AtomTypeEncodeToUInt(a.element);
					a.weight = WeightForAtomicNumber(a.atomicNumber);
					MoleculeCreateAnAtom(mol, &a, mol->natoms);
				} else {
					Atom *ap;
					if (i >= mol->natoms) {
						s_append_asprintf(errbuf, "Line %d: too many atoms", lineNumber);
						retval = 3;
						goto exit_loop;
					}
					if ((ap = ATOM_AT_INDEX(mol->atoms, i))->atomicNumber != dval[0]) {
						s_append_asprintf(errbuf, "Line %d: atomic number does not match", lineNumber);
						retval = 4;
						goto exit_loop;
					}
					vbuf[i].x = dval[1];
					vbuf[i].y = dval[2];
					vbuf[i].z = dval[3];
				}
				/*  Skip until a blank line is found  */
				/*  2013.6.11. Line including "PM3" is also recognized as the end of atom  */
				while ((status = sReadLineWithInterrupt(buf, sizeof buf, fp, &lineNumber)) > 0) {
					for (j = 0; buf[j] == ' '; j++);
					if (buf[j] == '\n' || strncmp(buf + j, "PM3", 3) == 0)
						break;
				}
				i++;
			}
			natoms = i;
			if (!newmol) {
				/*  Set atom positions  */
				IntGroup *ig;
				if (natoms < mol->natoms) {
					s_append_asprintf(errbuf, "Line %d: too few atoms", lineNumber);
					retval = 5;
					goto exit_loop;
				}
				ig = IntGroupNewWithPoints(0, natoms, -1);
				MolActionCreateAndPerform(mol, gMolActionSetAtomPositions, ig, natoms, vbuf);
				IntGroupRelease(ig);
			}
			if (vbuf == NULL)
				vbuf = (Vector *)calloc(sizeof(Vector), natoms);
			nframes = MoleculeGetNumberOfFrames(mol);
			if (status < 0)
				break;
			continue;
		} else if (strstr(buf, "DATA FROM NSERCH") != NULL || (strstr(buf, "RESULTS FROM SUCCESSFUL") != NULL && (n1 = 1))) {
			/*  Skip until the separator line is read (three or four lines)  */
			i = 0;
			do {
				if (i++ >= 4) {
					s_append_asprintf(errbuf, "Line %d: the separator line at the top of the coordinates is not found: bad format?", lineNumber);
					retval = 6;
					goto exit_loop;
				}
				ReadLine(buf, sizeof buf, fp, &lineNumber);
			} while (strstr(buf, "----------------------------") == NULL);
			for (i = 0; i < natoms; i++) {
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					s_append_asprintf(errbuf, "Unexpected end of file in reading NSERCH data");
					retval = 6;
					goto exit_loop;
				}
				if (sscanf(buf, "%12s %lf %lf %lf %lf", sval, &dval[0], &dval[1], &dval[2], &dval[3]) < 5) {
					s_append_asprintf(errbuf, "Line %d: bad format in NSERCH coordinate data", lineNumber);
					retval = 6;
					goto exit_loop;
				}
				vbuf[i].x = dval[1];
				vbuf[i].y = dval[2];
				vbuf[i].z = dval[3];
			}
			ig = IntGroupNewWithPoints(nframes, 1, -1);
			MolActionCreateAndPerform(mol, gMolActionInsertFrames, ig, natoms, vbuf, 0, NULL);
			IntGroupRelease(ig);
			nframes++;
			if (n1 == 0)
				optimizing = 1;  /*  Flag to skip reading the VEC group  */
			else
				optimizing = 0;
			continue;
		} else if (strstr(buf, "E(UHF)") != NULL || (strstr(buf, "E(RHF)") != NULL && (n1 = 1)) || (strstr(buf, "E(ROHF)") != NULL && (n1 = 2))) {
			if (mol->bset == NULL) {
				i = MoleculeSetMOInfo(mol, n1, 0, 0);
				if (i != 0) {
					s_append_asprintf(errbuf, "Line %d: cannot allocate basis set internal buffer", lineNumber);
					retval = 8;
					goto exit_loop;
				}
			}
		} else if (strncmp(buf, " $VEC", 5) == 0) {
			Double *coeffs;
			/*  Read the vec group  */
			if (mol->bset == NULL || mol->bset->ncomps == 0)
				continue;  /*  Just ignore  */
			if (optimizing)
				continue;  /*  Ignore VEC group during optimization  */
			coeffs = (Double *)calloc(sizeof(Double), mol->bset->ncomps);
			if (coeffs == NULL) {
				s_append_asprintf(errbuf, "Line %d: low memory during $VEC", lineNumber);
				retval = 9;
				goto exit_loop;
			}
			i = k = 0;
			while ((status = sReadLineWithInterrupt(buf, sizeof buf, fp, &lineNumber)) > 0) {
				len = strlen(buf);
				if (strncmp(buf, " $END", 5) == 0)
					break;
				while ((j = 5 + (k % 5) * 15) <= len && buf[j] != 0 && buf[j] != '\n') {
					strncpy(sval, buf + j, 15);
					sval[15] = 0;
					coeffs[k] = strtod(sval, NULL);
					k++;
					if ((k % 5) == 0)
						break;
				}
				if (k < mol->bset->ncomps)
					continue;
				j = MoleculeSetMOCoefficients(mol, i + 1, -1000000, k, coeffs);
				if (j != 0) {
					s_append_asprintf(errbuf, "Line %d: cannot set coefficients for MO %d", lineNumber, i + 1);
					free(coeffs);
					retval = 10;
					goto exit_loop;
				}
				i++;
				k = 0;
			}
			if (status < 0)
				break;
			continue;
		} else if ((strstr(buf, "ELECTRIC POTENTIAL") != NULL || strstr(buf, "ELECTROSTATIC POTENTIAL") != NULL) && strstr(buf, "ELPOTT") != NULL) {
			i = 0;
			while ((status = sReadLineWithInterrupt(buf, sizeof buf, fp, &lineNumber)) > 0) {
				Elpot *ep;
				if (strstr(buf, "TOTAL NUMBER OF GRID POINTS") != NULL)
					continue;
				if (sscanf(buf, "%d %lf %lf %lf %lf", &ival[0], &dval[0], &dval[1], &dval[2], &dval[3]) < 5)
					break;
				ep = AssignArray(&mol->elpots, &mol->nelpots, sizeof(Elpot), i, NULL);
				ep->pos.x = dval[0];
				ep->pos.y = dval[1];
				ep->pos.z = dval[2];
				ep->esp = dval[3];
				i++;
			}
			if (status > 0)
				goto redo;  /*  This section has no end line, so the last line should be processed again  */
			else break;    /*  End of file encountered or interrupted */
		}  /*  TODO: read MOLPLT info if present  */
	}
	if (status < 0) {
		s_append_asprintf(errbuf, "User interrupt at line %d", lineNumber);
		retval = 11;
	}
exit_loop:
	if (vbuf != NULL)
		free(vbuf);
	if (mol->natoms > 0)
		retval = 0;  /*  Return the partially constructed molecule  */
	if (newmol && mol->nbonds == 0) {
		/*  Guess bonds  */
		Int nbonds, *bonds;
		MoleculeGuessBonds(mol, 0.0, &nbonds, &bonds);
		if (nbonds > 0) {
			MolActionCreateAndPerform(mol, gMolActionAddBonds, nbonds * 2, bonds, NULL);
			free(bonds);
		}
	}
	return 0;
}

int
MoleculeReadCoordinatesFromFile(Molecule *mp, const char *fname, const char *ftype, char **errbuf)
{
	int retval;
	if (ftype == NULL || *ftype == 0) {
		const char *cp;
		cp = strrchr(fname, '.');
		if (cp != NULL)
			ftype = cp + 1;
		else {
			cp = guessMoleculeType(fname);
			if (strcmp(cp, "???") != 0)
				ftype = cp;
		}
	}
	if (strcasecmp(ftype, "pdb") == 0) {
		retval = MoleculeReadCoordinatesFromPdbFile(mp, fname, errbuf);
	}
	if (retval != 0) {
		/*  Try all formats once again  */
		retval = MoleculeReadCoordinatesFromPdbFile(mp, fname, errbuf);
	}
	return retval;
}

int
MoleculeReadCoordinatesFromPdbFile(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	char buf[1024];
	char *p;
	int lineNumber;
	int i, j, new_unit, retval;
	Atom *ap;
	IntGroup *ig;
	Vector *vp = NULL;
	Int ibuf[12];
	Int entries = 0;
	retval = 0;
	*errbuf = NULL;
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot open file");
		return -1;
	}
/*	flockfile(fp); */
	if (mp->natoms == 0)
		new_unit = 1;
	else {
		/*  Allocate buffer for undo-capable modification  */
		vp = (Vector *)calloc(sizeof(Vector), mp->natoms);
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			/*  Retain current position if the atom info is missing in the input file  */
			vp[i] = ap->r;
		}
		new_unit = 0;
	}
	lineNumber = 0;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		if (strncmp(buf, "END", 3) == 0)
			break;
		if (strncmp(buf, "HETATM", 6) == 0 || strncmp(buf, "ATOM", 4) == 0) {
			struct {
				Int serial, intCharge, resSeq;
				Vector r;
				Double occ, temp;
				char segName[5], resName[4], atomName[5], resSeqStr[5], atomType[3], element[3], occStr[6];
			} w;
			memset(&w, 0, sizeof(w));
			ReadFormat(buf, "x6 I5 x1 S4 x1 S3 x1 x1 S4 x1 x3 F8 F8 F8 S6 F6 x6 S4 S2 I2",
				&w.serial, w.atomName, w.resName, w.resSeqStr, &w.r.x, &w.r.y, &w.r.z,
				w.occStr, &w.temp, w.segName, w.element, &w.intCharge);
			if (w.atomName[0] == 0) {
				continue;  /*  Atom name is empty  */
			}
			/*  A workaround for residue number >= 10000 (XPLOR style)  */
			if (w.resSeqStr[0] >= 'A' && w.resSeqStr[0] <= 'Z') {
				w.resSeq = (w.resSeqStr[0] - 'A' + 10) * 1000 + atoi(w.resSeqStr + 1);
			} else {
				w.resSeq = atoi(w.resSeqStr);
			}
			if (w.element[0] == 0) {
				/*  $element = ($name =~ /([A-Za-z]{1,2})/); # in Perl  */
				for (p = w.atomName; *p != 0; p++) {
					if (isalpha(*p) && *p != '_') {
						w.element[0] = toupper(*p);
						if (isalpha(p[1]) && p[1] != '_') {
							w.element[1] = toupper(p[1]);
							w.element[2] = 0;
						} else {
							w.element[1] = 0;
						}
						break;
					}
				}
			}
			if (w.occStr[0] == 0)
				w.occ = 1.0;
			else
				w.occ = atof(w.occStr);
			if (w.serial <= 0) {
				s_append_asprintf(errbuf, "line %d: non-positive atom number %d", lineNumber, w.serial);
				retval = 1;
				goto abort;
			}
			w.serial--;  /*  The internal atom number is 0-based  */
			if (w.serial >= mp->natoms) {
				if (new_unit) {
					/*  Create a new atom entry  */
					ap = AssignArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, w.serial, NULL);
				} else {
					s_append_asprintf(errbuf, "line %d: the atom number %d does not exist in the structure file", lineNumber, w.serial+1);
					retval = 1;
					goto abort;
				}
			}
			if (new_unit) {
				ap = ATOM_AT_INDEX(mp->atoms, w.serial);
				ap->r = w.r;
				ap->occupancy = w.occ;
				ap->tempFactor = w.temp;
				if (w.segName[0] == 0)
					strncpy(w.segName, "MAIN", 4);
				strncpy(ap->segName, w.segName, 4);
				ap->resSeq = w.resSeq;
				strncpy(ap->resName, w.resName, 4);
				strncpy(ap->aname, w.atomName, 4);
				strncpy(ap->element, w.element, 2);
				ap->element[2] = 0;
				ap->atomicNumber = ElementToInt(ap->element);
				ap->type = AtomTypeEncodeToUInt(ap->element);
				ap->weight = WeightForAtomicNumber(ap->atomicNumber);
				ap->intCharge = w.intCharge;
				if (ap->resSeq > 0) {
					if (ap->resSeq < mp->nresidues) {
						/*  Update the resName according to residues[]  */
						strncpy(ap->resName, mp->residues[ap->resSeq], 4);
					} else {
						/*  Register the resName to residues[]  */
						AssignArray(&mp->residues, &mp->nresidues, 4, ap->resSeq, w.resName);
					}
				} else {
					ap->resSeq = 0;
					strcpy(ap->resName, "XXX");
					if (mp->nresidues == 0)
						AssignArray(&mp->residues, &mp->nresidues, 4, 0, ap->resName);
				}
				i = ElementToInt(ap->element);
				if (i >= 0)
					ap->weight = gElementParameters[i].weight;
			} else {
				/*  Not a new unit: only the atom position is updated  */
				vp[w.serial] = w.r;
			}
			entries++;
		} else if (strncmp(buf, "CONECT", 6) == 0 && new_unit) {
			i = ReadFormat(buf, "x6 I5I5I5I5I5I5I5I5I5I5I5I5",
				ibuf, ibuf + 1, ibuf + 2, ibuf + 3,
				ibuf + 4, ibuf + 5, ibuf + 6, ibuf + 7,
				ibuf + 8, ibuf + 9, ibuf + 10, ibuf + 11);
			if (i >= 2) {
				Int bbuf[25];
				int bi;
				for (j = 0; j < i; j++) {
					if (ibuf[j] < 0 || ibuf[j] > mp->natoms) {
						s_append_asprintf(errbuf, "line %d: The CONECT record contains non-existent atom %d", lineNumber, ibuf[j]);
						retval = 1;
						goto abort;
					} else if (ibuf[j] == 0)
						break;
				}
				i = j;
				if (i < 2)
					continue;
				for (j = 1, bi = 0; j < i; j++) {
					if (ibuf[0] < ibuf[j]) {
						if (MoleculeLookupBond(mp, ibuf[0], ibuf[j]) >= 0) {
							s_append_asprintf(errbuf, "line %d: warning: duplicate bond %d-%d\n", lineNumber, ibuf[0], ibuf[j]);
						} else {
							bbuf[bi * 2] = ibuf[0] - 1;
							bbuf[bi * 2 + 1] = ibuf[j] - 1;
							bi++;
						}
					}
				}
				if (bi == 0)
					continue;
				bbuf[bi * 2] = -1;
				retval = MoleculeAddBonds(mp, bi, bbuf, NULL, 1);
				if (retval < 0) {
					s_append_asprintf(errbuf, "line %d: bad bond specification", lineNumber);
					retval = 1;
					goto abort;
				}
			}
		}
	}
/*	funlockfile(fp); */
	fclose(fp);
	if (new_unit) {
		/*  Renumber atoms if some atom number is unoccupied  */
		int *old2new, oldidx, newidx;
		old2new = (int *)calloc(sizeof(int), mp->natoms);
		if (old2new == NULL) {
			s_append_asprintf(errbuf, "Out of memory");
			retval = 1;
			goto abort;
		}
		for (oldidx = newidx = 0; oldidx < mp->natoms; oldidx++) {
			ap = ATOM_AT_INDEX(mp->atoms, oldidx);
			if (ap->aname[0] != 0) {
				old2new[oldidx] = newidx;
				if (oldidx > newidx)
					memmove(ATOM_AT_INDEX(mp->atoms, newidx), ap, gSizeOfAtomRecord);
				newidx++;
			}
		}
		mp->natoms = newidx;
		if (oldidx > newidx) {
			/*  Renumber the connects and bonds  */
			Int *cp;
			for (i = 0; i < mp->natoms; i++) {
				ap = ATOM_AT_INDEX(mp->atoms, i);
				cp = AtomConnectData(&ap->connect);
				for (j = 0; j < ap->connect.count; j++) {
					cp[j] = old2new[cp[j]];
				}
			}
			for (i = 0; i < mp->nbonds * 2; i++) {
				mp->bonds[i] = old2new[mp->bonds[i]];
			}
		}
		retval = MoleculeRebuildTablesFromConnects(mp);
		if (retval != 0) {
			/*  This error may not happen  */
			s_append_asprintf(errbuf, "Cannot build angle/dihedral/improper tables");
			retval = 1;
			goto abort;
		}
		/*  Undo action: delete all atoms  */
		{
			MolAction *act;
			ig = IntGroupNewWithPoints(0, mp->natoms, -1);
			act = MolActionNew(gMolActionUnmergeMolecule, ig);
			act->frame = mp->cframe;
			MolActionCallback_registerUndo(mp, act);
			MolActionRelease(act);
			IntGroupRelease(ig);
		}
	} else {
		/*  Set the new atom positions  */
		ig = IntGroupNewWithPoints(0, mp->natoms, -1);
		MolActionCreateAndPerform(mp, gMolActionSetAtomPositions, ig, mp->natoms, vp);
		IntGroupRelease(ig);
		free(vp);
		vp = NULL;
	}
	mp->nframes = -1;  /*  Should be recalculated later  */
	if (entries == 0)
		return 1;  /*  No atoms  */
	return 0;
	abort:
	if (fp != NULL) {
	/*	funlockfile(fp); */
		fclose(fp);
	}
	if (vp != NULL)
		free(vp);
	if (entries == 0)
		return 1;  /*  Maybe different format?  */
	return retval;
}

int
MoleculeReadCoordinatesFromDcdFile(Molecule *mp, const char *fname, char **errbuf)
{
	DcdRecord dcd;
	SFloat32 *xp, *yp, *zp;
	Vector *vp, *cp;
	IntGroup *ig;
	int n, errcount = 0;
	*errbuf = NULL;
	if (mp == NULL || mp->natoms == 0) {
		s_append_asprintf(errbuf, "Molecule is empty");
		return 1;
	}
	n = DcdOpen(fname, &dcd);
	if (n != 0) {
		switch (n) {
			case -2: s_append_asprintf(errbuf, "Cannot open file"); break;
			case 1:  s_append_asprintf(errbuf, "Premature EOF encountered"); break;
			case 2:  s_append_asprintf(errbuf, "Bad block length of the first section"); break;
			case 3:  s_append_asprintf(errbuf, "\"CORD\" signature is missing"); break;
			case 4:  s_append_asprintf(errbuf, "Bad termination of the first section"); break;
			case 5:  s_append_asprintf(errbuf, "The title section is not correct"); break;
			case 6:  s_append_asprintf(errbuf, "The atom number section is not correct"); break;
			default: s_append_asprintf(errbuf, "Read error in dcd file"); break;
		}
		errcount++;
	} else {
		if (dcd.natoms == 0) {
			s_append_asprintf(errbuf, "No atoms were found in the dcd file");
			errcount++;
		} else if (dcd.nframes == 0) {
			s_append_asprintf(errbuf, "No frames were found in the dcd file");
			errcount++;
		}
	}
	if (errcount > 0) {
		if (n == 0)
			DcdClose(&dcd);
		return 1;
	}

	vp = (Vector *)calloc(sizeof(Vector), mp->natoms * dcd.nframes);
	if (dcd.nextra)
		cp = (Vector *)calloc(sizeof(Vector), dcd.nframes * 4);
	else cp = NULL;
	xp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	yp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	zp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	ig = IntGroupNewWithPoints(MoleculeGetNumberOfFrames(mp), dcd.nframes, -1);
	if (vp == NULL || xp == NULL || yp == NULL || zp == NULL || ig == NULL) {
		s_append_asprintf(errbuf, "Cannot allocate memory");
		if (vp) free(vp);
		if (cp) free(cp);
		if (xp) free(xp);
		if (yp) free(yp);
		if (zp) free(zp);
		if (ig) IntGroupRelease(ig);
		return 1;
	}
	for (n = 0; n < dcd.nframes; n++) {
		int i;
		Vector *vpp;
		SFloat32 dcdcell[6];
		if (DcdReadFrame(&dcd, n, xp, yp, zp, dcdcell)) {
			s_append_asprintf(errbuf, "Read error in dcd file");
			goto exit;
		}
		for (i = 0, vpp = &vp[n * mp->natoms]; i < dcd.natoms && i < mp->natoms; i++, vpp++) {
			vpp->x = xp[i];
			vpp->y = yp[i];
			vpp->z = zp[i];
		}
		if (cp != NULL) {
			Double sing;
			vpp = &cp[n * 4];
			/*  dcdcell = {a, gamma, b, beta, alpha, c} */
			/*  angles are described either in cosines (Charmm and NAMD > 2.5) or degrees (NAMD 2.5)  */
			if (dcdcell[1] < -1.0 || dcdcell[1] > 1.0 || dcdcell[3] < -1.0 || dcdcell[3] > 1.0 || dcdcell[4] < -1.0 || dcdcell[4] > 1.0) {
				dcdcell[4] = cos(dcdcell[4] * kDeg2Rad);  /*  cos(alpha)  */
				dcdcell[3] = cos(dcdcell[3] * kDeg2Rad);  /*  cos(beta)  */
				dcdcell[1] = cos(dcdcell[1] * kDeg2Rad);  /*  cos(gamma)  */
			}
			/*  a axis lies along the cartesian x axis  */
			sing = sqrt(1 - dcdcell[1] * dcdcell[1]);
			vpp[0].x = dcdcell[0];
			vpp[0].y = 0;
			vpp[0].z = 0;
			vpp[1].x = dcdcell[2] * dcdcell[1];
			vpp[1].y = dcdcell[2] * sing;
			vpp[1].z = 0;
			vpp[2].x = dcdcell[5] * dcdcell[3];
			vpp[2].y = dcdcell[5] * (dcdcell[4] - dcdcell[3] * dcdcell[1]) / sing;
			vpp[2].z = sqrt(dcdcell[5] * dcdcell[5] - vpp[2].x * vpp[2].x - vpp[2].y * vpp[2].y);
			vpp[3].x = vpp[3].y = vpp[3].z = 0.0;
			if (mp->cell == NULL) {
				/*  Create periodicity if not present  */
				MolActionCreateAndPerform(mp, gMolActionSetBox, &vpp[0], &vpp[1], &vpp[2], &vpp[3], 7, 0);
			}
		}
	}
	if (MolActionCreateAndPerform(mp, gMolActionInsertFrames, ig, mp->natoms * dcd.nframes, vp, (cp == NULL ? 0 : dcd.nframes * 4), cp) != 0)
		s_append_asprintf(errbuf, "Cannot insert frames");
	mp->startStep = dcd.nstart;
	mp->stepsPerFrame = dcd.ninterval;
	mp->psPerStep = dcd.delta;
exit:
	DcdClose(&dcd);
	if (cp != NULL)
		free(cp);
	free(vp);
	free(xp);
	free(yp);
	free(zp);
	IntGroupRelease(ig);
	if (errcount == 0)
		return 0;
	else return 1;
}

int
MoleculeReadExtendedInfo(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	char buf[1024];
	int lineNumber;
	int i, retval;
	Vector v[3], vv;
	double d[3];
	int n, flag;
	char flags[3];
	*errbuf = NULL;
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot open file");
		return -1;
	}
	errbuf[0] = 0;
	lineNumber = 0;
	retval = 0;
	flags[0] = flags[1] = flags[2] = 0;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		if (strncmp(buf, "Bounding box:", 13) == 0) {
			for (i = 0; i < 3; i++) {
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					s_append_asprintf(errbuf, "line %d: missing %d component of the bounding box", lineNumber, i + 1);
					retval = 1;
					goto abort;
				}
				n = sscanf(buf, "%lf %lf %lf %d", &d[0], &d[1], &d[2], &flag);
				if (n < 3) {
					vv.x = vv.y = vv.z = 0.0;
					switch (i) {
						case 0: vv.x = d[0]; break;
						case 1: vv.y = d[0]; break;
						case 2: vv.z = d[0]; break;
					}
					if (n == 1 || (n == 2 && d[1] != 0.0))
						flags[i] = 1;
				} else {
					vv.x = d[0];
					vv.y = d[1];
					vv.z = d[2];
					if (n == 4)
						flags[i] = (flag != 0);
					else
						flags[i] = (VecLength2(vv) != 0);
				}
				v[i] = vv;
			}
			if (mp->cell != NULL)
				vv = mp->cell->origin;
			else
				vv.x = vv.y = vv.z = 0.0;
			MoleculeSetPeriodicBox(mp, &v[0], &v[1], &v[2], &vv, flags, 0);
		} else if (strncmp(buf, "Bounding box origin:", 20) == 0) {
			if (mp->cell != NULL) {
				v[0] = mp->cell->axes[0];
				v[1] = mp->cell->axes[1];
				v[2] = mp->cell->axes[2];
				memmove(flags, mp->cell->flags, 3);
			} else {
				v[0].x = 1.0; v[0].y = v[0].z = 0.0;
				v[1].y = 1.0; v[1].x = v[1].z = 0.0;
				v[2].z = 1.0; v[2].x = v[2].y = 0.0;
				flags[0] = flags[1] = flags[2] = 1.0;
			}
			if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0 || (n = sscanf(buf, "%lf %lf %lf", &d[0], &d[1], &d[2]) < 3)) {
				s_append_asprintf(errbuf, "line %d: wrong format for the bounding box origin", lineNumber);
				retval = 1;
				goto abort;
			}
			vv.x = d[0];
			vv.y = d[1];
			vv.z = d[2];
			MoleculeSetPeriodicBox(mp, &v[0], &v[1], &v[2], &vv, flags, 0);
		}
	}
	fclose(fp);
	return 0;
abort:
	if (fp != NULL)
		fclose(fp);
	return retval;
}
			
int
MoleculeWriteToFile(Molecule *mp, const char *fname, const char *ftype, char **errbuf)
{
	int retval;
	*errbuf = NULL;
	if (ftype == NULL || *ftype == 0) {
		const char *cp;
		cp = strrchr(fname, '.');
		if (cp != NULL)
			ftype = cp + 1;
		else {
			cp = guessMoleculeType(fname);
			if (strcmp(cp, "???") != 0)
				ftype = cp;
		}
	}
	if (strcasecmp(ftype, "psf") == 0) {
		retval = MoleculeWriteToPsfFile(mp, fname, errbuf);
	} else if (strcasecmp(ftype, "pdb") == 0) {
		retval = MoleculeWriteToPdbFile(mp, fname, errbuf);
	} else if (strcasecmp(ftype, "tep") == 0) {
		retval = MoleculeWriteToTepFile(mp, fname, errbuf);
	} else {
		s_append_asprintf(errbuf, "The file format should be specified");
		retval = 1;
	}
	if (retval == 0)
		MoleculeSetPath(mp, fname);
	return retval;
}

int
MoleculeWriteToMbsfFile(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	Int i, j, k, n1, n2, n3, n_aniso, nframes, nanchors, n_uff;
	Atom *ap;
	char *p;
	char bufs[6][8];

	*errbuf = NULL;
	fp = fopen(fname, "wb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot write to file %s", fname);
		return 1;
	}
	errbuf[0] = 0;

	nframes = MoleculeFlushFrames(mp);

	fprintf(fp, "!:atoms\n");
	fprintf(fp, "! idx seg_name res_seq res_name name type charge weight element atomic_number occupancy temp_factor int_charge\n");
	n1 = n2 = n3 = n_aniso = nanchors = n_uff = 0;
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		strncpy(bufs[0], ap->segName, 4);
		bufs[0][4] = 0;
		strncpy(bufs[1], ap->resName, 4);
		bufs[1][4] = 0;
		strncpy(bufs[2], ap->aname, 4);
		bufs[2][4] = 0;
		AtomTypeDecodeToString(ap->type, bufs[3]);
		bufs[3][6] = 0;
		strncpy(bufs[4], ap->element, 4);
		bufs[4][2] = 0;
		for (j = 0; j < 5; j++) {
			if (bufs[j][0] == 0) {
				bufs[j][0] = '_';
				bufs[j][1] = 0;
			}
			for (k = 0; k < 6; k++) {
				if (bufs[j][k] == 0)
					break;
				if (bufs[j][k] > 0 && bufs[j][k] < ' ')
					bufs[j][k] = '_';
			}
		}
		if (SYMOP_ALIVE(ap->symop))
			n1++;
		if (ap->fix_force != 0)
			n2++;
		if (ap->mm_exclude || ap->periodic_exclude)
			n3++;
		if (ap->aniso != NULL)
			n_aniso++;
		if (ap->anchor != NULL)
			nanchors++;
		if (ap->uff_type[0] != 0)
			n_uff++;
		fprintf(fp, "%d %s %d %s %s %s %.5f %.5f %s %d %f %f %d\n", i, bufs[0], ap->resSeq, bufs[1], bufs[2], bufs[3], ap->charge, ap->weight, bufs[4], ap->atomicNumber, ap->occupancy, ap->tempFactor, ap->intCharge);
	}
	fprintf(fp, "\n");
	
	if (n_uff > 0) {
		fprintf(fp, "!:uff_type\n");
		fprintf(fp, "! idx uff_type\n");
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			fprintf(fp, "%d %.5s\n", i, ap->uff_type);
		}
		fprintf(fp, "\n");
	}
	
	if (n1 > 0) {
		fprintf(fp, "!:atoms_symop\n");
		fprintf(fp, "! idx symop symbase\n");
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			int n;
			n = ap->symop.sym * 1000000 + ap->symop.dx * 10000 + ap->symop.dy * 100 + ap->symop.dz;
			fprintf(fp, "%d %d %d\n", i, n, ap->symbase);
		}
		fprintf(fp, "\n");
	}
	
	if (n2 > 0) {
		fprintf(fp, "!:atoms_fix\n");
		fprintf(fp, "! idx fix_force fix_pos\n");
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			fprintf(fp, "%d %f %f %f %f\n", i, ap->fix_force, ap->fix_pos.x, ap->fix_pos.y, ap->fix_pos.z);
		}
		fprintf(fp, "\n");
	}
	
	if (n3 > 0) {
		fprintf(fp, "!:mm_exclude\n");
		fprintf(fp, "! idx mm_exclude periodic_exclude\n");
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			fprintf(fp, "%d %d %d\n", i, ap->mm_exclude, ap->periodic_exclude);
		}
		fprintf(fp, "\n");
	}
	
	if (nanchors > 0) {
		fprintf(fp, "!:pi_anchor\n");
		fprintf(fp, "! idx count; n1 weight1; n2 weight2; ...; nN weightN\n");
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			Int *ip;
			if (ap->anchor == NULL)
				continue;
			k = ap->anchor->connect.count;
			ip = AtomConnectData(&ap->anchor->connect);
			fprintf(fp, "%d %d\n", i, k);
			for (j = 0; j < k; j++) {
				fprintf(fp, "%d %f\n", ip[j], ap->anchor->coeffs[j]);
			}
		}
		fprintf(fp, "\n");
	}
				
	n1 = nframes;
	if (n1 > 0)
		n2 = mp->cframe;
	else
		n2 = 0;
	for (i = 0; (i == n2 || i < n1); i++) {
		fprintf(fp, "!:positions ; frame %d\n", i);
		fprintf(fp, "! idx x y z [sx sy sz]\n");
		for (j = 0, ap = mp->atoms; j < mp->natoms; j++, ap = ATOM_NEXT(ap)) {
			Vector *vp;
			Byte sig_flag = 0;
			if (i != n2 && i < ap->nframes)
				vp = ap->frames + i;
			else {
				vp = &(ap->r);
				if (ap->sigma.x != 0.0 || ap->sigma.y != 0.0 || ap->sigma.z != 0.0)
					sig_flag = 1;
			}
			fprintf(fp, "%d %.8f %.8f %.8f", j, vp->x, vp->y, vp->z);
			if (sig_flag) {
				fprintf(fp, " %.8f %.8f %.8f", ap->sigma.x, ap->sigma.y, ap->sigma.z);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
	}
	
	if (mp->nbonds > 0) {
		fprintf(fp, "!:bonds\n");
		fprintf(fp, "! from1 to1 from2 to2 from3 to3 from4 to4\n");
		for (i = 0; i < mp->nbonds; i++) {
			fprintf(fp, "%d %d%c", mp->bonds[i * 2], mp->bonds[i * 2 + 1], (i % 4 == 3 || i == mp->nbonds - 1 ? '\n' : ' '));
		}
		fprintf(fp, "\n");
	}

	if (mp->nbondOrders > 0) {
		fprintf(fp, "!:bond_orders\n");
		fprintf(fp, "! order1 order2 order3 order4\n");
		for (i = 0; i < mp->nbondOrders; i++) {
			fprintf(fp, "%.6f%c", mp->bondOrders[i], (i % 4 == 3 || i == mp->nbondOrders - 1 ? '\n' : ' '));
		}
		fprintf(fp, "\n");
	}
	
	if (mp->nangles > 0) {
		fprintf(fp, "!:angles\n");
		fprintf(fp, "! a1 b1 c1 a2 b2 c2 a3 b3 c3\n");
		for (i = 0; i < mp->nangles; i++) {
			fprintf(fp, "%d %d %d%c", mp->angles[i * 3], mp->angles[i * 3 + 1], mp->angles[i * 3 + 2], (i % 3 == 2 || i == mp->nangles - 1 ? '\n' : ' '));
		}
		fprintf(fp, "\n");
	}
	
	if (mp->ndihedrals > 0) {
		fprintf(fp, "!:dihedrals\n");
		fprintf(fp, "! a1 b1 c1 d1 a2 b2 c2 d2\n");
		for (i = 0; i < mp->ndihedrals; i++) {
			fprintf(fp, "%d %d %d %d%c", mp->dihedrals[i * 4], mp->dihedrals[i * 4 + 1], mp->dihedrals[i * 4 + 2], mp->dihedrals[i * 4 + 3], (i % 2 == 1 || i == mp->ndihedrals - 1 ? '\n' : ' '));
		}
		fprintf(fp, "\n");
	}
	
	if (mp->nimpropers > 0) {
		fprintf(fp, "!:impropers\n");
		fprintf(fp, "! a1 b1 c1 d1 a2 b2 c2 d2\n");
		for (i = 0; i < mp->nimpropers; i++) {
			fprintf(fp, "%d %d %d %d%c", mp->impropers[i * 4], mp->impropers[i * 4 + 1], mp->impropers[i * 4 + 2], mp->impropers[i * 4 + 3], (i % 2 == 1 || i == mp->nimpropers - 1 ? '\n' : ' '));
		}
		fprintf(fp, "\n");
	}
	
	if (mp->cell != NULL) {
		fprintf(fp, "!:xtalcell\n");
		fprintf(fp, "! a b c alpha beta gamma\n");
		fprintf(fp, "! This information is redundant and overridden by the following periodic_box info\n");
		fprintf(fp, "%f %f %f %f %f %f\n", mp->cell->cell[0], mp->cell->cell[1], mp->cell->cell[2], mp->cell->cell[3], mp->cell->cell[4], mp->cell->cell[5]);
		fprintf(fp, "\n");

		fprintf(fp, "!:periodic_box\n");
		fprintf(fp, "! ax ay az; bx by bz; cx cy cz; ox oy oz; fa fb fc [sigma; sa sb sc s_alpha s_beta s_gamma]\n");
		for (i = 0; i < 3; i++)
			fprintf(fp, "%15.8f %15.8f %15.8f\n", mp->cell->axes[i].x, mp->cell->axes[i].y, mp->cell->axes[i].z);
		fprintf(fp, "%15.8f %15.8f %15.8f\n", mp->cell->origin.x, mp->cell->origin.y, mp->cell->origin.z);
		fprintf(fp, "%d %d %d%s\n", mp->cell->flags[0], mp->cell->flags[1], mp->cell->flags[2], (mp->cell->has_sigma ? " 1" : ""));
		if (mp->cell->has_sigma) {
			fprintf(fp, "%f %f %f %f %f %f\n", mp->cell->cellsigma[0], mp->cell->cellsigma[1], mp->cell->cellsigma[2], mp->cell->cellsigma[3], mp->cell->cellsigma[4], mp->cell->cellsigma[5]);
		}
		fprintf(fp, "\n");
	}
	
	if (mp->nframe_cells > 0) {
		fprintf(fp, "!:frame_periodic_boxes\n");
		fprintf(fp, "! ax ay az; bx by bz; cx cy cz; ox oy oz\n");
		for (i = 0; i < mp->nframe_cells * 4; i++) {
			fprintf(fp, "%15.8f %15.8f %15.8f\n", mp->frame_cells[i].x, mp->frame_cells[i].y, mp->frame_cells[i].z);
		}
		fprintf(fp, "\n");
	}
	
	if (mp->nsyms > 0) {
		fprintf(fp, "!:symmetry_operations\n");
		fprintf(fp, "! a11 a12 a13; a21 a22 a23; a31 a32 a33; t1 t2 t3\n");
		for (i = 0; i < mp->nsyms; i++) {
			Transform *tp = mp->syms + i;
			const unsigned char s_index_order[12] = {0, 3, 6, 1, 4, 7, 2, 5, 8, 9, 10, 11};
			for (j = 0; j < 12; j++)
				fprintf(fp, "%11.6f%c", (*tp)[s_index_order[j]], (j % 3 == 2 ? '\n' : ' '));
		}
		fprintf(fp, "\n");
	}
	
	if (n_aniso > 0) {
		fprintf(fp, "!:anisotropic_thermal_parameters\n");
		fprintf(fp, "! b11 b22 b33 b12 b13 b23 [sigma; sb11 sb22 sb33 sb12 sb13 sb23]\n");
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ap->aniso != NULL) {
				Double *bp = ap->aniso->bij;
				fprintf(fp, "%14.8g %14.8g %14.8g %14.8g %14.8g %14.8g%s\n", bp[0], bp[1], bp[2], bp[3], bp[4], bp[5], (ap->aniso->has_bsig ? " 1" : ""));
				if (ap->aniso->has_bsig) {
					bp = ap->aniso->bsig;
					fprintf(fp, "%14.8g %14.8g %14.8g %14.8g %14.8g %14.8g\n", bp[0], bp[1], bp[2], bp[3], bp[4], bp[5]);
				}
			} else {
				fprintf(fp, "0 0 0 0 0 0\n");
			}
		}
		fprintf(fp, "\n");		
	}
	
	if (mp->arena != NULL) {
		MDArena *arena = mp->arena;
		fprintf(fp, "!:md_parameters\n");
		fprintf(fp, "log_file %s\n", arena->log_result_name);
		fprintf(fp, "coord_file %s\n", arena->coord_result_name);
		fprintf(fp, "vel_file %s\n", arena->vel_result_name);
		fprintf(fp, "force_file %s\n", arena->force_result_name);
		fprintf(fp, "debug_file %s\n", arena->debug_result_name);
		fprintf(fp, "debug_output_level %d\n", arena->debug_output_level);
		fprintf(fp, "step %d\n", arena->step);
		fprintf(fp, "coord_output_freq %d\n", arena->coord_output_freq);
		fprintf(fp, "energy_output_freq %d\n", arena->energy_output_freq);
		fprintf(fp, "coord_frame %d\n", arena->coord_result_frame);
		fprintf(fp, "timestep %g\n", arena->timestep);
		fprintf(fp, "cutoff %g\n", arena->cutoff);
		fprintf(fp, "electro_cutoff %g\n", arena->electro_cutoff);
		fprintf(fp, "pairlist_distance %g\n", arena->pairlist_distance);
		fprintf(fp, "switch_distance %g\n", arena->switch_distance);
		fprintf(fp, "temperature %g\n", arena->temperature);
		fprintf(fp, "andersen_freq %d\n", arena->andersen_thermo_freq);
		fprintf(fp, "andersen_coupling %g\n", arena->andersen_thermo_coupling);
		fprintf(fp, "random_seed %d\n", arena->random_seed);
		fprintf(fp, "dielectric %g\n", arena->dielectric);
		fprintf(fp, "gradient_convergence %g\n", arena->gradient_convergence);
		fprintf(fp, "coordinate_convergence %g\n", arena->coordinate_convergence);
		fprintf(fp, "use_xplor_shift %d\n", arena->use_xplor_shift);
		fprintf(fp, "scale14_vdw %g\n", arena->scale14_vdw);
		fprintf(fp, "scale14_elect %g\n", arena->scale14_elect);
		fprintf(fp, "relocate_center %d\n", arena->relocate_center);
		fprintf(fp, "surface_probe_radius %g\n", arena->probe_radius);
		fprintf(fp, "surface_tension %g\n", arena->surface_tension);
		fprintf(fp, "surface_potential_freq %d\n", arena->surface_potential_freq);
		fprintf(fp, "use_graphite %d\n", arena->use_graphite);
		fprintf(fp, "alchemical_lambda %g\n", arena->alchem_lambda);
		fprintf(fp, "alchemical_delta_lambda %g\n", arena->alchem_dlambda);
		if (arena->nalchem_flags > 0) {
			fprintf(fp, "alchem_flags %d", arena->nalchem_flags);
			for (i = 0; i < arena->nalchem_flags; i++) {
				if (i % 60 == 0)
					fputc('\n', fp);
				else if (i % 10 == 0)
					fputc(' ', fp);
				fputc('0' + arena->alchem_flags[i], fp);
			}
			fputc('\n', fp);
		}
		if (arena->pressure != NULL) {
			Double *dp;
			fprintf(fp, "pressure_freq %d\n", arena->pressure->freq);
			fprintf(fp, "pressure_coupling %g\n", arena->pressure->coupling);
			dp = arena->pressure->apply;
			fprintf(fp, "pressure %g %g %g %g %g %g %g %g %g\n", dp[0], dp[1], dp[2], dp[3], dp[4], dp[5], dp[6], dp[7], dp[8]);
			dp = arena->pressure->cell_flexibility;
			fprintf(fp, "pressure_cell_flexibility %g %g %g %g %g %g %g %g\n", dp[0], dp[1], dp[2], dp[3], dp[4], dp[5], dp[6], dp[7]);
			fprintf(fp, "pressure_fluctuate_cell_origin %g\n", arena->pressure->fluctuate_cell_origin);
			fprintf(fp, "pressure_fluctuate_cell_orientation %g\n", arena->pressure->fluctuate_cell_orientation);
		}
		fprintf(fp, "\n");

		if (mp->par != NULL) {
			Parameter *par = mp->par;
			fprintf(fp, "!:parameters\n");
			ParameterAppendToFile(par, fp);
			fprintf(fp, "\n");
		}
		
		fprintf(fp, "!:velocity\n");
		fprintf(fp, "! idx vx vy vz\n");
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			fprintf(fp, "%d %.8f %.8f %.8f\n", i, ap->v.x, ap->v.y, ap->v.z);
		}
		fprintf(fp, "\n");

		fprintf(fp, "!:force\n");
		fprintf(fp, "! idx fx fy fz\n");
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			fprintf(fp, "%d %.8f %.8f %.8f\n", i, ap->f.x, ap->f.y, ap->f.z);
		}
		fprintf(fp, "\n");
	}
	
	if (mp->mview != NULL) {
		double f[4];
		if (mp->mview->track != NULL) {
			fprintf(fp, "!:trackball\n");
			fprintf(fp, "! scale; trx try trz; theta_deg x y z\n");
			f[0] = TrackballGetScale(mp->mview->track);
			fprintf(fp, "%f\n", f[0]);
			TrackballGetTranslate(mp->mview->track, f);
			fprintf(fp, "%f %f %f\n", f[0], f[1], f[2]);
			TrackballGetRotate(mp->mview->track, f);
			fprintf(fp, "%f %f %f %f\n", f[0], f[1], f[2], f[3]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "!:view\n");
		fprintf(fp, "show_unit_cell %d\n", mp->mview->showUnitCell);
		fprintf(fp, "show_periodic_box %d\n", mp->mview->showPeriodicBox);
		fprintf(fp, "show_expanded_atoms %d\n", mp->mview->showExpandedAtoms);
		fprintf(fp, "show_ellipsoids %d\n", mp->mview->showEllipsoids);
		fprintf(fp, "show_hydrogens %d\n", mp->mview->showHydrogens);
		fprintf(fp, "show_dummy_atoms %d\n", mp->mview->showDummyAtoms);
		fprintf(fp, "show_rotation_center %d\n", mp->mview->showRotationCenter);
		fprintf(fp, "show_graphite_flag %d\n", mp->mview->showGraphiteFlag);
		fprintf(fp, "show_graphite %d\n", mp->mview->showGraphite);
		fprintf(fp, "show_periodic_image_flag %d\n", mp->mview->showPeriodicImageFlag);
		fprintf(fp, "show_periodic_image %d %d %d %d %d %d\n",
				mp->mview->showPeriodicImage[0], mp->mview->showPeriodicImage[1],
				mp->mview->showPeriodicImage[2], mp->mview->showPeriodicImage[3],
				mp->mview->showPeriodicImage[4], mp->mview->showPeriodicImage[5]);
		if (mp->mview->atomRadius != 0.2)
			fprintf(fp, "atom_radius %f\n", mp->mview->atomRadius);
		if (mp->mview->bondRadius != 0.1)
			fprintf(fp, "bond_radius %f\n", mp->mview->bondRadius);
		if (mp->mview->atomResolution != 12)
			fprintf(fp, "atom_resolution %d\n", mp->mview->atomResolution);
		if (mp->mview->bondResolution != 8)
			fprintf(fp, "bond_resolution %d\n", mp->mview->bondResolution);
		fprintf(fp, "\n");
	}

	if (mp->nmolprops > 0) {
		MolProp *prp;
		for (i = 0, prp = mp->molprops; i < mp->nmolprops; i++, prp++) {
			/*  Encode the property name if necessary  */
			char enc[1024];
			n1 = n2 = 0;
			for (p = prp->propname; *p != 0 && n1 < 900; p++) {
				if (*p > ' ' && *p != '%' && *p < 0x7f) {
					enc[n1++] = *p;
					n2 = n1;
				} else {
					sprintf(enc + n1, "%%%02x", *p);
					n1 += 3;
				}
			}
			if (*p == 0)
				enc[n1] = 0;
			else {
				enc[n2] = 0; /* Truncate after last ASCII character */
				n1 = n2;
			}
			if (n1 == 0) {
				sprintf(enc, "prop_%d", i + 1);
				n1 = strlen(enc);
			}
			fprintf(fp, "!:property ; %s\n", enc);
			for (j = 0; j < nframes; j++) {
				fprintf(fp, "%.18g\n", prp->propvals[j]);
			}
			fprintf(fp, "\n");
		}
	}
	
	if (mp->bset != NULL) {
		/*  Gaussian primitive info  */
		ShellInfo *sp;
		PrimInfo *pp;
		fprintf(fp, "!:gaussian_primitives\n");
		fprintf(fp, "! sym nprims a_idx; A C Csp\n");
		for (i = 0, sp = mp->bset->shells; i < mp->bset->nshells; i++, sp++) {
			switch (sp->sym) {
				case kGTOType_S:  p = "S";  break;
				case kGTOType_P:  p = "P";  break;
				case kGTOType_SP: p = "SP"; break;
				case kGTOType_D:  p = "D";  break;
				case kGTOType_D5: p = "D5"; break;
				case kGTOType_F:  p = "F";  break;
				case kGTOType_F7: p = "F7"; break;
				case kGTOType_G:  p = "G";  break;
				case kGTOType_G9: p = "G9"; break;
				default: snprintf(bufs[0], 8, "X%d", sp->sym); p = bufs[0]; break;
			}
			fprintf(fp, "%s %d %d\n", p, sp->nprim, sp->a_idx);
			pp = mp->bset->priminfos + sp->p_idx;
			for (j = 0; j < sp->nprim; j++, pp++) {
				fprintf(fp, "%.18g %.18g %.18g\n", pp->A, pp->C, pp->Csp);
			}
		}
		fprintf(fp, "\n");
		
		/*  MO info  */
		fprintf(fp, "!:mo_info\n");
		fprintf(fp, "! uhf|rhf|rohf ne_alpha ne_beta\n");
		switch (mp->bset->rflag) {
			case 0: p = "UHF"; break;
			case 1: p = "RHF"; break;
			case 2: p = "ROHF"; break;
			default: p = "(unknown)"; break;
		}
		fprintf(fp, "%s %d %d\n", p, mp->bset->ne_alpha, mp->bset->ne_beta);
		fprintf(fp, "\n");

		/*  MO coefficients  */
		fprintf(fp, "!:mo_coefficients\n");
		for (i = 0; i < mp->bset->nmos; i++) {
			fprintf(fp, "MO %d %.18g\n", i + 1, mp->bset->moenergies[i]);
			for (j = 0; j < mp->bset->ncomps; j++) {
				fprintf(fp, "%.18g%c", mp->bset->mo[i * mp->bset->ncomps + j], (j % 6 == 5 || j == mp->bset->ncomps - 1 ? '\n' : ' '));
			}
		}
		fprintf(fp, "\n");
	}

	if (mp->mview != NULL && mp->mview->ngraphics > 0) {
		MainViewGraphic *gp;
		fprintf(fp, "!:graphics\n");
		for (i = 0; i < mp->mview->ngraphics; i++) {
			gp = mp->mview->graphics + i;
			switch (gp->kind) {
				case kMainViewGraphicLine: fprintf(fp, "line\n"); break;
				case kMainViewGraphicPoly: fprintf(fp, "poly\n"); break;
				case kMainViewGraphicCylinder: fprintf(fp, "cylinder\n"); break;
				case kMainViewGraphicCone: fprintf(fp, "cone\n"); break;
				case kMainViewGraphicEllipsoid: fprintf(fp, "ellipsoid\n"); break;
				default: fprintf(fp, "unknown\n"); break;
			}
			fprintf(fp, "%d %d\n", gp->closed, gp->visible);
			fprintf(fp, "%.4f %.4f %.4f %.4f\n", gp->rgba[0], gp->rgba[1], gp->rgba[2], gp->rgba[3]);
			fprintf(fp, "%d\n", gp->npoints);
			for (j = 0; j < gp->npoints; j++)
				fprintf(fp, "%.6f %.6f %.6f\n", gp->points[j * 3], gp->points[j * 3 + 1], gp->points[j * 3 + 2]);
			fprintf(fp, "%d\n", gp->nnormals);
			for (j = 0; j < gp->nnormals; j++)
				fprintf(fp, "%.6f %.6f %.6f\n", gp->normals[j * 3], gp->normals[j * 3 + 1], gp->normals[j * 3 + 2]);
		}
		fprintf(fp, "\n");
	}
	
	/*  Plug-in in the Ruby world  */
	{
		char *outMessage;
		if (MolActionCreateAndPerform(mp, SCRIPT_ACTION(";s"),
									  "proc { savembsf_plugin rescue \"Plug-in error: #{$!.to_s}\" }", &outMessage) == 0) {
			if (outMessage[0] != 0) {
				if (strncmp(outMessage, "Plug-in", 7) == 0) {
					s_append_asprintf(errbuf, "%s", outMessage);
				} else {
					fprintf(fp, "%s\n", outMessage);
				}
			}
			free(outMessage);
		}
	}
	
	fclose(fp);
	return 0;
}

int
MoleculeWriteToPsfFile(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	int i;
	Atom *ap;
	*errbuf = NULL;
	fp = fopen(fname, "wb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot write to file %s", fname);
		return 1;
	}
	fprintf(fp, "PSF\n\n");
	fprintf(fp, "       1 !NTITLE\n");
	fprintf(fp, " REMARKS FILENAME=\n");
	fprintf(fp, "\n");
	
	/*  Atoms  */
	fprintf(fp, "%8d !NATOM\n", mp->natoms);
	for (i = 0; i < mp->natoms; i++) {
		const char *fmt;
		ap = ATOM_AT_INDEX(mp->atoms, i);
		fprintf(fp, "%8d ", i + 1);
		if (ap->resSeq >= 10000) {
			fmt = "%-3.3s %-5d ";
		} else {
			fmt = "%-4.4s %-4d ";
		}
		fprintf(fp, fmt, ap->segName, ap->resSeq);
		fprintf(fp, "%-3.3s  %-4.4s %-4.4s   %12.6f  %8.4f           0\n",
			ap->resName, ap->aname, AtomTypeDecodeToString(ap->type, NULL), ap->charge, ap->weight);
	}
	fprintf(fp, "\n");
	
	/*  Bonds  */
	fprintf(fp, "%8d !NBOND: bonds\n", mp->nbonds);
	for (i = 0; i < mp->nbonds * 2; i++) {
		fprintf(fp, "%8d", mp->bonds[i] + 1);
		if (i % 8 == 7)
			fprintf(fp, "\n");
	}
	if (i % 8 != 0)
		fprintf(fp, "\n");
	fprintf(fp, "\n");
	
	/*  Angles  */
	fprintf(fp, "%8d !NTHETA: angles\n", mp->nangles);
	for (i = 0; i < mp->nangles * 3; i++) {
		fprintf(fp, "%8d", mp->angles[i] + 1);
		if (i % 9 == 8)
			fprintf(fp, "\n");
	}
	if (i % 9 != 0)
		fprintf(fp, "\n");
	fprintf(fp, "\n");
	
	/*  Dihedrals  */
	fprintf(fp, "%8d !NPHI: dihedrals\n", mp->ndihedrals);
	for (i = 0; i < mp->ndihedrals * 4; i++) {
		fprintf(fp, "%8d", mp->dihedrals[i] + 1);
		if (i % 8 == 7)
			fprintf(fp, "\n");
	}
	if (i % 8 != 0)
		fprintf(fp, "\n");
	fprintf(fp, "\n");
	
	/*  Dihedrals  */
	fprintf(fp, "%8d !NIMPHI: impropers\n", mp->nimpropers);
	for (i = 0; i < mp->nimpropers * 4; i++) {
		fprintf(fp, "%8d", mp->impropers[i] + 1);
		if (i % 8 == 7)
			fprintf(fp, "\n");
	}
	if (i % 8 != 0)
		fprintf(fp, "\n");
	fprintf(fp, "\n");
	
	fprintf(fp, "%8d !NDON: donors\n\n", 0);
	fprintf(fp, "%8d !NACC: acceptors\n\n", 0);
	fprintf(fp, "%8d !NNB: non-bonding exclusions\n\n", 0);
	for (i = 0; i < mp->natoms; i++) {
		fprintf(fp, "%8d", 0);
		if (i % 8 == 7)
			fprintf(fp, "\n");
	}
	if (i % 8 != 0)
		fprintf(fp, "\n");
	fprintf(fp, "\n");
	fprintf(fp, "%8d !NGRP: groups\n", 1);
	fprintf(fp, "       0       0       0\n");
	fprintf(fp, "\n");
	
	i = strlen(fname);
	if (i > 5 && strcmp(fname + i - 5, ".psfx") == 0) {
		/*  Extended psf (with coordinates and other info)  */
		fprintf(fp, "%8d !COORD: coordinates\n", mp->natoms);
		for (i = 0; i < mp->natoms; i++) {
			Vector r;
			ap = ATOM_AT_INDEX(mp->atoms, i);
			r = ap->r;
			fprintf(fp, " %.8g %.8g %.8g ! %d,%.4s\n", r.x, r.y, r.z, i + 1, ap->aname);
		}
		fprintf(fp, "\n");
	}
		
	fclose(fp);
	return 0;
}

int
MoleculeWriteToPdbFile(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	int i, j;
	Atom *ap;
	*errbuf = NULL;
	fp = fopen(fname, "wb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot write to file %s", fname);
		return 1;
	}
	for (i = 0; i < mp->natoms; i++) {
		char buf[6];
		ap = ATOM_AT_INDEX(mp->atoms, i);
		if (ap->resSeq >= 10000) {
			snprintf(buf, sizeof buf, "%c%03d", 'A' + (ap->resSeq - 10000) / 1000, ap->resSeq % 1000);
		} else {
			snprintf(buf, sizeof buf, "%4d", ap->resSeq);
		}
		fprintf(fp, "ATOM  %5d %-4.4s%1.1s%-3.3s %1.1s%4.4s%1.1s   "
					"%8.3f%8.3f%8.3f %5.2f %5.2f      "
					"%-4.4s%-2.2s%-2d\n",
			i + 1, ap->aname, " ", ap->resName, " ", buf, " ",
			ap->r.x, ap->r.y, ap->r.z, ap->occupancy, ap->tempFactor,
			ap->segName, ap->element, ap->intCharge);
	}
	for (i = 0; i < mp->natoms; i++) {
		Int *cp;
		ap = ATOM_AT_INDEX(mp->atoms, i);
		cp = AtomConnectData(&ap->connect);
		for (j = 0; j < ap->connect.count; j++) {
			if (j % 4 == 0) {
				if (j > 0)
					fprintf(fp, "\n");
				fprintf(fp, "CONECT%5d", i + 1);
			}
			fprintf(fp, "%5d", cp[j] + 1);
		}
		if (j > 0)
			fprintf(fp, "\n");
	}
	fprintf(fp, "END\n");
	fclose(fp);
	return 0;
}

int
MoleculeWriteToDcdFile(Molecule *mp, const char *fname, char **errbuf)
{
	DcdRecord dcd;
	SFloat32 *xp, *yp, *zp;
	int n;
	*errbuf = NULL;
	if (mp == NULL || mp->natoms == 0) {
		s_append_asprintf(errbuf, "Molecule is empty");
		return 1;
	}
	memset(&dcd, 0, sizeof(dcd));
	dcd.natoms = mp->natoms;
	dcd.nframes = MoleculeGetNumberOfFrames(mp);
	if (dcd.nframes == 0) {
		s_append_asprintf(errbuf, "no frame is present");
		return 1;
	}
	dcd.nstart = mp->startStep;
	dcd.ninterval = mp->stepsPerFrame;
	if (dcd.ninterval == 0)
		dcd.ninterval = 1;
	dcd.nend = dcd.nstart + (dcd.nframes - 1) * dcd.ninterval;
	if (mp->cell != NULL)
		dcd.nextra = 1;
	dcd.delta = mp->psPerStep;
	if (dcd.delta == 0.0)
		dcd.delta = 1.0;
	dcd.ncharmver = 24;
	n = DcdCreate(fname, &dcd);
	if (n != 0) {
		if (n < 0)
			s_append_asprintf(errbuf, "Cannot create dcd file");
		else
			s_append_asprintf(errbuf, "Cannot write dcd header");
		DcdClose(&dcd);
		return 1;
	}
	
	xp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	yp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	zp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	if (xp == NULL || yp == NULL || zp == NULL) {
		s_append_asprintf(errbuf, "Cannot allocate memory");
		if (xp) free(xp);
		if (yp) free(yp);
		if (zp) free(zp);
		DcdClose(&dcd);
		return 1;
	}
	for (n = 0; n < dcd.nframes; n++) {
		int i;
		Atom *ap;
		for (i = 0, ap = mp->atoms; i < dcd.natoms && mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			Vector r;
			if (ap->frames == NULL || n >= ap->nframes)
				r = ap->r;
			else
				r = ap->frames[n];
			xp[i] = r.x;
			yp[i] = r.y;
			zp[i] = r.z;
		}
		if (i < dcd.natoms) {
			size_t sz = (dcd.natoms - i) * sizeof(SFloat32);
			memset(xp + i, 0, sz);
			memset(yp + i, 0, sz);
			memset(zp + i, 0, sz);
		}
		if (n < mp->nframe_cells && mp->frame_cells != NULL) {
			Vector *cp = &(mp->frame_cells[n * 4]);
			dcd.globalcell[0] = VecLength(cp[0]);
			dcd.globalcell[2] = VecLength(cp[1]);
			dcd.globalcell[5] = VecLength(cp[2]);
			dcd.globalcell[1] = VecDot(cp[0], cp[1]) / (dcd.globalcell[0] * dcd.globalcell[2]);
			dcd.globalcell[3] = VecDot(cp[0], cp[2]) / (dcd.globalcell[0] * dcd.globalcell[5]);
			dcd.globalcell[4] = VecDot(cp[1], cp[2]) / (dcd.globalcell[2] * dcd.globalcell[5]);			
		}			
		if (DcdWriteFrame(&dcd, n, xp, yp, zp, dcd.globalcell)) {
			s_append_asprintf(errbuf, "Write error in dcd file");
			goto exit;
		}
	}
	
exit:
	DcdClose(&dcd);
	free(xp);
	free(yp);
	free(zp);
	if (errbuf[0] == 0)
		return 0;
	else return 1;
}

int
MoleculeWriteExtendedInfo(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	int i;
	Vector v;
	*errbuf = NULL;
	fp = fopen(fname, "wb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot write to file %s", fname);
		return 1;
	}
	if (mp->cell != NULL) {
		fprintf(fp, "Bounding box:\n");
		for (i = 0; i < 3; i++) {
			v = mp->cell->axes[i];
			fprintf(fp, "%.3f %.3f %.3f %d\n", v.x, v.y, v.z, mp->cell->flags[i]);
		}
		fprintf(fp, "Bounding box origin:\n");
		v = mp->cell->origin;
		fprintf(fp, "%.3f %.3f %.3f\n", v.x, v.y, v.z);
	}
	fclose(fp);
	return 0;
}
		
 static int
sCompareByElement(const void *ap, const void *bp)
{
	return ((*(Atom **)bp)->atomicNumber - (*(Atom **)ap)->atomicNumber);
}

static int
sMakeAdc(int n, int base, Symop symop)
{
	int an, sym;
	if (SYMOP_ALIVE(symop)) {
		an = base;
		sym = (symop.dx + 5) * 10000 + (symop.dy + 5) * 1000 + (symop.dz + 5) * 100 + symop.sym + 1;
	} else {
		an = n;
		sym = 55501;
	}
	return (an + 1) * 100000 + sym;
}

static int
sCompareAdc(const void *ap, const void *bp)
{
	int n = *((Int *)ap) % 100000 - *((Int *)bp) % 100000;
	if (n == 0)
		n = *((Int *)ap) / 100000 - *((Int *)bp) / 100000;
	return n;
}

static void
sOutputAtomListInstructions(FILE *fp, int natoms, Atom *atoms)
{
	int i, j, k, an, sym;
	Atom *ap;
	Int *adc;
	adc = (Int *)malloc(sizeof(Int) * natoms);
	if (adc == NULL)
		return;
	for (i = 0, ap = atoms; i < natoms; i++, ap++) {
		if (ap->exflags & kAtomHiddenFlag)
			continue;
		adc[i] = sMakeAdc(i, ap->symbase, ap->symop);
	}
	mergesort(adc, natoms, sizeof(Int), sCompareAdc);
	
	/*  Create the atom list  */
	an = sym = -1;
	for (i = j = k = 0; i < natoms; i++) {
		int an1 = adc[i] / 100000;
		int sym1 = adc[i] % 100000;
		if (sym == sym1 && an1 == an + 1) {
			/*  Continuous  */
			an = an1;
			k++;
			continue;
		}
		if (k > 0)
			/*  Output the last atom with a minus sign  */
			adc[j++] = -(an * 100000 + sym);
		/*  Output this atom  */
		adc[j++] = adc[i];
		an = an1;
		sym = sym1;
		k = 0;
	}
	if (k > 0)
		adc[j++] = -(an * 100000 + sym);
	
	/*  Create the instruction cards  */
	for (i = k = 0; i < j; i++) {
		if (k == 0)
			fprintf(fp, "      401");
		fprintf(fp, "%9d", adc[i]);
		k++;
		if (i == j - 1 || k == 6 || (k == 5 && i < j - 2 && adc[i + 2] < 0)) {
			fprintf(fp, "\n");
			k = 0;
		}
	}
	free(adc);
}

static int
sEllipsoidType(int an)
{
	return (an >= 18 ? 3 : (an >= 2 && an != 6 ? 2 : (an > 0 ? 1 : 0)));
}

static void
sOutputAtomTypeInstructions(FILE *fp, int natoms, Atom *atoms)
{
	int i;
	Atom *ap;
	int etype, elast, istart, ilast, n1, n2;
	elast = istart = ilast = -1;
	for (i = 0, ap = atoms; i <= natoms; i++, ap++) {
		if (i < natoms) {
			if (SYMOP_ALIVE(ap->symop))
				continue;
			if (ap->exflags & kAtomHiddenFlag)
				continue;
			etype = sEllipsoidType(ap->atomicNumber);
			if (elast < 0) {
				istart = ilast = i;
				elast = etype;
				continue;
			} else if (elast == etype && ilast == i - 1) {
				ilast++;
				continue;
			}
		}
		/*  Output the instruction card for the 'last' block of atoms  */
		switch (etype) {
			case 2:
				n1 = 4; n2 = 0; break;
			case 3:
				n1 = 4; n2 = 5; break;
			default:
				n1 = 1; n2 = 0; break;
		}
		fprintf(fp, "  1   715 %8d        0 %8d        0    0.100    0.000    0.000\n", n1, n2);
		fprintf(fp, "                           %9d%9d\n", istart + 1, ilast + 1);
		elast = etype;
		ilast = istart = i;
	}
}

static int
sCompareBondType(const void *ap, const void *bp)
{
	/*  Descending order  */
	return *((int *)bp) - *((int *)ap);
}

static void
sOutputBondInstructions(FILE *fp, int natoms, Atom *atoms, int overlap_correction)
{
	Atom *ap, *ap2;
	char buf[96];
	int i, j, n[5], an, count, n1, n2, k;
	Int *cp;
	Int nexbonds;
	Int *exbonds;
	static const float sBondRad[4] = {0.060, 0.060, 0.060, 0.040};
	static const int sBondShade[4] = {5, 3, 1, 1};

	n[0] = n[1] = n[2] = n[3] = 0;  /*  Start index of 3rd row atoms (and higher), 2nd row, 1st row, and H */
	n[4] = natoms;
	for (i = natoms - 1, ap = atoms + i; i >= 0; i--, ap--) {
		an = ap->atomicNumber;
		if (an < 2)
			n[3] = i;
		if (an < 10)
			n[2] = i;
		if (an < 18)
			n[1] = i;
	}
	nexbonds = 0;
	exbonds = NULL;
	count = 0;

	if (overlap_correction)
		strcpy(buf, "  2  1001    0.000\n");
	else
		strcpy(buf, "  2   812\n");
	
	for (i = 0; i < 4; i++) {
		for (j = i; j < 4; j++) {
			/*  Examine bonds between "group i" and "group j"  */
			Vector dr;
			double d;
			double min_bond = 10000.0;     /*  Minimum distance between bound atoms  */
			double min_nonbond = 10000.0;  /*  Minimum distance between non-bound atoms  */
			double max_bond = -10000.0;    /*  Maximum distance between bound atoms  */
			int count_exbond = 0;          /*  Number of explicit bonds in this group  */
			for (n1 = n[i], ap = atoms + n1; n1 < n[i + 1]; n1++, ap++) {
				for (n2 = n[j], ap2 = atoms + n2; n2 < n[j + 1]; n2++, ap2++) {
					if (n1 == n2)
						continue;
					VecSub(dr, ap->r, ap2->r);
					d = VecLength(dr);
					cp = AtomConnectData(&ap->connect);
					for (k = ap->connect.count - 1; k >= 0; k--) {
						if (cp[k] == n2)
							break;
					}
					if (k >= 0) {
						/*  n1 and n2 are bound  */
						if (d < min_bond)
							min_bond = d;
						if (d > max_bond)
							max_bond = d;
					} else {
						/*  n1 and n2 are not bound  */
						if (d < min_nonbond)
							min_nonbond = d;
					}
				}
			}
			if (min_bond == 10000.0)
				continue;  /*  No bonds between these groups  */
			min_bond *= 0.9;
			if (max_bond + 0.002 < min_nonbond)
				max_bond += 0.002;
			else {
				max_bond = min_nonbond - 0.002;
				/*  Some bonds may be omitted, so scan all bonds again  */
				for (n1 = n[i], ap = ATOM_AT_INDEX(atoms, n1); n1 < n[i + 1]; n1++, ap = ATOM_NEXT(ap)) {
					cp = AtomConnectData(&ap->connect);
					for (k = ap->connect.count - 1; k >= 0; k--) {
						n2 = cp[k];
						if (n2 < n[j] || n2 >= n[j + 1])
							continue;
						ap2 = atoms + n2;
						VecSub(dr, ap->r, ap2->r);
						d = VecLength(dr);
						if (d > max_bond) {
							/*  This bond should be explicitly defined  */
							Int adc1, adc2;
							if (count_exbond == 0) {
								adc1 = -(i + 1);  /*  Bond type  */
								AssignArray(&exbonds, &nexbonds, sizeof(Int), nexbonds, &adc1);
							}
							adc1 = sMakeAdc(n1, ap->symbase, ap->symop);
							adc2 = sMakeAdc(n2, ap2->symbase, ap2->symop);
							AssignArray(&exbonds, &nexbonds, sizeof(Int), nexbonds, &adc1);
							AssignArray(&exbonds, &nexbonds, sizeof(Int), nexbonds, &adc2);
							count_exbond++;
						}
					}
				}
			}
			/*  Output the last instruction card  */
			fputs(buf, fp);
			/*  Make a new trailer card  */
			snprintf(buf, sizeof(buf), "  2      %3d%3d%3d%3d%3d%6.3f%6.3f%6.3f\n", n[i]+1, n[i+1], n[j]+1, n[j+1], sBondShade[i], min_bond, max_bond, sBondRad[i]);
			count++;
		}
	}
	if (count > 0) {
		/*  Output the last trailer card  */
		buf[2] = ' ';
		fputs(buf, fp);
	}
	if (nexbonds > 0) {
		if (count == 0 && overlap_correction) {
			/*  1001 card is not yet written, so write it  */
			buf[2] = ' ';
			fputs(buf, fp);
		}
		snprintf(buf, sizeof(buf), "  1   %3d", (overlap_correction ? 821 : 811));
		k = -exbonds[0] - 1;  /*  Bond type for the first block  */
		i = 1;  /*  Index for exbonds[]  */
		j = 0;  /*  Count in this block  */
		while (i <= nexbonds) {
			if (j >= 29 || i == nexbonds || exbonds[i] < 0) {
				/*  End of block  */
				buf[2] = '2';
				fputs(buf, fp);
				/*  The trailer card  */
				fprintf(fp, "                     %3d            %6.3f\n", sBondShade[k], sBondRad[k]);
				if (i == nexbonds)
					break;
				if (exbonds[i] < 0)
					k = -exbonds[i++] - 1;  /*  The new bond type  */
				j = 0;
			} else if (j > 0 && j % 3 == 0) {
				buf[2] = '1';
				fputs(buf, fp);
			}
			n1 = exbonds[i++];
			n2 = exbonds[i++];
			snprintf(buf + 9 + (j % 3) * 18, sizeof(buf) - 9 - (j % 3) * 18, "%9d%9d\n", n1, n2);
			j++;
		}
		free(exbonds);
	}
}

int
MoleculeWriteToTepFile(Molecule *mp, const char *fname, char **errbuf)
{
	FILE *fp;
	int i, j, natoms, *ip;
	Int *cp;
	Atom *ap, *atoms, **app;
	Double *dp;
	static Double sUnit[] = {1, 1, 1, 90, 90, 90};
	
	*errbuf = NULL;

	/*  Create sorted array of atoms  */
	natoms = mp->natoms;
	atoms = (Atom *)calloc(sizeof(Atom), natoms);
	app = (Atom **)calloc(sizeof(Atom *), natoms);
	ip = (int *)calloc(sizeof(int), natoms);
	if (atoms == NULL || app == NULL || ip == NULL) {
		s_append_asprintf(errbuf, "Cannot allocate memory");
		return 1;
	}
	/*  Sort the atom pointer by atomic number  */
	for (i = 0, ap = mp->atoms; i < natoms; i++, ap = ATOM_NEXT(ap))
		app[i] = ap;
	mergesort(app, natoms, sizeof(Atom *), sCompareByElement);
	for (i = 0; i < natoms; i++) {
		/*  ip[old_index] is new_index  */
		ip[app[i] - mp->atoms] = i;
	}
	/*  Copy the atom record to atoms[]  */
	/*  The 'v' member contains crystallographic coordinates  */
	/*  The connection table and symbase are renumbered  */
	/*  Hidden flags are modified to reflect the visibility in the MainView  */
	for (i = 0, ap = atoms; i < natoms; i++, ap++) {
		AtomDuplicateNoFrame(ap, app[i]);
	/*	memmove(ap, app[i], gSizeOfAtomRecord); */
		MoleculeCartesianToXtal(mp, &(ap->v), &(ap->r));
		cp = AtomConnectData(&ap->connect);
		for (j = ap->connect.count - 1; j >= 0; j--) {
			cp[j] = ip[cp[j]];
		}
		if (SYMOP_ALIVE(ap->symop))
			ap->symbase = ip[ap->symbase];
		if (MainView_isAtomHidden(mp->mview, i)) {
			ap->exflags |= kAtomHiddenFlag;
		} else {
			ap->exflags &= ~kAtomHiddenFlag;
		}
	}
	free(ip);
	free(app);
	
	fp = fopen(fname, "wb");
	if (fp == NULL) {
		s_append_asprintf(errbuf, "Cannot write to file %s", fname);
		return 1;
	}

	/*  Title line  */
	fprintf(fp, "Generated by Molby\n");
	
	/*  XtalCell  */
	if (mp->cell != NULL) {
		dp = mp->cell->cell;
	} else {
		dp = sUnit;
	}
	fprintf(fp, "%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f\n", dp[0], dp[1], dp[2], dp[3], dp[4], dp[5]);
	
	/*  Symmetry operations  */
	if (mp->nsyms > 0) {
		for (i = 0; i < mp->nsyms; i++) {
			dp = mp->syms[i];
			fprintf(fp, "%c%14g%3g%3g%3g%15g%3g%3g%3g%15g%3g%3g%3g\n", (i == mp->nsyms - 1 ? '1' : ' '), dp[9], dp[0], dp[1], dp[2], dp[10], dp[3], dp[4], dp[5], dp[11], dp[6], dp[7], dp[8]);
		}
	} else {
		fprintf(fp, "1             0  1  0  0              0  0  1  0              0  0  0  1\n");
	}
	
	/*  Atoms  */
	for (i = 0, ap = atoms; i < natoms; i++, ap++) {
		/*  The 'v' field contains crystallographic coordinates  */
		fprintf(fp, " %4.4s%22s%9.4f%9.4f%9.4f%9d\n", ap->aname, "", ap->v.x, ap->v.y, ap->v.z, 0);
		if (ap->aniso != NULL) {
			dp = ap->aniso->bij;
			fprintf(fp, " %8.5f%9.6f%9.6f%9.6f%9.6f%9.6f%9d\n", dp[0], dp[1], dp[2], dp[3], dp[4], dp[5], 0);
		} else {
			Double temp = ap->tempFactor;
			if (temp <= 0)
				temp = 1.2;
			fprintf(fp, " %8.3f%9g%9g%9g%9g%9g%9d\n", temp, 0.0, 0.0, 0.0, 0.0, 0.0, 6);
		}
	}
	/*  Special points  */
	{
		Vector camera, lookat, up, xvec, yvec, zvec;
		MainView_getCamera(mp->mview, &camera, &lookat, &up);
		VecSub(zvec, lookat, camera);
		VecCross(xvec, zvec, up);
		NormalizeVec(&xvec, &xvec);
		NormalizeVec(&yvec, &up);
		VecInc(xvec, lookat);
		VecInc(yvec, lookat);
		MoleculeCartesianToXtal(mp, &lookat, &lookat);
		MoleculeCartesianToXtal(mp, &xvec, &xvec);
		MoleculeCartesianToXtal(mp, &yvec, &yvec);
		fprintf(fp, " ORGN                      %9g%9g%9g        0\n", 0.0, 0.0, 0.0);
		fprintf(fp, " %8.3f%9g%9g%9g%9g%9g%9d\n", 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 6);
		fprintf(fp, " CNTR                      %9g%9g%9g        0\n", lookat.x, lookat.y, lookat.z);
		fprintf(fp, " %8.3f%9g%9g%9g%9g%9g%9d\n", 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 6);
		fprintf(fp, " X                         %9g%9g%9g        0\n", xvec.x, xvec.y, xvec.z);
		fprintf(fp, " %8.3f%9g%9g%9g%9g%9g%9d\n", 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 6);
		fprintf(fp, " Y                         %9g%9g%9g        0\n", yvec.x, yvec.y, yvec.z);
		fprintf(fp, "1%8.3f%9g%9g%9g%9g%9g%9d\n", 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 6);
	}
	
	/*  Instructions  */
	fprintf(fp, "      201\n");
	fprintf(fp, "      205       12\n");
	fprintf(fp, "      301      6.6      6.6        0      0.8\n");
	sOutputAtomListInstructions(fp, natoms, atoms);
	fprintf(fp, "      501%4d55501%4d55501%4d55501%4d55501%4d55501                 1\n", natoms + 2, natoms + 2, natoms + 3, natoms + 2, natoms + 4);
	fprintf(fp, "      502        1      0.0        2      0.0        3      0.0\n");
	fprintf(fp, "      604                               1.538\n");

	sOutputBondInstructions(fp, natoms, atoms, 1);
	sOutputAtomTypeInstructions(fp, natoms, atoms);
	sOutputBondInstructions(fp, natoms, atoms, 0);

	for (i = 0; i < natoms; i++) {
		AtomClean(atoms + i);
	}
	free(atoms);

	fprintf(fp, "      202\n");
	fprintf(fp, "  0    -1\n");
	fclose(fp);
	return 0;
}

void
MoleculeDump(Molecule *mol)
{
	int i, j;
	Int *cp;
	Atom *ap;
	for (i = 0; i < mol->natoms; i++) {
		char buf1[8];
		ap = ATOM_AT_INDEX(mol->atoms, i);
		snprintf(buf1, sizeof buf1, "%3.4s.%d", ap->resName, ap->resSeq);
		fprintf(stderr, "%4d %-7s %-4.6s %-4.6s %-2.2s %7.3f %7.3f %7.3f %6.3f [", i, buf1, ap->aname, AtomTypeDecodeToString(ap->type, NULL), ap->element, ap->r.x, ap->r.y, ap->r.z, ap->charge);
		cp = AtomConnectData(&ap->connect);
		for (j = 0; j < ap->connect.count; j++) {
			fprintf(stderr, "%s%d", (j > 0 ? "," : ""), cp[j]);
		}
		fprintf(stderr, "]\n");
	}
}

#pragma mark ====== MD support (including modification of Molecule) ======

/*  Call md_prepare for the MDArena. If MDArena has not been created, a new arena is created.
	If something goes wrong, returns 1 (for missing parameters) or -1 (more serious error).
    If retmsg is not NULL, a message describing the problem is returned there. This message
    must be free'd by the caller.  */
int
MoleculePrepareMDArena(Molecule *mol, int check_only, char **retmsg)
{
	const char *msg;
	Int nangles, *angles, ndihedrals, *dihedrals, nimpropers, *impropers;
	Int missing = 0;
	IntGroup *ig1, *ig2, *ig3;
	MDArena *arena = mol->arena;

	if (arena == NULL) {
		md_arena_new(mol);
		arena = mol->arena;
	} else if (arena->xmol != mol)
		md_arena_set_molecule(arena, mol);

	arena->is_initialized = 0;
	
	/*  Rebuild the tables  */
	ig1 = ig2 = ig3 = NULL;
	nangles = MoleculeFindMissingAngles(mol, &angles);
	ndihedrals = MoleculeFindMissingDihedrals(mol, &dihedrals);
	nimpropers = MoleculeFindMissingImpropers(mol, &impropers);
	if (nangles > 0) {
		ig1 = IntGroupNewWithPoints(mol->nangles, nangles, -1);
		MolActionCreateAndPerform(mol, gMolActionAddAngles, nangles * 3, angles, ig1);
		free(angles);
		IntGroupRelease(ig1);
	}
	if (ndihedrals > 0) {
		ig2 = IntGroupNewWithPoints(mol->ndihedrals, ndihedrals, -1);
		MolActionCreateAndPerform(mol, gMolActionAddDihedrals, ndihedrals * 4, dihedrals, ig2);
		free(dihedrals);
		IntGroupRelease(ig2);
	}
	if (nimpropers > 0) {
		ig3 = IntGroupNewWithPoints(mol->nimpropers, nimpropers, -1);
		MolActionCreateAndPerform(mol, gMolActionAddImpropers, nimpropers * 4, impropers, ig3);
		free(impropers);
		IntGroupRelease(ig3);
	}
	
	{
		/*  Update the path information of the molecule before MD setup  */
		char *buf = (char *)malloc(4096);
		MoleculeCallback_pathName(mol, buf, sizeof buf);
		MoleculeSetPath(mol, buf);
		free(buf);
	}
		
	/*  Prepare parameters and internal information  */
	msg = md_prepare(arena, check_only);
	
	/*  Some parameters are missing?  */
	if (msg != NULL) {
		if (strstr(msg, "parameter") != NULL && strstr(msg, "missing") != NULL)
			missing = 1;
		else {
			if (retmsg != NULL)
				asprintf(retmsg, "cannot initialize for MD: %s", msg);
			return -1;
		}
	}
	
	/*  The local parameter list is updated  */
	{
		Int parType, idx;
		if (mol->par == NULL)
			mol->par = ParameterNew();
		for (parType = kFirstParType; parType <= kLastParType; parType++) {
			/*  Delete global and undefined parameters  */
			UnionPar *up, *upbuf;
			Int nparams, count;
			ig1 = IntGroupNew();
			for (idx = 0; (up = ParameterGetUnionParFromTypeAndIndex(mol->par, parType, idx)) != NULL; idx++) {
				if (up->bond.src != 0)
					IntGroupAdd(ig1, idx, 1);
			}
			if (IntGroupGetCount(ig1) > 0)
				MolActionCreateAndPerform(mol, gMolActionDeleteParameters, parType, ig1);
			IntGroupRelease(ig1);
			/*  Copy global and undefined parameters from arena and insert to mol->par  */
			nparams = ParameterGetCountForType(arena->par, parType);
			if (nparams == 0)
				continue;
			upbuf = (UnionPar *)calloc(sizeof(UnionPar), nparams);
			ig1 = IntGroupNew();
			ig2 = IntGroupNew();
			for (idx = 0; (up = ParameterGetUnionParFromTypeAndIndex(arena->par, parType, idx)) != NULL; idx++) {
				if (up->bond.src > 0)
					IntGroupAdd(ig1, idx, 1); /* Global parameter */
				else if (up->bond.src < 0)
					IntGroupAdd(ig2, idx, 1); /* Undefined parameter */
			}
			if ((count = IntGroupGetCount(ig1)) > 0) {
				/*  Insert global parameters (at the top)  */
				ParameterCopy(arena->par, parType, upbuf, ig1);
				ig3 = IntGroupNewWithPoints(0, count, -1);
				MolActionCreateAndPerform(mol, gMolActionAddParameters, parType, ig3, count, upbuf);
				IntGroupRelease(ig3);
			}
			if ((count = IntGroupGetCount(ig2)) > 0) {
				/*  Insert undefined parameters (at the bottom)  */
				ParameterCopy(arena->par, parType, upbuf, ig2);
				idx = ParameterGetCountForType(mol->par, parType);
				ig3 = IntGroupNewWithPoints(idx, count, -1);
				MolActionCreateAndPerform(mol, gMolActionAddParameters, parType, ig3, count, upbuf);
				IntGroupRelease(ig3);
			}
			IntGroupRelease(ig2);
			IntGroupRelease(ig1);
			free(upbuf);
		}
		mol->needsMDRebuild = 0;  /*  We know the "modified" parameters are consistent with the MDArena  */
	}
	
	if (missing) {
		if (retmsg != NULL)
			*retmsg = strdup(msg);
		return 1;
	} else return 0;
}

#pragma mark ====== Serialize ======

Molecule *
MoleculeDeserialize(const char *data, Int length, Int *timep)
{
	Molecule *mp;
	Parameter *par;
	Atom *ap;
/*	int result; */

	mp = MoleculeNew();
	if (mp == NULL)
		goto out_of_memory;
	par = ParameterNew();
	if (par == NULL)
		goto out_of_memory;

	while (length >= 12) {
		const char *ptr = data + 8 + sizeof(Int);
		int len = *((const Int *)(data + 8));
		int i, j, n;
		if (strcmp(data, "ATOM") == 0) {
			n = len / gSizeOfAtomRecord;
			NewArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, n);
			memmove(mp->atoms, ptr, len);
		} else if (strcmp(data, "ANISO") == 0) {
			n = len / (sizeof(Int) + sizeof(Aniso));
			for (i = 0; i < n; i++) {
				j = *((const Int *)ptr);
				if (j < 0 || j >= mp->natoms)
					goto bad_format;
				ap = ATOM_AT_INDEX(mp->atoms, j);
				ap->aniso = (Aniso *)calloc(sizeof(Aniso), 1);
				if (ap->aniso == NULL)
					goto out_of_memory;
				*(ap->aniso) = *((Aniso *)(ptr + sizeof(Int)));
				ptr += sizeof(Int) + sizeof(Aniso);
			}
		} else if (strcmp(data, "FRAME") == 0) {
			for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
				if (ap->nframes == 0)
					continue;
				n = ap->nframes;
				ap->frames = NULL;
				ap->nframes = 0;
				NewArray(&ap->frames, &ap->nframes, sizeof(Vector), n);
				if (ap->frames == NULL)
					goto out_of_memory;
				memmove(ap->frames, ptr, sizeof(Vector) * ap->nframes);
				ptr += sizeof(Vector) * ap->nframes;
			}
		} else if (strcmp(data, "EXTCON") == 0) {
			for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
				if (ap->connect.count <= ATOM_CONNECT_LIMIT)
					continue;
				n = ap->connect.count;
				ap->connect.count = 0;
				ap->connect.u.ptr = NULL;
				NewArray(&(ap->connect.u.ptr), &(ap->connect.count), sizeof(Int), n);
				memmove(ap->connect.u.ptr, ptr, sizeof(Int) * n);
				ptr += sizeof(Int) * n;
			}
		} else if (strcmp(data, "BOND") == 0) {
			n = len / (sizeof(Int) * 2);
			NewArray(&mp->bonds, &mp->nbonds, sizeof(Int) * 2, n);
			memmove(mp->bonds, ptr, len);
		} else if (strcmp(data, "ANGLE") == 0) {
			n = len / (sizeof(Int) * 3);
			NewArray(&mp->angles, &mp->nangles, sizeof(Int) * 3, n);
			memmove(mp->angles, ptr, len);
		} else if (strcmp(data, "DIHED") == 0) {
			n = len / (sizeof(Int) * 4);
			NewArray(&mp->dihedrals, &mp->ndihedrals, sizeof(Int) * 4, n);
			memmove(mp->dihedrals, ptr, len);
		} else if (strcmp(data, "IMPROP") == 0) {
			n = len / (sizeof(Int) * 4);
			NewArray(&mp->impropers, &mp->nimpropers, sizeof(Int) * 4, n);
			memmove(mp->impropers, ptr, len);
		} else if (strcmp(data, "RESIDUE") == 0) {
			n = len / 4;
			NewArray(&mp->residues, &mp->nresidues, 4, n);
			memmove(mp->residues, ptr, len);
		} else if (strcmp(data, "CELL") == 0) {
			mp->cell = (XtalCell *)malloc(sizeof(XtalCell));
			if (mp->cell == NULL)
				goto out_of_memory;
			memmove(mp->cell, ptr, sizeof(XtalCell));
		} else if (strcmp(data, "SYMOP") == 0) {
			n = len / sizeof(Transform);
			NewArray(&mp->syms, &mp->nsyms, sizeof(Transform), n);
			memmove(mp->syms, ptr, len);
		} else if (strcmp(data, "ANCHOR") == 0) {
			const char *ptr2 = ptr + len;
			while (ptr < ptr2) {
				PiAnchor an;
				memset(&an, 0, sizeof(an));
				i = *((Int *)ptr);
				if (i >= 0 && i < mp->natoms) {
					n = *((Int *)(ptr + sizeof(Int)));
					AtomConnectResize(&(an.connect), n);
					memmove(AtomConnectData(&(an.connect)), ptr + sizeof(Int) * 2, sizeof(Int) * n);
					NewArray(&an.coeffs, &an.ncoeffs, sizeof(Double), n);
					memmove(an.coeffs, ptr + sizeof(Int) * (2 + n), sizeof(Double) * n);
					ap = ATOM_AT_INDEX(mp->atoms, i);
					ap->anchor = (PiAnchor *)malloc(sizeof(PiAnchor));
					memmove(ap->anchor, &an, sizeof(PiAnchor));
				}
				ptr += sizeof(Int) * (2 + n) + sizeof(Double) * n;
			}
		} else if (strcmp(data, "TIME") == 0) {
			if (timep != NULL)
				*timep = *((Int *)ptr);
		} else if (strcmp(data, "BONDPAR") == 0) {
			mp->par = par;
			n = len / sizeof(BondPar);
			NewArray(&par->bondPars, &par->nbondPars, sizeof(BondPar), n);
			memmove(par->bondPars, ptr, len);
		} else if (strcmp(data, "ANGPAR") == 0) {
			mp->par = par;
			n = len / sizeof(AnglePar);
			NewArray(&par->anglePars, &par->nanglePars, sizeof(AnglePar), n);
			memmove(par->anglePars, ptr, len);
		} else if (strcmp(data, "DIHEPAR") == 0) {
			mp->par = par;
			n = len / sizeof(TorsionPar);
			NewArray(&par->dihedralPars, &par->ndihedralPars, sizeof(TorsionPar), n);
			memmove(par->dihedralPars, ptr, len);
		} else if (strcmp(data, "IMPRPAR") == 0) {
			mp->par = par;
			n = len / sizeof(TorsionPar);
			NewArray(&par->improperPars, &par->nimproperPars, sizeof(TorsionPar), n);
			memmove(par->improperPars, ptr, len);
		} else if (strcmp(data, "VDWPAR") == 0) {
			mp->par = par;
			n = len / sizeof(VdwPar);
			NewArray(&par->vdwPars, &par->nvdwPars, sizeof(VdwPar), n);
			memmove(par->vdwPars, ptr, len);
		} else if (strcmp(data, "VDWPPAR") == 0) {
			mp->par = par;
			n = len / sizeof(VdwPairPar);
			NewArray(&par->vdwpPars, &par->nvdwpPars, sizeof(VdwPairPar), n);
			memmove(par->vdwpPars, ptr, len);
		} else if (strcmp(data, "VCUTPAR") == 0) {
			mp->par = par;
			n = len / sizeof(VdwCutoffPar);
			NewArray(&par->vdwCutoffPars, &par->nvdwCutoffPars, sizeof(VdwCutoffPar), n);
			memmove(par->vdwCutoffPars, ptr, len);
		}
		len += 8 + sizeof(Int);
		data += len;
		length -= len;
	}
	if (mp->par == NULL)
		ParameterRelease(par);
/*	result = MoleculeRebuildTablesFromConnects(mp);
	if (result != 0)
		goto bad_format; */
	return mp;
	
  out_of_memory:
	Panic("Low memory while deserializing molecule data");
	return NULL; /* Not reached */

  bad_format:
	Panic("internal error: bad format during deserializing molecule data");
	return NULL; /* Not reached */
}

char *
MoleculeSerialize(Molecule *mp, Int *outLength, Int *timep)
{
	char *ptr, *p;
	int len, len_all, i, naniso, nframes, nconnects, nanchors;
	Atom *ap;

	/*  Array of atoms  */
	len = 8 + sizeof(Int) + gSizeOfAtomRecord * mp->natoms;
	ptr = (char *)malloc(len);
	if (ptr == NULL)
		goto out_of_memory;
	memmove(ptr, "ATOM\0\0\0\0", 8);
	*((Int *)(ptr + 8)) = gSizeOfAtomRecord * mp->natoms;
	p = ptr + 8 + sizeof(Int);
	memmove(p, mp->atoms, gSizeOfAtomRecord * mp->natoms);
	naniso = nframes = nconnects = nanchors = 0;
	for (i = 0; i < mp->natoms; i++) {
		ap = ATOM_AT_INDEX(p, i);
		if (ap->aniso != NULL) {
			naniso++;
			ap->aniso = NULL;
		}
		if (ap->frames != NULL) {
			nframes += ap->nframes;
			ap->frames = NULL;
		}
		if (ap->connect.count > ATOM_CONNECT_LIMIT) {
			nconnects += ap->connect.count;
			ap->connect.u.ptr = NULL;
		}
		if (ap->anchor != NULL) {
			nanchors++;
			ap->anchor = NULL;
		}
	}
	len_all = len;

	/*  Array of aniso  */
	if (naniso > 0) {
		len = 8 + sizeof(Int) + (sizeof(Int) + sizeof(Aniso)) * naniso;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "ANISO\0\0\0", 8);
		*((Int *)(p + 8)) = (sizeof(Int) + sizeof(Aniso)) * naniso;
		p += 8 + sizeof(Int);
		for (i = 0; i < mp->natoms; i++) {
			ap = ATOM_AT_INDEX(mp->atoms, i);
			if (ap->aniso != NULL) {
				*((Int *)p) = i;
				*((Aniso *)(p + sizeof(Int))) = *(ap->aniso);
				p += sizeof(Int) + sizeof(Aniso);
			}
		}
		len_all += len;
	}
	
	/*  Array of frames  */
	if (nframes > 0) {
		len = 8 + sizeof(Int) + sizeof(Vector) * nframes;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "FRAME\0\0\0", 8);
		*((Int *)(p + 8)) = sizeof(Vector) * nframes;
		p += 8 + sizeof(Int);
		for (i = 0; i < mp->natoms; i++) {
			ap = ATOM_AT_INDEX(mp->atoms, i);
			if (ap->frames != NULL) {
				memmove(p, ap->frames, sizeof(Vector) * ap->nframes);
				p += sizeof(Vector) * ap->nframes;
			}
		}
		len_all += len;
	}
	
	/*  Array of connects  */
	if (nconnects > 0) {
		len = 8 + sizeof(Int) + sizeof(Int) * nconnects;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "EXTCON\0\0", 8);
		*((Int *)(p + 8)) = sizeof(Int) * nconnects;
		p += 8 + sizeof(Int);
		for (i = 0; i < mp->natoms; i++) {
			ap = ATOM_AT_INDEX(mp->atoms, i);
			if (ap->connect.count > ATOM_CONNECT_LIMIT) {
				memmove(p, ap->connect.u.ptr, sizeof(Int) * ap->connect.count);
				p += sizeof(Int) * ap->connect.count;
			}
		}
		len_all += len;
	}
	
	/*  Bonds, angles, dihedrals, impropers  */
	if (mp->nbonds > 0) {
		len = 8 + sizeof(Int) + sizeof(Int) * 2 * mp->nbonds;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "BOND\0\0\0\0", 8);
		*((Int *)(p + 8)) = sizeof(Int) * 2 * mp->nbonds;
		p += 8 + sizeof(Int);
		memmove(p, mp->bonds, sizeof(Int) * 2 * mp->nbonds);
		len_all += len;
	}
	if (mp->nangles > 0) {
		len = 8 + sizeof(Int) + sizeof(Int) * 3 * mp->nangles;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "ANGLE\0\0\0", 8);
		*((Int *)(p + 8)) = sizeof(Int) * 3 * mp->nangles;
		p += 8 + sizeof(Int);
		memmove(p, mp->angles, sizeof(Int) * 3 * mp->nangles);
		len_all += len;
	}
	if (mp->ndihedrals > 0) {
		len = 8 + sizeof(Int) + sizeof(Int) * 4 * mp->ndihedrals;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "DIHED\0\0\0", 8);
		*((Int *)(p + 8)) = sizeof(Int) * 4 * mp->ndihedrals;
		p += 8 + sizeof(Int);
		memmove(p, mp->dihedrals, sizeof(Int) * 4 * mp->ndihedrals);
		len_all += len;
	}
	if (mp->nimpropers > 0) {
		len = 8 + sizeof(Int) + sizeof(Int) * 4 * mp->nimpropers;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "IMPROP\0\0", 8);
		*((Int *)(p + 8)) = sizeof(Int) * 4 * mp->nimpropers;
		p += 8 + sizeof(Int);
		memmove(p, mp->impropers, sizeof(Int) * 4 * mp->nimpropers);
		len_all += len;
	}
	
	/*  Array of residues  */
	if (mp->nresidues > 0) {
		len = 8 + sizeof(Int) + 4 * mp->nresidues;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "RESIDUE\0", 8);
		*((Int *)(p + 8)) = 4 * mp->nresidues;
		p += 8 + sizeof(Int);
		memmove(p, mp->residues, 4 * mp->nresidues);
		len_all += len;
	}

	/*  Unit cell  */
	if (mp->cell != NULL) {
		len = 8 + sizeof(Int) + sizeof(XtalCell);
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "CELL\0\0\0\0", 8);
		*((Int *)(p + 8)) = sizeof(XtalCell);
		p += 8 + sizeof(Int);
		memmove(p, mp->cell, sizeof(XtalCell));
		len_all += len;
	}
	
	/*  Symmetry operations  */
	if (mp->nsyms > 0) {
		len = 8 + sizeof(Int) + sizeof(Transform) * mp->nsyms;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "SYMOP\0\0\0", 8);
		*((Int *)(p + 8)) = sizeof(Transform) * mp->nsyms;
		p += 8 + sizeof(Int);
		memmove(p, mp->syms, sizeof(Transform) * mp->nsyms);
		len_all += len;
	}
	
	/*  Pi-anchors  */
	if (nanchors > 0) {
		/*  Estimate the necessary storage first  */
		/*  One entry consists of { atom_index (Int), number_of_connects (Int), connects (Int's), weights (Double's) }  */
		len = 8 + sizeof(Int);
		for (i = 0; i < mp->natoms; i++) {
			ap = ATOM_AT_INDEX(mp->atoms, i);
			if (ap->anchor != NULL)
				len += sizeof(Int) * 2 + (sizeof(Int) + sizeof(Double)) * ap->anchor->connect.count;
		}
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "ANCHOR\0\0", 8);
		*((Int *)(p + 8)) = len - (8 + sizeof(Int));
		p += 8 + sizeof(Int);
		for (i = 0; i < mp->natoms; i++) {
			Int count, *ip;
			ap = ATOM_AT_INDEX(mp->atoms, i);
			if (ap->anchor != NULL) {
				count = ap->anchor->connect.count;
				*((Int *)p) = i;
				*((Int *)(p + sizeof(Int))) = count;
				p += sizeof(Int) * 2;
				ip = AtomConnectData(&(ap->anchor->connect));
				memmove(p, ip, sizeof(Int) * count);
				p += sizeof(Int) * count;
				memmove(p, ap->anchor->coeffs, sizeof(Double) * count);
				p += sizeof(Double) * count;
			}
		}
		len_all += len;
	}
	
	/*  Parameters  */
	if (mp->par != NULL) {
		int type;
		for (type = kFirstParType; type <= kLastParType; type++) {
			const char *parname;
			Int parsize, parcount;
			void *parptr;
			switch (type) {
				case kBondParType:
					parname = "BONDPAR\0";
					parsize = sizeof(BondPar);
					parcount = mp->par->nbondPars;
					parptr = mp->par->bondPars;
					break;
				case kAngleParType:
					parname = "ANGPAR\0\0";
					parsize = sizeof(AnglePar);
					parcount = mp->par->nanglePars;
					parptr = mp->par->anglePars;
					break;
				case kDihedralParType:
					parname = "DIHEPAR\0";
					parsize = sizeof(TorsionPar);
					parcount = mp->par->ndihedralPars;
					parptr = mp->par->dihedralPars;
					break;
				case kImproperParType:
					parname = "IMPRPAR\0";
					parsize = sizeof(TorsionPar);
					parcount = mp->par->nimproperPars;
					parptr = mp->par->improperPars;
					break;
				case kVdwParType:
					parname = "VDWPAR\0\0";
					parsize = sizeof(VdwPar);
					parcount = mp->par->nvdwPars;
					parptr = mp->par->vdwPars;
					break;
				case kVdwPairParType:
					parname = "VDWPPAR\0";
					parsize = sizeof(VdwPairPar);
					parcount = mp->par->nvdwpPars;
					parptr = mp->par->vdwpPars;
					break;
				case kVdwCutoffParType:
					parname = "VCUTPAR\0";
					parsize = sizeof(VdwCutoffPar);
					parcount = mp->par->nvdwCutoffPars;
					parptr = mp->par->vdwCutoffPars;
					break;
				default:
					continue;
			}
			if (parcount > 0) {
				len = 8 + sizeof(Int) + parsize * parcount;
				ptr = (char *)realloc(ptr, len_all + len);
				if (ptr == NULL)
					goto out_of_memory;
				p = ptr + len_all;
				memmove(p, parname, 8);
				*((Int *)(p + 8)) = parsize * parcount;
				p += 8 + sizeof(Int);
				memmove(p, parptr, parsize * parcount);
				len_all += len;
			}
		}
	}
	
	/*  Time stamp  */
	{
		time_t tm = time(NULL);
		len = 8 + sizeof(Int) + sizeof(Int);
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "TIME\0\0\0\0", 8);
		*((Int *)(p + 8)) = sizeof(Int);
		p += 8 + sizeof(Int);
		*((Int *)p) = (Int)tm;
		len_all += len;
		if (timep != NULL)
			*timep = (Int)tm;
	}
	
	if (outLength != NULL)
		*outLength = len_all;
	return ptr;

  out_of_memory:
    Panic("Low memory while serializing a molecule data");
	return NULL; /* Not reached */	
}

#pragma mark ====== Search for bonds, angles, dihedrals, impropers ======

static IntGroup *
sMoleculeSearchIncludingAtoms(int nitems, Int *items, int nsize, IntGroup *atomgroup, const char *msg)
{
	int i, j;
	Int *ip;
	IntGroup *gp = NULL;
	if (atomgroup == NULL)
		return NULL;
	for (i = 0, ip = items; i < nitems; i++, ip += nsize) {
		for (j = 0; j < nsize; j++) {
			if (IntGroupLookup(atomgroup, ip[j], NULL) != 0) {
				if (gp == NULL)
					gp = IntGroupNew();
				if (gp == NULL || IntGroupAdd(gp, i, 1) != 0)
					Panic("Low memory while searching %s", msg);
				break;
			}
		}
	}
	return gp;
}

IntGroup *
MoleculeSearchBondsIncludingAtoms(Molecule *mp, IntGroup *atomgroup)
{
	if (mp == NULL)
		return NULL;
	return sMoleculeSearchIncludingAtoms(mp->nbonds, mp->bonds, 2, atomgroup, "bonds");
}

IntGroup *
MoleculeSearchAnglesIncludingAtoms(Molecule *mp, IntGroup *atomgroup)
{
	if (mp == NULL)
		return NULL;
	return sMoleculeSearchIncludingAtoms(mp->nangles, mp->angles, 3, atomgroup, "angles");
}

IntGroup *
MoleculeSearchDihedralsIncludingAtoms(Molecule *mp, IntGroup *atomgroup)
{
	if (mp == NULL)
		return NULL;
	return sMoleculeSearchIncludingAtoms(mp->ndihedrals, mp->dihedrals, 4, atomgroup, "dihedrals");
}

IntGroup *
MoleculeSearchImpropersIncludingAtoms(Molecule *mp, IntGroup *atomgroup)
{
	if (mp == NULL)
		return NULL;
	return sMoleculeSearchIncludingAtoms(mp->nimpropers, mp->impropers, 4, atomgroup, "impropers");
}

static IntGroup *
sMoleculeSearchAcrossAtomGroup(int nitems, Int *items, int nsize, IntGroup *atomgroup, const char *msg)
{
	int i, j;
	Int *ip;
	IntGroup *gp = NULL;
	if (atomgroup == NULL)
		return NULL;
	for (i = 0, ip = items; i < nitems; i++, ip += nsize) {
		int k = -1;
		for (j = 0; j < nsize; j++) {
			int kk;
			kk = (IntGroupLookup(atomgroup, ip[j], NULL) != 0);
			if (k < 0)
				k = kk;
			else if (k != kk) {
				/*  This bond etc. crosses the atom group border  */
				if (gp == NULL)
					gp = IntGroupNew();
				if (gp == NULL || IntGroupAdd(gp, i, 1) != 0)
					Panic("Low memory while searching %s", msg);
				break;
			}
		}
	}
	return gp;
}

IntGroup *
MoleculeSearchBondsAcrossAtomGroup(Molecule *mp, IntGroup *atomgroup)
{
	if (mp == NULL)
		return NULL;
	return sMoleculeSearchAcrossAtomGroup(mp->nbonds, mp->bonds, 2, atomgroup, "bonds");
}

IntGroup *
MoleculeSearchAnglesAcrossAtomGroup(Molecule *mp, IntGroup *atomgroup)
{
	if (mp == NULL)
		return NULL;
	return sMoleculeSearchAcrossAtomGroup(mp->nangles, mp->angles, 3, atomgroup, "angles");
}

IntGroup *
MoleculeSearchDihedralsAcrossAtomGroup(Molecule *mp, IntGroup *atomgroup)
{
	if (mp == NULL)
		return NULL;
	return sMoleculeSearchAcrossAtomGroup(mp->ndihedrals, mp->dihedrals, 4, atomgroup, "dihedrals");
}

IntGroup *
MoleculeSearchImpropersAcrossAtomGroup(Molecule *mp, IntGroup *atomgroup)
{
	if (mp == NULL)
		return NULL;
	return sMoleculeSearchAcrossAtomGroup(mp->nimpropers, mp->impropers, 4, atomgroup, "impropers");
}

/*  Subroutine for MoleculeGuessBonds. It can be also used independently, but make sure that *outNbonds/*outBonds 
    _correctly_ represents an array of two integers (as in mp->nbonds/mp->bonds).  */
/*  Find atoms within the given "distance" from the given position.  */
/*  If limit is negative, its absolute value denotes the threshold distance in angstrom; otherwise,
 the threshold distance is given by the sum of van der Waals radii times limit, and radius is
 the van der Waals radius of the atom at the given position. */
/*  Index is the atom index of the given atom; it is only used in returning the "bond" array
 to the caller. If index is negative, then (-index) is the real atom index, and
 only atoms with lower indices than (-index) are looked for.  */
int
MoleculeFindCloseAtoms(Molecule *mp, const Vector *vp, Double radius, Double limit, Int *outNbonds, Int **outBonds, Int index)
{
	Int n2, j, nlim, newbond[2];
	Double a2, alim;
	Vector dr, r2;
	if (index < 0) {
		nlim = index = -index;
	} else {
		nlim = mp->natoms;
	}
	for (j = 0; j < nlim; j++) {
		Atom *bp = ATOM_AT_INDEX(mp->atoms, j);
		if (index == j)
			continue;
		n2 = bp->atomicNumber;
		if (n2 >= 0 && n2 < gCountElementParameters)
			a2 = gElementParameters[n2].radius;
		else a2 = gElementParameters[6].radius;
		r2 = bp->r;
		VecSub(dr, *vp, r2);
		if (limit < 0)
			alim = -limit;
		else
			alim = limit * (radius + a2);
		if (VecLength2(dr) < alim * alim) {
			newbond[0] = index;
			newbond[1] = j;
			/*	MoleculeAddBonds(mp, 1, newbonds); */
			AssignArray(outBonds, outNbonds, sizeof(Int) * 2, *outNbonds, newbond);
		}
	}
	return 0;
}

/*  Guess the bonds from the coordinates  */
/*  If limit is negative, its absolute value denotes the threshold distance in angstrom; otherwise,
    the threshold distance is given by the sum of van der Waals radii times limit.  */
int
MoleculeGuessBonds(Molecule *mp, Double limit, Int *outNbonds, Int **outBonds)
{
	Int nbonds, *bonds, i, newbond[2];
	Atom *ap;
	nbonds = 0;
	bonds = NULL;
	if (limit == 0.0)
		limit = 1.2;
	for (i = 1, ap = ATOM_NEXT(mp->atoms); i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		Vector r = ap->r;
		Int an = ap->atomicNumber;
		Double rad;
		if (an >= 0 && an < gCountElementParameters)
			rad = gElementParameters[an].radius;
		else rad = gElementParameters[6].radius;
		MoleculeFindCloseAtoms(mp, &r, rad, limit, &nbonds, &bonds, -i);
	}
	if (nbonds > 0) {
		newbond[0] = kInvalidIndex;
		newbond[1] = 0;
		AssignArray(&bonds, &nbonds, sizeof(Int) * 2, nbonds, newbond);
		nbonds--;
	}
	if (outNbonds != NULL)
		*outNbonds = nbonds;
	if (outBonds != NULL)
		*outBonds = bonds;
	return 0;
}

/*  Rebuild the bond/angle/dihedral/improper tables from atom.connects[] information  */
int
MoleculeRebuildTablesFromConnects(Molecule *mp)
{
	int i, j, k, retval;
	Atom *ap;
	Int ibuf[6], *cp;
	
	__MoleculeLock(mp);

	/*  Find bonds   */
	if (mp->nbonds == 0) {
		for (i = 0; i < mp->natoms; i++) {
			ap = ATOM_AT_INDEX(mp->atoms, i);
			cp = AtomConnectData(&ap->connect);
			for (j = 0; j < ap->connect.count; j++) {
				k = cp[j];
				if (i >= k)
					continue;
				ibuf[0] = i;
				ibuf[1] = k;
				/*  MoleculeAddBonds() should not be used, because it assumes connects[] and
				    bonds are already in sync  */
				AssignArray(&mp->bonds, &mp->nbonds, sizeof(Int) * 2, mp->nbonds, ibuf);
			/*	retval = MoleculeAddBonds(mp, 1, ibuf);
				if (retval != 0)
					goto abort; */
			}
		}
	}
	
	/*  Find angles  */
	if (mp->nangles == 0) {
		for (i = 0; i < mp->natoms; i++) {
			ap = ATOM_AT_INDEX(mp->atoms, i);
			cp = AtomConnectData(&ap->connect);
			for (j = 0; j < ap->connect.count; j++) {
				for (k = j + 1; k < ap->connect.count; k++) {
					ibuf[0] = cp[j];
					ibuf[1] = i;
					ibuf[2] = cp[k];
					ibuf[3] = -1;
					retval = MoleculeAddAngles(mp, ibuf, NULL);
					if (retval < 0)
						goto abort;
				}
			}
		}
	}
	
	/*  Find dihedrals  */
	if (mp->ndihedrals == 0) {
		for (i = 0; i < mp->natoms; i++) {
			ap = ATOM_AT_INDEX(mp->atoms, i);
			cp = AtomConnectData(&ap->connect);
			for (j = 0; j < ap->connect.count; j++) {
				int jj, kk, mm, m;
				Atom *apjj;
				Int *cpjj;
				jj = cp[j];
				if (i >= jj)
					continue;
				apjj = ATOM_AT_INDEX(mp->atoms, jj);
				cpjj = AtomConnectData(&apjj->connect);
				for (k = 0; k < ap->connect.count; k++) {
					if (k == j)
						continue;
					kk = cp[k];
					for (m = 0; m < apjj->connect.count; m++) {
						mm = cpjj[m];
						if (mm == i || mm == kk)
							continue;
						ibuf[0] = kk;
						ibuf[1] = i;
						ibuf[2] = jj;
						ibuf[3] = mm;
						ibuf[4] = -1;
						retval = MoleculeAddDihedrals(mp, ibuf, NULL);
						if (retval < 0)
							goto abort;
					}
				}
			}
		}
	}
	
	/*  Find impropers  */
	if (mp->nimpropers == 0) {
		for (i = 0; i < mp->natoms; i++) {
			int i1, i2, i4, n1, n2, n4;
			ap = ATOM_AT_INDEX(mp->atoms, i);
			cp = AtomConnectData(&ap->connect);
			for (i1 = 0; i1 < ap->connect.count; i1++) {
				n1 = cp[i1];
				for (i2 = i1 + 1; i2 < ap->connect.count; i2++) {
					n2 = cp[i2];
					for (i4 = i2 + 1; i4 < ap->connect.count; i4++) {
						n4 = cp[i4];
						ibuf[0] = n1;
						ibuf[1] = n2;
						ibuf[2] = i;
						ibuf[3] = n4;
						ibuf[4] = -1;
						retval = MoleculeAddImpropers(mp, ibuf, NULL);
						if (retval < 0)
							goto abort;
					}
				}
			}
		}
	}

	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);
	return 0;

  abort:
	__MoleculeUnlock(mp);
	return retval;
}

int
MoleculeAreAtomsConnected(Molecule *mol, int idx1, int idx2)
{
	Atom *ap1 = ATOM_AT_INDEX(mol->atoms, idx1);
	if (AtomConnectHasEntry(&ap1->connect, idx2))
		return 1;
	else if (ap1->anchor != NULL && AtomConnectHasEntry(&(ap1->anchor->connect), idx2))
		return 2;
	else return 0;
}

#pragma mark ====== Atom names ======

/*  Look for the n1-th atom in resno-th residue (n1 is 0-based)  */
int
MoleculeLookupAtomInResidue(Molecule *mp, int n1, int resno)
{
	int i, j, lasti;
	Atom *ap;
	if (mp == NULL || mp->natoms == 0)
		return -1;
	lasti = -1;
	for (i = j = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (ap->resSeq == resno) {
			lasti = i;
			if (j++ == n1)
				return i;
		}
	}
	if (n1 == -1)
		return lasti; /* max */
	return -1;
}

int
MoleculeAnalyzeAtomName(const char *s, char *resName, int *resSeq, char *atomName)
{
    int n;
    char *p;
	n = strtol(s, &p, 0);
	if (p > s) {
		while (isspace(*p))
			p++;
		if (*p == 0) {
		  resName[0] = 0;
		  *resSeq = -1;
		  atomName[0] = 0;
		  return n;
		}
	}

	if ((p = strchr(s, ':')) != NULL) {
		/*  Residue is specified  */
		char *pp;
		if ((pp = strchr(s, '.')) != NULL && pp < p) {
			/*  Residue number is also specified  */
			char *ppp;
			n = pp - s;
			*resSeq = strtol(pp + 1, &ppp, 0);
			if (ppp == pp + 1)
				return -2;  /*  Bad format  */
			while (isspace(*ppp))
				ppp++;
			if (ppp != p)
				return -2;  /*  Bad format  */
		} else {
			*resSeq = -1;
			/*  Check whether the "residue name" is an integer  */
			n = strtol(s, &pp, 0);
			if (pp > s) {
				while (isspace(*pp))
					pp++;
				if (*pp == 0 || *pp == ':') {
					*resSeq = n;
					if (*resSeq < 0)
						return -2;  /*  Bad format  */
				}
			}
			if (*resSeq >= 0)
				n = 0;
			else
				n = p - s;
		}
		if (n >= sizeof(resName))
			n = sizeof(resName) - 1;
		strncpy(resName, s, n);
		resName[n] = 0;
		p++;
	} else {
		resName[0] = 0;
		*resSeq = -1;
		p = (char *)s;
	}
	strncpy(atomName, p, 4);
	atomName[4] = 0;
	return 0;
}

/*  Convert a string to atom index, where string = "((\w+\.)?(\d+):)?(\w+)" or an integer  */
int
MoleculeAtomIndexFromString(Molecule *mp, const char *s)
{
	char resName[6];
	int resSeq, n;
	char atomName[6];
	/*	char *p; */

	n = MoleculeAnalyzeAtomName(s, resName, &resSeq, atomName);
	if (atomName[0] == 0) {
	  if (n >= mp->natoms)
	    n = -1;  /* Out of range */
	  return n;
	}
	for (n = 0; n < mp->natoms; n++) {
		Atom *ap = ATOM_AT_INDEX(mp->atoms, n);
		if ((resName[0] == 0 || strncmp(resName, ap->resName, 4) == 0)
			&& (resSeq < 0 || ap->resSeq == resSeq)
			&& strncmp(atomName, ap->aname, 4) == 0) {
			return n;
		}
	}
	return -1;  /*  Not found  */
}

void
MoleculeGetAtomName(Molecule *mp, int index, char *buf, int bufsize)
{
	Atom *ap;
	int n;
	if (mp == NULL || index < 0 || index >= mp->natoms) {
		buf[0] = 0;
		return;
	}
	ap = mp->atoms + index;
	if (ap->resSeq != 0) {
		n = snprintf(buf, bufsize, "%s%d:", ap->resName, ap->resSeq);
		buf += n;
		bufsize -= n;
	}
	snprintf(buf, bufsize, "%.4s", ap->aname);
}

#pragma mark ====== Selection ======

static void
sMoleculeNotifyChangeSelection(Molecule *mp)
{
	/*  TODO: Finer control of notification types may be necessary  */
	MoleculeCallback_notifyModification(mp, 0);
}

void
MoleculeSetSelection(Molecule *mp, IntGroup *select)
{
	if (mp == NULL)
		return;
	if (select != NULL)
		IntGroupRetain(select);
	if (mp->selection != NULL)
		IntGroupRelease(mp->selection);
	mp->selection = select;
	sMoleculeNotifyChangeSelection(mp);
}

IntGroup *
MoleculeGetSelection(Molecule *mp)
{
	if (mp == NULL)
		return NULL;
	else return mp->selection;
}

void
MoleculeSelectAtom(Molecule *mp, int n1, int extending)
{
	if (mp->selection == NULL)
		mp->selection = IntGroupNew();
	if (!extending)
		IntGroupClear(mp->selection);
	IntGroupAdd(mp->selection, n1, 1);
	sMoleculeNotifyChangeSelection(mp);
}

void
MoleculeUnselectAtom(Molecule *mp, int n1)
{
	if (mp->selection != NULL)
		IntGroupRemove(mp->selection, n1, 1);
	sMoleculeNotifyChangeSelection(mp);
}

void
MoleculeToggleSelectionOfAtom(Molecule *mp, int n1)
{
	if (mp->selection == NULL)
		mp->selection = IntGroupNew();
	IntGroupReverse(mp->selection, n1, 1);
	sMoleculeNotifyChangeSelection(mp);
}

int
MoleculeIsAtomSelected(Molecule *mp, int n1)
{
	if (mp != NULL && mp->selection != NULL && IntGroupLookup(mp->selection, n1, NULL))
		return 1;
	else return 0;
}

int
MoleculeIsBondSelected(Molecule *mp, int n1, int n2)
{
	if (mp != NULL && MoleculeAreAtomsConnected(mp, n1, n2) && mp->selection != NULL && IntGroupLookup(mp->selection, n1, NULL) && IntGroupLookup(mp->selection, n2, NULL))
		return 1;
	else return 0;
}

IntGroup *
MoleculeModifySelectionByRemovingAtoms(Molecule *mp, IntGroup *selection, IntGroup *remove)
{
	int status;
	IntGroup *remain, *ig1, *ig2;
	ig1 = ig2 = NULL;
	remain = IntGroupNewFromIntGroup(remove);
	if (remain == NULL)
		status = -1;
	else
		status = IntGroupReverse(remain, 0, mp->natoms);
	if (status == 0) {
		ig1 = IntGroupNew();
		if (ig1 == NULL)
			status = -1;
		else
			status = IntGroupDifference(selection, remove, ig1);
	}
	if (status == 0) {
		ig2 = IntGroupNew();
		if (ig2 == NULL)
			status = -1;
		else
			status = IntGroupDeconvolute(ig1, remain, ig2);
	}
	if (remain != NULL)
		IntGroupRelease(remain);
	if (ig1 != NULL)
		IntGroupRelease(ig1);
	if (status == 0)
		return ig2;
	else {
		if (ig2 != NULL)
			IntGroupRelease(ig2);
		return NULL;
	}
}

#pragma mark ====== Atom Equivalence ======

struct sEqList {
	int i[2];
	struct sEqList *next;
	struct sEqList *link;
};

static struct sEqList *sListBase = NULL;
static struct sEqList *sListFree = NULL;

static struct sEqList *
sAllocEqList(void)
{
	struct sEqList *lp;
	if (sListFree != NULL) {
		lp = sListFree;
		sListFree = lp->next;
		lp->i[0] = lp->i[1] = 0;
		lp->next = NULL;
		return lp;
	}
	lp = (struct sEqList *)calloc(sizeof(struct sEqList), 1);
	lp->link = sListBase;
	sListBase = lp;
	return lp;
}

static void
sFreeEqList(struct sEqList *list)
{
	list->next = sListFree;
	sListFree = list;
}

static void
sDeallocateEqLists(void)
{
	struct sEqList *lp, *lp_link;
	for (lp = sListBase; lp != NULL; lp = lp_link) {
		lp_link = lp->link;
		free(lp);
	}
	sListBase = NULL;
	sListFree = NULL;
}

static int
sExistInEqList(int i, int idx, struct sEqList *list)
{
	while (list != NULL) {
		if (list->i[idx] == i)
			return 1;
		list = list->next;
	}
	return 0;
}

static struct sEqList *
sMoleculeCheckEquivalence(Molecule *mol, int i, int j, struct sEqList *list, int **db, IntGroup *ig)
{
	Atom *api, *apj;
	struct sEqList *list1, *list2;
	Int ii, jj, ni, nj, *cpi, *cpj;
	api = ATOM_AT_INDEX(mol->atoms, i);
	apj = ATOM_AT_INDEX(mol->atoms, j);
	if (api->atomicNumber != apj->atomicNumber)
		return NULL;
	list1 = sAllocEqList();
	if (list1 == NULL)
		return NULL;
	list1->i[0] = i;
	list1->i[1] = j;
	list1->next = list;
	if (i == j || (db[i] != NULL && db[i] == db[j]))
		return list1;
	cpi = AtomConnectData(&api->connect);
	cpj = AtomConnectData(&apj->connect);
	for (ni = 0; ni < api->connect.count; ni++) {
		ii = cpi[ni];
		if (ig != NULL && IntGroupLookupPoint(ig, ii) < 0)
			continue;
		if (sExistInEqList(ii, 0, list1))
			continue;
		list2 = NULL;
		for (nj = 0; nj < apj->connect.count; nj++) {
			jj = cpj[nj];
			if (ig != NULL && IntGroupLookupPoint(ig, jj) < 0)
				continue;
			if (sExistInEqList(jj, 1, list1))
				continue;
			list2 = sMoleculeCheckEquivalence(mol, ii, jj, list1, db, ig);
			if (list2 != NULL)
				break;
		}
		if (list2 == NULL) {
			sFreeEqList(list1);
			return NULL;    /*  No equivalent to ii  */
		}
		list1 = list2;      /*  ii is OK, try next  */
	}
	return list1;
}

int
sDBInclude(Int *ip, int i)
{
	int j;
	if (ip == NULL)
		return -1;
	for (j = ip[0] - 1; j >= 0; j--) {
		if (ip[j] == i)
			return j;
	}
	return -1;
}

Int *
MoleculeSearchEquivalentAtoms(Molecule *mol, IntGroup *ig)
{
	Int **db;  /*  List of equivalents for each atom  */
	Int *ip, *result;
	Atom *api, *apj, *apk;
	Int *cpi, *cpj, *ibuf, nibuf;
	int i, j, k, ii, jj, kk;
	if (mol == NULL || mol->natoms == 0)
		return NULL;
	db = (Int **)calloc(sizeof(Int *), mol->natoms);
	ibuf = NULL;
	nibuf = 0;

	/*  Find the equivalent univalent atoms  */
	for (i = 0, api = mol->atoms; i < mol->natoms; i++, api = ATOM_NEXT(api)) {
		if (api->connect.count < 2)
			continue;
		cpi = AtomConnectData(&api->connect);
		for (j = 0; j < api->connect.count; j++) {
			Int n;
			n = 0;
			jj = cpi[j];
			if (ig != NULL && IntGroupLookupPoint(ig, jj) < 0)
				continue;
			AssignArray(&ibuf, &nibuf, sizeof(Int), n, &jj);
			n++;
			apj = ATOM_AT_INDEX(mol->atoms, jj);
			if (apj->connect.count != 1 || db[jj] != NULL)
				continue;
			cpj = AtomConnectData(&apj->connect);
			for (k = j + 1; k < api->connect.count; k++) {
				kk = cpj[k];
				if (ig != NULL && IntGroupLookupPoint(ig, kk) < 0)
					continue;
				apk = ATOM_AT_INDEX(mol->atoms, kk);
				if (apk->connect.count != 1 || db[kk] != NULL)
					continue;
				if (apj->atomicNumber == apk->atomicNumber) {
					AssignArray(&ibuf, &nibuf, sizeof(Int), n, &kk);
					n++;
				}
			}
			if (n > 1) {
				ip = (Int *)calloc(sizeof(Int), n + 1);
				if (ip == NULL)
					return NULL;
				ip[0] = n;
				memmove(ip + 1, ibuf, sizeof(Int) * n);
				for (k = 0; k < n; k++)
					db[ip[k + 1]] = ip;
			}
		}
	}
	if (ibuf != NULL) {
		free(ibuf);
		ibuf = NULL;
	}
	
	/*  Try matching (i,j) pair  */
	for (i = 0, api = mol->atoms; i < mol->natoms; i++, api = ATOM_NEXT(api)) {
		if (ig != NULL && IntGroupLookupPoint(ig, i) < 0)
			continue;
		for (j = i + 1, apj = ATOM_AT_INDEX(mol->atoms, j); j < mol->natoms; j++, apj = ATOM_NEXT(apj)) {
			struct sEqList *list;
			if (ig != NULL && IntGroupLookupPoint(ig, j) < 0)
				continue;
			if (api->atomicNumber != apj->atomicNumber)
				continue;  /*  Different elements do not match  */
			if (db[i] != NULL && db[i] == db[j])
				continue;  /*  Already equivalent  */
			list = sMoleculeCheckEquivalence(mol, i, j, NULL, db, ig);
			if (list == NULL)
				continue;  /*  (i,j) do not match  */
			while (list != NULL) {
				ii = list->i[0];
				jj = list->i[1];
				if (ii != jj && (db[ii] == NULL || db[ii] != db[jj])) {
					/*  Merge db[ii] and db[jj]  */
					k = (db[ii] == NULL ? 1 : db[ii][0]) + (db[jj] == NULL ? 1 : db[jj][0]);
					ip = (Int *)calloc(sizeof(Int), k + 1);
					if (ip == NULL)
						return NULL;  /*  Out of memory  */
					if (db[ii] == NULL) {
						ip[1] = ii;
						k = 2;
					} else {
						memmove(ip + 1, db[ii] + 1, db[ii][0] * sizeof(Int));
						k = db[ii][0] + 1;
					}
					if (db[jj] == NULL) {
						ip[k++] = jj;
					} else {
						memmove(ip + k, db[jj] + 1, db[jj][0] * sizeof(Int));
						k += db[jj][0];
					}
					ip[0] = k - 1;
					/*  Free old ones  */
					if (db[ii] != NULL)
						free(db[ii]);
					if (db[jj] != NULL)
						free(db[jj]);
					for (k = 0; k < ip[0]; k++)
						db[ip[k + 1]] = ip;
					if (0) {
						/*  For debug  */
						printf("(%d,%d) matched: ", ii, jj);
						for (k = 0; k < ip[0]; k++) {
							printf("%c%d", (k == 0 ? '[' : ','), ip[k + 1]);
						}
						printf("]\n");
					}
				}
				list = list->next;
			}
		}
	}
	
	/*  Record the equivalent atoms with the lowest index for each atom  */
	result = (Int *)calloc(sizeof(Int), mol->natoms);
	for (i = 0; i < mol->natoms; i++)
		result[i] = -1;
	for (i = 0; i < mol->natoms; i++) {
		if (result[i] >= 0 || (ip = db[i]) == NULL)
			continue;
		k = mol->natoms;
		for (j = 0; j < ip[0]; j++) {
			kk = ip[j + 1];
			if (kk < k)
				k = kk;
		}
		for (j = 0; j < ip[0]; j++) {
			result[ip[j + 1]] = k;
			db[ip[j + 1]] = NULL;
		}
		free(ip);
	}
	sDeallocateEqLists();
	return result;
}

#pragma mark ====== Symmetry expansion ======

int
MoleculeGetTransformForSymop(Molecule *mp, Symop symop, Transform *tf, int is_cartesian)
{
	Transform t;
	if (mp == NULL || mp->cell == NULL)
		return -1;
	if (symop.sym >= mp->nsyms && symop.sym != 0)
		return -2;
	memmove(*tf, SYMMETRY_AT_INDEX(mp->syms, symop.sym), sizeof(Transform));
	(*tf)[9] += symop.dx;
	(*tf)[10] += symop.dy;
	(*tf)[11] += symop.dz;
	if (is_cartesian) {
		TransformMul(t, *tf, mp->cell->rtr);
		TransformMul(*tf, mp->cell->tr, t);
	}
	return 0;
}

int
MoleculeGetSymopForTransform(Molecule *mp, const Transform tf, Symop *symop, int is_cartesian)
{
	Transform t;
	int i, j, n[3];
	if (mp == NULL || mp->cell == NULL)
		return -1;
	if (is_cartesian) {
		TransformMul(t, tf, mp->cell->tr);
		TransformMul(t, mp->cell->rtr, t);
	} else {
		memmove(t, tf, sizeof(Transform));
	}
	for (i = 0; i < mp->nsyms || i == 0; i++) {
		Transform *tp = &(SYMMETRY_AT_INDEX(mp->syms, i));
		for (j = 0; j < 9; j++) {
			if (fabs((*tp)[j] - t[j]) > 1e-4)
				break;
		}
		if (j == 9) {
			for (j = 9; j < 12; j++) {
				double f1 = t[j] - (*tp)[j];
				double f2 = floor(f1 + 0.5);
				if (fabs(f1 - f2) > 1e-4)
					break;
				n[j - 9] = f2;
			}
			if (j == 12) {
				/*  Found  */
				symop->sym = i;
				symop->dx = n[0];
				symop->dy = n[1];
				symop->dz = n[2];
				symop->alive = (SYMOP_ALIVE((*symop)) != 0);
				return 0;
			}
		}
	}
	return -3;  /*  Not found  */
}

int
MoleculeTransformBySymop(Molecule *mp, const Vector *vpin, Vector *vpout, Symop symop)
{
	if (mp == NULL)
		return 1;
	if (symop.sym >= mp->nsyms && symop.sym != 0)
		return 2;
	if (mp->cell != NULL /* && !mp->is_xtal_coord */) {
		TransformVec(vpout, mp->cell->rtr, vpin);
		TransformVec(vpout, SYMMETRY_AT_INDEX(mp->syms, symop.sym), vpout);
		vpout->x += symop.dx;
		vpout->y += symop.dy;
		vpout->z += symop.dz;
		TransformVec(vpout, mp->cell->tr, vpout);
	} else {
		TransformVec(vpout, SYMMETRY_AT_INDEX(mp->syms, symop.sym), vpin);
		vpout->x += symop.dx;
		vpout->y += symop.dy;
		vpout->z += symop.dz;
	}
	return 0;
}

/*  Add expanded atoms. Returns the number of newly created atoms.
	If indices is non-NULL, it should be an array of Int with at least 
	IntGroupGetCount(group) entries, and on return it contains the
    indices of the expanded atoms (may be existing atoms if the expanded
    atoms are already present)
    If allowOverlap is non-zero, then the new atom is created even when the
    coordinates coincide with the some other atom (special position) of the
    same element; otherwise, such atom will not be created and the existing
    atom is returned in indices[].  */
int
MoleculeAddExpandedAtoms(Molecule *mp, Symop symop, IntGroup *group, Int *indices, Int allowOverlap)
{
	int i, n, n0, n1, n2, base, count, *table;
	Atom *ap;
	IntGroupIterator iter;
	Transform tr, t1;
	Symop symop1;
	Atom *ap2;
	Vector nr, dr;
	
	if (mp == NULL || mp->natoms == 0 || group == NULL || (count = IntGroupGetCount(group)) == 0)
		return -1;
	if (symop.sym != 0 && symop.sym >= mp->nsyms)
		return -2;

	/*  Create atoms, with avoiding duplicates  */
	n0 = n1 = mp->natoms;
	table = (int *)malloc(sizeof(int) * n0);
	if (table == NULL)
		return -3;
	for (i = 0; i < n0; i++)
		table[i] = -1;
	IntGroupIteratorInit(group, &iter);
	MoleculeGetTransformForSymop(mp, symop, &tr, 0);
	__MoleculeLock(mp);
	for (i = 0; i < count; i++) {
		n = IntGroupIteratorNext(&iter);
		ap = ATOM_AT_INDEX(mp->atoms, n);
		if (SYMOP_ALIVE(ap->symop)) {
			/*  Calculate the cumulative symop  */
			Transform tr2;
			MoleculeGetTransformForSymop(mp, ap->symop, &t1, 0);
			TransformMul(tr2, tr, t1);
			if (MoleculeGetSymopForTransform(mp, tr2, &symop1, 0) != 0) {
				if (indices != NULL)
					indices[i] = -1;
				continue;  /*  Skip this atom  */
			}
			base = ap->symbase;
		} else {
			symop1 = symop;
			base = n;
		}

		/*  Calculate the expande position  */
		MoleculeTransformBySymop(mp, &(ap->r), &nr, symop);
		
		/*  Is this expansion already present?  */
		for (n2 = 0, ap2 = mp->atoms; n2 < n0; n2++, ap2 = ATOM_NEXT(ap2)) {
			/*  Symmetry operation and the base atom are the same  */
			if (ap2->symbase == base && SYMOP_EQUAL(symop1, ap2->symop))
				break;
			/*  Atomic number and the position are the same  */
			if (ap2->atomicNumber == ap->atomicNumber && allowOverlap == 0) {
				VecSub(dr, ap2->r, nr);
				if (VecLength2(dr) < 1e-6)
					break;
			}
		}
		if (n2 < n0) {
			/*  If yes, then skip it  */
			if (indices != NULL)
				indices[i] = n2;
			continue;
		} else {
			/*  Create a new atom  */
			Atom newAtom;
			AtomDuplicate(&newAtom, ap);
			MoleculeCreateAnAtom(mp, &newAtom, -1);
			AtomClean(&newAtom);
			ap2 = ATOM_AT_INDEX(mp->atoms, mp->natoms - 1);
			ap2->r = nr;
			ap2->symbase = base;
			ap2->symop = symop1;
			ap2->symop.alive = (symop1.dx != 0 || symop1.dy != 0 || symop1.dz != 0 || symop1.sym != 0);
			table[n] = n1;  /*  The index of the new atom  */
			MoleculeSetAnisoBySymop(mp, n1);  /*  Recalculate anisotropic parameters according to symop  */
			if (indices != NULL)
				indices[i] = n1;
			n1++;
		}
	}
	IntGroupIteratorRelease(&iter);

	/*  Create bonds  */
	for (i = n0; i < n1; i++) {
		Int b[2], j;
		ap = ATOM_AT_INDEX(mp->atoms, i);
		if (SYMOP_ALIVE(ap->symop) && MoleculeGetTransformForSymop(mp, ap->symop, &tr, 1) == 0) {
			/*  For each connected atom, look for the transformed atom  */
			Int *cp;
			ap2 = ATOM_AT_INDEX(mp->atoms, ap->symbase);
			cp = AtomConnectData(&ap2->connect);
			n2 = ap2->connect.count;
			for (n = 0; n < n2; n++) {
				Atom *apn = ATOM_AT_INDEX(mp->atoms, cp[n]);
				nr = apn->r;
				TransformVec(&nr, tr, &nr);
				/*  Look for the bonded atom transformed by ap->symop  */
				for (j = 0, ap2 = mp->atoms; j < mp->natoms; j++, ap2 = ATOM_NEXT(ap2)) {
					if (ap2->symbase == cp[n] && SYMOP_EQUAL(ap->symop, ap2->symop))
						break;
					VecSub(dr, nr, ap2->r);
					if (ap2->atomicNumber == apn->atomicNumber && VecLength2(dr) < 1e-6)
						break;
				}
				if (j < mp->natoms) {
					/*  Bond i-j is created  */
					b[0] = i;
					b[1] = j;
					if (MoleculeLookupBond(mp, b[0], b[1]) < 0)
						MoleculeAddBonds(mp, 1, b, NULL, 1);
				}
			}
		}
	}
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);
	free(table);
	return n1 - n0;  /*  The number of added atoms  */
}

/*  Recalculate the coordinates of symmetry expanded atoms.
    (Also recalculate the positions of pi-anchor atoms)
	Returns the number of affected atoms.
    If group is non-NULL, only the expanded atoms whose base atoms are in the
    given group are considered.
	If groupout and vpout are non-NULL, the indices of the affected atoms
	and the original positions are returned (for undo operation).
	The pointers returned in *groupout and *vpout must be released and 
	free()'ed by the caller  */
int
MoleculeAmendBySymmetry(Molecule *mp, IntGroup *group, IntGroup **groupout, Vector **vpout)
{
	int i, count;
	Atom *ap, *bp;
	Vector nr, dr;
	IntGroup *ig = NULL;
	Vector *vp = NULL;
	
	if (mp == NULL || mp->natoms == 0)
		return 0;

	__MoleculeLock(mp);
	count = 0;
	if (mp->nsyms != 0) {
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (!SYMOP_ALIVE(ap->symop))
				continue;
			if (group != NULL && IntGroupLookup(group, ap->symbase, NULL) == 0)
				continue;
			bp = ATOM_AT_INDEX(mp->atoms, ap->symbase);
			MoleculeTransformBySymop(mp, &(bp->r), &nr, ap->symop);
			VecSub(dr, nr, ap->r);
			if (VecLength2(dr) < 1e-20)
				continue;
			if (groupout != NULL) {
				if (ig == NULL) {
					ig = IntGroupNew();
					vp = (Vector *)calloc(sizeof(Vector), mp->natoms);
				}
				vp[count] = ap->r;
				IntGroupAdd(ig, i, 1);
			}
			ap->r = nr;
			count++;
		}
	}
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		Int *ip, j, n;
		if (ap->anchor == NULL)
			continue;
		if (group != NULL) {
			if (IntGroupLookup(group, i, NULL) == 0) {
				n = ap->anchor->connect.count;
				ip = AtomConnectData(&(ap->anchor->connect));
				for (j = 0; j < n; j++) {
					if (IntGroupLookup(group, ip[j], NULL) != 0)
						break;
				}
				if (j == n)
					continue;  /*  This pi-anchor should not be modified  */
			}
		}
		nr = ap->r;
		MoleculeCalculatePiAnchorPosition(mp, i);
		VecSub(dr, nr, ap->r);
		if (VecLength2(dr) < 1e-20) {
			ap->r = nr;  /*  No change  */
			continue;
		}
		if (groupout != NULL) {
			if (ig == NULL) {
				ig = IntGroupNew();
				vp = (Vector *)calloc(sizeof(Vector), mp->natoms);
			}
			vp[count] = nr;
			IntGroupAdd(ig, i, 1);
		}
		count++;
	}
	mp->needsMDCopyCoordinates = 1;
	__MoleculeUnlock(mp);

	if (count > 0) {
		if (groupout != NULL && vpout != NULL) {
			*groupout = ig;
			*vpout = (Vector *)realloc(vp, sizeof(Vector) * count);
		} else {
			IntGroupRelease(ig);
			free(vp);
		}
	} else {
		if (groupout != NULL && vpout != NULL) {
			*groupout = NULL;
			*vpout = NULL;
		}
	}
	return count;
}

#pragma mark ====== Show/hide atoms ======

static void
sMoleculeNotifyChangeAppearance(Molecule *mp)
{
	/*  TODO: Finer control of notification types may be necessary  */
	MoleculeCallback_notifyModification(mp, 0);
}


static void
sMoleculeUnselectHiddenAtoms(Molecule *mp)
{
	int i;
	if (mp == NULL || mp->selection == NULL)
		return;
	for (i = 0; i < mp->natoms; i++) {
		Atom *ap = ATOM_AT_INDEX(mp->atoms, i);
		if ((ap->exflags & kAtomHiddenFlag) && IntGroupLookupPoint(mp->selection, i) >= 0)
			IntGroupRemove(mp->selection, i, 1);
	}
	sMoleculeNotifyChangeAppearance(mp);
}

int
MoleculeShowAllAtoms(Molecule *mp)
{
	int i;
	if (mp == NULL)
		return 0;
	for (i = 0; i < mp->natoms; i++) {
		Atom *ap = ATOM_AT_INDEX(mp->atoms, i);
		ap->exflags &= ~kAtomHiddenFlag;
	}
	sMoleculeNotifyChangeAppearance(mp);
	return 1;
}

int
MoleculeShowReverse(Molecule *mp)
{
	int i;
	if (mp == NULL)
		return 0;
	for (i = 0; i < mp->natoms; i++) {
		Atom *ap = ATOM_AT_INDEX(mp->atoms, i);
		ap->exflags ^= kAtomHiddenFlag;
	}
	sMoleculeUnselectHiddenAtoms(mp);
	sMoleculeNotifyChangeAppearance(mp);
	return 1;
}

int
MoleculeHideAtoms(Molecule *mp, IntGroup *ig)
{
	int i;
	if (mp == NULL || ig == NULL)
		return 0;
	for (i = 0; i < mp->natoms; i++) {
		Atom *ap = ATOM_AT_INDEX(mp->atoms, i);
		if (ap->exflags & kAtomHiddenFlag)
			continue;  /*  Already hidden  */
		if (IntGroupLookupPoint(ig, i) >= 0)
			ap->exflags |= kAtomHiddenFlag;
	}
	sMoleculeUnselectHiddenAtoms(mp);
	sMoleculeNotifyChangeAppearance(mp);
	return 1;
}

#pragma mark ====== Reversible Editing ======

/*
static void
sMoleculeNotifyModification(Molecule *mp)
{
	**  TODO: Finer control of notification types may be necessary  **
	MoleculeCallback_notifyModification(mp, 0);
}
*/

/*  Insert new[0,1,2,...] to old[n0,n1,n2,...], where {n0,n1,n2,...} is the points in IntGroup  */
int
sInsertElementsToArrayAtPositions(void *objs, int nobjs, const void *newobjs, int nnewobjs, size_t size, IntGroup *where)
{
	int n1, n2, n3, i;
	if (where == NULL) {
		/*  Append the new objects at the end  */
		memmove((char *)objs + size * nobjs, (char *)newobjs, size * nnewobjs);
		return 0;
	}
	n1 = IntGroupGetCount(where);  /*  Position to get new object  */
	n2 = nobjs;                    /*  Position to get old object  */
	n3 = n1 + n2;                  /*  Position to place new/old object  */
	for (i = IntGroupGetIntervalCount(where) - 1; i >= 0; i--) {
		int start = IntGroupGetStartPoint(where, i);
		int end = IntGroupGetEndPoint(where, i);
		if (end < n3) {
			/*  old[end-(n3-n2)..n2-1] is moved to old[end..n3-1]  */
			memmove((char *)objs + size * end, (char *)objs + size * (end - (n3 - n2)), size * (n3 - end));
			n2 = end - (n3 - n2);
			n3 = end;
		}
		/*  new[n1-(end-start)..n1-1] is moved to old[n3-(end-start)..n3-1]  */
		memmove((char *)objs + size * (n3 - (end - start)), (char *)newobjs + size * (n1 - (end - start)), size * (end - start));
		n3 -= end - start;
		n1 -= end - start;
	}
	return 0;
}

/*  Move objs[n0,n1,n2,...] to clip[0,1,2,...], where {n0,n1,n2,...} is the points in IntGroup  */
int
sRemoveElementsFromArrayAtPositions(void *objs, int nobjs, void *clip, size_t size, IntGroup *where)
{
	int n1, n2, n3, start, end, i;
	if (where == NULL || IntGroupGetCount(where) == 0)
		return 0;  /*  No operation  */
	if (objs == NULL || nobjs == 0)
		return 1;  /*  Bad argument  */
	n1 = 0;  /*  Position to move remaining elements to */
	n2 = 0;  /*  Position to move remaining elements from  */
	n3 = 0;  /*  Position to move removed elements to  */
	for (i = 0; (start = IntGroupGetStartPoint(where, i)) >= 0; i++) {
		end = IntGroupGetEndPoint(where, i);
		if (n2 < start) {
			/*  Move (start - n2) elements from objs[n2] to objs[n1]  */
			if (n1 < n2)
				memmove((char *)objs + size * n1, (char *)objs + size * n2, size * (start - n2));
			n1 += start - n2;
			n2 = start;
		}
		/*  Move (end - start) elements from objs[n2] to clip[n3]  */
		if (clip != NULL)
			memmove((char *)clip + size * n3, (char *)objs + size * n2, size * (end - start));
		n3 += (end - start);
		n2 += (end - start);
	}
	/*  Move (nobjs - n2) elements from objs[n2] to objs[n1]  */
	if (nobjs > n2)
		memmove((char *)objs + size * n1, (char *)objs + size * n2, size * (nobjs - n2));
	return 0;
}

/*  Copy objs[n0,n1,n2,...] to clip[0,1,2,...], where {n0,n1,n2,...} is the points in IntGroup  */
int
sCopyElementsFromArrayAtPositions(void *objs, int nobjs, void *clip, size_t size, IntGroup *where)
{
	int n1, start, end, i;
	if (objs == NULL || where == NULL)
		return 1;  /*  Bad argument  */
	n1 = 0;  /*  Position to move removed elements to  */
	for (i = 0; (start = IntGroupGetStartPoint(where, i)) >= 0; i++) {
		end = IntGroupGetEndPoint(where, i);
		/*  Copy (end - start) elements from objs[start] to clip[n1]  */
		if (clip != NULL)
			memmove((char *)clip + size * n1, (char *)objs + size * start, size * (end - start));
		n1 += (end - start);
	}
	return 0;
}

/*  Create a new atom with no bonding information. ap must _not_ be inside the given molecule
   (Use AtomDuplicate() first) */
int
MoleculeCreateAnAtom(Molecule *mp, const Atom *ap, int pos)
{
    Atom *ap1, *api;
	int i;
	if (mp == NULL || ap == NULL || mp->noModifyTopology)
		return -1;
	__MoleculeLock(mp);
	if (pos < 0 || pos >= mp->natoms)
		pos = mp->natoms;
	ap1 = AssignArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, mp->natoms, NULL);
	if (ap1 == NULL)
		goto error;  /*  Out of memory  */
	ap1 = ATOM_AT_INDEX(mp->atoms, pos);
	if (pos < mp->natoms - 1) {
		memmove(ATOM_AT_INDEX(mp->atoms, pos + 1), ATOM_AT_INDEX(mp->atoms, pos), gSizeOfAtomRecord * (mp->natoms - 1 - pos));
	}
	if (AtomDuplicate(ap1, ap) == NULL) {
		/*  Cannot duplicate: restore the original state  */
		memmove(ATOM_AT_INDEX(mp->atoms, pos), ATOM_AT_INDEX(mp->atoms, pos + 1), gSizeOfAtomRecord * (mp->natoms - 1 - pos));
		mp->natoms--;
		goto error;
	}
	ap1->connect.count = 0;
	if (ap1->resSeq >= mp->nresidues)
		AssignArray(&mp->residues, &mp->nresidues, 4, ap1->resSeq, ap1->resName);
	if (ap1->resName[0] == 0)
	  strncpy(ap1->resName, mp->residues[ap1->resSeq], 4);
	if (ap1->segName[0] == 0)
	  strncpy(ap1->segName, "MAIN", 4);
	if (pos < mp->natoms - 1) {
		/*  Renumber the connect table, bonds, angles, etc. */
		for (i = 0, api = ATOM_AT_INDEX(mp->atoms, i); i < mp->natoms; i++, api = ATOM_NEXT(api)) {
			int j;
			Int *cp;
			cp = AtomConnectData(&api->connect);
			for (j = 0; j < api->connect.count; j++) {
				if (cp[j] >= pos)
					cp[j]++;
			}
			if (api->anchor != NULL) {
				cp = AtomConnectData(&api->anchor->connect);
				for (j = 0; j < api->anchor->connect.count; j++) {
					if (cp[j] >= pos)
						cp[j]++;
				}
			}
		}
		for (i = 0; i < mp->nbonds * 2; i++) {
			if (mp->bonds[i] >= pos)
				mp->bonds[i]++;
		}
		for (i = 0; i < mp->nangles * 3; i++) {
			if (mp->angles[i] >= pos)
				mp->angles[i]++;
		}
		for (i = 0; i < mp->ndihedrals * 4; i++) {
			if (mp->dihedrals[i] >= pos)
				mp->dihedrals[i]++;
		}
		for (i = 0; i < mp->nimpropers * 4; i++) {
			if (mp->impropers[i] >= pos)
				mp->impropers[i]++;
		}
	}
	mp->nframes = -1;  /*  Should be recalculated later  */
	MoleculeIncrementModifyCount(mp);
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);
	return pos;
error:
	__MoleculeUnlock(mp);
	return -1;
}

#if defined(DEBUG)

static int s_error_count;

static int
s_fprintf(FILE *fp, const char *fmt, ...)
{
	va_list va;
	va_start(va, fmt);
	s_error_count++;
	return vfprintf(fp, fmt, va);
}

int
MoleculeCheckSanity(Molecule *mol)
{
	const char *fail = "Sanity check failure";
	Int i, j, *ip, c[4];
	Atom *ap;
	s_error_count = 0;
	for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (ap->resSeq >= mol->nresidues)
			s_fprintf(stderr, "%s: atom %d residue %d but nresidues %d\n", fail, i, ap->resSeq, mol->nresidues);
		if (ap->type != 0 && ap->type < kAtomTypeMinimum)
			s_fprintf(stderr, "%s: atom %d atom type %d less than minimum\n", fail, i, ap->type);
		if (ap->atomicNumber < 0 || ap->atomicNumber > 113)
			s_fprintf(stderr, "%s: atom %d atomic number %d\n", fail, i, ap->atomicNumber);
		ip = AtomConnectData(&ap->connect);
		for (j = 0; j < ap->connect.count; j++) {
			if (ip[j] < 0 || ip[j] >= mol->natoms)
				s_fprintf(stderr, "%s: atom %d connect[%d] = %d out of range\n", fail, i, j, ip[j]);
			if (AtomConnectHasEntry(&(ATOM_AT_INDEX(mol->atoms, ip[j])->connect), i) == 0)
				s_fprintf(stderr, "%s: atom %d has connect %d but atom %d has no connect %d\n", fail, i, ip[j], ip[j], i);
		}
	}
	for (i = 0, ip = mol->bonds; i < mol->nbonds; i++, ip += 2) {
		if (ip[0] < 0 || ip[0] >= mol->natoms || ip[1] < 0 || ip[1] >= mol->natoms)
			s_fprintf(stderr, "%s: bond %d %d-%d out of range\n", fail, i, ip[0], ip[1]);
		if (AtomConnectHasEntry(&(ATOM_AT_INDEX(mol->atoms, ip[0])->connect), ip[1]) == 0)
			s_fprintf(stderr, "%s: bond %d %d-%d but atom %d has no connect %d\n", fail, i, ip[0], ip[1], ip[0], ip[1]);
	}
	for (i = 0, ip = mol->angles; i < mol->nangles; i++, ip += 3) {
		if (ip[0] < 0 || ip[0] >= mol->natoms || ip[1] < 0 || ip[1] >= mol->natoms || ip[2] < 0 || ip[2] >= mol->natoms)
			s_fprintf(stderr, "%s: angle %d %d-%d-%d out of range\n", fail, i, ip[0], ip[1], ip[2]);
		c[0] = MoleculeAreAtomsConnected(mol, ip[1], ip[0]);
		if (c[0] == 0)
			s_fprintf(stderr, "%s: angle %d %d-%d-%d but atom %d has no connect %d\n", fail, i, ip[0], ip[1], ip[2], ip[1], ip[0]);
		c[1] = MoleculeAreAtomsConnected(mol, ip[1], ip[2]);
		if (c[1] == 0)
			s_fprintf(stderr, "%s: angle %d %d-%d-%d but atom %d has no connect %d\n", fail, i, ip[0], ip[1], ip[2], ip[1], ip[2]);
		if (c[0] == 2 && c[1] == 2)
			s_fprintf(stderr, "%s: angle %d %d-%d-%d but bonds %d-%d and %d-%d are both virtual\n", fail, i, ip[0], ip[1], ip[2], ip[1], ip[0], ip[1], ip[2]);
	}
	for (i = 0, ip = mol->dihedrals; i < mol->ndihedrals; i++, ip += 4) {
		if (ip[0] < 0 || ip[0] >= mol->natoms || ip[1] < 0 || ip[1] >= mol->natoms || ip[2] < 0 || ip[2] >= mol->natoms || ip[3] < 0 || ip[3] >= mol->natoms)
			s_fprintf(stderr, "%s: dihedral %d %d-%d-%d%d out of range\n", fail, i, ip[0], ip[1], ip[2], ip[3]);
		c[0] = MoleculeAreAtomsConnected(mol, ip[1], ip[0]);
		c[1] = MoleculeAreAtomsConnected(mol, ip[1], ip[2]);
		c[2] = MoleculeAreAtomsConnected(mol, ip[2], ip[3]);
		if (c[0] == 0)
			s_fprintf(stderr, "%s: dihedral %d %d-%d-%d-%d but atom %d has no connect %d\n", fail, i, ip[0], ip[1], ip[2], ip[3], ip[1], ip[0]);
		if (c[1] == 0)
			s_fprintf(stderr, "%s: dihedral %d %d-%d-%d-%d but atom %d has no connect %d\n", fail, i, ip[0], ip[1], ip[2], ip[3], ip[1], ip[2]);
		if (c[2] == 0)
			s_fprintf(stderr, "%s: dihedral %d %d-%d-%d-%d but atom %d has no connect %d\n", fail, i, ip[0], ip[1], ip[2], ip[3], ip[2], ip[3]);
	}
	for (i = 0, ip = mol->impropers; i < mol->nimpropers; i++, ip += 4) {
		if (ip[0] < 0 || ip[0] >= mol->natoms || ip[1] < 0 || ip[1] >= mol->natoms || ip[2] < 0 || ip[2] >= mol->natoms || ip[3] < 0 || ip[3] >= mol->natoms)
			s_fprintf(stderr, "%s: improper %d %d-%d-%d%d out of range\n", fail, i, ip[0], ip[1], ip[2], ip[3]);
		c[0] = MoleculeAreAtomsConnected(mol, ip[2], ip[0]);
		c[1] = MoleculeAreAtomsConnected(mol, ip[2], ip[1]);
		c[2] = MoleculeAreAtomsConnected(mol, ip[2], ip[3]);
		if (c[0] == 0)
			s_fprintf(stderr, "%s: improper %d %d-%d-%d-%d but atom %d has no connect %d\n", fail, i, ip[0], ip[1], ip[2], ip[3], ip[2], ip[0]);
		if (c[1] == 0)
			s_fprintf(stderr, "%s: improper %d %d-%d-%d-%d but atom %d has no connect %d\n", fail, i, ip[0], ip[1], ip[2], ip[3], ip[2], ip[1]);
		if (c[2] == 0)
			s_fprintf(stderr, "%s: improper %d %d-%d-%d-%d but atom %d has no connect %d\n", fail, i, ip[0], ip[1], ip[2], ip[3], ip[2], ip[3]);
	}
	return s_error_count;
}
#endif

/*  Merge two molecules. We use this procedure for all add-atom operations.  */
/*  resSeqOffset is an offset to add to the (non-zero) residue numbers in src. */
/*  If nactions and actions are non-NULL, then the corresponding undo actions are created and returned. */
/*  If forUndo is non-zero, then only the atoms are inserted; other information should be inserted
    separately by other undo actions.  */
int
MoleculeMerge(Molecule *dst, Molecule *src, IntGroup *where, Int resSeqOffset, Int *nactions, MolAction ***actions, Int forUndo)
{
	Int nsrc, ndst;
	Int i, j, n1, n2, n3, n4, *cp;
	Int *new2old, *old2new;
	IntGroup *ig;
	Atom *ap;
	MolAction *act;
	
	if (dst == NULL || src == NULL || src->natoms == 0 || (where != NULL && IntGroupGetIntervalCount(where) == 0))
		return 0;  /*  Do nothing  */

	if (dst->noModifyTopology)
		return 1;  /*  Prohibited operation  */

	if (where != NULL && IntGroupGetCount(where) != src->natoms)
		return 1;  /*  Bad parameter  */

	if (nactions != NULL)
		*nactions = 0;
	if (actions != NULL)
		*actions = NULL;
	act = NULL;

	__MoleculeLock(dst);

	nsrc = src->natoms;
	ndst = dst->natoms;
	if (resSeqOffset < 0)
		resSeqOffset = 0;

	/*  Atom index table. For "old" index, 0..ndst-1 are for atoms in dst,
	    and ndst..ndst+nsrc-1 are for atoms in src.  */ 
	new2old = (Int *)calloc(sizeof(Int), (ndst + nsrc) * 2);
	if (new2old == NULL)
		goto panic;
	old2new = new2old + ndst + nsrc;
	n1 = 0;  /*  dst index  */
	n2 = 0;  /*  src index  */
	n3 = 0;  /*  "merged" index  */
	i = 0;
	while (n1 < ndst || n2 < nsrc) {
		if (where == NULL || (n4 = IntGroupGetStartPoint(where, i)) < 0)
			n4 = ndst - n1;
		else n4 -= n3;
		/*  n4 elements from dst[n1] will go to merged[n3]  */
		for (j = 0; j < n4; j++) {
			old2new[n1 + j] = n3 + j;
			new2old[n3 + j] = n1 + j;
		}
		n3 += n4;
		n1 += n4;
		if (where == NULL || (n4 = IntGroupGetInterval(where, i)) < 0)
			n4 = nsrc - n2;
		/*  n4 elements from src[n2] will go to merged[n3]  */
		for (j = 0; j < n4; j++) {
			old2new[ndst + n2 + j] = n3 + j;
			new2old[n3 + j] = ndst + n2 + j;
		}
		n3 += n4;
		n2 += n4;
		i++;
	}

	/*  Expand the destination array  */
	if (AssignArray(&(dst->atoms), &(dst->natoms), gSizeOfAtomRecord, ndst + nsrc - 1, NULL) == NULL)
		goto panic;

	/*  Move the atoms  */
	if (where == NULL) {
		/*  Duplicate atoms to the end of the destination array  */
		for (i = 0; i < nsrc; i++) {
			ap = ATOM_AT_INDEX(dst->atoms, ndst + i);
			if (AtomDuplicate(ap, ATOM_AT_INDEX(src->atoms, i)) == NULL)
				goto panic;
			if (forUndo)  /*  For undo action, all bonds come from another undo action, so connection info are cleared */
				AtomConnectResize(&ap->connect, 0);
		}
	} else {
		/*  Duplicate to a temporary storage and then insert  */
		Atom *tempatoms = (Atom *)malloc(gSizeOfAtomRecord * nsrc);
		if (tempatoms == NULL)
			goto panic;
		for (i = 0; i < nsrc; i++) {
			ap = ATOM_AT_INDEX(tempatoms, i);
			if (AtomDuplicate(ap, ATOM_AT_INDEX(src->atoms, i)) == NULL)
				goto panic;
			if (forUndo)  /*  See above  */
				AtomConnectResize(&ap->connect, 0);				
		}
		if (sInsertElementsToArrayAtPositions(dst->atoms, ndst, tempatoms, nsrc, gSizeOfAtomRecord, where) != 0)
			goto panic;
		free(tempatoms);
	}
	dst->natoms = ndst + nsrc;

	/*  Renumber the atom indices in connect[] and symbase, and modify the residue numbers  */
	for (i = 0, ap = dst->atoms; i < dst->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (new2old[i] < ndst) {
			/*  This atom is from dst  */
			n1 = 0;
		} else {
			/*  This atom is from src  */
			n1 = ndst;  /*  Offset to the internal number  */
			if (ap->resSeq != 0)
				ap->resSeq += resSeqOffset;  /*  Modify residue number  */
		}
		cp = AtomConnectData(&ap->connect);
		for (j = 0; j < ap->connect.count; j++)
			cp[j] = old2new[cp[j] + n1];
		if (SYMOP_ALIVE(ap->symop))
			ap->symbase = old2new[ap->symbase + n1];
		if (ap->anchor != NULL) {
			cp = AtomConnectData(&ap->anchor->connect);
			for (j = 0; j < ap->anchor->connect.count; j++)
				cp[j] = old2new[cp[j] + n1];
		}
	}
	
	/*  Move the bonds, angles, dihedrals, impropers  */
	for (i = 0; i < 4; i++) {
		Int *nitems, *nitems_src;
		Int **items, **items_src;
		Int nsize;  /*  Number of Ints in one element  */
		switch (i) {
			case 0:
				nitems = &dst->nbonds; items = &dst->bonds; nsize = 2; break;
			case 1:
				nitems = &dst->nangles; items = &dst->angles; nsize = 3; break;
			case 2:
				nitems = &dst->ndihedrals; items = &dst->dihedrals; nsize = 4; break;
			case 3:
				nitems = &dst->nimpropers; items = &dst->impropers; nsize = 4; break;
		}
		nitems_src = (Int *)((char *)src + ((char *)nitems - (char *)dst));
		items_src = (Int **)((char *)src + ((char *)items - (char *)dst));
		if (forUndo) {
			/*  During undo, no bonds etc. are copied from src; they will be taken care later
			    by undo actions  */
			n1 = *nitems;
			n2 = 0;
		} else {
			/*  Keep the old number of entries in dst, because it is updated by AssignArray()  */
			n1 = *nitems;
			/*  Also keep the old number of entries in src, in case src and dst point the same molecule  */
			n2 = *nitems_src;
			/*  Expand the array  */
			if (AssignArray(items, nitems, sizeof(Int) * nsize, *nitems + *nitems_src - 1, NULL) == NULL)
				goto panic;
			/*  Copy the items  */
			memmove(*items + n1 * nsize, *items_src, sizeof(Int) * nsize * n2);
			if (i == 0) {
				/*  Copy the bond order info if present */
				Int nn1 = dst->nbondOrders;
				if (dst->bondOrders != NULL || src->bondOrders != NULL) {
					if (AssignArray(&dst->bondOrders, &dst->nbondOrders, sizeof(Double), dst->nbonds - 1, NULL) == NULL)
						goto panic;
					memset(dst->bondOrders + nn1, 0, sizeof(Double) * (dst->nbonds - nn1));
					if (src->bondOrders != NULL)
						memmove(dst->bondOrders + n1, src->bondOrders, sizeof(Double) * n2);
				}
			}
		}
		/*  Renumber  */
		for (j = 0; j < n1 * nsize; j++)
			(*items)[j] = old2new[(*items)[j]];
		for (j = n1 * nsize; j < (n1 + n2) * nsize; j++)
			(*items)[j] = old2new[(*items)[j] + ndst];
		if (forUndo == 0 && actions != NULL) {
			ig = IntGroupNewWithPoints(n1, n2, -1);
			switch (i) {
				case 0: act = MolActionNew(gMolActionDeleteBonds, ig); break;
				case 1: act = MolActionNew(gMolActionDeleteAngles, ig); break;
				case 2: act = MolActionNew(gMolActionDeleteDihedrals, ig); break;
				case 3: act = MolActionNew(gMolActionDeleteImpropers, ig); break;
			}
			IntGroupRelease(ig);
			AssignArray(actions, nactions, sizeof(MolAction *), *nactions, &act);
			act = NULL;
		}
	}
	
	/*  Renumber existing parameters  */
	if (dst->par != NULL) {
		int type;
		for (type = kFirstParType; type <= kLastParType; type++) {
			UnionPar *up1;
			n1 = ParameterGetCountForType(dst->par, type);
			for (i = 0; i < n1; i++) {
				up1 = ParameterGetUnionParFromTypeAndIndex(dst->par, type, i);
				ParameterRenumberAtoms(type, up1, ndst, old2new);
			}
		}
	}

	/*  Merge parameters from src  */
	if (src->par != NULL && forUndo == 0) {
		UnionPar *up1, *up2;
		int type;
		if (dst->par == NULL)
			dst->par = ParameterNew();
		else {
			/*  Renumber existing parameters  */
			for (type = kFirstParType; type <= kLastParType; type++) {
				n1 = ParameterGetCountForType(dst->par, type);
				for (i = 0; i < n1; i++) {
					up1 = ParameterGetUnionParFromTypeAndIndex(dst->par, type, i);
					ParameterRenumberAtoms(type, up1, ndst, old2new);
				}
			}
		}
		ig = IntGroupNew();
		for (type = kFirstParType; type <= kLastParType; type++) {
			n1 = ParameterGetCountForType(src->par, type);
			n2 = ParameterGetCountForType(dst->par, type);
			if (n1 == 0)
				continue;
			/*  Determine which parameter should be copied from src to dst  */
			for (i = 0; i < n1; i++) {
				UInt types[4];
				up1 = ParameterGetUnionParFromTypeAndIndex(src->par, type, i);
				n3 = ParameterGetAtomTypes(type, up1, types);
				for (j = 0; j < n3; j++) {
					/*  If it includes explicit atom index, then it should be copied  */
					if (types[j] < kAtomTypeMinimum) {
						IntGroupAdd(ig, i, 1);
						break;
					}
				}
				if (j == n3) {
					for (j = 0; j < n2; j++) {
						up2 = ParameterGetUnionParFromTypeAndIndex(dst->par, type, j);
						if (ParameterCompare(up1, up2, type))
							break;
					}
					if (j >= n2)
						/*  This is an unknown parameter; should be copied  */
						IntGroupAdd(ig, i, 1);
				}
			}
			n1 = IntGroupGetCount(ig);
			if (n1 == 0)
				continue;
			up1 = (UnionPar *)calloc(sizeof(UnionPar), n1);
			if (up1 == NULL)
				goto panic;
			/*  Copy parameters and renumber indices if necessary  */
			for (i = j = 0; i < n1; i++) {
				up2 = ParameterGetUnionParFromTypeAndIndex(src->par, type, IntGroupGetNthPoint(ig, i));
				if (up2 == NULL)
					continue;
				up1[j] = *up2;
				ParameterRenumberAtoms(type, up1 + j, nsrc, old2new + ndst);
				j++;
			}
			/*  Merge parameters  */
			IntGroupClear(ig);
			IntGroupAdd(ig, n2, j);
			if (ParameterInsert(dst->par, type, up1, ig) < j)
				goto panic;
			if (actions != NULL) {
				act = MolActionNew(gMolActionDeleteParameters, type, ig);
				AssignArray(actions, nactions, sizeof(MolAction *), *nactions, &act);
				act = NULL;
			}
			IntGroupClear(ig);
			free(up1);
		}
		IntGroupRelease(ig);
	}
	
	/*  Copy the residues if necessary  */
	/*  src[1..src->nresidues-1] should become dst[1+resSeqOffset..src->nresidues+resSeqOffset-1];
	    However, 1+resSeqOffset should not overwrite the existing residue in dst;
		i.e. if 1+resSeqOffset is less than dst->nresidues, copy should start from src[dst->nresidues-resSeqOffset] instead of src[1].  */
	if (forUndo == 0) {
		n1 = dst->nresidues;
		if (1 + resSeqOffset < n1) {
			n2 = n1;
		} else n2 = 1 + resSeqOffset; /* n2 is the start index of residues from src[] */
		if (src->nresidues > 1 && n1 < src->nresidues + resSeqOffset) {
			if (AssignArray(&dst->residues, &dst->nresidues, sizeof(dst->residues[0]), src->nresidues + resSeqOffset - 1, NULL) == NULL)
				goto panic;
			memmove(dst->residues + n2, src->residues + n2 - resSeqOffset, sizeof(dst->residues[0]) * (src->nresidues - (n2 - resSeqOffset)));
			if (nactions != NULL) {
				act = MolActionNew(gMolActionChangeNumberOfResidues, n1);
				AssignArray(actions, nactions, sizeof(MolAction *), *nactions, &act);
				act = NULL;
			}
		}
	}

	MoleculeCleanUpResidueTable(dst);
	
	free(new2old);
	dst->nframes = -1;  /*  Should be recalculated later  */

	MoleculeIncrementModifyCount(dst);
	dst->needsMDRebuild = 1;
	__MoleculeUnlock(dst);
	return 0;

  panic:
	__MoleculeUnlock(dst);
    Panic("Low memory while adding atoms");
	return 1;  /*  Not reached  */
}

/*  Unmerge the molecule. If necessary, the undo actions are stored in nactions/actions array.
    (The nactions/actions array must be initialized by the caller)  */
static int
sMoleculeUnmergeSub(Molecule *src, Molecule **dstp, IntGroup *where, int resSeqOffset, int moveFlag, Int *nactions, MolAction ***actions, Int forUndo)
{
	Int nsrc, ndst, nsrcnew;
	Int i, j, n1, n2, n3, n4, *cp;
	Int *new2old, *old2new;
	IntGroup *move_g, *del_g, *remain_g, *dst_par_g, *remove_par_g;
	Molecule *dst;
	Atom *ap, *dst_ap;
	UnionPar *up;
	MolAction *act;

	if (src == NULL || src->natoms == 0 || where == NULL || IntGroupGetIntervalCount(where) == 0) {
		/*  Do nothing  */
		if (dstp != NULL)
			*dstp = NULL;
		return 0;
	}
	
	if (src->noModifyTopology && moveFlag)
		return 1;  /*  Prohibit editing  */

	if ((ndst = IntGroupGetCount(where)) > src->natoms)
		return 1;  /*  Bad parameter  */

	__MoleculeLock(src);
	
	act = NULL;
	
	nsrc = src->natoms;
	nsrcnew = nsrc - ndst;
	if (resSeqOffset < 0)
		resSeqOffset = 0;

	/*  Atom index table. For "new" index, 0..nsrcnew-1 are for atoms remaining in src,
	    and nsrcnew..nsrc-1 are for atoms moved into dst.  */ 
	new2old = (Int *)calloc(sizeof(Int), nsrc * 2);
	if (new2old == NULL)
		goto panic;
	old2new = new2old + nsrc;
	n1 = 0;  /*  src index  */
	n2 = 0;  /*  dst index  */
	n3 = 0;  /*  src index after "unmerge"  */
	i = 0;
	while (n1 < nsrc || n2 < ndst) {
		if ((n4 = IntGroupGetStartPoint(where, i)) < 0)
			n4 = nsrc - n1;
		else n4 -= n1;
		/*  n4 elements from src[n1] will go to unmerged[n3]  */
		for (j = 0; j < n4; j++) {
			old2new[n1 + j] = n3 + j;
			new2old[n3 + j] = n1 + j;
		}
		n3 += n4;
		n1 += n4;
		if ((n4 = IntGroupGetInterval(where, i)) < 0)
			n4 = nsrc - n1;
		/*  n4 elements from src[n1] will go to dst[n2]  */
		for (j = 0; j < n4; j++) {
			old2new[n1 + j] = nsrcnew + n2 + j;
			new2old[nsrcnew + n2 + j] = n1 + j;
		}
		n1 += n4;
		n2 += n4;
		i++;
	}

	/*  Atoms to remain in the source group  */
	if (moveFlag) {
		remain_g = IntGroupNewWithPoints(0, nsrc, -1);
		IntGroupRemoveIntGroup(remain_g, where);
	} else remain_g = NULL;
	
	/*  Find parameters to be moved to the dst (dst_par_g), and to be removed from the src (remove_par_g) */
	if (src->par != NULL) {
		dst_par_g = IntGroupNew();
		if (moveFlag)
			remove_par_g = IntGroupNew();
		else remove_par_g = NULL;
		for (n1 = kFirstParType; n1 <= kLastParType; n1++) {
			n2 = ParameterGetCountForType(src->par, n1);
			if (n2 == 0)
				continue;
			for (i = 0; i < n2; i++) {
				up = ParameterGetUnionParFromTypeAndIndex(src->par, n1, i);
				if (ParameterIsRelevantToAtomGroup(n1, up, src->atoms, where)) {
					/*  This parameter is to be copied to dst  */
					IntGroupAdd(dst_par_g, i + (n1 - kFirstParType) * kParameterIndexOffset, 1);
				}
				if (moveFlag && !ParameterIsRelevantToAtomGroup(n1, up, src->atoms, remain_g)) {
					/*  This parameter is to be removed  */
					IntGroupAdd(remove_par_g, i + (n1 - kFirstParType) * kParameterIndexOffset, 1);
				}
			}
		}
	} else dst_par_g = remove_par_g = NULL;
	
	/*  Pi anchors should be modified if the anchor and its component atoms become separated between
	    src anc dst  */
	if (moveFlag) {
		Int ibufsize, *ibuf, flag_i, flag_j;
		ibufsize = 8;
		ibuf = (Int *)malloc(sizeof(Int) * ibufsize);
		for (i = 0, ap = src->atoms; i < src->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ap->anchor == NULL)
				continue;
			flag_i = (old2new[i] < nsrcnew);
			cp = AtomConnectData(&ap->anchor->connect);
			for (j = n1 = 0; j < ap->anchor->connect.count; j++) {
				flag_j = (old2new[cp[j]] < nsrcnew);
				if (flag_i == flag_j) {
					if (n1 >= ibufsize) {
						ibufsize += 8;
						ibuf = (Int *)realloc(ibuf, sizeof(Int) * ibufsize);
					}
					ibuf[n1++] = cp[j];
				}
			}
			if (n1 < j) {
				/*  Need to modify the pi anchor list  */
				if (n1 <= 1)
					n1 = 0;
				MolActionCreateAndPerform(src, SCRIPT_ACTION("isI"), "set_atom_attr", i, "anchor_list", n1, ibuf);
			}
		}
	}
	
	/*  Make a new molecule  */
	if (dstp != NULL) {
		dst = MoleculeNew();
		if (dst == NULL)
			goto panic;
		/*  Expand the destination array  */
		if (AssignArray(&(dst->atoms), &(dst->natoms), gSizeOfAtomRecord, ndst - 1, NULL) == NULL)
			goto panic;
		dst_ap = dst->atoms;
	} else {
		dst = NULL;
		dst_ap = (Atom *)calloc(sizeof(Atom), ndst);
		if (dst_ap == NULL)
			goto panic;
	}
	
	/*  Move the atoms  */
	if (moveFlag) {
		if (sRemoveElementsFromArrayAtPositions(src->atoms, src->natoms, dst_ap, gSizeOfAtomRecord, where) != 0)
			goto panic;
		src->natoms = nsrcnew;
		if (dst == NULL) {
			/*  The atom record must be deallocated correctly  */
			for (i = 0; i < ndst; i++)
				AtomClean(ATOM_AT_INDEX(dst_ap, i));
		}
	} else {
		if (dst != NULL) {
			for (i = 0; (n1 = IntGroupGetNthPoint(where, i)) >= 0; i++)
				AtomDuplicate(ATOM_AT_INDEX(dst_ap, i), ATOM_AT_INDEX(src->atoms, n1));
		}
	}
	
	if (dst == NULL) {
		/*  The dummy destination array is no longer needed  */
		free(dst_ap);
		dst_ap = NULL;
	}
	
	/*  Renumber the atom indices in connect[] (src) */
	if (moveFlag) {
		for (i = 0, ap = src->atoms; i < src->natoms; i++, ap = ATOM_NEXT(ap)) {
			cp = AtomConnectData(&ap->connect);
			for (j = n1 = 0; j < ap->connect.count; j++) {
				n2 = old2new[cp[j]];
				if (n2 < nsrcnew)
					cp[n1++] = n2;
			}
			AtomConnectResize(&ap->connect, n1);
			if (ap->anchor != NULL) {
				cp = AtomConnectData(&ap->anchor->connect);
				for (j = n1 = 0; j < ap->anchor->connect.count; j++) {
					n2 = old2new[cp[j]];
					if (n2 < nsrcnew)
						cp[n1++] = n2;
				}
				if (n1 != ap->anchor->connect.count) {
					/*  This should not happen!!  */
					AtomConnectResize(&ap->anchor->connect, n1);
					fprintf(stderr, "Internal error in sMoleculeUnmergeSub (line %d)\n", __LINE__);
					if (n1 == 0) {
						free(ap->anchor->coeffs);
						free(ap->anchor);
						ap->anchor = NULL;
					}
				}
			}
		}
	}
	
	/*  Renumber the atom indices in connect[] (dst)  */
	if (dst != NULL) {
		for (i = 0, ap = dst->atoms; i < dst->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ap->resSeq != 0 && ap->resSeq - resSeqOffset >= 0)
				ap->resSeq -= resSeqOffset;
			else ap->resSeq = 0;
			cp = AtomConnectData(&ap->connect);
			for (j = n1 = 0; j < ap->connect.count; j++) {
				n2 = old2new[cp[j]] - nsrcnew;
				if (n2 >= 0)
					cp[n1++] = n2;
			}
			AtomConnectResize(&ap->connect, n1);
			if (ap->anchor != NULL) {
				cp = AtomConnectData(&ap->anchor->connect);
				for (j = n1 = 0; j < ap->anchor->connect.count; j++) {
					n2 = old2new[cp[j]] - nsrcnew;
					if (n2 >= 0)
						cp[n1++] = n2;
				}
				if (n1 != ap->anchor->connect.count) {
					/*  This can happen, and the anchor info is silently modified  */
					if (n1 <= 1) {
						AtomConnectResize(&ap->anchor->connect, 0);
						free(ap->anchor->coeffs);
						free(ap->anchor);
						ap->anchor = NULL;
					} else {
						Double d;
						AtomConnectResize(&ap->anchor->connect, n1);
						d = 0.0;
						for (j = 0; j < n1; j++)
							d += ap->anchor->coeffs[j];
						for (j = 0; j < n1; j++)
							ap->anchor->coeffs[j] /= d;
						MoleculeCalculatePiAnchorPosition(dst, i);
					}
				}
			}
		}
	}

	/*  Separate the bonds, angles, dihedrals, impropers  */
	/*  TODO: Improper torsions should also be copied!  */
	move_g = IntGroupNew();
	if (move_g == NULL)
		goto panic;
	for (i = 3; i >= 0; i--) {
		Int *nitems, *nitems_dst;
		Int **items, **items_dst;
		Int nsize;  /*  Number of Ints in one element  */
		unsigned char *counts;
		del_g = IntGroupNew();
		switch (i) {
			case 0:
				nitems = &src->nbonds; items = &src->bonds; nsize = 2; break;
			case 1:
				nitems = &src->nangles; items = &src->angles; nsize = 3; break;
			case 2:
				nitems = &src->ndihedrals; items = &src->dihedrals; nsize = 4; break;
			case 3:
				nitems = &src->nimpropers; items = &src->impropers; nsize = 4; break;
			default:
				nitems = NULL; items = NULL; nsize = 0; break;  /*  Not reached  */
		}
		if (dst != NULL) {
			nitems_dst = (Int *)((char *)dst + ((char *)nitems - (char *)src));
			items_dst = (Int **)((char *)dst + ((char *)items - (char *)src));
		} else {
			nitems_dst = NULL;
			items_dst = NULL;
		}
		counts = (unsigned char *)calloc(1, *nitems);
		/*  Find the entries that should be moved to dst  */
		n2 = 0;
		for (j = 0; j < *nitems * nsize; j++) {
			n1 = old2new[(*items)[j]];
			if (n1 >= nsrcnew)
				counts[j / nsize]++; /* Count the atom belonging to dst */ 
		}
		for (j = n2 = n3 = 0; j < *nitems; j++) {
			if (counts[j] > 0) {
				/*  Remove from src  */
				n2++;
				if (IntGroupAdd(del_g, j, 1) != 0)
					goto panic;
				if (counts[j] == nsize) {
					/*  Move to dst  */
					n3++;
					if (IntGroupAdd(move_g, j, 1) != 0)
						goto panic;
				}
			}
		}
		if (n2 > 0) {
			/*  Expand the destination array  */
			if (items_dst != NULL && n3 > 0) {
				if (AssignArray(items_dst, nitems_dst, sizeof(Int) * nsize, n3 - 1, NULL) == NULL)
					goto panic;
				if (sCopyElementsFromArrayAtPositions(*items, *nitems, *items_dst, sizeof(Int) * nsize, move_g) != 0)
					goto panic;
				if (i == 0 && src->bondOrders != NULL) {
					if (AssignArray(&dst->bondOrders, &dst->nbondOrders, sizeof(Double), n3 - 1, NULL) == NULL)
						goto panic;
					if (sCopyElementsFromArrayAtPositions(src->bondOrders, src->nbondOrders, dst->bondOrders, sizeof(Double), move_g) != 0)
						goto panic;
				}
			}
			/*  Remove from src  */
			if (moveFlag && forUndo == 0) {
				if (nactions != NULL) {
					Int k, *ip;
					Double *dp;
					ip = (Int *)malloc(sizeof(Int) * nsize * n2);
					for (j = 0; (k = IntGroupGetNthPoint(del_g, j)) >= 0; j++)
						memmove(ip + j * nsize, *items + k * nsize, sizeof(Int) * nsize);
					if (i == 0 && src->bondOrders != NULL) {
						dp = (Double *)malloc(sizeof(Double) * n2);
						for (j = 0; (k = IntGroupGetNthPoint(del_g, j)) >= 0; j++)
							dp[j] = src->bondOrders[k];
					} else dp = NULL;
					switch (i) {
						case 0:
							act = MolActionNew(gMolActionAddBondsForUndo, n2 * nsize, ip, del_g); break;
						case 1:
							act = MolActionNew(gMolActionAddAngles, n2 * nsize, ip, del_g); break;
						case 2:
							act = MolActionNew(gMolActionAddDihedrals, n2 * nsize, ip, del_g); break;
						case 3:
							act = MolActionNew(gMolActionAddImpropers, n2 * nsize, ip, del_g); break;
					}
					if (act != NULL) {
						AssignArray(actions, nactions, sizeof(MolAction *), *nactions, &act);
						act = NULL;
					}
					free(ip);
					if (dp != NULL) {
						act = MolActionNew(gMolActionAssignBondOrders, n2, dp, del_g);
						AssignArray(actions, nactions, sizeof(MolAction *), *nactions, &act);
						act = NULL;
						free(dp);
					}
				}
				if (sRemoveElementsFromArrayAtPositions(*items, *nitems, NULL, sizeof(Int) * nsize, del_g) != 0)
					goto panic;
				(*nitems) -= n2;
			}
		}
		/*  Renumber the entries  */
		if (moveFlag) {
			for (j = 0; j < *nitems * nsize; j++) {
				(*items)[j] = old2new[(*items)[j]];
			}
		}
		if (items_dst != NULL) {
			for (j = 0; j < *nitems_dst * nsize; j++) {
				(*items_dst)[j] = old2new[(*items_dst)[j]] - nsrcnew;
			}
		}
		free(counts);
		IntGroupClear(move_g);
		IntGroupRelease(del_g);
	}
	IntGroupRelease(move_g);
	
	/*  Copy the residues  */
	if (dst != NULL) {
		/*  src[i] will become dst[i - resSeqOffset] (src->nresidues > i >= 1 + resSeqOffset)  */
		n1 = src->nresidues - resSeqOffset;  /*  This will be dst->nresidues (if >0)  */
		if (AssignArray(&dst->residues, &dst->nresidues, sizeof(dst->residues[0]), (n1 > 0 ? n1 - 1: 0), NULL) == NULL)
			goto panic;
		if (n1 > 1) {
			memmove(dst->residues + 1, src->residues + resSeqOffset + 1, sizeof(dst->residues[0]) * (n1 - 1));
		}
	}

	/*  Copy the parameters to dst */
	if (dst != NULL && dst_par_g != NULL && (n2 = IntGroupGetCount(dst_par_g)) > 0) {
		IntGroup *dst_new_g = IntGroupNew();
		Int dst_par_count[kLastParType - kFirstParType + 1];
		if (dst_new_g == NULL)
			goto panic;
		for (i = 0; i <= kLastParType - kFirstParType; i++)
			dst_par_count[i] = 0;
		up = (UnionPar *)calloc(sizeof(UnionPar), n2);
		if (up == NULL)
			goto panic;
		if (ParameterCopy(src->par, kFirstParType, up, dst_par_g) < n2)
			goto panic;
		/*  Renumber the explicit atom indices  */
		for (i = 0; i < nsrc; i++)
			old2new[i] -= nsrcnew;  /*  new indices for atoms in dst; otherwise negative numbers  */
		for (i = 0; i < n2; i++) {
			/*  Renumber the indices, and count the number of parameters for each type  */
			n1 = kFirstParType + IntGroupGetNthPoint(dst_par_g, i) / kParameterIndexOffset;
			dst_par_count[n1 - kFirstParType]++;
			ParameterRenumberAtoms(n1, up + i, nsrc, old2new);
		}
		for (i = 0; i < nsrc; i++)
			old2new[i] += nsrcnew;
		if (dst->par == NULL)
			dst->par = ParameterNew();
		for (i = 0; i <= kLastParType - kFirstParType; i++) {
			if (dst_par_count[i] > 0)
				IntGroupAdd(dst_new_g, i * kParameterIndexOffset, dst_par_count[i]);
		}
		if (ParameterInsert(dst->par, kFirstParType, up, dst_new_g) < n2)
			goto panic;
		free(up);
		IntGroupRelease(dst_new_g);
	}
	IntGroupRelease(dst_par_g);

	/*  Remove the unused parameter. Note: the parameters that are in remove_par_g and not in 
	    dst_par_g will disappear. To support undo, these parameters should be taken care separately.  */
	if (forUndo == 0 && remove_par_g != NULL && (n2 = IntGroupGetCount(remove_par_g)) > 0) {
		UnionPar *up = (UnionPar *)malloc(sizeof(UnionPar) * n2);
		ParameterDelete(src->par, kFirstParType, up, remove_par_g);
		if (nactions != NULL) {
			act = MolActionNew(gMolActionAddParameters, kFirstParType, remove_par_g, n2, up);
			AssignArray(actions, nactions, sizeof(MolAction *), *nactions, &act);
			act = NULL;
		}
		free(up);
	}
	IntGroupRelease(remove_par_g);
	
	/*  Renumber the parameter records remaining in the src  */
	if (moveFlag) {
		for (n1 = kFirstParType; n1 <= kLastParType; n1++) {
			n2 = ParameterGetCountForType(src->par, n1);
			for (i = 0; i < n2; i++) {
				up = ParameterGetUnionParFromTypeAndIndex(src->par, n1, i);
				ParameterRenumberAtoms(n1, up, nsrc, old2new);
			}
		}
	}

	/*  Clean up  */
	IntGroupRelease(remain_g);
	MoleculeCleanUpResidueTable(src);
	if (dst != NULL)
		MoleculeCleanUpResidueTable(dst);
	free(new2old);

	src->nframes = -1;  /*  Should be recalculated later  */
	if (dst != NULL)
		dst->nframes = -1;  /*  Should be recalculated later  */

	
	if (dstp != NULL)
		*dstp = dst;

	MoleculeIncrementModifyCount(src);
	src->needsMDRebuild = 1;
	__MoleculeUnlock(src);
	
	return 0;

  panic:
	__MoleculeUnlock(src);
/*    Panic("Low memory while removing atoms"); */
	return -1;
}

/*  Separate molecule into two parts. The atoms specified by 'where' are moved
    from src to a new molecule, which is returned as *dstp. Dstp can be NULL, 
	in which case the moved atoms are discarded.  */
int
MoleculeUnmerge(Molecule *src, Molecule **dstp, IntGroup *where, int resSeqOffset, Int *nactions, MolAction ***actions, Int forUndo)
{
	return sMoleculeUnmergeSub(src, dstp, where, resSeqOffset, 1, nactions, actions, forUndo);
}

/*  Extract atoms from a given molecule into two parts. The atoms specified by 
	'where' are copied from src to a new molecule, which is returned as *dstp.
    If dummyFlag is non-zero, then the atoms that are not included in the group 
	but are connected to any atoms in the group are converted to "dummy" atoms 
	(i.e. with element "Du" and names beginning with an underscore) and included 
	in the new molecule object.  */
int
MoleculeExtract(Molecule *src, Molecule **dstp, IntGroup *where, int dummyFlag)
{
	int retval;

	/*  Extract the fragment  */
	retval = sMoleculeUnmergeSub(src, dstp, where, 0, 0, NULL, NULL, 0);
	if (retval != 0)
		return retval;

	if (dummyFlag) {

		/*  Search bonds crossing the molecule border  */
		IntGroup *ig = MoleculeSearchBondsAcrossAtomGroup(src, where);
		if (ig != NULL) {
			IntGroupIterator iter;
			Int i, idx;
			idx = 1;
			IntGroupIteratorInit(ig, &iter);
			while ((i = IntGroupIteratorNext(&iter)) >= 0) {
				/*  The atoms at the border  */
				Int n1, n2, nn[3];
				Atom a, *ap;
				n1 = src->bonds[i*2];
				n2 = src->bonds[i*2+1];
				if ((nn[0] = IntGroupLookupPoint(where, n1)) < 0) {
					int w = n1;
					n1 = n2;
					n2 = w;
					if ((nn[0] = IntGroupLookupPoint(where, n1)) < 0)
						continue;  /*  Actually this is an internal error  */
				}
				/*  n1 is in *where, n2 is not; nn[0] is the index of atom n1 in the new molecule  */
				/*  Create a new dummy atom with the same segment/residue info with n1
				    and the same position as n2  */
				ap = ATOM_AT_INDEX(src->atoms, n1);
				memset(&a, 0, gSizeOfAtomRecord);
				a.segSeq = ap->segSeq;
				memmove(a.segName, ap->segName, 4);
				a.resSeq = ap->resSeq;
				memmove(a.resName, ap->resName, 4);
				ElementToString(0, a.element);  /*  "Du"  */
				snprintf(a.aname, 4, "_%d", idx++);
				a.r = ATOM_AT_INDEX(src->atoms, n2)->r;
				/*  Add the dummy atom to the new molecule; nn[1] is the index
				    of the new dummy atom in the new molecule  */
				nn[1] = MoleculeCreateAnAtom(*dstp, &a, -1);
				/*  Connect nn1 and nn2  */
				nn[2] = kInvalidIndex;
				MoleculeAddBonds(*dstp, 1, nn, NULL, 1);
			}
			IntGroupIteratorRelease(&iter);
			IntGroupRelease(ig);
		}
	}
	
	return 0;
}

int
MoleculeAddBonds(Molecule *mp, Int nbonds, const Int *bonds, IntGroup *where, Int autoGenerate)
{
	Int nangles, ndihedrals;
	Int *angles, *dihedrals;
	Int i, j, k, kk, n1, n2, cn1, cn2;
	Int *cp1, *cp2;
	Int temp[4];
	Atom *ap1, *ap2, *ap3;
	
	if (mp == NULL || bonds == NULL || nbonds <= 0)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */

	/*  Note: Duplicates and validity are not checked (the caller must do that)  */

	__MoleculeLock(mp);

	n1 = mp->nbonds;
	if (AssignArray(&(mp->bonds), &(mp->nbonds), sizeof(Int) * 2, n1 + nbonds - 1, NULL) == NULL
		|| sInsertElementsToArrayAtPositions(mp->bonds, n1, bonds, nbonds, sizeof(Int) * 2, where) != 0) {
		__MoleculeUnlock(mp);
		return -4;  /*  Out of memory  */
	}
	if (mp->bondOrders != NULL) {
		/*  Expand the bond order info (all new entries are zero)  */
		Double *dp = (Double *)calloc(sizeof(Double), nbonds);
		if (dp == NULL)
			return -4;
		if (AssignArray(&(mp->bondOrders), &(mp->nbondOrders), sizeof(Double), n1 + nbonds - 1, NULL) == NULL
			|| sInsertElementsToArrayAtPositions(mp->bondOrders, n1, dp, nbonds, sizeof(Double), where) != 0) {
			__MoleculeUnlock(mp);
			free(dp);
			return -4;
		}
		free(dp);
	}
	
	angles = dihedrals = NULL;
	nangles = ndihedrals = 0;
	
	/*  Add connects[], and angles/dihedrals (if autoGenerate is true)  */
	for (i = 0; i < nbonds; i++) {
		
		/*  One entry at time  */
		/*  (Otherwise, duplicate entries of angles and dihedrals result)  */
		n1 = bonds[i * 2];
		n2 = bonds[i * 2 + 1];
		
		ap1 = ATOM_AT_INDEX(mp->atoms, n1);
		AtomConnectInsertEntry(&ap1->connect, -1, n2);
		ap2 = ATOM_AT_INDEX(mp->atoms, n2);
		AtomConnectInsertEntry(&ap2->connect, -1, n1);
	
		/*  Add angles and dihedrals  */
		if (autoGenerate) {
			AtomConnect *ac1, *ac2;
			if (ap1->anchor == NULL || ap2->anchor == NULL) {
				/*  N1-N2-{XY} or N2-N1-{XY} angles (X: connected atom, Y: constitute atom of pi-anchor)  */
				for (j = 0; j < 4; j++) {
					switch (j) {
						case 0: temp[0] = n1; temp[1] = n2; ac1 = &ap2->connect; break;  /* N1-N2-X */
						case 1: if (ap2->anchor == NULL) continue; else ac1 = &ap2->anchor->connect; break; /* N1-N2-Y */
						case 2: temp[0] = n2; temp[1] = n1; ac1 = &ap1->connect; break;  /* N2-N1-X */
						case 3: if (ap1->anchor == NULL) continue; else ac1 = &ap1->anchor->connect; break; /* N2-N1-Y */
					}
					cp1 = AtomConnectData(ac1);
					cn1 = ac1->count;
					for (k = 0; k < cn1; k++) {
						temp[2] = cp1[k];
						if (temp[2] == temp[0])
							continue;
						ap3 = ATOM_AT_INDEX(mp->atoms, temp[2]);
						if (ap3->anchor != NULL) {
							/*  Avoid X-anchor-anchor angle (anchor-X-anchor is allowed)  */
							if ((j < 2 && ap2->anchor != NULL) || (j >= 2 && ap1->anchor != NULL))
								continue;
						}
						if (AssignArray(&angles, &nangles, sizeof(Int) * 3, nangles, temp) == NULL)
							goto panic;
						/*  Dihedrals N1-N2-X-{XY} or N2-N1-X-{XY}  */
						if (j == 1 || j == 3)
							continue;
						cp2 = AtomConnectData(&ap3->connect);
						for (kk = 0; kk < ap3->connect.count; kk++) {
							temp[3] = cp2[kk];
							if (temp[3] == temp[0] || temp[3] == temp[1])
								continue;
							if (AssignArray(&dihedrals, &ndihedrals, sizeof(Int) * 4, ndihedrals, temp) == NULL)
								goto panic;
						}
						if (ap3->anchor != NULL) {
							/*  N1-N2-X-Y or N2-N1-X-Y  */
							/*  for Y, only the first constitute atom is considered  */
							cp2 = AtomConnectData(&ap3->anchor->connect);
							temp[3] = cp2[0];
							if (temp[3] == temp[0] || temp[3] == temp[1])
								continue;
							if (AssignArray(&dihedrals, &ndihedrals, sizeof(Int) * 4, ndihedrals, temp) == NULL)
								goto panic;
						}
					}
				}
			}
			/*  X-N1-N2-X dihedrals  */
			/*  Y-N1-N2-anchor is allowed, but the force may be zero if the angle N1-N2-anchor is */
			/*  close to 180 deg (e.g. in ferrocene, C-anchor-Fe-anchor dihedral should be k=0)  */
			if (ap1->anchor == NULL) {
				ac1 = &ap1->connect;
				cn1 = ac1->count;
			} else {
				ac1 = &ap1->anchor->connect;
				cn1 = 1;  /*  Only the first constitute atom of pi-anchor is considered  */
			}
			if (ap2->anchor == NULL) {
				ac2 = &ap2->connect;
				cn2 = ac2->count;
			} else {
				ac2 = &ap2->anchor->connect;
				cn2 = 1;  /*  Only the first constitute atom of pi-anchor is considered  */
			}
			temp[1] = n1;
			temp[2] = n2;
			cp1 = AtomConnectData(ac1);
			cp2 = AtomConnectData(ac2);
			for (j = 0; j < cn1; j++) {
				temp[0] = cp1[j];
				if (temp[0] == temp[2])
					continue;
				for (k = 0; k < cn2; k++) {
					temp[3] = cp2[k];
					if (temp[3] == temp[0] || temp[3] == temp[1])
						continue;
					if (AssignArray(&dihedrals, &ndihedrals, sizeof(Int) * 4, ndihedrals, temp) == NULL)
						goto panic;
				}
			}
		}
	}
	
	if (angles != NULL) {
		temp[0] = kInvalidIndex;
		if (AssignArray(&angles, &nangles, sizeof(Int) * 3, nangles, temp) == NULL)
			goto panic;
		MoleculeAddAngles(mp, angles, NULL);
		free(angles);
	}
	if (dihedrals != NULL) {
		temp[0] = kInvalidIndex;
		if (AssignArray(&dihedrals, &ndihedrals, sizeof(Int) * 4, ndihedrals, temp) == NULL)
			goto panic;
		MoleculeAddDihedrals(mp, dihedrals, NULL);
		free(dihedrals);
	}

	MoleculeIncrementModifyCount(mp);
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);

	return nbonds;

  panic:
	__MoleculeUnlock(mp);
	Panic("Low memory while adding bonds");
	return -1;  /*  Not reached  */
}

/*  Delete bonds  */
/*  The deleted angles and dihedrals are stored in outRemoval.  */
/*  (*outRemoval) is an array of integers, containing:
      [0..na*3-1]: the angle indices
      [na*3..na*3+nd*4-1]: the dihedral indices
	  [na*3+nd*4..na*3+nd*4+ni*4-1]: the improper indices
    *outRemovedPos is an intgroup denoting the positions of the removed angles/dihedrals/impropers.
	  the angle indices are included as they are,
      the dihedral indices are offset by ATOMS_MAX_NUMBER,
      the improper indices are offset by ATOMS_MAX_NUMBER*2.
    Note: the removed bond indices are not returned, because the caller should already know them.  */
int
MoleculeDeleteBonds(Molecule *mp, Int *bonds, IntGroup *where, Int **outRemoved, IntGroup **outRemovedPos)
{
	Int i, j, n1, n2, nw;
	Int *ip, *jp, na, nd, ni;
	IntGroup *ag, *dg, *ig;
	Atom *ap;
	IntGroupIterator iter;

	if (mp == NULL)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */

	__MoleculeLock(mp);

	/*  Update connects[]  */
	IntGroupIteratorInit(where, &iter);
	while ((i = IntGroupIteratorNext(&iter)) >= 0) {
		n1 = mp->bonds[i * 2];
		n2 = mp->bonds[i * 2 + 1];
		ap = ATOM_AT_INDEX(mp->atoms, n1);
		ip = AtomConnectData(&ap->connect);
		for (j = 0; j < ap->connect.count; j++) {
			if (ip[j] == n2) {
				AtomConnectDeleteEntry(&ap->connect, j);
				break;
			}
		}
		ap = ATOM_AT_INDEX(mp->atoms, n2);
		ip = AtomConnectData(&ap->connect);
		for (j = 0; j < ap->connect.count; j++) {
			if (ip[j] == n1) {
				AtomConnectDeleteEntry(&ap->connect, j);
				break;
			}
		}
	}
	
	/*  Remove bonds, angles, dihedrals, impropers  */
	ag = IntGroupNew();
	dg = ig = NULL;
	na = nd = ni = 0;
	
	nw = IntGroupGetCount(where);
	jp = (Int *)malloc(sizeof(Int) * nw * 2);
	j = 0;
	IntGroupIteratorReset(&iter);
	while ((i = IntGroupIteratorNext(&iter)) >= 0) {
		jp[j++] = mp->bonds[i * 2];
		jp[j++] = mp->bonds[i * 2 + 1];
	}
	IntGroupIteratorRelease(&iter);

	for (i = 0, ip = mp->angles; i < mp->nangles; i++, ip += 3) {
		for (j = 0; j < nw; j++) {
			n1 = jp[j * 2];
			n2 = jp[j * 2 + 1];
			if ((ip[0] == n1 && ip[1] == n2)
				|| (ip[1] == n1 && ip[0] == n2)
				|| (ip[1] == n1 && ip[2] == n2)
				|| (ip[2] == n1 && ip[1] == n2)) {
				if (IntGroupAdd(ag, i, 1) != 0)
					goto panic;
				na++;
				break;
			}
		}
	}
	for (i = 0, ip = mp->dihedrals; i < mp->ndihedrals; i++, ip += 4) {
		for (j = 0; j < nw; j++) {
			n1 = jp[j * 2];
			n2 = jp[j * 2 + 1];
			if ((ip[0] == n1 && ip[1] == n2)
			 || (ip[1] == n1 && ip[0] == n2)
			 || (ip[1] == n1 && ip[2] == n2)
			 || (ip[2] == n1 && ip[1] == n2)
			 || (ip[2] == n1 && ip[3] == n2)
			 || (ip[3] == n1 && ip[2] == n2)) {
				if (dg == NULL)
					dg = IntGroupNew();
				if (IntGroupAdd(dg, i, 1) != 0)
					goto panic;
				nd++;
				break;
			}
		}
	}
	for (i = 0, ip = mp->impropers; i < mp->nimpropers; i++, ip += 4) {
		for (j = 0; j < nw; j++) {
			n1 = jp[j * 2];
			n2 = jp[j * 2 + 1];
			if ((ip[0] == n1 && ip[2] == n2)
			 || (ip[1] == n1 && ip[2] == n2)
			 || (ip[3] == n1 && ip[2] == n2)
			 || (ip[0] == n2 && ip[2] == n1)
			 || (ip[1] == n2 && ip[2] == n1)
			 || (ip[3] == n2 && ip[2] == n1)) {
				if (ig == NULL)
					ig = IntGroupNew();
				if (IntGroupAdd(ig, i, 1) != 0)
					goto panic;
				ni++;
				break;
			}
		}
	}
	free(jp);
	
	if (sRemoveElementsFromArrayAtPositions(mp->bonds, mp->nbonds, bonds, sizeof(Int) * 2, where) != 0)
		goto panic;
	mp->nbonds -= IntGroupGetCount(where);
	if (mp->nbonds == 0) {
		free(mp->bonds);
		mp->bonds = NULL;
	}
	if (mp->bondOrders != NULL) {
		if (sRemoveElementsFromArrayAtPositions(mp->bondOrders, mp->nbondOrders, NULL, sizeof(Double), where) != 0)
			goto panic;
		mp->nbondOrders -= IntGroupGetCount(where);
		if (mp->nbondOrders == 0) {
			free(mp->bondOrders);
			mp->bondOrders = NULL;
		}
	}
	if (na == 0 && nd == 0 && ni == 0)
		ip = NULL;
	else
		ip = (Int *)malloc(sizeof(Int) * (na * 3 + nd * 4 + ni * 4));
	if (na > 0)
		MoleculeDeleteAngles(mp, ip, ag);
	if (nd > 0)
		MoleculeDeleteDihedrals(mp, ip + na * 3, dg);
	if (ni > 0)
		MoleculeDeleteImpropers(mp, ip + na * 3 + nd * 4, ig);
	if (ip != NULL) {
		IntGroupOffset(dg, ATOMS_MAX_NUMBER);
		IntGroupOffset(ig, ATOMS_MAX_NUMBER * 2);
		IntGroupAddIntGroup(ag, dg);
		IntGroupAddIntGroup(ag, ig);
		IntGroupRelease(dg);
		IntGroupRelease(ig);
	}

	if (IntGroupGetCount(ag) == 0) {
		IntGroupRelease(ag);
		ag = NULL;
	}
	
	*outRemoved = ip;
	*outRemovedPos = ag;

	MoleculeIncrementModifyCount(mp);
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);

	return na * 3 + nd * 4 + ni * 4;

  panic:
	__MoleculeUnlock(mp);
	Panic("Low memory while removing bonds");
	return -1;  /*  Not reached  */
}

int
MoleculeAssignBondOrders(Molecule *mp, const Double *orders, IntGroup *where)
{
	Int i, j;
	IntGroupIterator iter;
	if (mp == NULL || orders == NULL || mp->nbonds == 0)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */
	if (mp->bondOrders == NULL) {
		AssignArray(&mp->bondOrders, &mp->nbondOrders, sizeof(Double), mp->nbonds - 1, NULL);
		memset(mp->bondOrders, 0, sizeof(Double) * mp->nbondOrders);
	}
	IntGroupIteratorInit(where, &iter);
	j = 0;
	while ((i = IntGroupIteratorNext(&iter)) >= 0) {
		if (i >= mp->nbondOrders)
			break;
		mp->bondOrders[i] = orders[j++];
	}
	IntGroupIteratorRelease(&iter);
	return 0;
}

int
MoleculeGetBondOrders(Molecule *mp, Double *outOrders, IntGroup *where)
{
	Int i, j;
	IntGroupIterator iter;
	if (mp == NULL || mp->nbonds == 0)
		return 0;
	if (mp->bondOrders == NULL) {
		/*  Returns all zero  */
		i = IntGroupGetCount(where);
		for (j = 0; j < i; j++)
			outOrders[j] = 0.0;
	} else {
		IntGroupIteratorInit(where, &iter);
		j = 0;
		while ((i = IntGroupIteratorNext(&iter)) >= 0) {
			if (i < mp->nbondOrders)
				outOrders[j] = mp->bondOrders[i];
			else outOrders[j] = 0.0;
			j++;
		}
	}
	return 0;
}

int
MoleculeAddAngles(Molecule *mp, const Int *angles, IntGroup *where)
{
	int n1, nc;
	if (mp == NULL || angles == NULL)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */

	__MoleculeLock(mp);
	if (where != NULL)
		nc = IntGroupGetCount(where);
	else {
		for (n1 = 0; angles[n1 * 3] >= 0; n1++)
			;
		nc = n1;
	}
	if (nc > 0) {
		n1 = mp->nangles;
		if (AssignArray(&(mp->angles), &(mp->nangles), sizeof(Int) * 3, n1 + nc - 1, NULL) == NULL
			|| sInsertElementsToArrayAtPositions(mp->angles, n1, angles, nc, sizeof(Int) * 3, where) != 0) {
			__MoleculeUnlock(mp);
			Panic("Low memory while adding angles");
		}
	}
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);
	return nc;
}

int
MoleculeDeleteAngles(Molecule *mp, Int *angles, IntGroup *where)
{
	int nc;
	if (mp == NULL || where == NULL)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */
	__MoleculeLock(mp);
	if (sRemoveElementsFromArrayAtPositions(mp->angles, mp->nangles, angles, sizeof(Int) * 3, where) != 0) {
		__MoleculeUnlock(mp);
		Panic("Bad argument while deleting angles");
	}
	mp->nangles -= (nc = IntGroupGetCount(where));
	if (mp->nangles == 0) {
		free(mp->angles);
		mp->angles = NULL;
	}
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);
	return nc;
}

int
MoleculeAddDihedrals(Molecule *mp, const Int *dihedrals, IntGroup *where)
{
	int n1, nc;
	if (mp == NULL || dihedrals == NULL)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */
	if (where != NULL)
		nc = IntGroupGetCount(where);
	else {
		for (n1 = 0; dihedrals[n1 * 4] >= 0; n1++)
			;
		nc = n1;
	}
	if (nc <= 0)
		return 0;
	n1 = mp->ndihedrals;
	__MoleculeLock(mp);
	if (AssignArray(&(mp->dihedrals), &(mp->ndihedrals), sizeof(Int) * 4, n1 + nc - 1, NULL) == NULL
	|| sInsertElementsToArrayAtPositions(mp->dihedrals, n1, dihedrals, nc, sizeof(Int) * 4, where) != 0) {
		__MoleculeUnlock(mp);
		Panic("Low memory while adding dihedrals");
	}
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);
	return nc;
}

int
MoleculeDeleteDihedrals(Molecule *mp, Int *dihedrals, IntGroup *where)
{	
	int nc;
	if (mp == NULL || where == NULL)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */
	__MoleculeLock(mp);
	if (sRemoveElementsFromArrayAtPositions(mp->dihedrals, mp->ndihedrals, dihedrals, sizeof(Int) * 4, where) != 0) {
		__MoleculeUnlock(mp);
		Panic("Internal error: bad argument while deleting dihedrals");
	}
	mp->ndihedrals -= (nc = IntGroupGetCount(where));
	if (mp->ndihedrals == 0) {
		free(mp->dihedrals);
		mp->dihedrals = NULL;
	}
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);
	return nc;
}

int
MoleculeAddImpropers(Molecule *mp, const Int *impropers, IntGroup *where)
{
	int n1, nc;
	if (mp == NULL || impropers == NULL)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */
	if (where != NULL)
		nc = IntGroupGetCount(where);
	else {
		for (n1 = 0; impropers[n1 * 4] >= 0; n1++)
			;
		nc = n1;
	}
	if (nc <= 0)
		return 0;
	n1 = mp->nimpropers;
	__MoleculeLock(mp);
	if (AssignArray(&(mp->impropers), &(mp->nimpropers), sizeof(Int) * 4, n1 + nc - 1, NULL) == NULL
	|| sInsertElementsToArrayAtPositions(mp->impropers, n1, impropers, nc, sizeof(Int) * 4, where) != 0) {
		__MoleculeUnlock(mp);
		Panic("Low memory while adding impropers");
	}
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);
	return nc;
}

int
MoleculeDeleteImpropers(Molecule *mp, Int *impropers, IntGroup *where)
{
	int nc;
	if (mp == NULL || where == NULL)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */
	__MoleculeLock(mp);
	if (sRemoveElementsFromArrayAtPositions(mp->impropers, mp->nimpropers, impropers, sizeof(Int) * 4, where) != 0) {
		__MoleculeUnlock(mp);
		Panic("Internal error: bad argument while deleting impropers");
	}
	mp->nimpropers -= (nc = IntGroupGetCount(where));
	if (mp->impropers == NULL) {
		free(mp->impropers);
		mp->impropers = NULL;
	}
	__MoleculeUnlock(mp);
	return nc;
}

int
MoleculeLookupBond(Molecule *mp, Int n1, Int n2)
{
	Int i, *ip;
	if (mp == NULL || mp->bonds == NULL)
		return -1;
	for (i = 0, ip = mp->bonds; i < mp->nbonds; i++, ip += 2) {
		if ((n1 == ip[0] && n2 == ip[1]) || (n1 == ip[1] && n2 == ip[0]))
			return i;
	}
	return -1;
}

int
MoleculeLookupAngle(Molecule *mp, Int n1, Int n2, Int n3)
{
	Int i, *ip;
	if (mp == NULL || mp->angles == NULL)
		return -1;
	for (i = 0, ip = mp->angles; i < mp->nangles; i++, ip += 3) {
		if ((n1 == ip[0] && n2 == ip[1] && n3 == ip[2]) ||
			(n1 == ip[2] && n2 == ip[1] && n3 == ip[0]))
			return i;
	}
	return -1;
}

int
MoleculeLookupDihedral(Molecule *mp, Int n1, Int n2, Int n3, Int n4)
{
	Int i, *ip;
	if (mp == NULL || mp->dihedrals == NULL)
		return -1;
	for (i = 0, ip = mp->dihedrals; i < mp->ndihedrals; i++, ip += 4) {
		if ((n1 == ip[0] && n2 == ip[1] && n3 == ip[2] && n4 == ip[3]) ||
			(n1 == ip[3] && n2 == ip[2] && n3 == ip[1] && n4 == ip[0]))
			return i;
	}
	return -1;
}

int
MoleculeLookupImproper(Molecule *mp, Int n1, Int n2, Int n3, Int n4)
{
	Int i, *ip;
	if (mp == NULL || mp->impropers == NULL)
		return -1;
	for (i = 0, ip = mp->impropers; i < mp->nimpropers; i++, ip += 4) {
		if (n3 != ip[2])
			continue;
		if ((n1 == ip[0] && ((n2 == ip[1] && n4 == ip[3]) || (n2 == ip[3] && n4 == ip[1]))) ||
			(n1 == ip[1] && ((n2 == ip[0] && n4 == ip[3]) || (n2 == ip[3] && n4 == ip[0]))) ||
			(n1 == ip[3] && ((n2 == ip[0] && n4 == ip[1]) || (n2 == ip[1] && n4 == ip[0]))))
			return i;
	}
	return -1;
}

/*  Remove the bond at bondIndex and create two dummy atoms instead.
    The dummy atoms are placed at the end of atoms[], and the residue
	numbers are the same as the root atoms (i.e. the atoms to which
	the dummy atoms are connected). The indices are returned in
	dummyIndices[0,1].  */
int
MoleculeConvertBondToDummies(Molecule *mp, Int bondIndex, Int *dummyIndices)
{
	Int roots[3], newBonds[5];
	Vector dr;
	Atom *rootp[2];
	Atom na[2], *nap;
	int i, natoms;
	IntGroup *ig;
	if (mp == NULL || mp->noModifyTopology)
		return 0;
	if (bondIndex < 0 || bondIndex >= mp->nbonds)
		return -1;
	roots[0] = mp->bonds[bondIndex * 2];
	roots[1] = mp->bonds[bondIndex * 2 + 1];
	roots[2] = kInvalidIndex;
	rootp[0] = ATOM_AT_INDEX(mp->atoms, roots[0]);
	rootp[1] = ATOM_AT_INDEX(mp->atoms, roots[1]);
	VecSub(dr, rootp[0]->r, rootp[1]->r);
	for (i = 0; i < 2; i++) {
		float w;
		nap = &na[i];
		memmove(nap, rootp[i], sizeof(na));
		nap->aname[0] = '*';
		strcpy(nap->element, "Du");
		nap->type = 0;
		nap->charge = nap->weight = 0.0;
		nap->atomicNumber = 0;
		nap->connect.count = 0;
		w = (i == 0 ? 0.4 : -0.4);
		VecScaleInc(nap->r, dr, w);
		VecZero(nap->v);
		VecZero(nap->f);
		nap->intCharge = 0;
		nap->exflags = 0;
	}

	/*  Expand atoms array and append the dummy atoms at the end  */
	__MoleculeLock(mp);
	natoms = mp->natoms;
	if (AssignArray(&(mp->atoms), &(mp->natoms), gSizeOfAtomRecord, natoms + 1, NULL) == NULL)
		goto panic;
	memmove(&mp->atoms[natoms], na, gSizeOfAtomRecord * 2);
	dummyIndices[0] = natoms;
	dummyIndices[1] = natoms + 1;

	/*  Remove the old bond and create new bonds  */
	ig = IntGroupNewWithPoints(bondIndex, 1, -1);
	if (ig == NULL)
		goto panic;
	MoleculeDeleteBonds(mp, NULL, ig, NULL, NULL);
	IntGroupRelease(ig);
	newBonds[0] = roots[0];
	newBonds[1] = dummyIndices[0];
	newBonds[2] = roots[1];
	newBonds[3] = dummyIndices[1];
	newBonds[4] = kInvalidIndex;
	
	i = (MoleculeAddBonds(mp, 2, newBonds, NULL, 1) < 0 ? -1 : 0);
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);
	return i;

panic:
	__MoleculeUnlock(mp);
	Panic("Low memory during creating dummy atoms");
	return 1;
}

/*  Remove two dummy atoms at dummyIndices[0], dummyIndices[1] and create
    a bond between the two root atoms. The value bondIndex is used as a
	hint where to store the new bond; if 0 <= bondIndex <= nbonds, then
	the new bond is stored as the bondIndex'th bond; otherwise, bondIndex
	is ignored and the new bond is stored at the end of bonds[].  */
int
MoleculeConvertDummiesToBond(Molecule *mp, Int bondIndex, Int *dummyIndices)
{
	return 0;
}

/*
Int
MoleculeReplaceAllAngles(Molecule *mol, Int nangles, const Int *angles, Int **outAngles)
{
	Int n1, *np1;
	if (mol == NULL || mol->noModifyTopology)
		return -1;
	n1 = mol->nangles;
	np1 = mol->angles;
	mol->nangles = 0;
	mol->angles = NULL;
	if (nangles > 0) {
		__MoleculeLock(mol);
		NewArray(&mol->angles, &mol->nangles, sizeof(Int) * 3, nangles);
		memmove(mol->angles, angles, sizeof(Int) * 3 * nangles);
		mol->needsMDRebuild = 1;
		__MoleculeUnlock(mol);
	}
	*outAngles = np1;
	return n1;
}
						
Int
MoleculeReplaceAllDihedrals(Molecule *mol, Int ndihedrals, const Int *dihedrals, Int **outDihedrals)
{
	Int n1, *np1;
	if (mol == NULL || mol->noModifyTopology)
		return -1;
	n1 = mol->ndihedrals;
	np1 = mol->dihedrals;
	mol->ndihedrals = 0;
	mol->dihedrals = NULL;
	if (ndihedrals > 0) {
		__MoleculeLock(mol);
		NewArray(&mol->dihedrals, &mol->ndihedrals, sizeof(Int) * 4, ndihedrals);
		memmove(mol->dihedrals, dihedrals, sizeof(Int) * 4 * ndihedrals);
		mol->needsMDRebuild = 1;
		__MoleculeUnlock(mol);
	}
	*outDihedrals = np1;
	return n1;
}

Int
MoleculeReplaceAllImpropers(Molecule *mol, Int nimpropers, const Int *impropers, Int **outImpropers)
{
	Int n1, *np1;
	if (mol == NULL || mol->noModifyTopology)
		return -1;
	n1 = mol->nimpropers;
	np1 = mol->impropers;
	mol->nimpropers = 0;
	mol->impropers = NULL;
	if (nimpropers > 0) {
		__MoleculeLock(mol);
		NewArray(&mol->impropers, &mol->nimpropers, sizeof(Int) * 4, nimpropers);
		memmove(mol->impropers, impropers, sizeof(Int) * 4 * nimpropers);
		mol->needsMDRebuild = 1;
		__MoleculeUnlock(mol);
	}
	*outImpropers = np1;
	return n1;
}
*/

Int
MoleculeFindMissingAngles(Molecule *mol, Int **outAngles)
{
	Int i, j, k, *ip;
	Atom *ap;
	Int nangles;
	Int *angles;
	
	if (mol == NULL || mol->natoms == 0 || mol->atoms == NULL)
		return 0;  /*  molecule is empty  */
	if (mol->noModifyTopology)
		return -1;
	nangles = 0;
	angles = NULL;
	for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
		Int *cp = AtomConnectData(&ap->connect);
		if (ap->anchor != NULL)
			continue;
		for (j = 0; j < ap->connect.count; j++) {
			Int j0 = cp[j];
			if (ATOM_AT_INDEX(mol->atoms, j0)->anchor != NULL)
				continue;
			for (k = j + 1; k < ap->connect.count; k++) {
				Int k0 = cp[k];
				if (ATOM_AT_INDEX(mol->atoms, k0)->anchor != NULL)
					continue;
				if (MoleculeLookupAngle(mol, j0, i, k0) < 0) {
					ip = (Int *)AssignArray(&angles, &nangles, sizeof(Int) * 3, nangles, NULL);
					ip[0] = j0;
					ip[1] = i;
					ip[2] = k0;
				}
			}
		}
	}
	if (nangles > 0) {
		ip = (Int *)AssignArray(&angles, &nangles, sizeof(Int) * 3, nangles, NULL);
		ip[0] = -1;
		nangles--;
	}
	if (outAngles != NULL)
		*outAngles = angles;
	return nangles;
}

Int
MoleculeFindMissingDihedrals(Molecule *mol, Int **outDihedrals)
{
	Int n1, n2, n3, n4, *ip, *cp2, *cp3;
	Atom *ap2, *ap3;
	Int ndihedrals;
	Int *dihedrals;
	
	if (mol == NULL || mol->natoms == 0 || mol->atoms == NULL)
		return 0;  /*  molecule is empty  */
	ndihedrals = 0;
	dihedrals = NULL;
	for (n2 = 0, ap2 = mol->atoms; n2 < mol->natoms; n2++, ap2 = ATOM_NEXT(ap2)) {
		Int i1, i3, i4, *ip;
		if (ap2->anchor != NULL)
			continue;
		cp2 = AtomConnectData(&ap2->connect);
		for (i3 = 0; i3 < ap2->connect.count; i3++) {
			n3 = cp2[i3];
			if (n2 > n3)
				continue;
			ap3 = ATOM_AT_INDEX(mol->atoms, n3);
			if (ap3->anchor != NULL)
				continue;
			cp3 = AtomConnectData(&ap3->connect);
			for (i1 = 0; i1 < ap2->connect.count; i1++) {
				n1 = cp2[i1];
				if (n1 == n3)
					continue;
				if (ATOM_AT_INDEX(mol->atoms, n1)->anchor != NULL)
					continue;
				for (i4 = 0; i4 < ap3->connect.count; i4++) {
					n4 = cp3[i4];
					if (n2 == n4 || n1 == n4)
						continue;
					if (ATOM_AT_INDEX(mol->atoms, n4)->anchor != NULL)
						continue;
					if (MoleculeLookupDihedral(mol, n1, n2, n3, n4) < 0) {
						ip = (Int *)AssignArray(&dihedrals, &ndihedrals, sizeof(Int) * 4, ndihedrals, NULL);
						ip[0] = n1;
						ip[1] = n2;
						ip[2] = n3;
						ip[3] = n4;
					}
				}
			}
		}
	}
	if (ndihedrals > 0) {
		ip = (Int *)AssignArray(&dihedrals, &ndihedrals, sizeof(Int) * 4, ndihedrals, NULL);
		ip[0] = -1;
		ndihedrals--;
	}
	if (outDihedrals != NULL)
		*outDihedrals = dihedrals;
	return ndihedrals;
}

Int
MoleculeFindMissingImpropers(Molecule *mol, Int **outImpropers)
{
	Int n1, n2, n3, n4, t1, t2, t3, t4, *ip, *cp;
	Parameter *par = mol->par;
	Atom *ap, *ap3;
	Int nimpropers;
	Int *impropers;
	
	if (mol == NULL || mol->natoms == 0 || mol->atoms == NULL)
		return 0;  /*  molecule is empty  */
	if ((par == NULL || par->nimproperPars == 0) && (gBuiltinParameters == NULL || gBuiltinParameters->nimproperPars == 0))
		return 0;  /*  No improper parameters are defined  */
	nimpropers = 0;
	impropers = NULL;
	ap = mol->atoms;
	for (n3 = 0, ap3 = ap; n3 < mol->natoms; n3++, ap3 = ATOM_NEXT(ap3)) {
		Int i1, i2, i4, found, *ip;
		t3 = ap3->type;
		cp = AtomConnectData(&ap3->connect);
		for (i1 = 0; i1 < ap3->connect.count; i1++) {
			n1 = cp[i1];
			t1 = ATOM_AT_INDEX(ap, n1)->type;
			for (i2 = i1 + 1; i2 < ap3->connect.count; i2++) {
				n2 = cp[i2];
				t2 = ATOM_AT_INDEX(ap, n2)->type;
				for (i4 = i2 + 1; i4 < ap3->connect.count; i4++) {
					n4 = cp[i4];
					t4 = ATOM_AT_INDEX(ap, n4)->type;
					found = 0;
					if (ParameterLookupImproperPar(par, t1, t2, t3, t4, n1, n2, n3, n4, 0) != NULL)
						found = 1;
					else if (ParameterLookupImproperPar(gBuiltinParameters, t1, t2, t3, t4, -1, -1, -1, -1, 0) != NULL)
						found = 1;
					if (found && MoleculeLookupImproper(mol, n1, n2, n3, n4) < 0) {
						ip = (Int *)AssignArray(&impropers, &nimpropers, sizeof(Int) * 4, nimpropers, NULL);
						ip[0] = n1;
						ip[1] = n2;
						ip[2] = n3;
						ip[3] = n4;
					}
				}
			}
		}
	}
	if (nimpropers > 0) {
		ip = (Int *)AssignArray(&impropers, &nimpropers, sizeof(Int) * 4, nimpropers, NULL);
		ip[0] = -1;
		nimpropers--;
	}
	if (outImpropers != NULL)
		*outImpropers = impropers;
	return nimpropers;
}

#pragma mark ====== Residues ======

void
MoleculeCleanUpResidueTable(Molecule *mp)
{
	int i, maxres;
	Atom *ap;
	if (mp == NULL || mp->natoms == 0)
		return;
	maxres = 0;
	__MoleculeLock(mp);
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (ap->resSeq >= maxres)
			maxres = ap->resSeq + 1;
		if (ap->resSeq < mp->nresidues) {
			if (strncmp(ap->resName, mp->residues[ap->resSeq], 4) != 0)
				strncpy(ap->resName, mp->residues[ap->resSeq], 4);
		} else {
			AssignArray(&mp->residues, &mp->nresidues, 4, ap->resSeq, ap->resName);
		}
	}
	if (maxres < mp->nresidues)
		mp->nresidues = maxres;
	__MoleculeUnlock(mp);
}

/*  Change the number of residues. If nresidues is greater than the current value,
    then the array mp->residues is expanded with null names. If nresidues is smaller
	than the current value, mp->nresidues is set to the smallest possible value
	that is no smaller than nresidues and larger than any of the resSeq values.  */
int
MoleculeChangeNumberOfResidues(Molecule *mp, int nresidues)
{
	int n;
	if (mp == NULL)
		return 0;
	if (mp->nresidues == nresidues)
		return nresidues;
	else if (mp->nresidues < nresidues) {
		__MoleculeLock(mp);
		n = mp->nresidues;
		AssignArray(&(mp->residues), &(mp->nresidues), 4, nresidues - 1, NULL);
		while (n < nresidues)
			mp->residues[n++][0] = 0;
		__MoleculeUnlock(mp);
		return nresidues;
	} else {
		int i;
		Atom *ap;
		n = nresidues;
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ap->resSeq >= n)
				n = ap->resSeq + 1;
		}
		mp->nresidues = n;
		return n;
	}
}

int
MoleculeChangeResidueNumberWithArray(Molecule *mp, IntGroup *group, Int *resSeqs)
{
	IntGroupIterator iter;
	int withArray, resSeq, maxSeq;
	int i, j;
	Atom *ap;
	
	/*  If LSB of resSeqs is 1, then a constant value is used for all specified atoms  */
	if (((int)resSeqs & 1) == 0) {
		withArray = 1;
		resSeq = 0;
	} else {
		withArray = 0;
		resSeq = ((int)resSeqs - 1) / 2;
	}
	
	IntGroupIteratorInit(group, &iter);

	/*  Change resSeqs  */
	maxSeq = 0;
	j = 0;
	__MoleculeLock(mp);
	while ((i = IntGroupIteratorNext(&iter)) >= 0) {
		ap = ATOM_AT_INDEX(mp->atoms, i);
		if (withArray)
			resSeq = resSeqs[j++];
		if (resSeq > maxSeq)
			maxSeq = resSeq;
		ap->resSeq = resSeq;
	}
	__MoleculeUnlock(mp);

	/*  Expand array if necessary  */
	if (maxSeq >= mp->nresidues)
		MoleculeChangeNumberOfResidues(mp, maxSeq + 1);

	/*  Synchronize resName and residues[]  */
	__MoleculeLock(mp);
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (group != NULL && IntGroupLookup(group, i, NULL) == 0)
			continue;
		strncpy(ap->resName, mp->residues[ap->resSeq], 4);
	}
	IntGroupIteratorRelease(&iter);
	__MoleculeUnlock(mp);
	
	MoleculeIncrementModifyCount(mp);
	
	return 0;
}

int
MoleculeChangeResidueNumber(Molecule *mp, IntGroup *group, int resSeq)
{
	return MoleculeChangeResidueNumberWithArray(mp, group, (Int *)(resSeq * 2 + 1));
}

/*  Offset the residue numbers by a certain amount. The argument nresidues, if non-negative,
    specifies the mp->nresidues after modifying the residue numbers.
	If all atoms are modified, then the table of residue names is also shifted. Otherwise,
	the table of residue names is not touched. */
int
MoleculeOffsetResidueNumbers(Molecule *mp, IntGroup *group, int offset, int nresidues)
{
	int i, maxSeq, nmodatoms;
	Atom *ap;
	IntGroupIterator iter;
	IntGroupIteratorInit(group, &iter);
	maxSeq = 0;
	if (nresidues < 0)
		nresidues = mp->nresidues;
	nmodatoms = 0;
	__MoleculeLock(mp);
	while ((i = IntGroupIteratorNext(&iter)) >= 0) {
		ap = ATOM_AT_INDEX(mp->atoms, i);
		ap->resSeq += offset;
		if (ap->resSeq < 0) {
			/*  Bad argument; undo change and returns this index + 1  */
			int bad_index = i;
			ap->resSeq -= offset;
			while ((i = IntGroupIteratorLast(&iter)) >= 0) {
				ap = ATOM_AT_INDEX(mp->atoms, i);
				ap->resSeq -= offset;
			}
			IntGroupIteratorRelease(&iter);
			return bad_index + 1;
		}
		if (ap->resSeq > maxSeq)
			maxSeq = ap->resSeq;
		nmodatoms++;
	}
	if (maxSeq >= nresidues)
		nresidues = maxSeq + 1;
	if (offset < 0 && nmodatoms == mp->natoms) {
		/*  Shift the residue names downward  */
		memmove(mp->residues, mp->residues - offset, 4 * (mp->nresidues + offset));
	}
	__MoleculeUnlock(mp);
	MoleculeChangeNumberOfResidues(mp, nresidues);
	if (offset > 0 && nmodatoms == mp->natoms) {
		/*  Shift the residue names upward  */
		__MoleculeLock(mp);
		memmove(mp->residues + offset, mp->residues, 4 * (mp->nresidues - offset));
		__MoleculeUnlock(mp);
	}
	IntGroupIteratorRelease(&iter);

	MoleculeIncrementModifyCount(mp);
	
	return 0;
}

/*  Change residue names for the specified residue numbers. Names is an array of
    chars containing argc*4 characters, and every 4 characters represent a
	residue name; characters '\x01'-'\x1f' are converted to '\0', which allow 
	names to be handled as a C string.  */
int
MoleculeChangeResidueNames(Molecule *mp, int argc, Int *resSeqs, char *names)
{
	int i, maxSeq;
	Atom *ap;
	maxSeq = 0;
	for (i = 0; i < argc; i++) {
		if (maxSeq < resSeqs[i])
			maxSeq = resSeqs[i];
	}
	if (maxSeq >= mp->nresidues)
		MoleculeChangeNumberOfResidues(mp, maxSeq + 1);
	__MoleculeLock(mp);
	for (i = 0; i < argc; i++) {
		char *p = mp->residues[resSeqs[i]];
		int j;
		strncpy(p, names + i * 4, 4);
		for (j = 0; j < 4; j++) {
			if (p[j] >= 0 && p[j] < 0x20)
				p[j] = 0;
		}
	}
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		strncpy(ap->resName, mp->residues[ap->resSeq], 4);
	}
	__MoleculeUnlock(mp);

	MoleculeIncrementModifyCount(mp);
	
	return 0;
}

/*  Returns the maximum residue number actually used  */
int
MoleculeMaximumResidueNumber(Molecule *mp, IntGroup *group)
{
	int i, maxSeq;
	Atom *ap;
	maxSeq = -1;
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (group != NULL && IntGroupLookup(group, i, NULL) == 0)
			continue;
		if (ap->resSeq > maxSeq)
			maxSeq = ap->resSeq;
	}
	return maxSeq;
}

/*  Returns the minimum residue number actually used  */
int
MoleculeMinimumResidueNumber(Molecule *mp, IntGroup *group)
{
	int i, minSeq;
	Atom *ap;
	minSeq = ATOMS_MAX_NUMBER;
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (group != NULL && IntGroupLookup(group, i, NULL) == 0)
			continue;
		if (ap->resSeq < minSeq)
			minSeq = ap->resSeq;
	}
	return (minSeq == ATOMS_MAX_NUMBER ? -1 : minSeq);
}

#pragma mark ====== Sort by Residues ======

static int
sAtomSortComparator(const void *a, const void *b)
{
	const Atom *ap, *bp;
	ap = *((const Atom **)a);
	bp = *((const Atom **)b);
	if (ap->resSeq == bp->resSeq) {
		/*  Retain the original order (i.e. atom with larger pointer address is larger)  */
		if (ap < bp)
			return -1;
		else if (ap > bp)
			return 1;
		else return 0;
	} else {
		/*  Compare the residue sequence. However, residue sequence 0 is always larger.  */
		if (ap->resSeq == 0)
			return 1;
		else if (bp->resSeq == 0)
			return -1;
		else if (ap->resSeq < bp->resSeq)
			return -1;
		else if (ap->resSeq > bp->resSeq)
			return 1;
		else return 0;
	}
}

static void
sMoleculeReorder(Molecule *mp)
{
	int i, res, prevRes;
	Atom **apArray;
	Int *old2new;
	Atom *newAtoms;
	if (mp == NULL || mp->natoms <= 1)
		return;

	/*  Sort the atoms, bonds, etc. */
	apArray = (Atom **)calloc(sizeof(Atom *), mp->natoms);
	old2new = (Int *)calloc(sizeof(Int), mp->natoms);
	newAtoms = (Atom *)calloc(gSizeOfAtomRecord, mp->natoms);
	if (apArray == NULL || old2new == NULL || newAtoms == NULL)
		Panic("Low memory during reordering atoms");
	for (i = 0; i < mp->natoms; i++)
		apArray[i] = ATOM_AT_INDEX(mp->atoms, i);

	/*  Sort the atoms. Note: apArray is an array of "Pointer to Atom"  */
	qsort(apArray, mp->natoms, sizeof(Atom *), sAtomSortComparator);
	
	/*  Make a table of 'which atom becomes which'  */
	for (i = 0; i < mp->natoms; i++) {
		int j = ((char *)(apArray[i]) - (char *)(mp->atoms)) / gSizeOfAtomRecord;
		old2new[j] = i;
	}
	
	/*  Renumber the bonds, etc.  */
	for (i = 0; i < mp->nbonds * 2; i++) {
		mp->bonds[i] = old2new[mp->bonds[i]];
	}
	for (i = 0; i < mp->nangles * 3; i++) {
		mp->angles[i] = old2new[mp->angles[i]];
	}
	for (i = 0; i < mp->ndihedrals * 4; i++) {
		mp->dihedrals[i] = old2new[mp->dihedrals[i]];
	}
	for (i = 0; i < mp->nimpropers * 4; i++) {
		mp->impropers[i] = old2new[mp->impropers[i]];
	}
	for (i = 0; i < mp->natoms; i++) {
		Int *ip, j;
		ip = AtomConnectData(&(apArray[i]->connect));
		for (j = 0; j < apArray[i]->connect.count; j++, ip++)
			*ip = old2new[*ip];
	}
	
	/*  Renumber the residues so that the residue numbers are contiguous  */
	res = prevRes = 0;
	for (i = 0; i < mp->natoms; i++) {
		if (apArray[i]->resSeq == 0)
			break;
		if (apArray[i]->resSeq != prevRes) {
			res++;
			prevRes = apArray[i]->resSeq;
			if (prevRes != res) {
				strncpy(mp->residues[res], mp->residues[prevRes], 4);
			}
		}
		apArray[i]->resSeq = res;
	}
	mp->nresidues = res + 1;

	/*  Sort the atoms and copy back to atoms[] */
	for (i = 0; i < mp->natoms; i++) {
		memmove(ATOM_AT_INDEX(newAtoms, i), apArray[i], gSizeOfAtomRecord);
	}
	memmove(mp->atoms, apArray, gSizeOfAtomRecord * mp->natoms);
	
	/*  Free the locally allocated storage  */
	free(newAtoms);
	free(old2new);
	free(apArray);
}

/*  Renumber atoms  */
int
MoleculeRenumberAtoms(Molecule *mp, const Int *new2old, Int *old2new_out, Int isize)
{
	Int *old2new, i, j, retval;
	Atom *saveAtoms;
	if (mp == NULL)
		return 0;
	if (mp->noModifyTopology)
		return -1;
	if (old2new_out != NULL)
		old2new = old2new_out;
	else
		old2new = (Int *)calloc(sizeof(Int), mp->natoms);
	saveAtoms = (Atom *)calloc(gSizeOfAtomRecord, mp->natoms);
	if (old2new == NULL || saveAtoms == NULL)
		Panic("Low memory during reordering atoms");
	memmove(saveAtoms, mp->atoms, gSizeOfAtomRecord * mp->natoms);
	__MoleculeLock(mp);
	for (i = 0; i < mp->natoms; i++)
		old2new[i] = -1;
	for (i = 0; i < isize && i < mp->natoms; i++) {
		j = new2old[i];
		if (j < 0 || j >= mp->natoms) {
			retval = 1; /* Out of range */
			goto end;
		}
		if (old2new[j] != -1) {
			retval = 2;  /*  Duplicate entry  */
			goto end;
		}
		old2new[j] = i;
	}
	if (i < mp->natoms) {
		for (j = 0; j < mp->natoms; j++) {
			if (old2new[j] != -1)
				continue;
			old2new[j] = i++;
		}
	}
	if (i != mp->natoms) {
		retval = 3;  /*  Internal inconsistency  */
		goto end;
	}

	/*  Renumber the bonds, etc.  */
	for (i = 0; i < mp->nbonds * 2; i++) {
		mp->bonds[i] = old2new[mp->bonds[i]];
	}
	for (i = 0; i < mp->nangles * 3; i++) {
		mp->angles[i] = old2new[mp->angles[i]];
	}
	for (i = 0; i < mp->ndihedrals * 4; i++) {
		mp->dihedrals[i] = old2new[mp->dihedrals[i]];
	}
	for (i = 0; i < mp->nimpropers * 4; i++) {
		mp->impropers[i] = old2new[mp->impropers[i]];
	}
	/*  Renumber the connection table and pi anchor table  */
	for (i = 0; i < mp->natoms; i++) {
		Atom *ap = ATOM_AT_INDEX(saveAtoms, i);
		Int *ip = AtomConnectData(&ap->connect);
		for (j = 0; j < ap->connect.count; j++, ip++)
			*ip = old2new[*ip];
		if (ap->anchor != NULL) {
			ip = AtomConnectData(&ap->anchor->connect);
			for (j = 0; j < ap->anchor->connect.count; j++, ip++)
				*ip = old2new[*ip];
		}
	}
	
	if (mp->par != NULL) {
		/*  Renumber the parameters  */
		int n;
		for (j = kFirstParType; j <= kLastParType; j++) {
			n = ParameterGetCountForType(mp->par, j);
			for (i = 0; i < n; i++) {
				UnionPar *up = ParameterGetUnionParFromTypeAndIndex(mp->par, j, i);
				if (up != NULL)
					ParameterRenumberAtoms(j, up, mp->natoms, old2new);
			}
		}
	}
	
	/*  Renumber the atoms  */
	for (i = 0; i < mp->natoms; i++)
		memmove(ATOM_AT_INDEX(mp->atoms, old2new[i]), ATOM_AT_INDEX(saveAtoms, i), gSizeOfAtomRecord);
	retval = 0;
	
	MoleculeIncrementModifyCount(mp);
	mp->needsMDRebuild = 1;

  end:
	__MoleculeUnlock(mp);
	free(saveAtoms);
	if (old2new_out == NULL)
		free(old2new);
	return retval;
}

#pragma mark ====== Coordinate Transform ======

void
MoleculeTransform(Molecule *mp, Transform tr, IntGroup *group)
{
	int i;
	Atom *ap;
	Symop new_symop;
	Transform rtr, symtr;
	if (mp == NULL || tr == NULL)
		return;
	TransformInvert(rtr, tr);
	__MoleculeLock(mp);
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (group == NULL || IntGroupLookup(group, i, NULL) != 0) {
			TransformVec(&ap->r, tr, &ap->r);
			if (!SYMOP_ALIVE(ap->symop))
				continue;
			/*  Transform symop  */
			if (MoleculeGetTransformForSymop(mp, ap->symop, &symtr, 1) != 0)
				continue;
			TransformMul(symtr, tr, symtr);
			if (group == NULL || IntGroupLookup(group, ap->symbase, NULL) != 0)
				TransformMul(symtr, symtr, rtr);
		} else {
			if (!SYMOP_ALIVE(ap->symop))
				continue;
			/*  Transform symop if the base atom is transformed  */
			if (group != NULL && IntGroupLookup(group, ap->symbase, NULL) == 0)
				continue;
			if (MoleculeGetTransformForSymop(mp, ap->symop, &symtr, 1) != 0)
				continue;
			TransformMul(symtr, symtr, rtr);
		}
		if (MoleculeGetSymopForTransform(mp, symtr, &new_symop, 1) != 0)
			continue;
		ap->symop = new_symop;
	}
	mp->needsMDCopyCoordinates = 1;
	__MoleculeUnlock(mp);
	sMoleculeNotifyChangeAppearance(mp);
}

/*
void
MoleculeMove(Molecule *mp, Transform tr, IntGroup *group)
{
	int i;
	Atom *ap;
	if (mp == NULL || tr == NULL)
		return;
	__MoleculeLock(mp);
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (group != NULL && IntGroupLookup(group, i, NULL) == 0)
			continue;
		TransformVec(&ap->r, tr, &ap->r);
	}
	mp->needsMDCopyCoordinates = 1;
	__MoleculeUnlock(mp);
	sMoleculeNotifyChangeAppearance(mp);
}
*/

void
MoleculeTranslate(Molecule *mp, const Vector *vp, IntGroup *group)
{
	Transform tr;
	if (mp == NULL || vp == NULL)
		return;
	memset(tr, 0, sizeof(tr));
	tr[0] = tr[4] = tr[8] = 1.0;
	tr[9] = vp->x;
	tr[10] = vp->y;
	tr[11] = vp->z;
	MoleculeTransform(mp, tr, group);
}

void
MoleculeRotate(Molecule *mp, const Vector *axis, Double angle, const Vector *center, IntGroup *group)
{
	Transform tr;
	TransformForRotation(tr, axis, angle, center);
	MoleculeTransform(mp, tr, group);
}

int
MoleculeCenterOfMass(Molecule *mp, Vector *center, IntGroup *group)
{
	int i;
	Atom *ap;
	Double w;
	if (mp == NULL || center == NULL)
		return 1;
	if (mp->natoms == 0 || (group != NULL && IntGroupGetCount(group) == 0))
		return 2;   /*  Empty molecule  */
	w = 0.0;
	center->x = center->y = center->z = 0.0;
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (group != NULL && IntGroupLookup(group, i, NULL) == 0)
			continue;
		VecScaleInc(*center, ap->r, ap->weight);
		w += ap->weight;
	}
	if (w < 1e-7)
		return 3;  /*  Atomic weights are not defined?  */
	w = 1.0 / w;
	VecScaleSelf(*center, w);
	return 0;
}

int
MoleculeBounds(Molecule *mp, Vector *min, Vector *max, IntGroup *group)
{
	Vector vmin, vmax;
	int i;
	Atom *ap;
	if (mp == NULL)
		return 1;
	if (mp->natoms == 0 || (group != NULL && IntGroupGetCount(group) == 0))
		return 2;   /*  Empty molecule  */
	vmin.x = vmin.y = vmin.z = 1e50;
	vmax.x = vmax.y = vmax.z = -1e50;
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (group != NULL && IntGroupLookup(group, i, NULL) == 0)
			continue;
		if (vmin.x > ap->r.x)
			vmin.x = ap->r.x;
		if (vmin.y > ap->r.y)
			vmin.y = ap->r.y;
		if (vmin.z > ap->r.z)
			vmin.z = ap->r.z;
		if (vmax.x < ap->r.x)
			vmax.x = ap->r.x;
		if (vmax.y < ap->r.y)
			vmax.y = ap->r.y;
		if (vmax.z < ap->r.z)
			vmax.z = ap->r.z;
	}
	if (min != NULL)
		*min = vmin;
	if (max != NULL)
		*max = vmax;
	return 0;	
}

#pragma mark ====== Measurements ======

Double
MoleculeMeasureBond(Molecule *mp, const Vector *vp1, const Vector *vp2)
{
	Vector r1, r2;
/*	if (mp->is_xtal_coord) {
		TransformVec(&r1, mp->cell->tr, vp1);
		TransformVec(&r2, mp->cell->tr, vp2);
	} else */ {
		r1 = *vp1;
		r2 = *vp2;
	}
	VecDec(r1, r2);
	return VecLength(r1);
}

Double
MoleculeMeasureAngle(Molecule *mp, const Vector *vp1, const Vector *vp2, const Vector *vp3)
{
	Vector r1, r2, r3;
	double w;
/*	if (mp->is_xtal_coord) {
		TransformVec(&r1, mp->cell->tr, vp1);
		TransformVec(&r2, mp->cell->tr, vp2);
		TransformVec(&r3, mp->cell->tr, vp3);
	} else */ {
		r1 = *vp1;
		r2 = *vp2;
		r3 = *vp3;
	}
	VecDec(r1, r2);
	VecDec(r3, r2);
	w = VecLength(r1) * VecLength(r3);
	if (w < 1e-20)
		return NAN;
	return acos(VecDot(r1, r3) / w) * kRad2Deg;
}

Double
MoleculeMeasureDihedral(Molecule *mp, const Vector *vp1, const Vector *vp2, const Vector *vp3, const Vector *vp4)
{
	Vector r1, r2, r3, r4, r21, r32, r43, v1, v2, v3;
	double w1, w2, w3;
/*	if (mp->is_xtal_coord) {
		TransformVec(&r1, mp->cell->tr, vp1);
		TransformVec(&r2, mp->cell->tr, vp2);
		TransformVec(&r3, mp->cell->tr, vp3);
		TransformVec(&r4, mp->cell->tr, vp4);
	} else */ {
		r1 = *vp1;
		r2 = *vp2;
		r3 = *vp3;
		r4 = *vp4;
	}
	VecSub(r21, r1, r2);
	VecSub(r32, r2, r3);
	VecSub(r43, r3, r4);
	VecCross(v1, r21, r32);
	VecCross(v2, r32, r43);
	VecCross(v3, r32, v1);
	w1 = VecLength(v1);
	w2 = VecLength(v2);
	w3 = VecLength(v3);
	if (w1 < 1e-10 || w2 < 1e-10 || w3 < 1e-10) {
		return NAN;
	} else {
		w1 = 1.0 / w1;
		w2 = 1.0 / w2;
		w3 = 1.0 / w3;
		VecScaleSelf(v1, w1);
		VecScaleSelf(v2, w2);
		VecScaleSelf(v3, w3);
		return -atan2(VecDot(v3, v2), VecDot(v1, v2)) * kRad2Deg;
	}
}

#pragma mark ====== XtalCell Parameters ======

void
MoleculeXtalToCartesian(Molecule *mp, Vector *dst, const Vector *src)
{
	if (mp->cell != NULL) {
		TransformVec(dst, mp->cell->tr, src);
	} else *dst = *src;
}

void
MoleculeCartesianToXtal(Molecule *mp, Vector *dst, const Vector *src)
{
	if (mp->cell != NULL) {
		TransformVec(dst, mp->cell->rtr, src);
	} else *dst = *src;
}

int
MoleculeCalculateCellFromAxes(XtalCell *cp, int calc_abc)
{
	static Transform identityTransform = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
	int n1, n2, n3;
	Vector *vp1, *vp2, *vp3;
	Vector v1, v2;

	if (cp == NULL)
		return 0;
	for (n1 = 0; n1 < 3; n1++) {
		if (cp->flags[n1] != 0)
			break;
	}
	if (n1 == 3) {
		/*  All directions are non-periodic  */
		memmove(&(cp->tr), &identityTransform, sizeof(Transform));
		memmove(&(cp->rtr), &identityTransform, sizeof(Transform));
	} else {
		n2 = (n1 + 1) % 3;
		n3 = (n1 + 2) % 3;
		vp1 = &(cp->axes[n1]);
		vp2 = &(cp->axes[n2]);
		vp3 = &(cp->axes[n3]);
		cp->tr[n1*3] = vp1->x;
		cp->tr[n1*3+1] = vp1->y;
		cp->tr[n1*3+2] = vp1->z;
		cp->tr[9] = cp->origin.x;
		cp->tr[10] = cp->origin.y;
		cp->tr[11] = cp->origin.z;
		if (cp->flags[n2] == 0 || cp->flags[n3] == 0) {
			/*  1-dimensional or 2-dimensional system  */
			/*  Create "dummy" axes, so that transforms between internal and cartesian coordinates are
			 possible with a single matrix  */
			if (cp->flags[n2] == 0 && cp->flags[n3] == 0) {
				/*  1-dimensional  */
				static Vector xvec = {1, 0, 0}, yvec = {0, 1, 0};
				VecCross(v1, *vp1, xvec);
				VecCross(v2, *vp1, yvec);
				if (VecLength2(v1) < VecLength2(v2))
					v1 = v2;
				VecCross(v2, *vp1, v1);
				if (NormalizeVec(&v1, &v1) || NormalizeVec(&v2, &v2))
					return -1;   /*  Non-regular transform  */
			} else if (cp->flags[n2] == 0) {
				v2 = *vp3;
				VecCross(v1, v2, *vp1);
				if (NormalizeVec(&v1, &v1))
					return -1;  /*  Non-regular transform  */
			} else {
				v1 = *vp2;
				VecCross(v2, *vp1, v1);
				if (NormalizeVec(&v2, &v2))
					return -1;  /*  Non-regular transform  */
			}
			cp->tr[n2*3] = v1.x;
			cp->tr[n2*3+1] = v1.y;
			cp->tr[n2*3+2] = v1.z;
			cp->tr[n3*3] = v2.x;
			cp->tr[n3*3+1] = v2.y;
			cp->tr[n3*3+2] = v2.z;
		} else {
			VecCross(v1, *vp1, *vp2);
			if (fabs(VecDot(v1, *vp3)) < 1e-7)
				return -1;  /*  Non-regular transform  */
			cp->tr[n2*3] = vp2->x;
			cp->tr[n2*3+1] = vp2->y;
			cp->tr[n2*3+2] = vp2->z;
			cp->tr[n3*3] = vp3->x;
			cp->tr[n3*3+1] = vp3->y;
			cp->tr[n3*3+2] = vp3->z;
		}
	}
	if (TransformInvert(cp->rtr, cp->tr))
		return -1;  /*  Non-regular transform  */

	/*  Calculate the reciprocal cell parameters  */
	cp->rcell[0] = sqrt(cp->rtr[0] * cp->rtr[0] + cp->rtr[3] * cp->rtr[3] + cp->rtr[6] * cp->rtr[6]);
	cp->rcell[1] = sqrt(cp->rtr[1] * cp->rtr[1] + cp->rtr[4] * cp->rtr[4] + cp->rtr[7] * cp->rtr[7]);
	cp->rcell[2] = sqrt(cp->rtr[2] * cp->rtr[2] + cp->rtr[5] * cp->rtr[5] + cp->rtr[8] * cp->rtr[8]);
	cp->rcell[3] = acos((cp->rtr[1] * cp->rtr[2] + cp->rtr[4] * cp->rtr[5] + cp->rtr[7] * cp->rtr[8]) / (cp->rcell[1] * cp->rcell[2])) * kRad2Deg;
	cp->rcell[4] = acos((cp->rtr[2] * cp->rtr[0] + cp->rtr[5] * cp->rtr[3] + cp->rtr[8] * cp->rtr[6]) / (cp->rcell[2] * cp->rcell[0])) * kRad2Deg;
	cp->rcell[5] = acos((cp->rtr[0] * cp->rtr[1] + cp->rtr[3] * cp->rtr[4] + cp->rtr[6] * cp->rtr[7]) / (cp->rcell[0] * cp->rcell[1])) * kRad2Deg;
	
	if (calc_abc) {
		/*  Calculate a, b, c, alpha, beta, gamma  */
		cp->cell[0] = sqrt(cp->tr[0] * cp->tr[0] + cp->tr[1] * cp->tr[1] + cp->tr[2] * cp->tr[2]);
		cp->cell[1] = sqrt(cp->tr[3] * cp->tr[3] + cp->tr[4] * cp->tr[4] + cp->tr[5] * cp->tr[5]);
		cp->cell[2] = sqrt(cp->tr[6] * cp->tr[6] + cp->tr[7] * cp->tr[7] + cp->tr[8] * cp->tr[8]);
		cp->cell[3] = acos((cp->tr[3] * cp->tr[6] + cp->tr[4] * cp->tr[7] + cp->tr[5] * cp->tr[8]) / (cp->cell[1] * cp->cell[2])) * kRad2Deg;
		cp->cell[4] = acos((cp->tr[6] * cp->tr[0] + cp->tr[7] * cp->tr[1] + cp->tr[8] * cp->tr[2]) / (cp->cell[2] * cp->cell[0])) * kRad2Deg;
		cp->cell[5] = acos((cp->tr[0] * cp->tr[3] + cp->tr[1] * cp->tr[4] + cp->tr[2] * cp->tr[5]) / (cp->cell[0] * cp->cell[1])) * kRad2Deg;
	}
	return 0;
}

void
MoleculeSetCell(Molecule *mp, Double a, Double b, Double c, Double alpha, Double beta, Double gamma, int convertCoordinates)
{
	XtalCell *cp;
	int i;
	Atom *ap;
	Transform cmat;
	if (mp == NULL)
		return;
	__MoleculeLock(mp);
	memset(&cmat, 0, sizeof(Transform));
	if (mp->cell != NULL)
		memmove(&cmat, &(mp->cell->rtr), sizeof(Transform));
	else
		memmove(&cmat, &gIdentityTransform, sizeof(Transform));
	if (a == 0.0) {
		if (mp->cell != NULL) {
			free(mp->cell);
			mp->needsMDRebuild = 1;
		}
		mp->cell = NULL;
	} else {
		cp = mp->cell;
		if (cp == NULL) {
			cp = (XtalCell *)calloc(sizeof(XtalCell), 1);
			if (cp == NULL)
				Panic("Low memory during setting cell parameters");
			mp->cell = cp;
			mp->needsMDRebuild = 1;
		}
		/*  alpha, beta, gamma are in degree  */
		cp->cell[0] = a;
		cp->cell[1] = b;
		cp->cell[2] = c;
		cp->cell[3] = alpha;
		cp->cell[4] = beta;
		cp->cell[5] = gamma;
		if (fabs(alpha - 90) < 0.0001 && fabs(beta - 90) < 0.0001 && fabs(gamma - 90) > 0.0001) {
			/*  c unique (hexagonal etc.)  */
			Double cosa, cosb, sinb, cosg;
			cosa = cos(alpha * kDeg2Rad);
			cosb = cos(beta * kDeg2Rad);
			sinb = sin(beta * kDeg2Rad);
			cosg = cos(gamma * kDeg2Rad);
			cp->axes[0].x = a * sinb;
			cp->axes[0].y = 0;
			cp->axes[0].z = a * cosb;
			cp->axes[1].x = b * (cosg - cosa * cosb) / sinb;
			cp->axes[1].z = b * cosa;
			cp->axes[1].y = sqrt(b * b - cp->axes[1].x * cp->axes[1].x - cp->axes[1].z * cp->axes[1].z);
			cp->axes[2].x = 0;
			cp->axes[2].y = 0;
			cp->axes[2].z = c;
		} else {
			/*  b unique  */
			Double cosg, sing, cosa, cosb;
			cosa = cos(alpha * kDeg2Rad);
			cosb = cos(beta * kDeg2Rad);
			cosg = cos(gamma * kDeg2Rad);
			sing = sin(gamma * kDeg2Rad);
			cp->axes[0].x = a * sing;
			cp->axes[0].y = a * cosg;
			cp->axes[0].z = 0;
			cp->axes[1].x = 0;
			cp->axes[1].y = b;
			cp->axes[1].z = 0;
			cp->axes[2].x = c * (cosb - cosa * cosg) / sing;
			cp->axes[2].y = c * cosa;
			cp->axes[2].z = sqrt(c * c - cp->axes[2].x * cp->axes[2].x - cp->axes[2].y * cp->axes[2].y);
		}
		cp->origin.x = cp->origin.y = cp->origin.z = 0.0;
		cp->flags[0] = cp->flags[1] = cp->flags[2] = 1;
		MoleculeCalculateCellFromAxes(cp, 0);
		TransformMul(cmat, cp->tr, cmat);
	}
	
	/*  Update the coordinates (if requested)  */
	if (convertCoordinates) {
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			TransformVec(&(ap->r), cmat, &(ap->r));
		}
	}
	
	/*  Update the anisotropic parameters  */
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		Aniso *anp = ap->aniso;
		if (anp != NULL) {
			MoleculeSetAniso(mp, i, 0, anp->bij[0], anp->bij[1], anp->bij[2], anp->bij[3], anp->bij[4], anp->bij[5], anp->bsig);
		}
	}
	__MoleculeUnlock(mp);
	sMoleculeNotifyChangeAppearance(mp);
}

void
MoleculeSetAniso(Molecule *mp, int n1, int type, Double x11, Double x22, Double x33, Double x12, Double x13, Double x23, const Double *sigmaptr)
{
	Double d, dx;
	int u = 0;
	const Double log2 = 0.693147180559945;
	const Double pi22 = 19.7392088021787;  /* 2*pi**2 */
	Transform m1, m2;
	Aniso *anp;
	XtalCell *cp;
	Vector axis[3];
	Double val[3];
	if (mp == NULL || n1 < 0 || n1 >= mp->natoms)
		return;
	anp = mp->atoms[n1].aniso;
	__MoleculeLock(mp);
	if (anp == NULL) {
		anp = (Aniso *)calloc(sizeof(Aniso), 1);
		if (anp == NULL) {
			__MoleculeUnlock(mp);
			Panic("Low memory during setting anisotropic atom parameters");
		}
		mp->atoms[n1].aniso = anp;
	}
	switch (type) {
		case 1: d = 1; dx = 0.5; break;
		case 2: d = log2; dx = log2; break;
		case 3: d = log2; dx = log2 * 0.5; break;
		case 4: u = 1; d = 0.25; dx = 0.25; break;
		case 5: u = 1; d = 0.25; dx = 0.125; break;
		case 8: u = 1; d = pi22; dx = pi22; break;
		case 9: u = 1; d = pi22; dx = pi22 * 0.5; break;
		case 10: d = pi22; dx = pi22; break;
		default: d = dx = 1; break;
	}
	anp->bij[0] = x11 * d;
	anp->bij[1] = x22 * d;
	anp->bij[2] = x33 * d;
	anp->bij[3] = x12 * dx;
	anp->bij[4] = x13 * dx;
	anp->bij[5] = x23 * dx;
	if (sigmaptr != NULL) {
		anp->has_bsig = 1;
		anp->bsig[0] = sigmaptr[0] * d;
		anp->bsig[1] = sigmaptr[1] * d;
		anp->bsig[2] = sigmaptr[2] * d;
		anp->bsig[3] = sigmaptr[3] * dx;
		anp->bsig[4] = sigmaptr[4] * dx;
		anp->bsig[5] = sigmaptr[5] * dx;
	} else {
		anp->has_bsig = 0;
		anp->bsig[0] = anp->bsig[1] = anp->bsig[2] = anp->bsig[3] = anp->bsig[4] = anp->bsig[5] = 0.0;
	}
	cp = mp->cell;
	if (cp != NULL && u == 1) {
		anp->bij[0] *= cp->rcell[0] * cp->rcell[0];
		anp->bij[1] *= cp->rcell[1] * cp->rcell[1];
		anp->bij[2] *= cp->rcell[2] * cp->rcell[2];
		anp->bij[3] *= cp->rcell[0] * cp->rcell[1]; /* * cos(cp->rcell[5] * kDeg2Rad); */
		anp->bij[4] *= cp->rcell[2] * cp->rcell[0]; /* * cos(cp->rcell[3] * kDeg2Rad); */
		anp->bij[5] *= cp->rcell[1] * cp->rcell[2]; /* * cos(cp->rcell[4] * kDeg2Rad); */
		if (sigmaptr != NULL) {
			anp->bsig[0] *= cp->rcell[0] * cp->rcell[0];
			anp->bsig[1] *= cp->rcell[1] * cp->rcell[1];
			anp->bsig[2] *= cp->rcell[2] * cp->rcell[2];
			anp->bsig[3] *= cp->rcell[0] * cp->rcell[1];
			anp->bsig[4] *= cp->rcell[2] * cp->rcell[0];
			anp->bsig[5] *= cp->rcell[1] * cp->rcell[2];
		}
	}
	
	/*  Calculate the principal axes (in Cartesian coordinates)  */
	/*  The principal axes are the eigenvectors of matrix At(B^-1)A, where
		B is (bij) and A is the reciprocal conversion matrix, i.e. x = Az
		in which x and z are the crystal-space and cartesian coordinates. */
	m1[0] = anp->bij[0] / pi22;
	m1[4] = anp->bij[1] / pi22;
	m1[8] = anp->bij[2] / pi22;
	m1[1] = m1[3] = anp->bij[3] / pi22;
	m1[2] = m1[6] = anp->bij[4] / pi22;
	m1[5] = m1[7] = anp->bij[5] / pi22;
	MatrixInvert(m1, m1);
	if (cp != NULL) {
		memmove(m2, cp->rtr, sizeof(Mat33));
		MatrixMul(m1, m1, m2);
		MatrixTranspose(m2, m2);
		MatrixMul(m1, m2, m1);
	}
	MatrixSymDiagonalize(m1, val, axis);
	for (u = 0; u < 3; u++) {
		if (val[u] < 0) {
			fprintf(stderr, "Non-positive definite thermal parameters for atom %.4s\n", mp->atoms[n1].aname);
			val[u] = 0.001;
		} else {
			val[u] = 1 / sqrt(val[u]);
		}
		anp->pmat[u*3] = axis[u].x * val[u];
		anp->pmat[u*3+1] = axis[u].y * val[u];
		anp->pmat[u*3+2] = axis[u].z * val[u];
	}
	__MoleculeUnlock(mp);
}

/*  Set the anisotropic parameter for atom idx according to the symop. If symop is not alive, nothing is done. */
void
MoleculeSetAnisoBySymop(Molecule *mp, int idx)
{
	Atom *ap, *ap2;
	Transform t1, t2;
	if (mp == NULL || idx < 0 || idx >= mp->natoms)
		return;
	ap = ATOM_AT_INDEX(mp->atoms, idx);
	if (!SYMOP_ALIVE(ap->symop))
		return;
	ap2 = ATOM_AT_INDEX(mp->atoms, ap->symbase);
	if (ap2->aniso == NULL) {
		if (ap->aniso != NULL) {
			free(ap->aniso);
			ap->aniso = NULL;
		}
		return;
	}
	if (ap->aniso == NULL)
		ap->aniso = (Aniso *)calloc(sizeof(Aniso), 1);
	if (ap->symop.sym == 0 || ap->symop.sym >= mp->nsyms) {
		/*  Just copy the aniso parameters  */
		memmove(ap->aniso, ap2->aniso, sizeof(Aniso));
		return;
	}
	memmove(t1, SYMMETRY_AT_INDEX(mp->syms, ap->symop.sym), sizeof(Transform));
	t1[9] = t1[10] = t1[11] = 0.0;
	memset(t2, 0, sizeof(Transform));
	t2[0] = ap2->aniso->bij[0];
	t2[4] = ap2->aniso->bij[1];
	t2[8] = ap2->aniso->bij[2];
	t2[1] = t2[3] = ap2->aniso->bij[3];
	t2[2] = t2[6] = ap2->aniso->bij[4];
	t2[5] = t2[7] = ap2->aniso->bij[5];
	TransformMul(t2, t1, t2);
	TransformInvert(t1, t1);
	TransformMul(t2, t2, t1);
	MoleculeSetAniso(mp, idx, 0, t2[0], t2[4], t2[8], t2[1], t2[2], t2[5], (ap2->aniso->has_bsig ? ap2->aniso->bsig : NULL));
}

int
MoleculeSetPeriodicBox(Molecule *mp, const Vector *ax, const Vector *ay, const Vector *az, const Vector *ao, const char *periodic, int convertCoordinates)
{
	static Vector zeroVec = {0, 0, 0};
	XtalCell b;
	Transform cmat;
	int i, n;
	Atom *ap;
	if (mp == NULL)
		return 0;
	if (mp->cell != NULL)
		memmove(&cmat, &(mp->cell->rtr), sizeof(Transform));
	else
		memmove(&cmat, &gIdentityTransform, sizeof(Transform));
	if (ax == NULL) {
		if (mp->cell != NULL) {
			free(mp->cell);
			mp->needsMDRebuild = 1;
		}
		mp->cell = NULL;
		return 0;
	}	
	memset(&b, 0, sizeof(b));
	b.axes[0] = (ax != NULL ? *ax : zeroVec);
	b.axes[1] = (ay != NULL ? *ay : zeroVec);
	b.axes[2] = (az != NULL ? *az : zeroVec);
	b.origin = *ao;
	memmove(b.flags, periodic, 3);
	if (MoleculeCalculateCellFromAxes(&b, 1) < 0)
		return -1;
	__MoleculeLock(mp);
	if (mp->cell == NULL) {
		mp->needsMDRebuild = 1;
	} else {
		if (mp->cell->has_sigma) {
			/*  Keep the sigma  */
			b.has_sigma = 1;
			memmove(b.cellsigma, mp->cell->cellsigma, sizeof(mp->cell->cellsigma));
		}
		if ((b.flags[0] != mp->cell->flags[0]) || (b.flags[1] != mp->cell->flags[1]) || (b.flags[2] != mp->cell->flags[2])) {
			mp->needsMDRebuild = 1;
		}
		free(mp->cell);
	}
	mp->cell = (XtalCell *)calloc(sizeof(XtalCell), 1);
	if (mp->cell != NULL) {
		memmove(mp->cell, &b, sizeof(XtalCell));
		TransformMul(cmat, b.tr, cmat);
		/*  Update the coordinates (if requested)  */
		if (convertCoordinates) {
			for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
				TransformVec(&(ap->r), cmat, &(ap->r));
			}
		}
		
		/*  Update the anisotropic parameters  */
		for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			Aniso *anp = ap->aniso;
			if (anp != NULL) {
				MoleculeSetAniso(mp, i, 0, anp->bij[0], anp->bij[1], anp->bij[2], anp->bij[3], anp->bij[4], anp->bij[5], anp->bsig);
			}
		}
		n = 0;
	} else n = -2;  /*  Out of memory  */
	__MoleculeUnlock(mp);
	sMoleculeNotifyChangeAppearance(mp);
	return n;
}

#pragma mark ====== Fragment manipulation ======

static void
sMoleculeFragmentSub(Molecule *mp, int idx, IntGroup *result, IntGroup *exatoms)
{
	Atom *ap;
	Int i, *cp, idx2;
	if (exatoms != NULL && IntGroupLookup(exatoms, idx, NULL))
		return;
	IntGroupAdd(result, idx, 1);
	ap = ATOM_AT_INDEX(mp->atoms, idx);
	cp = AtomConnectData(&ap->connect);
	for (i = 0; i < ap->connect.count; i++) {
		idx2 = cp[i];
		if (IntGroupLookup(result, idx2, NULL))
			continue;
		if (ap->anchor != NULL && ATOM_AT_INDEX(mp->atoms, idx2)->anchor != NULL)
			continue;  /*  bond between two pi_anchors is ignored  */
		sMoleculeFragmentSub(mp, idx2, result, exatoms);
	}
	if (ap->anchor != NULL) {
		cp = AtomConnectData(&ap->anchor->connect);
		for (i = 0; i < ap->anchor->connect.count; i++) {
			idx2 = cp[i];
			if (IntGroupLookup(result, idx2, NULL))
				continue;
			sMoleculeFragmentSub(mp, idx2, result, exatoms);
		}
	}
}

/*  The molecular fragment (= interconnected atoms) containing the atom n1 and
    not containing the atoms in exatoms  */
IntGroup *
MoleculeFragmentExcludingAtomGroup(Molecule *mp, int n1, IntGroup *exatoms)
{
	IntGroup *result;
	if (mp == NULL || mp->natoms == 0 || n1 < 0 || n1 >= mp->natoms)
		return NULL;
	result = IntGroupNew();
	sMoleculeFragmentSub(mp, n1, result, exatoms);
	return result;
}

/*  The molecular fragment (= interconnected atoms) containing the atom n1 and
    not containing the atoms n2, n3, ... (terminated by -1)  */
IntGroup *
MoleculeFragmentExcludingAtoms(Molecule *mp, int n1, int argc, int *argv)
{
	int i;
	IntGroup *exatoms, *result;
	if (mp == NULL || mp->natoms == 0 || n1 < 0 || n1 >= mp->natoms)
		return NULL;
	exatoms = IntGroupNew();
	for (i = 0; i < argc; i++)
		IntGroupAdd(exatoms, argv[i], 1);
	result = IntGroupNew();
	sMoleculeFragmentSub(mp, n1, result, exatoms);
	IntGroupRelease(exatoms);
	return result;
}

/*  The molecular fragment (= interconnected atoms) containing the atoms in inatoms and
    not containing the atoms in exatoms  */
IntGroup *
MoleculeFragmentWithAtomGroups(Molecule *mp, IntGroup *inatoms, IntGroup *exatoms)
{
	IntGroupIterator iter;
	IntGroup *result;
	int i;
	if (mp == NULL || mp->natoms == 0 || inatoms == NULL || IntGroupGetCount(inatoms) == 0)
		return NULL;
	IntGroupIteratorInit(inatoms, &iter);
	result = IntGroupNew();
	while ((i = IntGroupIteratorNext(&iter)) >= 0) {
		sMoleculeFragmentSub(mp, i, result, exatoms);
	}
	IntGroupIteratorRelease(&iter);
	return result;
}

/*  Returns non-zero if the given group is 'detachable' in the molecule, i.e. the
    group is bound to the rest of the molecule via only one bond.
	If the result is true, then the atoms belonging to the (only) bond are returned
	in *n1 and *n2, *n1 being the atom belonging to the fragment. The pointers n1
	and n2 can be NULL, if those informations are not needed.  */
int
MoleculeIsFragmentDetachable(Molecule *mp, IntGroup *group, int *n1, int *n2)
{
	Int i, i1, i2, j, k, bond_count, nval1, nval2, *cp;
	Atom *ap;
	if (mp == NULL || mp->natoms == 0 || group == NULL)
		return 0;  /*  Invalid arguments  */
	bond_count = 0;
	for (i = 0; (i1 = IntGroupGetStartPoint(group, i)) >= 0; i++) {
		i2 = IntGroupGetEndPoint(group, i);
		for (j = i1; j < i2; j++) {
			if (j < 0 || j >= mp->natoms)
				return 0;  /*  Invalid atom group  */
			ap = ATOM_AT_INDEX(mp->atoms, j);
			cp = AtomConnectData(&ap->connect);
			for (k = 0; k < ap->connect.count; k++) {
				if (ap->anchor != NULL && ATOM_AT_INDEX(mp->atoms, cp[k])->anchor != NULL)
					continue;  /*  Ignore bond between two pi_anchors  */
				if (IntGroupLookup(group, cp[k], NULL) == 0) {
					bond_count++;
					nval1 = j;
					nval2 = cp[k];
					if (bond_count > 1)
						return 0;  /*  Too many bonds  */
				}
			}
			if (ap->anchor != NULL) {
				cp = AtomConnectData(&ap->anchor->connect);
				for (k = 0; k < ap->anchor->connect.count; k++) {
					if (IntGroupLookup(group, cp[k], NULL) == 0) {
						bond_count++;
						nval1 = j;
						nval2 = cp[k];
						if (bond_count > 1)
							return 0;  /*  Too many bonds  */
					}
				}					
			}
		}
	}
	if (bond_count == 1) {
		if (n1 != NULL)
			*n1 = nval1;
		if (n2 != NULL)
			*n2 = nval2;
		return 1;
	} else {
		return 0;
	}	
}

/*  Returns non-zero if the given group is 'rotatable' in the molecule. The group
    is said to be 'rotatable' when either of the following conditions are met; (1)
	the group is detachable, or (2) the group consists of two bonded atoms that define
	a detachable fragment. If it is rotatable, the group to rotate is returned to rotGroup
	(either a new IntGroup or 'group' with incremented reference count; thus the caller
	is responsible for releasing the returned value).  */
int
MoleculeIsFragmentRotatable(Molecule *mp, IntGroup *group, int *n1, int *n2, IntGroup **rotGroup)
{
	int i1, i2;
	if (MoleculeIsFragmentDetachable(mp, group, n1, n2)) {
		if (rotGroup != NULL) {
			IntGroupRetain(group);
			*rotGroup = group;
		}
		return 1;
	}
	if (group != NULL && IntGroupGetCount(group) == 2) {
		i1 = IntGroupGetNthPoint(group, 0);
		i2 = IntGroupGetNthPoint(group, 1);
		if (MoleculeAreAtomsConnected(mp, i1, i2)) {
			IntGroup *frag = MoleculeFragmentExcludingAtoms(mp, i2, 1, &i1);
			if (frag == NULL)
				return 0;
			i1 = MoleculeIsFragmentDetachable(mp, frag, n1, n2);
			if (i1 == 0) {
				IntGroupRelease(frag);
				if (rotGroup != NULL)
					*rotGroup = NULL;
				return 0;
			}
			if (rotGroup != NULL)
				*rotGroup = frag;
			else if (frag != NULL)
				IntGroupRelease(frag);
			return i1;
		}
	}
	return 0;
}

#pragma mark ====== Multiple frame ======

int
MoleculeGetNumberOfFrames(Molecule *mp)
{
	if (mp == NULL)
		return 0;
	if (mp->nframes <= 0) {
		/*  Recalculate  */
		int i, n;
		Atom *ap;
		for (i = n = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ap->nframes > n)
				n = ap->nframes;
		}
		if (n == 0)
			n = 1;
		mp->nframes = n;
	}
	return mp->nframes;
}

int
MoleculeInsertFrames(Molecule *mp, IntGroup *group, const Vector *inFrame, const Vector *inFrameCell)
{
	int i, j, count, n_new, n_old, natoms, exframes, last_inserted;
	Vector *tempv, *vp;
	Atom *ap;
	MolProp *prp;
	Double *dp;
	
	if (mp == NULL || (natoms = mp->natoms) == 0 || (count = IntGroupGetCount(group)) <= 0)
		return -1;

	n_old = MoleculeGetNumberOfFrames(mp);
	n_new = n_old + count;
	last_inserted = IntGroupGetNthPoint(group, count - 1);
	if (n_new <= last_inserted) {
		exframes = last_inserted - n_new + 1;  /*  number of extra frames that will be silently inserted  */
		n_new += exframes;
	} else exframes = 0;

	tempv = (Vector *)malloc(sizeof(Vector) * n_new * 4);  /*  "*4" for handling cells  */
	if (tempv == NULL)
		return -1;

	__MoleculeLock(mp);

	/*  Copy back the current coordinates  */
	/*  No change in the current coordinates, but the frame buffer is updated  */
	MoleculeSelectFrame(mp, mp->cframe, 1); 
	
	/*  Expand ap->frames for all atoms  */
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		Int n = ap->nframes;
		AssignArray(&ap->frames, &ap->nframes, sizeof(Vector), n_new - 1, NULL);
		if (ap->frames == NULL) {
			__MoleculeUnlock(mp);
			return -1;
		}
		for (j = n; j < n_new; j++)
			ap->frames[j] = ap->r;
	}
	if (mp->cell != NULL) {
		j = mp->nframe_cells;
		AssignArray(&mp->frame_cells, &mp->nframe_cells, sizeof(Vector) * 4, n_new - 1, NULL);
		for (i = j; i < n_new; i++) {
			/*  Set the current cell parameters to the expanded frames  */
			mp->frame_cells[i * 4] = mp->cell->axes[0];
			mp->frame_cells[i * 4 + 1] = mp->cell->axes[1];
			mp->frame_cells[i * 4 + 2] = mp->cell->axes[2];
			mp->frame_cells[i * 4 + 3] = mp->cell->origin;
		}
	}
	
	/*  Expand propvals for all properties  */
	for (i = 0, prp = mp->molprops; i < mp->nmolprops; i++, prp++) {
		dp = (Double *)realloc(prp->propvals, sizeof(Double) * n_new);
		if (dp == NULL) {
			__MoleculeUnlock(mp);
			return -1;
		}
		for (j = n_old; j < n_new; j++)
			dp[j] = 0.0;
		prp->propvals = dp;
	}
	
	/*  group = [n0..n1-1, n2..n3-1, ...]  */
	/*  s = t = 0,  */
	/*  tempv[0..n0-1] <- ap[0..n0-1], s += n0,
	    tempv[n0..n1-1] <- inFrame[0..(n1-n0-1)], t += n1-n0,
		tempv[n1..n2-1] <- ap[s..s+(n2-n1-1)], s += n2-n1,
		tempv[n2..n3-1] <- inFrame[t..t+(n3-n2-1)], t += n3-n2,
		...
		tempv[nl..n_new-1] <- ap[s..s+(n_new-nl-1)], s += n_new-nl
		At last, s will become n_old and t will become count.  */
	for (i = 0, ap = mp->atoms, prp = mp->molprops; i <= mp->natoms + mp->nmolprops; i++) {
		int s, t, ns, ne, mult;
		Vector cr;
		ne = s = t = 0;
		if (i == mp->natoms) {
			if (mp->cell == NULL || mp->frame_cells == NULL)
				continue;
			vp = mp->frame_cells;
			mult = 4;
		} else if (i < mp->natoms) {
			cr = ap->r;
			vp = ap->frames;
			mult = 1;
		} else {
			dp = prp->propvals;
		}
		for (j = 0; (ns = IntGroupGetStartPoint(group, j)) >= 0; j++) {
			if (ns > ne) {
				if (i <= mp->natoms)
					memmove(tempv + ne * mult, vp + s * mult, sizeof(Vector) * mult * (ns - ne));
				else
					memmove((Double *)tempv + ne, dp + s, sizeof(Double) * (ns - ne)); 
				s += ns - ne;
			}
			ne = IntGroupGetEndPoint(group, j);
			while (ns < ne) {
				if (i == mp->natoms) {
					if (inFrameCell != NULL) {
						tempv[ns * 4] = inFrameCell[t * 4];
						tempv[ns * 4 + 1] = inFrameCell[t * 4 + 1];
						tempv[ns * 4 + 2] = inFrameCell[t * 4 + 2];
						tempv[ns * 4 + 3] = inFrameCell[t * 4 + 3];
					} else {
						tempv[ns * 4] = mp->cell->axes[0];
						tempv[ns * 4 + 1] = mp->cell->axes[1];
						tempv[ns * 4 + 2] = mp->cell->axes[2];
						tempv[ns * 4 + 3] = mp->cell->origin;
					}
				} else if (i < mp->natoms) {
					if (inFrame != NULL)
						tempv[ns] = inFrame[natoms * t + i];
					else
						tempv[ns] = cr;
				} else {
					((Double *)tempv)[ns] = 0.0;
				}
				t++;
				ns++;
			}
		}
		if (n_new > ne) {
			if (i <= mp->natoms)
				memmove(tempv + ne * mult, vp + s * mult, sizeof(Vector) * mult * (n_new - ne));
			else
				memmove((Double *)tempv + ne, dp + s, sizeof(Double) * (n_new - ne));
			s += n_new - ne;
		}
		if (i < mp->natoms)
			ap->nframes = n_new;
		if (i <= mp->natoms) {
			memmove(vp, tempv, sizeof(Vector) * mult * n_new);
			if (i < mp->natoms) {
				ap->nframes = n_new;
				ap = ATOM_NEXT(ap);
			}
		} else {
			memmove(dp, (Double *)tempv, sizeof(Double) * n_new);
			prp++;
		}
	}
	free(tempv);
	mp->nframes = n_new;
	MoleculeSelectFrame(mp, last_inserted, 0);
	MoleculeIncrementModifyCount(mp);
	__MoleculeUnlock(mp);
	return count;
}

int
MoleculeRemoveFrames(Molecule *mp, IntGroup *inGroup, Vector *outFrame, Vector *outFrameCell)
{
	int i, count, n_new, n_old, natoms, nframes, old_count, new_cframe;
	Vector *tempv, *vp;
	Atom *ap;
	MolProp *prp;
	IntGroup *group, *group2;

	if (mp == NULL || (natoms = mp->natoms) == 0 || (count = IntGroupGetCount(inGroup)) <= 0)
		return -1;

	/*  outFrame[] should have enough size for Vector * natoms * group.count  */
	memset(outFrame, 0, sizeof(Vector) * natoms * count);
	if (mp->cell != NULL && mp->frame_cells != NULL)
		memset(outFrameCell, 0, sizeof(Vector) * 4 * count);

	n_old = MoleculeGetNumberOfFrames(mp);
	if (n_old == 1)
		return -2;  /*  Cannot delete last frame  */

	group = IntGroupNew();
	group2 = IntGroupNewWithPoints(0, n_old, -1);
	IntGroupIntersect(inGroup, group2, group);
	IntGroupRelease(group2);
	count = IntGroupGetCount(group);
	n_new = n_old - count;
	if (n_new < 1) {
		IntGroupRelease(group);
		return -2;  /*  Trying to delete too many frames  */
	}
	tempv = (Vector *)malloc(sizeof(Vector) * n_old * 4);  /*  "*4" for handling cells  */
	if (tempv == NULL) {
		IntGroupRelease(group);
		return -1;
	}

	__MoleculeLock(mp);

	/*  Copy back the current coordinates  */
	/*  No change in the current coordinates, but the frame buffer is updated  */
	MoleculeSelectFrame(mp, mp->cframe, 1); 

	/*  Determine which frame should be selected after removal is completed  */
	{
		int n1;
		if (IntGroupLookup(group, mp->cframe, &i)) {
			/*  cframe will be removed  */
			n1 = IntGroupGetStartPoint(group, i) - 1;
			if (n1 < 0)
				n1 = IntGroupGetEndPoint(group, i);
		} else n1 = mp->cframe;
		/*  Change to that frame  */
		MoleculeSelectFrame(mp, n1, 0);
		group2 = IntGroupNewFromIntGroup(group);
		IntGroupReverse(group2, 0, n_old);
		new_cframe = IntGroupLookupPoint(group2, n1);
		if (new_cframe < 0)
			return -3;  /*  This cannot happen  */
		IntGroupRelease(group2);
	}

	/*  group = [n0..n1-1, n2..n3-1, ...]  */
	/*  s = t = 0, */
	/*  tempv[0..n0-1] -> ap[0..n0-1], s += n0,
	    tempv[n0..n1-1] -> outFrame[0..(n1-n0-1)], t += n1-n0,
		tempv[n1..n2-1] -> ap[s..s+(n2-n1-1)], s += n2-n1,
		tempv[n2..n3-1] -> outFrame[t..t+(n3-n2-1)], t += n3-n2,
		...
		tempv[nl..n_old-1] -> ap[s..s+(n_old-nl-1)], s += n_old-nl
		At last, s will become n_new and t will become count.  */
	nframes = 0;
	for (i = 0, ap = mp->atoms, prp = mp->molprops; i <= mp->natoms + mp->nmolprops; i++) {
		int s, t, j, ns, ne;
		int mult;
		/*  if i == mp->natoms, mp->frame_cells is handled  */
		if (i == mp->natoms) {
			if (mp->cell == NULL || mp->frame_cells == NULL)
				continue;
			mult = 4 * sizeof(Vector);
			vp = mp->frame_cells;
			old_count = n_old;
		} else if (i < mp->natoms) {
			mult = sizeof(Vector);
			vp = ap->frames;
			if (vp == NULL) {
				NewArray(&ap->frames, &ap->nframes, sizeof(Vector), n_old);
				if (ap->frames == NULL) {
					__MoleculeUnlock(mp);
					return -1;
				}
				vp = ap->frames;
			}
			old_count = ap->nframes;
		} else {
			mult = sizeof(Double);
			vp = (Vector *)prp->propvals;
			old_count = n_old;
		}

		/*  Copy vp to tempv  */
		memset(tempv, 0, mult * n_old);
		memmove(tempv, vp, mult * (old_count > n_old ? n_old : old_count));
		ne = ns = s = t = 0;
		for (j = 0; ns < n_old && (ns = IntGroupGetStartPoint(group, j)) >= 0; j++) {
			if (ns > n_old)
				ns = n_old;
			if (ns > ne) {
				memmove((char *)vp + s * mult, (char *)tempv + ne * mult, mult * (ns - ne));
				s += ns - ne;
			}
			ne = IntGroupGetEndPoint(group, j);
			if (ne > n_old)
				ne = n_old;
			while (ns < ne) {
				if (i < mp->natoms)
					outFrame[natoms * t + i] = tempv[ns];
				else if (i == mp->natoms) {
					if (outFrameCell != NULL) {
						outFrameCell[t * 4] = tempv[ns * 4];
						outFrameCell[t * 4 + 1] = tempv[ns * 4 + 1];
						outFrameCell[t * 4 + 2] = tempv[ns * 4 + 2];
						outFrameCell[t * 4 + 3] = tempv[ns * 4 + 3];
					}
				}
				t++;
				ns++;
			}
		}
		if (n_old > ne) {
			memmove((char *)vp + s * mult, (char *)tempv + ne * mult, mult * (n_old - ne));
			s += n_old - ne;
		}
		if (i < mp->natoms)
			ap->nframes = s;
		if (nframes < s)
			nframes = s;
		if (s <= 1) {
			if (i < mp->natoms) {
				free(ap->frames);
				ap->frames = NULL;
				ap->nframes = 0;
			} else if (i == mp->natoms) {
				free(mp->frame_cells);
				mp->frame_cells = NULL;
				mp->nframe_cells = 0;
			} else {
				prp->propvals = (Double *)realloc(prp->propvals, sizeof(Double));
			}
		} else {
			if (i < mp->natoms) {
				AssignArray(&ap->frames, &ap->nframes, sizeof(Vector), s - 1, NULL);
				ap->nframes = s;
			} else if (i == mp->natoms) {
				AssignArray(&mp->frame_cells, &mp->nframe_cells, sizeof(Vector) * 4, s - 1, NULL);
				mp->nframe_cells = s;
			} else {
				prp->propvals = (Double *)realloc(prp->propvals, sizeof(Double) * s);
			}
		}
		if (i < mp->natoms) {
			ap = ATOM_NEXT(ap);
		} else if (i > mp->natoms) {
			prp++;
		}
	}
	free(tempv);
	mp->nframes = nframes;
	
	/*  Select the "last" frame; do not "copy back" the coordinates to the frame table  */
/* 	i = (mp->cframe >= nframes ? nframes - 1 : mp->cframe); */
	MoleculeSelectFrame(mp, new_cframe, 0);

	IntGroupRelease(group);

	MoleculeIncrementModifyCount(mp);
	__MoleculeUnlock(mp);
	return count;
}

int
MoleculeSelectFrame(Molecule *mp, int frame, int copyback)
{
	int i, cframe, nframes, modified;
	Atom *ap;
	cframe = mp->cframe;
	nframes = MoleculeGetNumberOfFrames(mp);
	if (frame == -1)
		frame = mp->cframe;
	if (mp == NULL || mp->natoms == 0 || frame < 0 || frame >= nframes)
		return -1;
	modified = 0;
	__MoleculeLock(mp);
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (copyback && cframe >= 0 && cframe < ap->nframes) {
			/*  Write the current coordinate back to the frame array  */
			ap->frames[cframe] = ap->r;
		}
		if ((frame != cframe || copyback == 0) && frame >= 0 && frame < ap->nframes) {
			/*  Read the coordinate from the frame array  */
			ap->r = ap->frames[frame];
			modified = 1;
		}
	}

	if (mp->cell != NULL && mp->frame_cells != NULL) {
		/*  Write the current cell back to the frame_cells array  */
		if (copyback && cframe >= 0) {
			Vector *vp = (Vector *)AssignArray(&mp->frame_cells, &mp->nframe_cells, sizeof(Vector) * 4, cframe, NULL);
			vp[0] = mp->cell->axes[0];
			vp[1] = mp->cell->axes[1];
			vp[2] = mp->cell->axes[2];
			vp[3] = mp->cell->origin;
		}
		/*  Set the cell from the frame array  */
		if ((frame != cframe || copyback == 0) && frame >= 0 && frame < mp->nframe_cells) {
			MoleculeSetPeriodicBox(mp, &mp->frame_cells[frame * 4], &mp->frame_cells[frame * 4 + 1], &mp->frame_cells[frame * 4 + 2], &mp->frame_cells[frame * 4 + 3], mp->cell->flags, 0);
			modified = 1;
			MoleculeAmendBySymmetry(mp, NULL, NULL, NULL);
		}
	}
	mp->cframe = frame;
	if (modified)
		mp->needsMDCopyCoordinates = 1;
	__MoleculeUnlock(mp);
	sMoleculeNotifyChangeAppearance(mp);
	return frame;
}

/*  If molecule is multi-frame, then flush the current information to the frame buffer.
    Returns the number of frames.  */
int
MoleculeFlushFrames(Molecule *mp)
{
	int nframes = MoleculeGetNumberOfFrames(mp);
	if (nframes > 1)
		MoleculeSelectFrame(mp, mp->cframe, 1);
	return nframes;
}

int
MoleculeReorderFrames(Molecule *mp, const Int *old_idx)
{
	Int *ip, i, j, n, nframes;
	Double *dp;
	Atom *ap;
	MolProp *prp;
	if (mp == NULL || old_idx == NULL)
		return 0;
	nframes = MoleculeGetNumberOfFrames(mp);
	MoleculeFlushFrames(mp);
	ip = (Int *)malloc(sizeof(Int) * nframes);
	if (ip == NULL)
		return -1;  /*  Out of memory  */
	memset(ip, 0, sizeof(Int) * nframes);
	/*  Check the argument  */
	for (i = 0; i < nframes; i++) {
		j = old_idx[i];
		if (j < 0 || j >= nframes || ip[j] != 0) {
			free(ip);
			return -2;  /*  Bad argument  */
		}
		ip[j] = 1;
	}
	free(ip);
	dp = (Double *)malloc(sizeof(Double) * nframes * 12);
	for (i = 0, ap = mp->atoms, prp = mp->molprops; i <= mp->natoms + mp->nmolprops; i++) {
		for (j = 0; j < nframes; j++) {
			n = old_idx[j];
			if (i < mp->natoms) {
				((Vector *)dp)[j] = (n < ap->nframes ? ap->frames[n] : ap->r);
			} else if (i == mp->natoms) {
				if (mp->cell != NULL) {
					if (n < mp->nframe_cells && mp->frame_cells != NULL)
						memmove(dp + j * 12, mp->frame_cells + n * 4, sizeof(Vector) * 4);
					else {
						((Vector *)dp)[j * 4] = mp->cell->axes[0];
						((Vector *)dp)[j * 4] = mp->cell->axes[1];
						((Vector *)dp)[j * 4] = mp->cell->axes[2];
						((Vector *)dp)[j * 4] = mp->cell->origin;
					}
				}
			} else {
				dp[j] = prp->propvals[n];
			}
		}
		for (j = 0; j < nframes; j++) {
			if (i < mp->natoms) {
				if (ap->nframes <= j)
					AssignArray(&ap->frames, &ap->nframes, sizeof(Vector), nframes - 1, NULL);
				ap->frames[j] = ((Vector *)dp)[j];
			} else if (i == mp->natoms) {
				if (mp->cell != NULL) {
					AssignArray(&mp->frame_cells, &mp->nframe_cells, sizeof(Vector) * 4, nframes - 1, NULL);
					memmove(mp->frame_cells + j * 4, dp + j * 12, sizeof(Vector) * 4);
				}
			} else {
				prp->propvals[j] = dp[j];
			}
		}
		if (i < mp->natoms)
			ap = ATOM_NEXT(ap);
		else if (i > mp->natoms)
			prp++;
	}
	free(dp);
	MoleculeSelectFrame(mp, mp->cframe, 0);
	return 0;
}

#pragma mark ====== Molecule Propeties ======

int
MoleculeCreateProperty(Molecule *mp, const char *name)
{
	int i;
	MolProp *prp;
	for (i = 0, prp = mp->molprops; i < mp->nmolprops; i++, prp++) {
		if (strcmp(prp->propname, name) == 0)
			return -(i + 1);
	}
	prp = (MolProp *)calloc(sizeof(MolProp), 1);
	if (prp == NULL)
		return -10000;
	prp->propname = strdup(name);
	if (prp->propname == NULL)
		return -10000;
	i = MoleculeGetNumberOfFrames(mp);
	prp->propvals = (Double *)calloc(sizeof(Double), i);
	if (prp->propvals == NULL)
		return -10000;
	AssignArray(&mp->molprops, &mp->nmolprops, sizeof(MolProp), mp->nmolprops, prp);
	free(prp);
	return mp->nmolprops - 1;
}

int
MoleculeLookUpProperty(Molecule *mp, const char *name)
{
	int i;
	MolProp *prp;
	for (i = 0, prp = mp->molprops; i < mp->nmolprops; i++, prp++) {
		if (strcmp(prp->propname, name) == 0)
			return i;
	}
	return -1;
}

int
MoleculeDeletePropertyAtIndex(Molecule *mp, int idx)
{
	if (idx >= 0 && idx < mp->nmolprops) {
		free(mp->molprops[idx].propname);
		free(mp->molprops[idx].propvals);
		DeleteArray(&mp->molprops, &mp->nmolprops, sizeof(MolProp), idx, 1, NULL);
		return idx;
	}
	return -1;
}

int
MoleculeSetProperty(Molecule *mp, int idx, IntGroup *ig, const Double *values)
{
	IntGroupIterator iter;
	int i, n, nframes;
	if (idx < 0 || idx >= mp->nmolprops)
		return -1;
	IntGroupIteratorInit(ig, &iter);
	nframes = MoleculeGetNumberOfFrames(mp);
	n = 0;
	while ((i = IntGroupIteratorNext(&iter)) >= 0) {
		if (i >= nframes)
			break;
		mp->molprops[idx].propvals[i] = values[n];
		n++;
	}
	IntGroupIteratorRelease(&iter);
	return n;
}

int
MoleculeGetProperty(Molecule *mp, int idx, IntGroup *ig, Double *outValues)
{
	IntGroupIterator iter;
	int i, n, nframes;
	if (idx < 0 || idx >= mp->nmolprops)
		return -1;
	IntGroupIteratorInit(ig, &iter);
	nframes = MoleculeGetNumberOfFrames(mp);
	n = 0;
	while ((i = IntGroupIteratorNext(&iter)) >= 0) {
		if (i >= nframes)
			break;
		outValues[n] = mp->molprops[idx].propvals[i];
		n++;
	}
	IntGroupIteratorRelease(&iter);
	return n;
}

#pragma mark ====== Pi Atoms ======

static inline void
sMoleculeCalculatePiAnchorPosition(Atom *ap, Atom *atoms)
{
	Int *cp, j, n;
	Atom *ap2;
	cp = AtomConnectData(&ap->anchor->connect);
	n = ap->anchor->connect.count;
	VecZero(ap->r);
	for (j = 0; j < n; j++) {
		Double w = ap->anchor->coeffs[j];
		ap2 = ATOM_AT_INDEX(atoms, cp[j]);
		VecScaleInc(ap->r, ap2->r, w);
	}	
}

void
MoleculeUpdatePiAnchorPositions(Molecule *mol)
{
	Int i;
	Atom *ap;
	for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (ap->anchor == NULL)
			continue;
		sMoleculeCalculatePiAnchorPosition(ap, mol->atoms);
	}
}

void
MoleculeCalculatePiAnchorPosition(Molecule *mol, int idx)
{
	Atom *ap;
	if (mol == NULL || idx < 0 || idx >= mol->natoms)
		return;
	ap = ATOM_AT_INDEX(mol->atoms, idx);
	if (ap->anchor == NULL)
		return;
	sMoleculeCalculatePiAnchorPosition(ap, mol->atoms);
}

int
MoleculeSetPiAnchorList(Molecule *mol, Int idx, Int nentries, Int *entries, Double *weights, Int *nUndoActions, struct MolAction ***undoActions)
{
	Atom *ap;
	Int *ip, i, j, n, *np;
	Double d;
	if (mol == NULL || idx < 0 || idx >= mol->natoms || nentries <= 1)
		return -1;  /*  Invalid argument  */
	if (weights != NULL) {
		d = 0.0;
		for (i = 0; i < nentries; i++) {
			if (weights[i] <= 0.0) {
				return 10;  /*  Weights must be positive  */
			}
			d += weights[i];
		}
		d = 1.0 / d;
	} else d = 1.0 / nentries;
	ap = ATOM_AT_INDEX(mol->atoms, idx);
	if (ap->anchor != NULL) {
		/*  Already an anchor: check if bonds/angles/dihedrals have entries related to this anchor  */
		IntGroup *bg, *ag, *dg, *ig;
		Int *ibuf, ibufsize;
		MolAction *act;
		bg = ag = dg = ig = NULL;
		ip = AtomConnectData(&ap->anchor->connect);
		for (i = 0; i < ap->anchor->connect.count; i++) {
			n = ip[i];
			for (j = 0; j < nentries; j++) {
				if (n == entries[j])
					break;
			}
			if (j == nentries) {
				/*  This entry will disappear: if any bond/angle/dihedral has idx-n pair, that should be removed.  */
				for (j = 0, np = mol->bonds; j < mol->nbonds; j++, np += 2) {
					if ((idx == np[0] && n == np[1]) || (idx == np[1] && n == np[0])) {
						if (bg == NULL)
							bg = IntGroupNew();
						IntGroupAdd(bg, j, 1);
					}
				}
				for (j = 0, np = mol->angles; j < mol->nangles; j++, np += 3) {
					if ((idx == np[0] && n == np[1]) || (idx == np[1] && n == np[2]) ||
						(idx == np[1] && n == np[0]) || (idx == np[2] && n == np[1])) {
						if (ag == NULL)
							ag = IntGroupNew();
						IntGroupAdd(ag, j, 1);
					}
				}
				for (j = 0, np = mol->dihedrals; j < mol->ndihedrals; j++, np += 4) {
					if ((idx == np[0] && n == np[1]) || (idx == np[1] && n == np[2]) || (idx == np[2] && n == np[3]) ||
						(idx == np[1] && n == np[0]) || (idx == np[2] && n == np[1]) || (idx == np[3] && n == np[2])) {
						if (dg == NULL)
							dg = IntGroupNew();
						IntGroupAdd(dg, j, 1);
					}
				}
				for (j = 0, np = mol->impropers; j < mol->nimpropers; j++, np += 4) {
					if ((idx == np[0] && n == np[2]) || (idx == np[1] && n == np[2]) || (idx == np[3] && n == np[2]) ||
						(idx == np[2] && n == np[0]) || (idx == np[2] && n == np[1]) || (idx == np[2] && n == np[3])) {
						if (ig == NULL)
							ig = IntGroupNew();
						IntGroupAdd(ig, j, 1);
					}
				}
			}
		}
		ibuf = NULL;
		ibufsize = 0;
		if (ig != NULL) {
			/*  Delete impropers (with undo info) */
			i = IntGroupGetCount(ig);
			AssignArray(&ibuf, &ibufsize, sizeof(Int), i * 4 - 1, NULL);
			MoleculeDeleteImpropers(mol, ibuf, ig);
			if (nUndoActions != NULL && undoActions != NULL) {
				act = MolActionNew(gMolActionAddImpropers, i * 4, ibuf, ig);
				AssignArray(undoActions, nUndoActions, sizeof(MolAction *), *nUndoActions, &act);
			}
			IntGroupRelease(ig);
		}
		if (dg != NULL) {
			/*  Delete dihedrals (with undo info)  */
			i = IntGroupGetCount(dg);
			AssignArray(&ibuf, &ibufsize, sizeof(Int), i * 4 - 1, NULL);
			MoleculeDeleteDihedrals(mol, ibuf, dg);
			if (nUndoActions != NULL && undoActions != NULL) {
				act = MolActionNew(gMolActionAddDihedrals, i * 4, ibuf, dg);
				AssignArray(undoActions, nUndoActions, sizeof(MolAction *), *nUndoActions, &act);
			}
			IntGroupRelease(dg);
		}
		if (ag != NULL) {
			/*  Delete angles (with undo info) */
			i = IntGroupGetCount(ag);
			AssignArray(&ibuf, &ibufsize, sizeof(Int), i * 3 - 1, NULL);
			MoleculeDeleteAngles(mol, ibuf, ag);
			if (nUndoActions != NULL && undoActions != NULL) {
				act = MolActionNew(gMolActionAddAngles, i * 3, ibuf, ag);
				AssignArray(undoActions, nUndoActions, sizeof(MolAction *), *nUndoActions, &act);
			}
			IntGroupRelease(ag);
		}
		if (bg != NULL) {
			/*  Delete bonds (with undo info) */
			i = IntGroupGetCount(bg);
			AssignArray(&ibuf, &ibufsize, sizeof(Int), i * 2 - 1, NULL);
			MoleculeDeleteBonds(mol, ibuf, bg, NULL, NULL);
			if (nUndoActions != NULL && undoActions != NULL) {
				act = MolActionNew(gMolActionAddBondsForUndo, i * 2, ibuf, bg);
				AssignArray(undoActions, nUndoActions, sizeof(MolAction *), *nUndoActions, &act);
			}
			IntGroupRelease(bg);
		}
	} else {
		ap->anchor = (PiAnchor *)calloc(sizeof(PiAnchor), 1);
	}
	AtomConnectResize(&ap->anchor->connect, nentries);
	memmove(AtomConnectData(&ap->anchor->connect), entries, sizeof(Int) * nentries);
	AssignArray(&ap->anchor->coeffs, &ap->anchor->ncoeffs, sizeof(Double), nentries - 1, NULL);
	if (weights != NULL) {
		memmove(ap->anchor->coeffs, weights, sizeof(Double) * nentries);
		for (i = 0; i < nentries; i++)
			ap->anchor->coeffs[i] *= d;   /*  Normalize weight  */
	} else {
		for (i = 0; i < nentries; i++)
			ap->anchor->coeffs[i] = d;
	}
	MoleculeCalculatePiAnchorPosition(mol, idx);
	return 0;
}

#pragma mark ====== MO calculation ======

/*  Calculate an MO value for a single point.  */
/*  Index is the MO number (1-based); 0 denotes "arbitrary vector"  */
/*  tmp is an array of (natoms * 4) atoms, and used to store dr and |dr|^2 for each atom.  */
static Double
sCalcMOPoint(Molecule *mp, const BasisSet *bset, Int index, const Vector *vp, Double *tmp)
{
	ShellInfo *sp;
	PrimInfo *pp;
	Double val, tval, *cnp, *tmpp, *mobasep, *mop;
	Int i, j;
	/*  Cache dr and |dr|^2  */
	if (index == 0)
		index = bset->nmos + 1;
	for (i = 0; i < mp->natoms; i++) {
		Vector r;
		r = ATOM_AT_INDEX(mp->atoms, i)->r;
		tmp[i * 4] = r.x = (vp->x - r.x) * kAngstrom2Bohr;
		tmp[i * 4 + 1] = r.y = (vp->y - r.y) * kAngstrom2Bohr;
		tmp[i * 4 + 2] = r.z = (vp->z - r.z) * kAngstrom2Bohr;
		tmp[i * 4 + 3] = r.x * r.x + r.y * r.y + r.z * r.z;
	}
	/*  Iterate over all shells  */
	val = 0.0;
	mobasep = bset->mo + (index - 1) * bset->ncomps;
	for (i = 0, sp = bset->shells; i < bset->nshells; i++, sp++) {
		pp = bset->priminfos + sp->p_idx;
		cnp = bset->cns + sp->cn_idx;
		if (sp->a_idx >= mp->natoms)
			return 0.0; /*  This may happen when molecule is edited after setting up MO info  */
		tmpp = tmp + sp->a_idx * 4;
		mop = mobasep + sp->m_idx;
		switch (sp->sym) {
			case kGTOType_S: {
				tval = 0;
				for (j = 0; j < sp->nprim; j++) {
					tval += *cnp++ * exp(-pp->A * tmpp[3]);
					pp++;
				}
				val += mop[0] * tval;
				break;
			}
			case kGTOType_P: {
				Double x, y, z;
				x = y = z = 0;
				for (j = 0; j < sp->nprim; j++) {
					tval = exp(-pp->A * tmpp[3]);
					x += *cnp++ * tval;
					y += *cnp++ * tval;
					z += *cnp++ * tval;
					pp++;
				}
				x *= mop[0] * tmpp[0];
				y *= mop[1] * tmpp[1];
				z *= mop[2] * tmpp[2];
				val += x + y + z;
				break;
			}
			case kGTOType_SP: {
				Double t, x, y, z;
				t = x = y = z = 0;
				for (j = 0; j < sp->nprim; j++) {
					tval = exp(-pp->A * tmpp[3]);
					t += *cnp++ * tval;
					x += *cnp++ * tval;
					y += *cnp++ * tval;
					z += *cnp++ * tval;
					pp++;
				}
				t *= mop[0];
				x *= mop[1] * tmpp[0];
				y *= mop[2] * tmpp[1];
				z *= mop[3] * tmpp[2];
				val += t + x + y + z;
				break;
			}
			case kGTOType_D: {
				Double xx, yy, zz, xy, xz, yz;
				xx = yy = zz = xy = xz = yz = 0;
				for (j = 0; j < sp->nprim; j++) {
					tval = exp(-pp->A * tmpp[3]);
					xx += *cnp++ * tval;
					yy += *cnp++ * tval;
					zz += *cnp++ * tval;
					xy += *cnp++ * tval;
					xz += *cnp++ * tval;
					yz += *cnp++ * tval;
					pp++;
				}
				xx *= mop[0] * tmpp[0] * tmpp[0];
				yy *= mop[1] * tmpp[1] * tmpp[1];
				zz *= mop[2] * tmpp[2] * tmpp[2];
				xy *= mop[3] * tmpp[0] * tmpp[1];
				xz *= mop[4] * tmpp[0] * tmpp[2];
				yz *= mop[5] * tmpp[1] * tmpp[2];
				val += xx + yy + zz + xy + xz + yz;
				break;
			}
			case kGTOType_D5: {
				Double d0, d1p, d1n, d2p, d2n;
				d0 = d1p = d1n = d2p = d2n = 0;
				for (j = 0; j < sp->nprim; j++) {
					tval = exp(-pp->A * tmpp[3]);
					d0 += *cnp++ * tval;
					d1p += *cnp++ * tval;
					d1n += *cnp++ * tval;
					d2p += *cnp++ * tval;
					d2n += *cnp++ * tval;
					pp++;
				}
				d0 *= mop[0] * (3 * tmpp[2] * tmpp[2] - tmpp[3]);
				d1p *= mop[1] * tmpp[0] * tmpp[2];
				d1n *= mop[2] * tmpp[1] * tmpp[2];
				d2p *= mop[3] * (tmpp[0] * tmpp[0] - tmpp[1] * tmpp[1]);
				d2n *= mop[4] * tmpp[0] * tmpp[1];
				val += d0 + d1p + d1n + d2p + d2n;
				break;
			}
			/*  TODO: Support F/F7 and G/G9 type orbitals  */
		}
	}
	return val;
}

/*  Calculate one MO. The input vectors are angstrom unit (changed from bohr unit: 20140520)  */
/*  mono is the MO number (1-based); 0 denotes "arbitrary vector" */
int
MoleculeCalcMO(Molecule *mp, Int mono, const Vector *op, const Vector *dxp, const Vector *dyp, const Vector *dzp, Int nx, Int ny, Int nz, int (*callback)(double progress, void *ref), void *ref)
{
	int ix, iy, iz, n, nn;
	Cube *cp;
	Double *tmp;
	if (mp == NULL || mp->bset == NULL)
		return -1;
	if (mp->bset->cns == NULL) {
		if (sSetupGaussianCoefficients(mp->bset) != 0)
			return -1;
	}
	if (mp->bset->natoms_bs > mp->natoms)
		return -3;  /*  Number of atoms is smaller than expected (internal error)  */
	
	cp = (Cube *)calloc(sizeof(Cube), 1);
	if (cp == NULL) {
		return -1;
	}
	cp->dp = (Double *)calloc(sizeof(Double), nx * ny * nz);
	if (cp->dp == NULL) {
		free(cp);
		return -1;
	}
	cp->idn = mono;
	cp->origin = *op;
	cp->dx = *dxp;
	cp->dy = *dyp;
	cp->dz = *dzp;
	cp->nx = nx;
	cp->ny = ny;
	cp->nz = nz;
	
	/*  TODO: use multithread  */
	tmp = (Double *)calloc(sizeof(Double), mp->bset->natoms_bs * 4);
	if (tmp == NULL) {
		free(cp->dp);
		free(cp);
		return -1;
	}
	n = nn = 0;
	for (ix = 0; ix < nx; ix++) {
		Vector p;
		for (iy = 0; iy < ny; iy++) {
			for (iz = 0; iz < nz; iz++) {
				p.x = op->x + dxp->x * ix + dyp->x * iy + dzp->x * iz;
				p.y = op->y + dxp->y * ix + dyp->y * iy + dzp->y * iz;
				p.z = op->z + dxp->z * ix + dyp->z * iy + dzp->z * iz;
				cp->dp[n++] = sCalcMOPoint(mp, mp->bset, mono, &p, tmp);
			}
			if (callback != NULL && n - nn > 100) {
				nn = n;
				if ((*callback)((double)n / ((double)nx * ny * nz), ref) != 0) {
					free(cp->dp);
					free(cp);
					free(tmp);
					return -2;  /*  User interrupt  */
				}
			}
		}
	}
	free(tmp);

	AssignArray(&(mp->bset->cubes), &(mp->bset->ncubes), sizeof(Cube *), mp->bset->ncubes, &cp);
	return mp->bset->ncubes - 1;
}

/*  Output values are in angstrom unit (changed from bohr unit: 20140520)  */
int
MoleculeGetDefaultMOGrid(Molecule *mp, Int npoints, Vector *op, Vector *xp, Vector *yp, Vector *zp, Int *nx, Int *ny, Int *nz)
{
	int i;
	Vector rmin, rmax, r;
	Double dr, dx, dy, dz;
	Atom *ap;
	if (mp == NULL || mp->bset == NULL)
		return -1;
	if (npoints <= 0)
		npoints = 1000000;
	rmin.x = rmin.y = rmin.z = 1e10;
	rmax.x = rmax.y = rmax.z = -1e10;
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		dr = RadiusForAtomicNumber(ap->atomicNumber);
		r = ap->r;
		if (dr == 0.0)
			dr = 1.0;
		dr = dr * 3.0 + 2.0;
		if (rmin.x > r.x - dr)
			rmin.x = r.x - dr;
		if (rmin.y > r.y - dr)
			rmin.y = r.y - dr;
		if (rmin.z > r.z - dr)
			rmin.z = r.z - dr;
		if (rmax.x < r.x + dr)
			rmax.x = r.x + dr;
		if (rmax.y < r.y + dr)
			rmax.y = r.y + dr;
		if (rmax.z < r.z + dr)
			rmax.z = r.z + dr;
	}
	dx = rmax.x - rmin.x;
	dy = rmax.y - rmin.y;
	dz = rmax.z - rmin.z;
	dr = pow(dx * dy * dz / npoints, 1.0/3.0);
	*nx = floor(dx / dr + 0.5);
	*ny = floor(dy / dr + 0.5);
	*nz = floor(dz / dr + 0.5);
	if (*nx == 0)
		*nx = 1;
	if (*ny == 0)
		*ny = 1;
	if (*nz == 0)
		*nz = 1;
	*op = rmin;
	xp->x = yp->y = zp->z = dr;
	xp->y = xp->z = yp->x = yp->z = zp->x = zp->y = 0.0;
	return 0;
}

const Cube *
MoleculeGetCubeAtIndex(Molecule *mp, Int index)
{
	if (mp == NULL || mp->bset == NULL || index < 0 || index >= mp->bset->ncubes)
		return NULL;
	return mp->bset->cubes[index];
}

int
MoleculeLookUpCubeWithMONumber(Molecule *mp, Int mono)
{
	int i;
	if (mp == NULL || mp->bset == NULL)
		return -1;
	for (i = 0; i < mp->bset->ncubes; i++) {
		if (mp->bset->cubes[i]->idn == mono)
			return i;
	}
	return -1;
}

int
MoleculeClearCubeAtIndex(Molecule *mp, Int index)
{
	int n;
	if (mp == NULL || mp->bset == NULL || index < 0 || index >= (n = mp->bset->ncubes))
		return -1;
	CubeRelease(mp->bset->cubes[index]);
	if (index < n - 1)
		memmove(mp->bset->cubes + index, mp->bset->cubes + index + 1, sizeof(Cube *) * (n - index - 1));
	if (--(mp->bset->ncubes) == 0) {
		free(mp->bset->cubes);
		mp->bset->cubes = NULL;
	}
	return mp->bset->ncubes;
}

int
MoleculeOutputCube(Molecule *mp, Int index, const char *fname, const char *comment)
{
	const Cube *cp;
	int i, j, k, n;
	FILE *fp;
	if (mp == NULL || mp->bset == NULL)
		return -1;  /*  Molecule or the basis set information is empty  */
	cp = MoleculeGetCubeAtIndex(mp, index);
	if (cp == NULL)
		return -2;  /*  MO not yet calculated  */
	fp = fopen(fname, "wb");
	if (fp == NULL)
		return -3;  /*  Cannot create file  */

	/*  Comment lines  */
	fprintf(fp, "%s MO=%d\n", comment, cp->idn);
	fprintf(fp, " MO coefficients\n");
	
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", -(mp->bset->natoms_bs),
			cp->origin.x * kAngstrom2Bohr, cp->origin.y * kAngstrom2Bohr, cp->origin.z * kAngstrom2Bohr);
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", cp->nx,
			cp->dx.x * kAngstrom2Bohr, cp->dx.y * kAngstrom2Bohr, cp->dx.z * kAngstrom2Bohr);
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", cp->ny,
			cp->dy.x * kAngstrom2Bohr, cp->dy.y * kAngstrom2Bohr, cp->dy.z * kAngstrom2Bohr);
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", cp->nz,
			cp->dz.x * kAngstrom2Bohr, cp->dz.y * kAngstrom2Bohr, cp->dz.z * kAngstrom2Bohr);
	
	/*  Atomic information  */
	for (i = 0; i < mp->natoms; i++) {
		Atom *ap = ATOM_AT_INDEX(mp->atoms, i);
		/*  The second number should actually be the effective charge  */
		fprintf(fp, "%5d %11.6f %11.6f %11.6f %11.6f\n", ap->atomicNumber, (double)ap->atomicNumber,
				ap->r.x * kAngstrom2Bohr, ap->r.y * kAngstrom2Bohr, ap->r.z * kAngstrom2Bohr);
	}
	fprintf(fp, "%5d%5d\n", 1, 1);
	
	/*  3D data  */
	for (i = n = 0; i < cp->nx; i++) {
		for (j = 0; j < cp->ny; j++) {
			for (k = 0; k < cp->nz; k++) {
				/*  On Windows, the "%e" format writes the exponent in 3 digits, but
				    this is not standard. So we avoid using %e  */
				Double d = cp->dp[n++];
				int exponent = (int)floor(log10(fabs(d)));
				Double base = d * pow(10, -1.0 * exponent);
				fprintf(fp, " %8.5fe%+03d", base, exponent);
			/*	fprintf(fp, " %12.5e", d); */
				if (k == cp->nz - 1 || k % 6 == 5)
					fprintf(fp, "\n");
			}
		}
	}
	fclose(fp);
	return 0;
}

#pragma mark ====== Marching Cube (for isosurface) ======

MCube *
MoleculeClearMCube(Molecule *mol, Int nx, Int ny, Int nz, const Vector *origin, Double dx, Double dy, Double dz)
{
	MCube *mc = mol->mcube;
	int i;
	float rgba[8] = { 1, 1, 1, 0.6, 0, 0, 1, 0.6 };
	if (mc != NULL) {
		free(mc->dp);
		free(mc->radii);
		free(mc->c[0].fp);
		free(mc->c[0].cubepoints);
		free(mc->c[0].triangles);
		free(mc->c[1].fp);
		free(mc->c[1].cubepoints);
		free(mc->c[1].triangles);
		memmove(rgba, mc->c[0].rgba, sizeof(float) * 4);
		memmove(rgba + 4, mc->c[1].rgba, sizeof(float) * 4);
		free(mc);
		mol->mcube = NULL;
	}
	if (nx > 0 && ny > 0 && nz > 0) {
		mc = (MCube *)calloc(sizeof(MCube), 1);
		mc->idn = -1;
		/*  round up to nearest 4N+1 integer  */
		dx *= nx;
		dy *= ny;
		dz *= nz;
		mc->nx = (nx + 2) / 4 * 4 + 1;
		mc->ny = (ny + 2) / 4 * 4 + 1;
		mc->nz = (nz + 2) / 4 * 4 + 1;
		mc->dx = dx / mc->nx;
		mc->dy = dy / mc->ny;
		mc->dz = dz / mc->nz;
		mc->origin = *origin;
		mc->dp = (Double *)malloc(sizeof(Double) * mc->nx * mc->ny * mc->nz);
		if (mc->dp == NULL) {
			free(mc);
			return NULL;
		}
		mc->radii = (Double *)calloc(sizeof(Double), mol->natoms);
		if (mc->radii == NULL) {
			free(mc->dp);
			free(mc);
			return NULL;
		}
		mc->nradii = mol->natoms;
		mc->c[0].fp = (unsigned char *)calloc(sizeof(unsigned char), mc->nx * mc->ny * mc->nz);
		mc->c[1].fp = (unsigned char *)calloc(sizeof(unsigned char), mc->nx * mc->ny * mc->nz);
		if (mc->c[0].fp == NULL || mc->c[1].fp == NULL) {
			free(mc->c[0].fp);
			free(mc->c[1].fp);
			free(mc->dp);
			free(mc->radii);
			free(mc);
			return NULL;
		}
		for (i = 0; i < mc->nx * mc->ny * mc->nz; i++) {
			mc->dp[i] = DBL_MAX;
		}
		memmove(mc->c[0].rgba, rgba, sizeof(float) * 4);
		memmove(mc->c[1].rgba, rgba + 4, sizeof(float) * 4);
		mol->mcube = mc;
	}
	MoleculeCallback_notifyModification(mol, 0);
	return mol->mcube;
}

static int sMarchingCubeTable[256][16] = {
	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
	{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
	{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
	{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
	{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
	{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
	{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
	{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
	{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
	{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
	{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
	{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
	{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
	{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
	{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
	{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
	{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
	{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
	{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
	{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
	{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
	{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
	{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
	{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
	{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
	{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
	{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
	{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
	{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
	{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
	{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
	{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
	{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
	{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
	{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
	{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
	{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
	{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
	{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
	{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
	{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
	{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
	{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
	{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
	{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
	{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
	{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
	{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
	{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
	{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
	{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
	{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
	{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
	{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
	{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
	{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
	{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
	{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
	{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
	{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
	{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
	{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
	{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
	{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
	{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
	{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
	{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
	{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
	{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
	{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
	{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
	{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
	{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
	{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
	{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
	{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
	{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
	{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
	{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
	{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
	{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
	{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
	{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
	{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
	{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
	{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
	{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
	{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
	{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
	{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
	{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
	{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
	{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
	{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
	{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
	{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
	{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
	{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
	{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
	{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
	{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
	{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
	{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
	{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
	{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
	{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
	{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
	{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
	{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
	{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
	{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
	{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
	{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
	{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
	{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
	{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
	{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
	{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
	{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
	{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
	{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
	{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
	{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
	{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
	{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
	{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
	{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
	{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};

/*  Recalculate the MCube  */
/*  If idn < 0, then the current grid settings and values are unchanged, and */
/*  only the marching cubes are regenerated.  */
int
MoleculeUpdateMCube(Molecule *mol, int idn)
{
	Int retval, step, sn;
	Int n, ix, iy, iz, nx, ny, nz;
	Int nn, iix, iiy, iiz;
	Int ncubepoints, c1, c2, c3;
	Int *ip;
	Double thres, *tmp, dd;
	Vector p;
	MCube *mc;
	MCubePoint *mcp;
	Atom *ap;

	if (mol == NULL || mol->bset == NULL || mol->mcube == NULL)
		return -1;
	if (mol->bset->cns == NULL) {
		if (sSetupGaussianCoefficients(mol->bset) != 0)
			return -1;
	}
	if (mol->bset->natoms_bs > mol->natoms)
		return -1;  /*  Number of atoms is smaller than expected  */

	mc = mol->mcube;
	if (idn >= 0) {
		ShellInfo *sp;
		Double *mobasep, *mop, mopmax;
		Double xmin, xmax, ymin, ymax, zmin, zmax;
		/*  Clear mcube values  */
		for (ix = 0; ix < mc->nx * mc->ny * mc->nz; ix++) {
			mc->dp[ix] = DBL_MAX;
			mc->c[0].fp[ix] = 0;
			mc->c[1].fp[ix] = 0;
		}
		mc->idn = idn;
		/*  Estimate the orbital sizes  */
		mc->radii = (Double *)realloc(mc->radii, sizeof(Double) * mol->natoms);
		if (mc->radii == NULL)
			return -2;  /*  Out of memory  */
		mc->nradii = mol->natoms;
		if (mc->idn == mol->bset->nmos + 1) {
			/*  Total electron density  */
			for (ix = 0; ix < mol->natoms; ix++)
				mc->radii[ix] = 1.0;
			mopmax = 1.0;
		} else {
			memset(mc->radii, 0, sizeof(Double) * mc->nradii);
			mobasep = mol->bset->mo + (mc->idn == 0 ? mol->bset->nmos : mc->idn - 1) * mol->bset->ncomps;
			mopmax = 0.0;
			for (ix = 0, sp = mol->bset->shells; ix < mol->bset->nshells; ix++, sp++) {
				if (sp->a_idx >= mol->natoms)
					continue;  /*  This may happen when molecule is edited after setting up MO info  */
				mop = mobasep + sp->m_idx;
				for (iy = 0; iy < sp->ncomp; iy++) {
					dd = fabs(mop[iy]);
					if (dd > mc->radii[sp->a_idx])
						mc->radii[sp->a_idx] = dd;
					if (dd > mopmax)
						mopmax = dd;
				}
			}
		}
		xmin = ymin = zmin = 1e10;
		xmax = ymax = zmax = -1e10;
		for (ix = 0, ap = mol->atoms; ix < mol->natoms; ix++, ap = ATOM_NEXT(ap)) {
			dd = RadiusForAtomicNumber(ap->atomicNumber);
			dd = (dd * 2.0 + 1.0) * (mc->radii[ix] / mopmax) * (mc->expand > 0.0 ? mc->expand : 1.0);
			mc->radii[ix] = dd;
			p = ap->r;
			dd += 0.1;
			if (p.x - dd < xmin)
				xmin = p.x - dd;
			if (p.y - dd < ymin)
				ymin = p.y - dd;
			if (p.z - dd < zmin)
				zmin = p.z - dd;
			if (p.x + dd > xmax)
				xmax = p.x + dd;
			if (p.y + dd > ymax)
				ymax = p.y + dd;
			if (p.z + dd > zmax)
				zmax = p.z + dd;
		}
		mc->origin.x = xmin;
		mc->origin.y = ymin;
		mc->origin.z = zmin;
		mc->dx = (xmax - xmin) / mc->nx;
		mc->dy = (ymax - ymin) / mc->ny;
		mc->dz = (zmax - zmin) / mc->nz;
	}
	
	/*  Temporary work area  */
	tmp = (Double *)calloc(sizeof(Double), mol->bset->natoms_bs * 4);
	if (tmp == NULL)
		return -2;
	
	/*  TODO: use multithread  */
	nx = mc->nx;
	ny = mc->ny;
	nz = mc->nz;
	step = 4;
	
#if 1
	/*  Calculate points within certain distances from atoms  */
	for (nn = 0, ap = mol->atoms; nn < mol->natoms; nn++, ap = ATOM_NEXT(ap)) {
	/*	dd = RadiusForAtomicNumber(ap->atomicNumber);
		if (dd == 0.0)
			dd = 1.0;
		dd = dd * 1.5 + 1.0; */
		dd = mc->radii[nn];
		p.x = ap->r.x - dd - mc->origin.x;
		p.y = ap->r.y - dd - mc->origin.y;
		p.z = ap->r.z - dd - mc->origin.z;
		c1 = p.x / mc->dx;
		c2 = p.y / mc->dy;
		c3 = p.z / mc->dz;
		iix = c1 + ceil(dd * 2.0 / mc->dx);
		iiy = c2 + ceil(dd * 2.0 / mc->dy);
		iiz = c3 + ceil(dd * 2.0 / mc->dz);
		if (c1 < 0)
			c1 = 0;
		if (c2 < 0)
			c2 = 0;
		if (c3 < 0)
			c3 = 0;
		if (iix >= nx)
			iix = nx - 1;
		if (iiy >= ny)
			iiy = ny - 1;
		if (iiz >= nz)
			iiz = nz - 1;
		for (ix = c1; ix <= iix; ix++) {
			p.x = mc->origin.x + mc->dx * ix;
			for (iy = c2; iy <= iiy; iy++) {
				p.y = mc->origin.y + mc->dy * iy;
				for (iz = c3; iz <= iiz; iz++) {
					n = (ix * ny + iy) * nz + iz;
					if (mc->dp[n] == DBL_MAX) {
						p.z = mc->origin.z + mc->dz * iz;
						if (mc->idn == mol->bset->nmos + 1) {
							/*  Total electron density  */
							Int ne_alpha, ne_beta;
							mc->dp[n] = 0.0;
							ne_alpha = mol->bset->ne_alpha;
							ne_beta = mol->bset->ne_beta;
							if (mol->bset->rflag == 2 && ne_alpha < ne_beta) {
								/*  ROHF case: ensure ne_alpha >= ne_beta  */
								ne_beta = ne_alpha;
								ne_alpha = mol->bset->ne_beta;
							}
							for (sn = 1; sn <= ne_alpha; sn++) {
								dd = sCalcMOPoint(mol, mol->bset, sn, &p, tmp);
								dd = dd * dd;
								if (mol->bset->rflag != 0 && sn <= ne_beta)
									dd *= 2;
								mc->dp[n] += dd;
							}
							if (mol->bset->rflag == 0) {
								for (sn = 1; sn <= ne_beta; sn++) {
									dd = sCalcMOPoint(mol, mol->bset, sn + mol->bset->ncomps, &p, tmp);
									mc->dp[n] += dd * dd;
								}
							}
						} else {
							mc->dp[n] = sCalcMOPoint(mol, mol->bset, mc->idn, &p, tmp);
						}
					}
				}
			}
		}
	}
	
#else
	/*  (i * step, j * step, k * step)  */
	for (ix = 0; ix < nx; ix += step) {
		for (iy = 0; iy < ny; iy += step) {
			for (iz = 0; iz < nz; iz += step) {
				n = (ix * ny + iy) * nz + iz;
				if (mc->dp[n] == DBL_MAX) {
					p.x = mc->origin.x + mc->dx * ix;
					p.y = mc->origin.y + mc->dy * iy;
					p.z = mc->origin.z + mc->dz * iz;
					mc->dp[n] = sCalcMOPoint(mol, mol->bset, mc->idn, &p, tmp);
				}
				n += step;
			}
		}
	}
	
	/*  Intermediate points  */
	for (step = 4; step > 1; step /= 2) {
		hstep = step / 2;
		for (sn = 0; sn <= 1; sn++) {
			n = 0;
			for (ix = 0; ix < nx - 1; ix += step) {
				for (iy = 0; iy < ny - 1; iy += step) {
					for (iz = 0; iz < nz - 1; iz += step) {
						flags = 0;
						thres = mc->thres * (sn == 0 ? 1 : -1);
						n = (ix * ny + iy) * nz + iz;
						if (mc->dp[n] == DBL_MAX || mc->dp[n + step * (nz * (ny + 1) + 1)] == DBL_MAX)
							continue;
						/*  (ix, iy, iz)  */
						if (mc->dp[n] >= thres)
							flags |= 1;
						/*  (ix + step, iy, iz)  */
						if (mc->dp[n + step * ny * nz] >= thres)
							flags |= 2;
						/*  (ix, iy + step, iz)  */
						if (mc->dp[n + step * nz] >= thres)
							flags |= 4;
						/*  (ix + 4, iy + step, iz)  */
						if (mc->dp[n + step * nz * (ny + 1)] >= thres)
							flags |= 8;
						/*  (ix, iy, iz + step)  */
						if (mc->dp[n + step] >= thres)
							flags |= 16;
						if (mc->dp[n + step * (ny * nz + 1)] >= thres)
							flags |= 32;
						/*  (ix, iy + step, iz + step)  */
						if (mc->dp[n + step * (nz + 1)] >= thres)
							flags |= 64;
						/*  (ix + step, iy + step, iz + step)  */
						if (mc->dp[n + step * (nz * (ny + 1) + 1)] >= thres)
							flags |= 128;
						if (flags != 0 && flags != 255) {
							/*  Calc the intermediate points  */
							for (iix = 0; iix <= step; iix += hstep) {
								for (iiy = 0; iiy <= step; iiy += hstep) {
									for (iiz = 0; iiz <= step; iiz += hstep) {
										if (iix % step == 0 && iiy % step == 0 && iiz % step == 0)
											continue;
										nn = n + (iix * ny + iiy) * nz + iiz;
										if (mc->dp[nn] == DBL_MAX) {
											p.x = mc->origin.x + mc->dx * (ix + iix);
											p.y = mc->origin.y + mc->dy * (iy + iiy);
											p.z = mc->origin.z + mc->dz * (iz + iiz);
											mc->dp[nn] = sCalcMOPoint(mol, mol->bset, mc->idn, &p, tmp);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
#endif

	free(tmp);
	
	/*  Calculate vertex positions and normal vectors  */
	for (sn = 0; sn <= 1; sn++) {
		n = 0;
		thres = mc->thres * (sn == 0 ? 1 : -1);
		VecZero(p);
		for (ix = 0; ix < nx - 1; ix++) {
			for (iy = 0; iy < ny - 1; iy++) {
				for (iz = 0; iz < nz - 1; iz++) {
					Double dd0, dd1;
					nn = (ix * ny + iy) * nz + iz;
					dd0 = mc->dp[nn];
					if (dd0 == DBL_MAX)
						continue;
					if (0) {
						dd1 = mc->dp[nn + ny * nz];
						if (dd1 != DBL_MAX)
							p.x = (dd1 - dd0) / mc->dx;
						else if (ix > 0 && (dd1 = mc->dp[nn - ny * nz]) != DBL_MAX)
							p.x = (dd0 - dd1) / mc->dx;
						else continue;  /*  Cannot define gradient  */
						dd1 = mc->dp[nn + nz];
						if (dd1 != DBL_MAX)
							p.y = (dd1 - dd0) / mc->dy;
						else if (iy > 0 && (dd1 = mc->dp[nn - nz]) != DBL_MAX)
							p.y = (dd0 - dd1) / mc->dy;
						else continue;
						dd1 = mc->dp[nn + 1];
						if (dd1 != DBL_MAX)
							p.z = (dd1 - dd0) / mc->dz;
						else if (iz > 0 && (dd1 = mc->dp[nn - 1]) != DBL_MAX)
							p.z = (dd0 - dd1) / mc->dz;
						else continue;
						NormalizeVec(&p, &p);
					}
					if (n + 3 >= mc->c[sn].ncubepoints) {
						/*  Expand cubepoints[] array  */
						mc->c[sn].cubepoints = (MCubePoint *)realloc(mc->c[sn].cubepoints, sizeof(MCubePoint) * (mc->c[sn].ncubepoints + 8192));
						if (mc->c[sn].cubepoints == NULL) {
							mc->c[sn].ncubepoints = 0;
							retval = -3;
							goto end;
						}
						mc->c[sn].ncubepoints += 8192;
					}
					mcp = mc->c[sn].cubepoints + n;
					iix = (dd0 >= thres ? 1 : -1);
					/*  (x, y, z)->(x + 1, y, z)  */
					dd1 = mc->dp[nn + ny * nz];
					if (dd1 != DBL_MAX) {
						iiy = (dd1 >= thres ? 1 : -1);
						if (iix != iiy) {
							/*  Register  */
							mcp->key = nn * 3;
							mcp->d = (thres - dd0) / (dd1 - dd0);
							mcp->pos[0] = mc->origin.x + mc->dx * (ix + mcp->d);
							mcp->pos[1] = mc->origin.y + mc->dy * iy;
							mcp->pos[2] = mc->origin.z + mc->dz * iz;
							mcp->grad[0] = p.x;
							mcp->grad[1] = p.y;
							mcp->grad[2] = p.z;
							mcp++;
							n++;
						}
					}
					/*  (x, y, z)->(x, y + 1, z)  */
					dd1 = mc->dp[nn + nz];
					if (dd1 != DBL_MAX) {
						iiy = (dd1 >= thres ? 1 : -1);
						if (iix != iiy) {
							/*  Register  */
							mcp->key = nn * 3 + 1;
							mcp->d = (thres - dd0) / (dd1 - dd0);
							mcp->pos[0] = mc->origin.x + mc->dx * ix;
							mcp->pos[1] = mc->origin.y + mc->dy * (iy + mcp->d);
							mcp->pos[2] = mc->origin.z + mc->dz * iz;
							mcp->grad[0] = p.x;
							mcp->grad[1] = p.y;
							mcp->grad[2] = p.z;
							mcp++;
							n++;
						}
					}
					/*  (x, y, z)->(x, y, z + 1)  */
					dd1 = mc->dp[nn + 1];
					if (dd1 != DBL_MAX) {
						iiy = (dd1 >= thres ? 1 : -1);
						if (iix != iiy) {
							/*  Register  */
							mcp->key = nn * 3 + 2;
							mcp->d = (thres - dd0) / (dd1 - dd0);
							mcp->pos[0] = mc->origin.x + mc->dx * ix;
							mcp->pos[1] = mc->origin.y + mc->dy * iy;
							mcp->pos[2] = mc->origin.z + mc->dz * (iz + mcp->d);
							mcp->grad[0] = p.x;
							mcp->grad[1] = p.y;
							mcp->grad[2] = p.z;
							mcp++;
							n++;
						}
					}
				}
			}
		}
		if (n < mc->c[sn].ncubepoints)
			mc->c[sn].cubepoints[n].key = -1;  /*  End mark  */
		ncubepoints = n;
		if (ncubepoints < 3) {
			/*  Less than 3 points: no triangles  */
			if (mc->c[sn].ntriangles > 0)
				mc->c[sn].triangles[0] = -1;  /*  End mark  */
			continue;
		}
		
		/*  Create triangle table  */
		n = 0;
		for (ix = 0; ix < nx - 1; ix++) {
			for (iy = 0; iy < ny - 1; iy++) {
				for (iz = 0; iz < nz - 1; iz++) {
					nn = (ix * ny + iy) * nz + iz;
					iix = 0;
					if ((dd = mc->dp[nn]) == DBL_MAX)
						continue;
					else if (dd >= thres)
						iix |= 1;
					if ((dd = mc->dp[nn + ny * nz]) == DBL_MAX)
						continue;
					else if (dd >= thres)
						iix |= 2;
					if ((dd = mc->dp[nn + ny * nz + nz]) == DBL_MAX)
						continue;
					else if (dd >= thres)
						iix |= 4;
					if ((dd = mc->dp[nn + nz]) == DBL_MAX)
						continue;
					else if (dd >= thres)
						iix |= 8;
					if ((dd = mc->dp[nn + 1]) == DBL_MAX)
						continue;
					else if (dd >= thres)
						iix |= 16;
					if ((dd = mc->dp[nn + ny * nz + 1]) == DBL_MAX)
						continue;
					else if (dd >= thres)
						iix |= 32;
					if ((dd = mc->dp[nn + ny * nz + nz + 1]) == DBL_MAX)
						continue;
					else if (dd >= thres)
						iix |= 64;
					if ((dd = mc->dp[nn + nz + 1]) == DBL_MAX)
						continue;
					else if (dd >= thres)
						iix |= 128;
					for (iiy = 0; iiy < 15; iiy++) {
						nn = sMarchingCubeTable[iix][iiy];
						if (nn < 0)
							break;
						/*  key index for edges 0-11  */
						switch (nn) {
							case 0:  iiz = (( ix      * ny + iy    ) * nz + iz    ) * 3;     break;
							case 1:  iiz = (((ix + 1) * ny + iy    ) * nz + iz    ) * 3 + 1; break;
							case 2:  iiz = (( ix      * ny + iy + 1) * nz + iz    ) * 3;     break;
							case 3:  iiz = (( ix      * ny + iy    ) * nz + iz    ) * 3 + 1; break;
							case 4:  iiz = (( ix      * ny + iy    ) * nz + iz + 1) * 3;     break;
							case 5:  iiz = (((ix + 1) * ny + iy    ) * nz + iz + 1) * 3 + 1; break;
							case 6:  iiz = (( ix      * ny + iy + 1) * nz + iz + 1) * 3;     break;
							case 7:  iiz = (( ix      * ny + iy    ) * nz + iz + 1) * 3 + 1; break;
							case 8:  iiz = (( ix      * ny + iy    ) * nz + iz    ) * 3 + 2; break;
							case 9:  iiz = (((ix + 1) * ny + iy    ) * nz + iz    ) * 3 + 2; break;
							case 10: iiz = (((ix + 1) * ny + iy + 1) * nz + iz    ) * 3 + 2; break;
							case 11: iiz = (( ix      * ny + iy + 1) * nz + iz    ) * 3 + 2; break;
							default:
								/*  Skip this triangle  */
								iiy = (iiy - iiy % 3) + 2;
								n = n - n % 3;
								continue;
						}
						/*  Look for the key index in cubepoints  */
						c1 = 0;
						c3 = ncubepoints - 1;
						mcp = mc->c[sn].cubepoints;
						while (1) {
							int w;
							/*  c1 is always less than c3  */
							if (c1 + 1 == c3) {
								/*  end of search  */
								if (mcp[c1].key == iiz) {
									c2 = c1;
								} else if (mcp[c3].key == iiz) {
									c2 = c3;
								} else {
									c2 = -1;
								}
								break;
							}
							c2 = (c1 + c3) / 2;
							w = mcp[c2].key - iiz;
							if (w == 0)
								break;
							if (w < 0) {
								c1 = c2;
							} else {
								c3 = c2;
							}
						}
						if (c2 < 0) {
							/*  Not found: skip this triangle  */
							iiy = (iiy - iiy % 3) + 2;
							n = n - n % 3;
							continue;
						}
						if (n + 1 >= mc->c[sn].ntriangles) {
							/*  Expand triangles[] array  */
							mc->c[sn].triangles = (Int *)realloc(mc->c[sn].triangles, sizeof(Int) * (mc->c[sn].ntriangles + 8192));
							if (mc->c[sn].triangles == NULL) {
								mc->c[sn].ntriangles = 0;
								retval = -4;
								goto end;
							}
							mc->c[sn].ntriangles += 8192;
						}
						mc->c[sn].triangles[n] = c2;
						n++;
					}
				}
			}
		}
		if (n < mc->c[sn].ntriangles)
			mc->c[sn].triangles[n] = -1;  /*  End mark  */
		
		/*  Estimate the normal vector  */
		for (n = 0, ip = mc->c[sn].triangles; ip[n] >= 0; n += 3) {
			Vector v[3];
			for (ix = 0; ix < 3; ix++) {
				mcp = &(mc->c[sn].cubepoints[ip[n + ix]]);
				v[ix].x = mcp->pos[0];
				v[ix].y = mcp->pos[1];
				v[ix].z = mcp->pos[2];
			}
			VecDec(v[2], v[0]);
			VecDec(v[1], v[0]);
			VecCross(v[0], v[1], v[2]);
			NormalizeVec(v, v);
			for (ix = 0; ix < 3; ix++) {
				mcp = &(mc->c[sn].cubepoints[ip[n + ix]]);
				mcp->grad[0] += v[0].x;
				mcp->grad[1] += v[0].y;
				mcp->grad[2] += v[0].z;
			}
		}
		for (n = 0, mcp = mc->c[sn].cubepoints; mcp->key >= 0; mcp++) {
			if (mcp->grad[0] != 0.0 || mcp->grad[1] != 0.0 || mcp->grad[2] != 0.0) {
				dd = 1.0 / sqrt(mcp->grad[0] * mcp->grad[0] + mcp->grad[1] * mcp->grad[1] + mcp->grad[2] * mcp->grad[2]);
				if (mc->thres < 0.0)
					dd = -dd;
				mcp->grad[0] *= dd;
				mcp->grad[1] *= dd;
				mcp->grad[2] *= dd;
			}
		}
	}
	retval = 0;
	MoleculeCallback_notifyModification(mol, 0);
end:
	/*  For debug  */
	if (0) {
		char *MyAppCallback_getDocumentHomeDir(void);
		FILE *fp;
		char *s;
		Double dmax, dmin;
		asprintf(&s, "%s/%s", MyAppCallback_getDocumentHomeDir(), "mcube_log.txt");
		fp = fopen(s, "w");
		dmax = -1e8;
		dmin = 1e8;
		for (n = 0; n < mc->nx * mc->ny * mc->nz; n++) {
			if (mc->dp[n] == DBL_MAX)
				continue;
			if (dmax < mc->dp[n])
				dmax = mc->dp[n];
			if (dmin > mc->dp[n])
				dmin = mc->dp[n];
		}
		dmax = fabs(dmax);
		dmin = fabs(dmin);
		if (dmax < dmin)
			dmax = dmin;
		dmax = 1.001 * dmax;
		fprintf(fp, "thres = %g = 100\n", mc->thres);
		for (iz = 0; iz < mc->nz; iz++) {
			fprintf(fp, "z = %d\n", iz);
			for (iy = 0; iy < mc->ny; iy++) {
				for (ix = 0; ix < mc->nx; ix++) {
					n = (ix * ny + iy) * nz + iz;
					dd = mc->dp[n];
					if (dd == DBL_MAX)
						fprintf(fp, " XXX ");
					else {
						dd = dd * 100 / mc->thres;
						if (dd > 999.0)
							dd = 999.0;
						else if (dd < -999.0)
							dd = -999.0;
						fprintf(fp, "%4d ", (int)(dd));
					}
				}
				fprintf(fp, "\n");
			}
			fprintf(fp, "\n");
		}
		
		for (sn = 0; sn <= 1; sn++) {
			for (n = 0; n < mc->c[sn].ncubepoints; n++) {
				MCubePoint *mcp = mc->c[sn].cubepoints + n;
				nn = mcp->key;
				if (nn == -1)
					break;
				iix = nn % 3;
				iz = nn / 3 % mc->nz;
				iy = nn / (3 * mc->nz) % mc->ny;
				ix = nn / (3 * mc->nz * mc->ny);
				fprintf(fp, "%c%d:[%d,%d,%d,%d] (%g,[%g,%g,%g],[%g,%g,%g])\n", (sn == 0 ? 'p' : 'P'),
						n, ix, iy, iz, iix,
						mcp->d, mcp->pos[0], mcp->pos[1], mcp->pos[2], mcp->grad[0], mcp->grad[1], mcp->grad[2]);
			}
			for (n = 0; n < mc->c[sn].ntriangles; n += 3) {
				if (mc->c[sn].triangles[n] < 0)
					break;
				fprintf(fp, "%c%d:(%d,%d,%d)\n", (sn == 0 ? 't' : 'T'), n / 3,
						mc->c[sn].triangles[n], mc->c[sn].triangles[n + 1], mc->c[sn].triangles[n + 2]);
			}
		}
		fclose(fp);
	}
	
	return retval;
}

void
MoleculeDeallocateMCube(MCube *mcube)
{
	free(mcube->dp);
	free(mcube->radii);
	free(mcube->c[0].cubepoints);
	free(mcube->c[0].triangles);
	free(mcube->c[1].cubepoints);
	free(mcube->c[1].triangles);
	free(mcube);
}
