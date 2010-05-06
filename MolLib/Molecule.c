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

#include "Missing.h"
#include "Dcd.h"
#include "MD/MDCore.h"
#include "MD/MDPressure.h"

static Molecule *sMoleculeRoot = NULL;
static int sMoleculeUntitledCount = 0;

Int gSizeOfAtomRecord = sizeof(Atom);

#pragma mark ====== Utility function ======

int
strlen_limit(const char *s, int limit)
{
	int len;
	for (len = 0; *s != 0 && (limit < 0 || len < limit); s++, len++);
	return len;
}

#pragma mark ======  Atom handling  ======

Atom *
AtomDuplicate(Atom *dst, const Atom *src)
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
	if (src->frames != NULL) {
		dst->frames = (Vector *)malloc(sizeof(Vector) * src->nframes);
		if (dst->frames != NULL) {
			memmove(dst->frames, src->frames, sizeof(Vector) * src->nframes);
			dst->nframes = src->nframes;
		} else {
			dst->nframes = 0;
		}
	}
	return dst;
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
	if (bset->pos != NULL)
		free(bset->pos);
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
MoleculeInitWithMolecule(Molecule *mp2, const Molecule *mp)
{
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

Molecule *
MoleculeRetain(Molecule *mp)
{
	ObjectIncrRefCount((Object *)mp);
	return mp;
}

void
MoleculeRelease(Molecule *mp)
{
	if (mp == NULL)
		return;
	if (ObjectDecrRefCount((Object *)mp) == 0) {
		if (mp->atoms != NULL) {
			int i;
			for (i = 0; i < mp->natoms; i++)
				AtomClean(mp->atoms + i);
			free(mp->atoms);
		}
		if (mp->bonds != NULL)
			free(mp->bonds);
		if (mp->angles != NULL)
			free(mp->angles);
		if (mp->dihedrals != NULL)
			free(mp->dihedrals);
		if (mp->impropers != NULL)
			free(mp->impropers);
		if (mp->residues != NULL)
			free(mp->residues);
		if (mp->bset != NULL)
			BasisSetRelease(mp->bset);
		if (mp->par != NULL)
			ParameterRelease(mp->par);
		if (mp->elpots != NULL)
			free(mp->elpots);
		ObjectDealloc((Object *)mp, (Object **)&sMoleculeRoot);
	}
}

void
MoleculeExchange(Molecule *mp1, Molecule *mp2)
{
	Molecule mp_temp;
	/*  'natoms' is the first member to be copied  */
	int ofs = offsetof(Molecule, natoms);
	memmove((char *)(&mp_temp) + ofs, (char *)mp1 + ofs, sizeof(Molecule) - ofs);
	memmove((char *)mp1 + ofs, (char *)mp2 + ofs, sizeof(Molecule) - ofs);
	memmove((char *)mp2 + ofs, (char *)(&mp_temp) + ofs, sizeof(Molecule) - ofs);
	if (mp1->arena != NULL && mp1->arena->mol == mp2)
		mp1->arena->mol = mp1;
	if (mp1->arena != NULL && mp1->arena->xmol == mp2)
		mp1->arena->xmol = mp1;
	if (mp2->arena != NULL && mp2->arena->mol == mp1)
		mp2->arena->mol = mp2;
	if (mp2->arena != NULL && mp2->arena->xmol == mp1)
		mp2->arena->xmol = mp2;
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

int
MoleculeLoadFile(Molecule *mp, const char *fname, const char *ftype, char *errbuf, int errbufsize)
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
		retval = MoleculeLoadPsfFile(mp, fname, errbuf, errbufsize);
	} else if (strcasecmp(ftype, "pdb") == 0) {
		retval = MoleculeReadCoordinatesFromPdbFile(mp, fname, errbuf, errbufsize);
	} else if (strcasecmp(ftype, "tep") == 0) {
		retval = MoleculeLoadTepFile(mp, fname, errbuf, errbufsize);
	} else if (strcasecmp(ftype, "res") == 0 || strcasecmp(ftype, "ins") == 0) {
		retval = MoleculeLoadShelxFile(mp, fname, errbuf, errbufsize);
	} else if (strcasecmp(ftype, "fchk") == 0) {
		retval = MoleculeLoadGaussianFchkFile(mp, fname, errbuf, errbufsize);
	} else {
		snprintf(errbuf, errbufsize, "Unknown format %s", ftype);
		return 1;
	}
/*	if (retval != 0) {
		retval = MoleculeLoadPsfFile(mp, fname, errbuf, errbufsize);
	} */
	return retval;
}

int
MoleculeLoadMbsfFile(Molecule *mol, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	Molecule *mp;
	char buf[1024];
	int i, j, k, n, err, fn, nframes;
	int lineNumber;
	int ibuf[12];
	Int iibuf[4];
	double dbuf[12];
	char cbuf[12][6];
	const char **pp;
	char *bufp, *valp, *comp;
	Int *ip;
	Double *dp;
	Vector v;
	Atom *ap;
	err = 0;
	if (errbuf == NULL) {
		errbuf = buf;
		errbufsize = 1024;
	}
	errbuf[0] = 0;
	mp = MoleculeNew();
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot open file");
		return 1;
	}
	/*	flockfile(fp); */
	lineNumber = 0;
	fn = 0;
	nframes = 0;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		if (strstr(buf, "!:atoms") == buf) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx seg_name res_seq res_name name type charge weight element atomic_number occupancy temp_factor int_charge */
				if (sscanf(buf, "%d %4s %d %4s %4s %4s %lf %lf %4s %d %lf %lf %d", &ibuf[0], cbuf[0], &ibuf[1], cbuf[1], cbuf[2], cbuf[3], &dbuf[0], &dbuf[1], cbuf[4], &ibuf[2], &dbuf[2], &dbuf[3], &ibuf[3]) < 13) {
					snprintf(errbuf, errbufsize, "line %d: coordinates cannot be read for atom %d", lineNumber, mp->natoms + 1);
					goto exit;
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
		} else if (strstr(buf, "!:atoms_symop") == buf) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx symop symbase */
				if (sscanf(buf, "%d %d %d", &ibuf[0], &ibuf[1], &ibuf[2]) < 3) {
					snprintf(errbuf, errbufsize, "line %d: symmetry operations cannot be read for atom %d", lineNumber, i + 1);
					goto exit;
				}
				if (i >= mp->natoms) {
					snprintf(errbuf, errbufsize, "line %d: too many atomic symmetry info\n", lineNumber);
					goto exit;
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
		} else if (strstr(buf, "!:atoms_fix") == buf) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx fix_force fix_pos */
				if (sscanf(buf, "%d %lf %lf %lf %lf", &ibuf[0], &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3]) < 5) {
					snprintf(errbuf, errbufsize, "line %d: fix atom info cannot be read for atom %d", lineNumber, i + 1);
					goto exit;
				}
				if (i >= mp->natoms) {
					snprintf(errbuf, errbufsize, "line %d: too many fix atom info\n", lineNumber);
					goto exit;
				}
				ap = ATOM_AT_INDEX(mp->atoms, i);
				ap->fix_force = dbuf[0];
				ap->fix_pos.x = dbuf[1];
				ap->fix_pos.y = dbuf[2];
				ap->fix_pos.z = dbuf[3];
				i++;
			}
			continue;
		} else if (strstr(buf, "!:positions") == buf) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx x y z */
				if (sscanf(buf, "%d %lf %lf %lf", &ibuf[0], &dbuf[0], &dbuf[1], &dbuf[2]) < 4) {
					snprintf(errbuf, errbufsize, "line %d: atom position cannot be read for atom %d frame %d", lineNumber, i + 1, nframes);
					goto exit;
				}
				if (i >= mp->natoms) {
					snprintf(errbuf, errbufsize, "line %d: too many atom position records\n", lineNumber);
					goto exit;
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
		} else if (strstr(buf, "!:bonds") == buf) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* from1 to1 from2 to2 from3 to3 from4 to4 */ 
				i = sscanf(buf, "%d %d %d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5], &ibuf[6], &ibuf[7]);
				if (i < 2 || i % 2 != 0) {
					snprintf(errbuf, errbufsize, "line %d: bad bond format", lineNumber);
					goto exit;
				}
				for (j = 0; j < i; j += 2) {
					iibuf[0] = ibuf[j];
					iibuf[1] = ibuf[j + 1];
					if (iibuf[0] < 0 || iibuf[0] >= mp->natoms || iibuf[1] < 0 || iibuf[1] >= mp->natoms || iibuf[0] == iibuf[1]) {
						snprintf(errbuf, errbufsize, "line %d: bad bond format", lineNumber);
						goto exit;
					}
					AssignArray(&mp->bonds, &mp->nbonds, sizeof(Int) * 2, mp->nbonds, iibuf);
					for (k = 0; k < 2; k++) {
						ap = ATOM_AT_INDEX(mp->atoms, iibuf[k]);
						for (n = 0; n < ap->nconnects; n++) {
							if (ap->connects[n] == iibuf[1 - k])
								break;
						}
						if (n >= ap->nconnects) {
							if (ap->nconnects >= ATOMS_MAX_CONNECTS - 1) {
								snprintf(errbuf, errbufsize, "line %d: too many bonds on atom %d", lineNumber, iibuf[k]);
								goto exit;
							}
							ap->connects[ap->nconnects++] = iibuf[1 - k];
						}
					}
				}
			}
			continue;
		} else if (strstr(buf, "!:angles") == buf) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a1 b1 c1 a2 b2 c2 a3 b3 c3 */ 
				i = sscanf(buf, "%d %d %d %d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5], &ibuf[6], &ibuf[7], &ibuf[8]);
				if (i == 0 || i % 3 != 0) {
					snprintf(errbuf, errbufsize, "line %d: bad angle format", lineNumber);
					goto exit;
				}
				for (j = 0; j < i; j += 3) {
					iibuf[0] = ibuf[j];
					iibuf[1] = ibuf[j + 1];
					iibuf[2] = ibuf[j + 2];
					if (iibuf[0] < 0 || iibuf[0] >= mp->natoms || iibuf[1] < 0 || iibuf[1] >= mp->natoms || iibuf[2] < 0 || iibuf[2] >= mp->natoms || iibuf[0] == iibuf[1] || iibuf[1] == iibuf[2]) {
						snprintf(errbuf, errbufsize, "line %d: bad angle format", lineNumber);
						goto exit;
					}
					AssignArray(&mp->angles, &mp->nangles, sizeof(Int) * 3, mp->nangles, iibuf);
				}
			}
			continue;
		} else if (strstr(buf, "!:dihedrals") == buf) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a1 b1 c1 d1 a2 b2 c2 d2 */ 
				i = sscanf(buf, "%d %d %d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5], &ibuf[6], &ibuf[7]);
				if (i == 0 || i % 4 != 0) {
					snprintf(errbuf, errbufsize, "line %d: bad dihedral format", lineNumber);
					goto exit;
				}
				for (j = 0; j < i; j += 4) {
					iibuf[0] = ibuf[j];
					iibuf[1] = ibuf[j + 1];
					iibuf[2] = ibuf[j + 2];
					iibuf[3] = ibuf[j + 3];
					if (iibuf[0] < 0 || iibuf[0] >= mp->natoms || iibuf[1] < 0 || iibuf[1] >= mp->natoms || iibuf[2] < 0 || iibuf[2] >= mp->natoms || iibuf[3] < 0 || iibuf[3] >= mp->natoms || iibuf[0] == iibuf[1] || iibuf[1] == iibuf[2] || iibuf[2] == iibuf[3] || iibuf[0] == iibuf[2] || iibuf[1] == iibuf[3] || iibuf[0] == iibuf[3]) {
						snprintf(errbuf, errbufsize, "line %d: bad dihedral format", lineNumber);
						goto exit;
					}
					AssignArray(&mp->dihedrals, &mp->ndihedrals, sizeof(Int) * 4, mp->ndihedrals, iibuf);
				}
			}
			continue;
		} else if (strstr(buf, "!:impropers") == buf) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a1 b1 c1 d1 a2 b2 c2 d2 */ 
				i = sscanf(buf, "%d %d %d %d %d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3], &ibuf[4], &ibuf[5], &ibuf[6], &ibuf[7]);
				if (i == 0 || i % 4 != 0) {
					snprintf(errbuf, errbufsize, "line %d: bad improper format", lineNumber);
					goto exit;
				}
				for (j = 0; j < i; j += 4) {
					iibuf[0] = ibuf[j];
					iibuf[1] = ibuf[j + 1];
					iibuf[2] = ibuf[j + 2];
					iibuf[3] = ibuf[j + 3];
					if (iibuf[0] < 0 || iibuf[0] >= mp->natoms || iibuf[1] < 0 || iibuf[1] >= mp->natoms || iibuf[2] < 0 || iibuf[2] >= mp->natoms || iibuf[3] < 0 || iibuf[3] >= mp->natoms || iibuf[0] == iibuf[1] || iibuf[1] == iibuf[2] || iibuf[2] == iibuf[3] || iibuf[0] == iibuf[2] || iibuf[1] == iibuf[3] || iibuf[0] == iibuf[3]) {
						snprintf(errbuf, errbufsize, "line %d: bad improper format", lineNumber);
						goto exit;
					}
					AssignArray(&mp->impropers, &mp->nimpropers, sizeof(Int) * 4, mp->nimpropers, iibuf);
				}
			}
			continue;
		} else if (strstr(buf, "!:xtalcell") == buf && mp->cell == NULL) {
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a b c alpha beta gamma */ 
				if (sscanf(buf, "%lf %lf %lf %lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5]) < 6) {
					snprintf(errbuf, errbufsize, "line %d: bad xtalcell format", lineNumber);
					goto exit;
				}
				MoleculeSetCell(mp, dbuf[0], dbuf[1], dbuf[2], dbuf[3], dbuf[4], dbuf[5], 0);
			}
			continue;
		} else if (strstr(buf, "!:symmetry_operations") == buf) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				Transform tr;
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a11 a12 a13; a21 a22 a23; a31 a32 a33; t1 t2 t3 */
				if (sscanf(buf, "%lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2]) < 3) {
					snprintf(errbuf, errbufsize, "line %d: bad symmetry_operation format", lineNumber);
					goto exit;
				}
				tr[i * 3] = dbuf[0];
				tr[i * 3 + 1] = dbuf[1];
				tr[i * 3 + 2] = dbuf[2];
				i++;
				if (i == 4) {
					AssignArray(&mp->syms, &mp->nsyms, sizeof(Transform), mp->nsyms, tr);
					i = 0;
				}
			}
			continue;
		} else if (strstr(buf, "!:periodic_box") == buf) {
			Vector vs[5];
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* a11 a12 a13; a21 a22 a23; a31 a32 a33; t1 t2 t3 */
				if (sscanf(buf, "%lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2]) < 3) {
					snprintf(errbuf, errbufsize, "line %d: bad symmetry_operation format", lineNumber);
					goto exit;
				}
				vs[i].x = dbuf[0];
				vs[i].y = dbuf[1];
				vs[i].z = dbuf[2];
				i++;
				if (i == 5) {
				/*	j = sscanf(buf, "%d %d %d %d", &ibuf[0], &ibuf[1], &ibuf[2], &ibuf[3]); */
					cbuf[0][0] = dbuf[0];
					cbuf[0][1] = dbuf[1];
					cbuf[0][2] = dbuf[2];
					MoleculeSetPeriodicBox(mp, vs, vs + 1, vs + 2, vs + 3, cbuf[0]);
				/*	if (j == 4)
						mp->is_xtal_coord = (ibuf[3] != 0); */
				}
			}
			continue;
		} else if (strstr(buf, "!:md_parameters") == buf) {
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
						   || (strcmp(comp, "temperature") == 0 && (dp = &arena->temperature) != NULL)
						   || (strcmp(comp, "andersen_coupling") == 0 && (dp = &arena->andersen_thermo_coupling) != NULL)
						   || (strcmp(comp, "dielectric") == 0 && (dp = &arena->dielectric) != NULL)
						   || (strcmp(comp, "gradient_convergence") == 0 && (dp = &arena->gradient_convergence) != NULL)
						   || (strcmp(comp, "coordinate_convergence") == 0 && (dp = &arena->coordinate_convergence) != NULL)
						   || (strcmp(comp, "scale14_vdw") == 0 && (dp = &arena->scale14_vdw) != NULL)
						   || (strcmp(comp, "scale14_elect") == 0 && (dp = &arena->scale14_elect) != NULL)
						   || (strcmp(comp, "surface_probe_radius") == 0 && (dp = &arena->probe_radius) != NULL)
						   || (strcmp(comp, "surface_tension") == 0 && (dp = &arena->surface_tension) != NULL)) {
					*dp = (valp == NULL ? 0.0 : strtod(valp, NULL));
				}
			}
			continue;
		} else if (strstr(buf, "!:pressure_control_parameters") == buf) {
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
						snprintf(errbuf, errbufsize, "line %d: bad format", lineNumber);
						goto exit;
					}
					for (i = 0; i < 9; i++)
						pressure->apply[i] = dbuf[i];
				} else if (strcmp(comp, "pressure_cell_flexibility") == 0) {
					if (sscanf(valp, "%lf %lf %lf %lf %lf %lf %lf %lf", &dbuf[0], &dbuf[1], &dbuf[2], &dbuf[3], &dbuf[4], &dbuf[5], &dbuf[6], &dbuf[7]) < 8) {
						snprintf(errbuf, errbufsize, "line %d: bad format", lineNumber);
						goto exit;
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
		} else if (strstr(buf, "!:velocity") == buf) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx vx vy vz */
				if (sscanf(buf, "%d %lf %lf %lf", &ibuf[0], &dbuf[0], &dbuf[1], &dbuf[2]) < 4) {
					snprintf(errbuf, errbufsize, "line %d: atom velocity cannot be read for atom %d", lineNumber, i + 1);
					goto exit;
				}
				if (i >= mp->natoms) {
					snprintf(errbuf, errbufsize, "line %d: too many atom velocity records\n", lineNumber);
					goto exit;
				}
				ap = ATOM_AT_INDEX(mp->atoms, i);
				ap->v.x = dbuf[0];
				ap->v.y = dbuf[1];
				ap->v.z = dbuf[2];
				i++;
			}
			continue;
		} else if (strstr(buf, "!:force") == buf) {
			i = 0;
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (buf[0] == '!')
					continue;
				if (buf[0] == '\n')
					break;
				/* idx fx fy fz */
				if (sscanf(buf, "%d %lf %lf %lf", &ibuf[0], &dbuf[0], &dbuf[1], &dbuf[2]) < 4) {
					snprintf(errbuf, errbufsize, "line %d: atom force cannot be read for atom %d", lineNumber, i + 1);
					goto exit;
				}
				if (i >= mp->natoms) {
					snprintf(errbuf, errbufsize, "line %d: too many atom force records\n", lineNumber);
					goto exit;
				}
				ap = ATOM_AT_INDEX(mp->atoms, i);
				ap->f.x = dbuf[0];
				ap->f.y = dbuf[1];
				ap->f.z = dbuf[2];
				i++;
			}
			continue;
		} else if (strstr(buf, "!:parameter") == buf) {
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
					snprintf(errbuf, errbufsize, "%s", bufp);
					goto exit;
				}
				i += j;
			}
			if (bufp != NULL) {
				MyAppCallback_setConsoleColor(1);
				MyAppCallback_showScriptMessage("%s", bufp);
				MyAppCallback_setConsoleColor(0);
				free(bufp);
			}
			continue;
		}
		/*  Unknown sections are silently ignored  */
	}

	MoleculeCleanUpResidueTable(mp);
	if (mp->arena != NULL)
		md_arena_set_molecule(mp->arena, mp);

exit:
	fclose(fp);
	if (errbuf[0] != 0) {
		free(mp);
		return -1;
	} else {
		MoleculeExchange(mp, mol);
		MoleculeRelease(mp);
	}
	return 0;
	
}

int
MoleculeLoadPsfFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
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
	if (errbuf == NULL) {
		errbuf = buf;
		errbufsize = 1024;
	}
	errbuf[0] = 0;
	if (mp == NULL)
		mp = MoleculeNew();
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot open file");
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
			#if 0
				if (fn == 1) {
					/*  Copy the coordinates of the first frame  */
					for (i = 0; i < mp->natoms; i++) {
						ap = ATOM_AT_INDEX(mp->atoms, i);
						frames[i] = ap->r;
					}
				}
				/*  Copy the coordinates of the last frame to the newly created frame  */
				memmove(frames + sizeof(Vector) * mp->natoms * fn, frames + sizeof(Vector) * mp->natoms * (fn - 1), sizeof(Vector) * mp->natoms);
			#endif
			}
			/*  Read coordinates  */
			for (i = 0; i < mp->natoms; i++) {
				double dval[3];
				Vector r;
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					err = 1;
					snprintf(errbuf, errbufsize, "line %d: premature end of file while reading coordinates (frame %d)", lineNumber, fn);
					goto exit;
				}
				if (sscanf(buf, "%lg %lg %lg", dval, dval + 1, dval + 2) != 3) {
					err = 1;
					snprintf(errbuf, errbufsize, "line %d: coordinates cannot be read for atom %d", lineNumber, i + 1);
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
			if (mp->atoms != NULL)
				free(mp->atoms);
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
					snprintf(errbuf, errbufsize, "line %d: premature end of file while reading atoms", lineNumber);
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
			if (mp->bonds != NULL)
				free(mp->bonds);
			if (NewArray(&mp->bonds, &mp->nbonds, sizeof(Int) * 2, nbonds) == NULL)
				goto panic;
			bp = mp->bonds;
			for (i = 0; i < nbonds; i += 4) {
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					snprintf(errbuf, errbufsize, "line %d: premature end of file while reading bonds", lineNumber);
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
						snprintf(errbuf, errbufsize, "line %d: The bond %d-%d includes non-existent atom", lineNumber, b1+1, b2+1);
						err = 1;
						goto exit;
					}
					*bp++ = b1;
					*bp++ = b2;
					ap = ATOM_AT_INDEX(mp->atoms, b1);
					if (ap->nconnects < ATOMS_MAX_CONNECTS)
						ap->connects[ap->nconnects++] = b2;
					else {
						snprintf(errbuf, errbufsize, "line %d: The atom %d has more than %d bonds", lineNumber, b1+1, ATOMS_MAX_CONNECTS);
						err = 1;
						goto exit;
					}
					ap = ATOM_AT_INDEX(mp->atoms, b2);
					if (ap->nconnects < ATOMS_MAX_CONNECTS)
						ap->connects[ap->nconnects++] = b1;
					else {
						snprintf(errbuf, errbufsize, "line %d: The atom %d has more than %d bonds", lineNumber, b2+1, ATOMS_MAX_CONNECTS);
						err = 1;
						goto exit;
					}
				}
			}
			continue;
		} else if (section == 4) {
			/*  Angles  */
			Int nangles;
			Int *gp;
			ReadFormat(buf, "I8", &nangles);
			if (mp->angles != NULL)
				free(mp->angles);
			if (NewArray(&mp->angles, &mp->nangles, sizeof(Int) * 3, nangles) == NULL)
				goto panic;
			gp = mp->angles;
			for (i = 0; i < nangles; i += 3) {
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					snprintf(errbuf, errbufsize, "line %d: premature end of file while reading angles", lineNumber);
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
						snprintf(errbuf, errbufsize, "line %d: The angle %d-%d-%d includes non-existent atom", lineNumber, a1+1, a2+1, a3+1);
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
					snprintf(errbuf, errbufsize, "line %d: premature end of file while reading %s", lineNumber, (section == 5 ? "dihedral" : "improper"));
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
						snprintf(errbuf, errbufsize, "line %d: The %s %d-%d-%d-%d angle includes non-existent atom", lineNumber, (section == 5 ? "dihedral" : "improper"), d1+1, d2+1, d3+1, d4+1);
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
			ap->frames = (Vector *)malloc(sizeof(Vector) * fn);
			if (ap->frames == NULL)
				goto panic;
			ap->nframes = fn;
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
				tr[i * 3] = sn;
				symop++;
			} else if (*symop == 'y' || *symop == 'Y') {
				tr[i * 3 + 1] = sn;
				symop++;
			} else if (*symop == 'z' || *symop == 'Z') {
				tr[i * 3 + 2] = sn;
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
MoleculeLoadTepFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	char buf[1024];
	int section = -1;
	int lineNumber;
	int cellType;
	Int ibuf[12];
	Double fbuf[12];
	Int *bonds, nbonds;
	if (errbuf == NULL) {
		errbuf = buf;
		errbufsize = 1024;
	}
	errbuf[0] = 0;
	if (mp == NULL)
		mp = MoleculeNew();
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot open file");
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
				tr[1] = fbuf[2];
				tr[2] = fbuf[3];
				tr[3] = fbuf[5];
				tr[4] = fbuf[6];
				tr[5] = fbuf[7];
				tr[6] = fbuf[9];
				tr[7] = fbuf[10];
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
					snprintf(errbuf, errbufsize, "line %d: bad symmetry specification", lineNumber);
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
				snprintf(errbuf, errbufsize, "unexpected end of file");
				return 1;
			}
			ReadFormat(buf, "I1F8F9F9F9F9F9F9", ibuf, fbuf, fbuf+1, fbuf+2, fbuf+3, fbuf+4, fbuf+5, fbuf+6);
			atomType = fbuf[6];
			if ((atomType >= 0 && atomType <= 5) || (atomType >= 8 && atomType <= 10)) { 
				/*  Anisotropic thermal parameters  */
				MoleculeSetAniso(mp, atomIndex, atomType, fbuf[0], fbuf[1], fbuf[2], fbuf[3], fbuf[5], fbuf[4]);
			}
			if (ibuf[0] != 0)
				section = 3;
			continue;
		}
	}
	fclose(fp);
	MoleculeGuessBonds(mp, 1.2, &nbonds, &bonds);
	if (nbonds > 0) {
		MoleculeAddBonds(mp, nbonds, bonds);
		free(bonds);
	}
	mp->nframes = -1;  /*  Should be recalculated later  */
	return 0;
  panic:
	Panic("low memory while reading structure file %s", fname);
	return -1; /* not reached */
}

int
MoleculeLoadShelxFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
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
	Int nsfacs = 0;
	Int nbonds, *bonds;
	char (*sfacs)[4] = NULL;

	if (errbuf == NULL) {
		errbuf = buf;
		errbufsize = 1024;
	}
	errbuf[0] = 0;
	if (mp == NULL)
		mp = MoleculeNew();
	currentResName[0] = 0;
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot open file");
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
				snprintf(errbuf, errbufsize, "line %d: bad symmetry specification", lineNumber);
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
						snprintf(errbuf, errbufsize, "line %d: unexpected end of file within the atom cards", lineNumber);
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
				if (n == 5)
					ap->tempFactor = fbuf[4] * 78.9568352087147; /* 8*pi*pi */
				else
					MoleculeSetAniso(mp, atomIndex, 8, fbuf[4], fbuf[5], fbuf[6], fbuf[9], fbuf[7], fbuf[8]);
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
		static Transform tr_r1 = {0, -1, 0, 1, -1, 0, 0, 0, 1, 0, 0, 0};
		static Transform tr_r2 = {-1, 1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0};
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
	
	MoleculeGuessBonds(mp, 1.2, &nbonds, &bonds);
	if (nbonds > 0) {
		MoleculeAddBonds(mp, nbonds, bonds);
		free(bonds);
	}
	mp->nframes = -1;  /*  Should be recalculated later  */
	return 0;
  panic:
	Panic("low memory while reading structure file %s", fname);
	return -1; /* not reached */
}

static void
sSeparateTokens(char *inString, char **outPtr, int size)
{
	char *p;
	int i;
	for (i = 0; i < size; i++) {
		p = strtok((i == 0 ? inString : NULL), " \n");
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
			}
		}
	}
	return 0;
}

int
MoleculeLoadGaussianFchkFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	char buf[1024];
	int lineNumber;
	int natoms, nbasis, i, j, k, n, mxbond, retval, nmos, nprims;
	BasisSet *bset;
	ShellInfo *sp;
	PrimInfo *pp;
	Int nary;
	Int *iary;
	Double *dary;
	Atom *ap;
	Vector *vp;
	Double w;

	if (errbuf == NULL) {
		errbuf = buf;
		errbufsize = 1024;
	}
	errbuf[0] = 0;
	if (mp == NULL)
		mp = MoleculeNew();
	bset = (BasisSet *)calloc(sizeof(BasisSet), 1);
	if (bset == NULL)
		goto panic;
	mp->bset = bset;
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot open file");
		return 1;
	}
	lineNumber = 0;
	natoms = nbasis = -1;
	mxbond = 0;
	nmos = 0;
	nprims = 0;
	nary = 0;
	iary = NULL;
	dary = NULL;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		char *tokens[16];
		char *p = buf + 41;
		while (p > buf && *p == ' ')
			p--;
		p[1] = 0;
		sSeparateTokens(buf + 42, tokens, 16);
		if (strcmp(buf, "Number of atoms") == 0) {
			if (tokens[1] == NULL || (natoms = atoi(tokens[1])) <= 0) {
				snprintf(errbuf, errbufsize, "Line %d: strange number of atoms: %s", lineNumber, tokens[1]);
				retval = 2;
				goto cleanup;
			}
			/*  Allocate atom records (all are empty for now)  */
			AssignArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, natoms - 1, NULL);
			/*  Also allocate atom position array for MO calculations  */
			AssignArray(&bset->pos, &bset->natoms, sizeof(Vector), natoms - 1, NULL);
			/*  Also allocate nuclear charge array  */
			bset->nuccharges = (Double *)calloc(sizeof(Double), natoms);
		} else if (strcmp(buf, "Number of electrons") == 0) {
			if (tokens[1] == NULL || (bset->nelectrons = atoi(tokens[1])) <= 0) {
				snprintf(errbuf, errbufsize, "Line %d: strange number of electrons: %s", lineNumber, tokens[1]);
				retval = 2;
				goto cleanup;
			}
		} else if (strcmp(buf, "Number of basis functions") == 0) {
			if (tokens[1] == NULL || (nbasis = atoi(tokens[1])) <= 0) {
				snprintf(errbuf, errbufsize, "Line %d: strange number of electrons: %s", lineNumber, tokens[1]);
				retval = 2;
				goto cleanup;
			}
		} else if (strcmp(buf, "Atomic numbers") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != natoms) {
				snprintf(errbuf, errbufsize, "Line %d: wrong or inconsistent number of atoms: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&iary, &nary, sizeof(Int), natoms, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read atomic numbers", lineNumber);
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
				snprintf(errbuf, errbufsize, "Line %d: wrong or inconsistent number of atoms: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), natoms, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read nuclear charges", lineNumber);
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
				snprintf(errbuf, errbufsize, "Line %d: wrong or inconsistent number of cartesian coordinates: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), natoms * 3, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read cartesian coordinates", lineNumber);
				retval = 2;
				goto cleanup;
			}
			for (i = 0, ap = mp->atoms, vp = bset->pos; i < natoms; i++, ap = ATOM_NEXT(ap), vp++) {
				vp->x = dary[i * 3];
				vp->y = dary[i * 3 + 1];
				vp->z = dary[i * 3 + 2];
				ap->r.x = vp->x * kBohr2Angstrom;
				ap->r.y = vp->y * kBohr2Angstrom;
				ap->r.z = vp->z * kBohr2Angstrom;
			}
			free(dary);
			dary = NULL;
		} else if (strcmp(buf, "MxBond") == 0) {
			if (tokens[1] == NULL || (mxbond = atoi(tokens[1])) <= 0) {
				snprintf(errbuf, errbufsize, "Line %d: strange number of bonds per atom: %s", lineNumber, tokens[1]);
				retval = 2;
				goto cleanup;
			}
		} else if (strcmp(buf, "IBond") == 0) {
			Int bonds[ATOMS_MAX_CONNECTS * 2 + 1], lim;
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != natoms * mxbond) {
				snprintf(errbuf, errbufsize, "Line %d: wrong or inconsistent number of bonds: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&iary, &nary, sizeof(Int), natoms * mxbond, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read bond information", lineNumber);
				retval = 2;
				goto cleanup;
			}
			lim = (mxbond > ATOMS_MAX_CONNECTS ? ATOMS_MAX_CONNECTS : mxbond);
			for (i = 0; i < natoms; i++) {
				for (j = k = 0; j < lim; j++) {
					n = iary[i * mxbond + j] - 1;
					if (n > i) {
						/*  Connect atom i and atom n  */
						bonds[k++] = i;
						bonds[k++] = n;
					}
				}
				if (k > 0) {
					bonds[k] = kInvalidIndex;
					MoleculeAddBonds(mp, k / 2, bonds);
				}
			}
			free(iary);
			iary = NULL;
		} else if (strcmp(buf, "Shell types") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0) {
				snprintf(errbuf, errbufsize, "Line %d: wrong number of shell types: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&iary, &nary, sizeof(Int), i, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read shell types", lineNumber);
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
						/*  TODO: Support F/F7 type orbitals  */
						/*	case 3: sp->sym = kGTOtype_F;  sp->ncomp = 10; break;
						 case -3: sp->sym = kGTOType_F7; sp->ncomp = 7; break; */
					default:
						snprintf(errbuf, errbufsize, "Line %d: unsupported shell type %d", lineNumber, iary[i]);
						retval = 2;
						goto cleanup;
				}
				sp->m_idx = n;
				n += sp->ncomp;
			}
			nmos = n;
			free(iary);
			iary = NULL;
		} else if (strcmp(buf, "Number of primitives per shell") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != bset->nshells) {
				snprintf(errbuf, errbufsize, "Line %d: wrong size of the primitive table: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&iary, &nary, sizeof(Int), i, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read primitive table", lineNumber);
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
				snprintf(errbuf, errbufsize, "Line %d: wrong size of the shell-to-atom map: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&iary, &nary, sizeof(Int), i, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read shell-to-atom table", lineNumber);
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
				snprintf(errbuf, errbufsize, "Line %d: wrong number of primitive exponents: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read primitive exponents", lineNumber);
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
				snprintf(errbuf, errbufsize, "Line %d: wrong number of contraction coefficients: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read contraction coefficients", lineNumber);
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
				snprintf(errbuf, errbufsize, "Line %d: wrong number of P(S=P) contraction coefficients: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&dary, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read P(S=P) contraction coefficients", lineNumber);
				retval = 2;
				goto cleanup;
			}
			for (i = 0, pp = bset->priminfos; i < bset->npriminfos; i++, pp++) {
				pp->Csp = dary[i];
			}
			free(dary);
			dary = NULL;
		} else if (strcmp(buf, "Alpha Orbital Energies") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != nmos) {
				snprintf(errbuf, errbufsize, "Line %d: wrong or inconsistent number of alpha orbitals: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&bset->moenergies, &bset->nmos, sizeof(Double), i, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read alpha orbital energies", lineNumber);
				retval = 2;
				goto cleanup;
			}
		} else if (strcmp(buf, "Alpha MO coefficients") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != bset->nmos * bset->nmos) {
				snprintf(errbuf, errbufsize, "Line %d: wrong or inconsistent number of MO coefficients: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&bset->mo, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read MO coefficients", lineNumber);
				retval = 2;
				goto cleanup;
			}
		} else if (strcmp(buf, "Total SCF Density") == 0) {
			if (tokens[2] == NULL || (i = atoi(tokens[2])) <= 0 || i != bset->nmos * (bset->nmos + 1) / 2) {
				snprintf(errbuf, errbufsize, "Line %d: wrong or inconsistent number of SCF densities: %s", lineNumber, tokens[2]);
				retval = 2;
				goto cleanup;
			}
			if (sReadNumberArray(&bset->scfdensities, &nary, sizeof(Double), i, fp, &lineNumber) != 0) {
				snprintf(errbuf, errbufsize, "Line %d: cannot read SCF densities", lineNumber);
				retval = 2;
				goto cleanup;
			}
		}
	}
	if (mp->natoms == 0) {
		snprintf(errbuf, errbufsize, "Atom information is missing");
		retval = 2;
		goto cleanup;
	}
	if (bset->shells == NULL || bset->priminfos == NULL) {
		snprintf(errbuf, errbufsize, "Gaussian primitive information is missing");
		retval = 2;
		goto cleanup;
	}
	if (bset->mo == NULL) {
		snprintf(errbuf, errbufsize, "MO coefficients were not found");
		retval = 2;
		goto cleanup;
	}
	if (sSetupGaussianCoefficients(bset) != 0) {
		snprintf(errbuf, errbufsize, "Internal error during setup MO calculation");
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
MoleculeLoadGamessDatFile(Molecule *mol, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	int newmol = 0;
	char buf[1024];
	int lineNumber, i, natoms = 0;
	int nframes = 0;
	int n1;
	int ival[8];
	double dval[8];
	char sval[16];
	Vector *vbuf = NULL;
	IntGroup *ig;
	if (errbuf == NULL) {
		errbuf = buf;
		errbufsize = 1024;
	}
	errbuf[0] = 0;
	if (mol == NULL) {
		mol = MoleculeNew();
	}
	if (mol->natoms == 0)
		newmol = 1;

	fp = fopen(fname, "rb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot open file");
		return 1;
	}
	
	/*  ESP is cleared (not undoable!)  */
	if (mol->elpots != NULL) {
		free(mol->elpots);
		mol->elpots = NULL;
		mol->nelpots = 0;
	}
	
	lineNumber = 0;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
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
			while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
				if (strncmp(buf, " $END", 5) == 0)
					break;
				if (sscanf(buf, "%12s %lf %lf %lf %lf", sval, &dval[0], &dval[1], &dval[2], &dval[3]) < 5) {
					snprintf(errbuf, errbufsize, "Line %d: bad format in $DATA section", lineNumber);
					return 2;
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
						snprintf(errbuf, errbufsize, "Line %d: too many atoms", lineNumber);
						return 3;
					}
					if ((ap = ATOM_AT_INDEX(mol->atoms, i))->atomicNumber != dval[0]) {
						snprintf(errbuf, errbufsize, "Line %d: atomic number does not match", lineNumber);
						return 4;
					}
					vbuf[i].x = dval[1];
					vbuf[i].y = dval[2];
					vbuf[i].z = dval[3];
				}
				/*  Skip until a blank line is found  */
				while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
					int j;
					for (j = 0; buf[j] == ' '; j++);
					if (buf[j] == '\n')
						break;
				}
				i++;
			}
			natoms = i;
			if (!newmol) {
				/*  Set atom positions  */
				IntGroup *ig;
				if (natoms < mol->natoms) {
					snprintf(errbuf, errbufsize, "Line %d: too few atoms", lineNumber);
					return 5;
				}
				ig = IntGroupNewWithPoints(0, natoms, -1);
				MolActionCreateAndPerform(mol, gMolActionSetAtomPositions, ig, natoms, vbuf);
				IntGroupRelease(ig);
			}
			if (vbuf == NULL)
				vbuf = (Vector *)calloc(sizeof(Vector), natoms);
			nframes = MoleculeGetNumberOfFrames(mol);
			continue;
		} else if (strstr(buf, "DATA FROM NSERCH") != NULL || (strstr(buf, "RESULTS FROM SUCCESSFUL") != NULL && (n1 = 1))) {
			/*  Skip until the separator line is read (three or four lines)  */
			i = 0;
			do {
				if (i++ >= 4) {
					snprintf(errbuf, errbufsize, "Line %d: the separator line at the top of the coordinates is not found: bad format?", lineNumber);
					return 6;
				}
				ReadLine(buf, sizeof buf, fp, &lineNumber);
			} while (strstr(buf, "----------------------------") == NULL);
			for (i = 0; i < natoms; i++) {
				if (ReadLine(buf, sizeof buf, fp, &lineNumber) <= 0) {
					snprintf(errbuf, errbufsize, "Unexpected end of file in reading NSERCH data");
					return 6;
				}
				if (sscanf(buf, "%12s %lf %lf %lf %lf", sval, &dval[0], &dval[1], &dval[2], &dval[3]) < 5) {
					snprintf(errbuf, errbufsize, "Line %d: bad format in NSERCH coordinate data", lineNumber);
					return 7;
				}
				vbuf[i].x = dval[1];
				vbuf[i].y = dval[2];
				vbuf[i].z = dval[3];
			}
			ig = IntGroupNewWithPoints(nframes, 1, -1);
			MolActionCreateAndPerform(mol, gMolActionInsertFrames, ig, natoms, vbuf, 0, NULL);
			IntGroupRelease(ig);
			nframes++;
			if (n1) {
				/*  TODO: read the VEC group  */
			}
			continue;
		} else if ((strstr(buf, "ELECTRIC POTENTIAL") != NULL || strstr(buf, "ELECTROSTATIC POTENTIAL") != NULL) && strstr(buf, "ELPOTT") != NULL) {
			i = 0;
			while ((n1 = ReadLine(buf, sizeof buf, fp, &lineNumber)) > 0) {
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
			if (n1 > 0)
				goto redo;  /*  This section has no end line, so the last line should be processed again  */
			else break;    /*  End of file encountered  */
		}  /*  TODO: read MOLPLT info if present  */
	}
	if (vbuf != NULL)
		free(vbuf);
	if (newmol && mol->nbonds == 0) {
		/*  Guess bonds  */
		Int nbonds, *bonds;
		MoleculeGuessBonds(mol, 1.2, &nbonds, &bonds);
		if (nbonds > 0) {
			MolActionCreateAndPerform(mol, gMolActionAddBonds, nbonds * 2, bonds);
			free(bonds);
		}
	}
	return 0;
}

int
MoleculeReadCoordinatesFromFile(Molecule *mp, const char *fname, const char *ftype, char *errbuf, int errbufsize)
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
		retval = MoleculeReadCoordinatesFromPdbFile(mp, fname, errbuf, errbufsize);
	}
	if (retval != 0) {
		/*  Try all formats once again  */
		retval = MoleculeReadCoordinatesFromPdbFile(mp, fname, errbuf, errbufsize);
	}
	return retval;
}

int
MoleculeReadCoordinatesFromPdbFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	char buf[1024];
	char *p;
	int lineNumber;
	int i, j, new_unit, retval;
	Atom *ap;
	Int ibuf[12];
	Int entries = 0;
	retval = 0;
	if (errbuf == NULL) {
		errbuf = buf;
		errbufsize = 1024;
	}
	errbuf[0] = 0;
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot open file");
		return -1;
	}
/*	flockfile(fp); */
	if (mp->natoms == 0)
		new_unit = 1;
	else new_unit = 0;
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
				snprintf(errbuf, errbufsize, "line %d: non-positive atom number %d", lineNumber, w.serial);
				retval = 1;
				goto abort;
			}
			w.serial--;  /*  The internal atom number is 0-based  */
			if (w.serial >= mp->natoms) {
				if (new_unit) {
					ap = AssignArray(&mp->atoms, &mp->natoms, gSizeOfAtomRecord, w.serial, NULL);
				} else {
					snprintf(errbuf, errbufsize, "line %d: the atom number %d does not exist in the structure file", lineNumber, w.serial+1);
					retval = 1;
					goto abort;
				}
			}
			ap = ATOM_AT_INDEX(mp->atoms, w.serial);
			ap->r = w.r;
			ap->occupancy = w.occ;
			ap->tempFactor = w.temp;
			if (new_unit) {
				if (w.segName[0] == 0)
					strncpy(w.segName, "MAIN", 4);
				strncpy(ap->segName, w.segName, 4);
				ap->resSeq = w.resSeq;
				strncpy(ap->resName, w.resName, 4);
				strncpy(ap->aname, w.atomName, 4);
				strncpy(ap->element, w.element, 2);
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
						snprintf(errbuf, errbufsize, "line %d: The CONECT record contains non-existent atom %d", lineNumber, ibuf[j]);
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
						bbuf[bi * 2] = ibuf[0] - 1;
						bbuf[bi * 2 + 1] = ibuf[j] - 1;
						bi++;
					}
				}
				if (bi == 0)
					continue;
				bbuf[bi * 2] = -1;
				retval = MoleculeAddBonds(mp, bi, bbuf);
				if (retval < 0) {
					snprintf(errbuf, errbufsize, "line %d: bad bond specification", lineNumber);
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
			snprintf(errbuf, errbufsize, "Out of memory");
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
			for (i = 0; i < mp->natoms; i++) {
				ap = ATOM_AT_INDEX(mp->atoms, i);
				for (j = 0; j < ap->nconnects; j++) {
					ap->connects[j] = old2new[ap->connects[j]];
				}
			}
			for (i = 0; i < mp->nbonds * 2; i++) {
				mp->bonds[i] = old2new[mp->bonds[i]];
			}
		}
		retval = MoleculeRebuildTablesFromConnects(mp);
		if (retval != 0) {
			/*  This error may not happen  */
			snprintf(errbuf, errbufsize, "Cannot build angle/dihedral/improper tables");
			retval = 1;
			goto abort;
		}

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
	if (entries == 0)
		return 1;  /*  Maybe different format?  */
	return retval;
}

int
MoleculeReadCoordinatesFromDcdFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	DcdRecord dcd;
	SFloat32 *xp, *yp, *zp;
	Vector *vp;
	IntGroup *ig;
	int n;
	errbuf[0] = 0;
	if (mp == NULL || mp->natoms == 0) {
		snprintf(errbuf, errbufsize, "Molecule is empty");
		return 1;
	}
	n = DcdOpen(fname, &dcd);
	if (n != 0) {
		switch (n) {
			case -2: snprintf(errbuf, errbufsize, "Cannot open file"); break;
			case 1:  snprintf(errbuf, errbufsize, "Premature EOF encountered"); break;
			case 2:  snprintf(errbuf, errbufsize, "Bad block length of the first section"); break;
			case 3:  snprintf(errbuf, errbufsize, "\"CORD\" signature is missing"); break;
			case 4:  snprintf(errbuf, errbufsize, "Bad termination of the first section"); break;
			case 5:  snprintf(errbuf, errbufsize, "The title section is not correct"); break;
			case 6:  snprintf(errbuf, errbufsize, "The atom number section is not correct"); break;
			default: snprintf(errbuf, errbufsize, "Read error in dcd file"); break;
		}
	} else {
		if (dcd.natoms == 0)
			snprintf(errbuf, errbufsize, "No atoms were found in the dcd file");
		else if (dcd.nframes == 0)
			snprintf(errbuf, errbufsize, "No frames were found in the dcd file");
	}
	if (errbuf[0] != 0) {
		if (n == 0)
			DcdClose(&dcd);
		return 1;
	}

	vp = (Vector *)calloc(sizeof(Vector), mp->natoms * dcd.nframes);
	xp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	yp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	zp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	ig = IntGroupNewWithPoints(MoleculeGetNumberOfFrames(mp), dcd.nframes, -1);
	if (vp == NULL || xp == NULL || yp == NULL || zp == NULL || ig == NULL) {
		snprintf(errbuf, errbufsize, "Cannot allocate memory");
		if (vp) free(vp);
		if (xp) free(xp);
		if (yp) free(yp);
		if (zp) free(zp);
		if (ig) IntGroupRelease(ig);
		return 1;
	}
	for (n = 0; n < dcd.nframes; n++) {
		int i;
		Vector *vpp;
		if (DcdReadFrame(&dcd, n, xp, yp, zp, dcd.globalcell)) {
			snprintf(errbuf, errbufsize, "Read error in dcd file");
			goto exit;
		}
		for (i = 0, vpp = &vp[n * mp->natoms]; i < dcd.natoms && i < mp->natoms; i++, vpp++) {
			vpp->x = xp[i];
			vpp->y = yp[i];
			vpp->z = zp[i];
		}
	}
	/*  TODO: implement frame-specific cells  */
	if (MoleculeInsertFrames(mp, ig, vp, NULL) < 0)
		snprintf(errbuf, errbufsize, "Cannot insert frames");
	if (dcd.with_unitcell) {
		Vector ax, ay, az, orig;
		char flags[3] = {1, 1, 1};
		ax.x = dcd.globalcell[0]; ax.y = ax.z = 0;
		ay.y = dcd.globalcell[1]; ay.x = ay.z = 0;
		az.z = dcd.globalcell[2]; az.x = az.y = 0;
		orig.x = dcd.globalcell[3];
		orig.y = dcd.globalcell[4];
		orig.z = dcd.globalcell[5];
		if (MoleculeSetPeriodicBox(mp, &ax, &ay, &az, &orig, flags) != 0) {
			snprintf(errbuf, errbufsize, "Cannot set unit cell");
			goto exit;
		}
	}
	mp->startStep = dcd.nstart;
	mp->stepsPerFrame = dcd.ninterval;
	mp->psPerStep = dcd.delta;
exit:
	DcdClose(&dcd);
	free(vp);
	free(xp);
	free(yp);
	free(zp);
	IntGroupRelease(ig);
	if (errbuf[0] == 0)
		return 0;
	else return 1;
}

int
MoleculeReadExtendedInfo(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	char buf[1024];
	int lineNumber;
	int i, retval;
	Vector v[3], vv;
	double d[3];
	int n, flag;
	char flags[3];
	fp = fopen(fname, "rb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot open file");
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
					snprintf(errbuf, errbufsize, "line %d: missing %d component of the bounding box", lineNumber, i + 1);
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
			MoleculeSetPeriodicBox(mp, &v[0], &v[1], &v[2], &vv, flags);
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
				snprintf(errbuf, errbufsize, "line %d: wrong format for the bounding box origin", lineNumber);
				retval = 1;
				goto abort;
			}
			vv.x = d[0];
			vv.y = d[1];
			vv.z = d[2];
			MoleculeSetPeriodicBox(mp, &v[0], &v[1], &v[2], &vv, flags);
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
MoleculeWriteToFile(Molecule *mp, const char *fname, const char *ftype, char *errbuf, int errbufsize)
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
		retval = MoleculeWriteToPsfFile(mp, fname, errbuf, errbufsize);
	} else if (strcasecmp(ftype, "pdb") == 0) {
		retval = MoleculeWriteToPdbFile(mp, fname, errbuf, errbufsize);
	} else if (strcasecmp(ftype, "tep") == 0) {
		retval = MoleculeWriteToTepFile(mp, fname, errbuf, errbufsize);
	} else {
		snprintf(errbuf, errbufsize, "The file format should be specified");
		retval = 1;
	}

	return retval;
}

int
MoleculeWriteToMbsfFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	int i, j, k, n1, n2;
	Atom *ap;
	char bufs[6][8];

	fp = fopen(fname, "wb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot write to file %s", fname);
		return 1;
	}
	errbuf[0] = 0;

	fprintf(fp, "!:atoms\n");
	fprintf(fp, "! idx seg_name res_seq res_name name type weight element atomic_number occupancy temp_factor int_charge\n");
	n1 = n2 = 0;
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		strncpy(bufs[0], ap->segName, 4);
		strncpy(bufs[1], ap->resName, 4);
		strncpy(bufs[2], ap->aname, 4);
		AtomTypeDecodeToString(ap->type, bufs[3]);
		strncpy(bufs[4], ap->element, 4);
		for (j = 0; j < 5; j++) {
			bufs[j][4] = 0;
			if (bufs[j][0] == 0) {
				bufs[j][0] = '_';
				bufs[j][1] = 0;
			}
			for (k = 0; k < 4; k++) {
				if (bufs[j][k] > 0 && bufs[j][k] < ' ')
					bufs[j][k] = '_';
			}
		}
		if (SYMOP_ALIVE(ap->symop))
			n1++;
		if (ap->fix_force != 0)
			n2++;
		fprintf(fp, "%d %s %d %s %s %s %.5f %.5f %s %d %f %f %d\n", i, bufs[0], ap->resSeq, bufs[1], bufs[2], bufs[3], ap->charge, ap->weight, bufs[4], ap->atomicNumber, ap->occupancy, ap->tempFactor, ap->intCharge);
	}
	fprintf(fp, "\n");
	
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
	
	if ((n1 = MoleculeGetNumberOfFrames(mp)) > 0)
		n2 = mp->cframe;
	else
		n2 = 0;
	for (i = 0; (i == n2 || i < n1); i++) {
		fprintf(fp, "!:positions ; frame %d\n", i);
		fprintf(fp, "! idx x y z\n");
		for (j = 0, ap = mp->atoms; j < mp->natoms; j++, ap = ATOM_NEXT(ap)) {
			Vector *vp;
			if (i != n2 && i < ap->nframes)
				vp = ap->frames + i;
			else
				vp = &(ap->r);
			fprintf(fp, "%d %.8f %.8f %.8f\n", j, vp->x, vp->y, vp->z);
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
		fprintf(fp, "!:periodic_box\n");
		fprintf(fp, "! ax ay az; bx by bz; cx cy cz; ox oy oz; fa fb fc\n");
		for (i = 0; i < 3; i++)
			fprintf(fp, "%15.8f %15.8f %15.8f\n", mp->cell->axes[i].x, mp->cell->axes[i].y, mp->cell->axes[i].z);
		fprintf(fp, "%15.8f %15.8f %15.8f\n", mp->cell->origin.x, mp->cell->origin.y, mp->cell->origin.z);
		fprintf(fp, "%d %d %d\n", mp->cell->flags[0], mp->cell->flags[1], mp->cell->flags[2]);
		fprintf(fp, "\n");

		fprintf(fp, "!:xtalcell\n");
		fprintf(fp, "! a b c alpha beta gamma; this info is redundant, periodic_box is used instead\n");
		fprintf(fp, "%f %f %f %f %f %f\n", mp->cell->cell[0], mp->cell->cell[1], mp->cell->cell[2], mp->cell->cell[3], mp->cell->cell[4], mp->cell->cell[5]);
		fprintf(fp, "\n");
	}
	
	if (mp->nsyms > 0) {
		fprintf(fp, "!:symmetry_operations\n");
		fprintf(fp, "! a11 a12 a13; a21 a22 a23; a31 a32 a33; t1 t2 t3\n");
		for (i = 0; i < mp->nsyms; i++) {
			Transform *tp = mp->syms + i;
			for (j = 0; j < 12; j++)
				fprintf(fp, "%11.6f%c", (*tp)[j], (j % 3 == 2 ? '\n' : ' '));
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

	fclose(fp);
	return 0;
}

int
MoleculeWriteToPsfFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	int i;
	Atom *ap;
	fp = fopen(fname, "wb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot write to file %s", fname);
		return 1;
	}
	errbuf[0] = 0;
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
#if 0
		if (mp->nframes > 0) {
			int fn;  /*  Frame number  */
			for (fn = 0; fn < ap->nframes; fn++) {
				fprintf(fp, "%8d !COORD: coordinates for frame %d\n", mp->natoms, fn);
				for (i = 0; i < mp->natoms; i++) {
					Vector r;
					ap = ATOM_AT_INDEX(mp->atoms, i);
					if (ap->frames == NULL || fn >= ap->nframes)
						r = ap->r;
					else
						r = ap->frames[fn];
					fprintf(fp, " %.8g %.8g %.8g ! %d,%.4s\n", r.x, r.y, r.z, i + 1, ap->name);
				}
				fprintf(fp, "\n");
			}
		}
#endif
	}
		
	fclose(fp);
	return 0;
}

int
MoleculeWriteToPdbFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	int i, j;
	Atom *ap;
	fp = fopen(fname, "wb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot write to file %s", fname);
		return 1;
	}
	errbuf[0] = 0;
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
		ap = ATOM_AT_INDEX(mp->atoms, i);
		for (j = 0; j < ap->nconnects; j++) {
			if (j % 4 == 0) {
				if (j > 0)
					fprintf(fp, "\n");
				fprintf(fp, "CONECT%5d", i + 1);
			}
			fprintf(fp, "%5d", ap->connects[j] + 1);
		}
		if (j > 0)
			fprintf(fp, "\n");
	}
	fprintf(fp, "END\n");
	fclose(fp);
	return 0;
}

int
MoleculeWriteToDcdFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	DcdRecord dcd;
	SFloat32 *xp, *yp, *zp;
	int n;
	errbuf[0] = 0;
	if (mp == NULL || mp->natoms == 0) {
		snprintf(errbuf, errbufsize, "Molecule is empty");
		return 1;
	}
	memset(&dcd, 0, sizeof(dcd));
	dcd.natoms = mp->natoms;
	dcd.nframes = MoleculeGetNumberOfFrames(mp);
	if (dcd.nframes == 0) {
		snprintf(errbuf, errbufsize, "no frame is present");
		return 1;
	}
	dcd.nstart = mp->startStep;
	dcd.ninterval = mp->stepsPerFrame;
	if (dcd.ninterval == 0)
		dcd.ninterval = 1;
	dcd.nend = dcd.nstart + (dcd.nframes - 1) * dcd.ninterval;
	if (mp->cell != NULL) {
		dcd.with_unitcell = 1;
		dcd.globalcell[0] = VecLength(mp->cell->axes[0]);
		dcd.globalcell[1] = VecLength(mp->cell->axes[1]);
		dcd.globalcell[2] = VecLength(mp->cell->axes[2]);
		dcd.globalcell[3] = mp->cell->origin.x;
		dcd.globalcell[4] = mp->cell->origin.y;
		dcd.globalcell[5] = mp->cell->origin.z;
	}
	dcd.delta = mp->psPerStep;
	if (dcd.delta == 0.0)
		dcd.delta = 1.0;
	n = DcdCreate(fname, &dcd);
	if (n != 0) {
		if (n < 0)
			snprintf(errbuf, errbufsize, "Cannot create dcd file");
		else
			snprintf(errbuf, errbufsize, "Cannot write dcd header");
		DcdClose(&dcd);
		return 1;
	}
	
	xp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	yp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	zp = (SFloat32 *)malloc(sizeof(SFloat32) * dcd.natoms);
	if (xp == NULL || yp == NULL || zp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot allocate memory");
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
		if (DcdWriteFrame(&dcd, n, xp, yp, zp, dcd.globalcell)) {
			snprintf(errbuf, errbufsize, "Write error in dcd file");
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
MoleculeWriteExtendedInfo(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	int i;
	Vector v;
	errbuf[0] = 0;
	fp = fopen(fname, "wb");
	if (fp == NULL) {
		snprintf(errbuf, errbufsize, "Cannot write to file %s", fname);
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
					for (k = ap->nconnects - 1; k >= 0; k--) {
						if (ap->connects[k] == n2)
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
				for (n1 = n[i], ap = atoms + n1; n1 < n[i + 1]; n1++, ap++) {
					for (k = ap->nconnects - 1; k >= 0; k--) {
						n2 = ap->connects[k];
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
	
#if 0
{
	/*  Explicit bond table, sorted by bond type  */
	for (i = j = 0; i < mp->nbonds; i++) {
		n1 = mp->bonds[i * 2];
		n2 = mp->bonds[i * 2 + 1];
		ap1 = ATOM_AT_INDEX(mp->atoms, n1);
		ap2 = ATOM_AT_INDEX(mp->atoms, n2);
		if ((ap1->exflags & kAtomHiddenFlag) || (ap2->exflags & kAtomHiddenFlag))
			continue;
		if (ap1->atomicNumber > 18 || ap2->atomicNumber > 18) {
			type = 3;
		} else if (ap1->atomicNumber > 1 && ap1->atomicNumber > 1) {
			type = 2;
		} else {
			type = 1;
		}
		ip[j * 3] = type;
		ip[j * 3 + 1] = sMakeAdc(n1, ap1->symbase, ap1->symop);
		ip[j * 3 + 2] = sMakeAdc(n2, ap2->symbase, ap2->symop);
		j++;
	}
	mergesort(ip, j, sizeof(int) * 3, sCompareBondType);
	
	/*  Output instruction cards  */
	strcpy(buf, "  1   811");
	for (i = n1 = 0; i < j; i++) {
		n2 = (n1 % 3) * 18 + 9;
		snprintf(buf + n2, 80 - n2, "%9d%9d\n", ip[i * 3 + 1], ip[i * 3 + 2]);
		if (i == j - 1 || n1 >= 29 || ip[i * 3] != ip[i * 3 + 3]) {
			/*  End of this instruction  */
			buf[2] = '2';
			fputs(buf, fp);
			switch (ip[i * 3]) {
				case 3: rad = 0.06; nshades = 5; break;
				case 2: rad = 0.06; nshades = 1; break;
				default: rad = 0.04; nshades = 1; break;
			}
			fprintf(fp, "                     %3d            %6.3f\n", nshades, rad);
			strcpy(buf, "  1   811");
			n1 = 0;
			continue;
		} else if (n1 % 3 == 2) {
			fputs(buf, fp);
			strcpy(buf, "  1      ");
		}
		n1++;
	}
	free(ip);
}
#endif

int
MoleculeWriteToTepFile(Molecule *mp, const char *fname, char *errbuf, int errbufsize)
{
	FILE *fp;
	int i, j, natoms, *ip;
	Atom *ap, *atoms, **app;
	Double *cp;
	
	errbuf[0] = 0;

	/*  Create sorted array of atoms  */
	natoms = mp->natoms;
	atoms = (Atom *)calloc(sizeof(Atom), natoms);
	app = (Atom **)calloc(sizeof(Atom *), natoms);
	ip = (int *)calloc(sizeof(int), natoms);
	if (atoms == NULL || app == NULL || ip == NULL) {
		snprintf(errbuf, errbufsize, "Cannot allocate memory");
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
	/*  Move the atom record to atoms[]  */
	/*  aniso and frames share the memory with the original  */
	/*  The 'v' member contains crystallographic coordinates  */
	/*  The connection table and symbase are renumbered  */
	/*  Hidden flags are modified to reflect the visibility in the MainView  */
	for (i = 0, ap = atoms; i < natoms; i++, ap++) {
		memmove(ap, app[i], gSizeOfAtomRecord);
		MoleculeCartesianToXtal(mp, &(ap->v), &(ap->r));
		for (j = ap->nconnects - 1; j >= 0; j--) {
			ap->connects[j] = ip[ap->connects[j]];
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
		snprintf(errbuf, errbufsize, "Cannot write to file %s", fname);
		return 1;
	}
	errbuf[0] = 0;

	/*  Title line  */
	fprintf(fp, "Generated by MyMolView\n");
	
	/*  XtalCell  */
	if (mp->cell != NULL) {
		cp = mp->cell->cell;
	} else {
		static Double sUnit[] = {1, 1, 1, 90, 90, 90};
		cp = sUnit;
	}
	fprintf(fp, "%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f\n", cp[0], cp[1], cp[2], cp[3], cp[4], cp[5]);
	
	/*  Symmetry operations  */
	if (mp->nsyms > 0) {
		for (i = 0; i < mp->nsyms; i++) {
			cp = mp->syms[i];
			fprintf(fp, "%c%14g%3g%3g%3g%15g%3g%3g%3g%15g%3g%3g%3g\n", (i == mp->nsyms - 1 ? '1' : ' '), cp[9], cp[0], cp[1], cp[2], cp[10], cp[3], cp[4], cp[5], cp[11], cp[6], cp[7], cp[8]);
		}
	} else {
		fprintf(fp, "1             0  1  0  0              0  0  1  0              0  0  0  1\n");
	}
	
	/*  Atoms  */
	for (i = 0, ap = atoms; i < natoms; i++, ap++) {
		/*  The 'v' field contains crystallographic coordinates  */
		fprintf(fp, " %4.4s%22s%9.4f%9.4f%9.4f%9d\n", ap->aname, "", ap->v.x, ap->v.y, ap->v.z, 0);
		if (ap->aniso != NULL) {
			cp = ap->aniso->bij;
			fprintf(fp, " %8.5f%9.6f%9.6f%9.6f%9.6f%9.6f%9d\n", cp[0], cp[1], cp[2], cp[3], cp[5], cp[4], 0);
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

	fprintf(fp, "      202\n");
	fprintf(fp, "  0    -1\n");
	fclose(fp);
	return 0;
}

void
MoleculeDump(Molecule *mol)
{
	int i, j;
	Atom *ap;
	for (i = 0; i < mol->natoms; i++) {
		char buf1[8];
		ap = ATOM_AT_INDEX(mol->atoms, i);
		snprintf(buf1, sizeof buf1, "%3.4s.%d", ap->resName, ap->resSeq);
		fprintf(stderr, "%4d %-7s %-4.6s %-4.6s %-2.2s %7.3f %7.3f %7.3f %6.3f [", i, buf1, ap->aname, AtomTypeDecodeToString(ap->type, NULL), ap->element, ap->r.x, ap->r.y, ap->r.z, ap->charge);
		for (j = 0; j < ap->nconnects; j++) {
			fprintf(stderr, "%s%d", (j > 0 ? "," : ""), ap->connects[j]);
		}
		fprintf(stderr, "]\n");
	}
}

#pragma mark ====== Serialize ======

Molecule *
MoleculeDeserialize(const char *data, Int length, Int *timep)
{
	Molecule *mp;
/*	int result; */

	mp = MoleculeNew();
	if (mp == NULL)
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
			Atom *ap;
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
			Atom *ap;
			for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
				if (ap->nframes == 0)
					continue;
				ap->frames = (Vector *)malloc(sizeof(Vector) * ap->nframes);
				if (ap->frames == NULL)
					goto out_of_memory;
				memmove(ap->frames, ptr, sizeof(Vector) * ap->nframes);
				ptr += sizeof(Vector) * ap->nframes;
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
		} else if (strcmp(data, "EXATOM") == 0) {
			n = len / sizeof(ExAtom);
			NewArray(&mp->exatoms, &mp->nexatoms, sizeof(ExAtom), n);
			memmove(mp->exatoms, ptr, len);
		} else if (strcmp(data, "EXBOND") == 0) {
			n = len / (sizeof(Int) * 2);
			NewArray(&mp->exbonds, &mp->nexbonds, sizeof(Int) * 2, n);
			memmove(mp->exbonds, ptr, len);
		} else if (strcmp(data, "TIME") == 0) {
			if (timep != NULL)
				*timep = *((Int *)ptr);
		}
		len += 8 + sizeof(Int);
		data += len;
		length -= len;
	}
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
	int len, len_all, i, naniso, nframes;

	/*  Array of atoms  */
	len = 8 + sizeof(Int) + gSizeOfAtomRecord * mp->natoms;
	ptr = (char *)malloc(len);
	if (ptr == NULL)
		goto out_of_memory;
	memmove(ptr, "ATOM\0\0\0\0", 8);
	*((Int *)(ptr + 8)) = gSizeOfAtomRecord * mp->natoms;
	p = ptr + 8 + sizeof(Int);
	memmove(p, mp->atoms, gSizeOfAtomRecord * mp->natoms);
	naniso = nframes = 0;
	for (i = 0; i < mp->natoms; i++) {
		Atom *ap = ATOM_AT_INDEX(p, i);
		if (ap->aniso != NULL) {
			naniso++;
			ap->aniso = NULL;
		}
		if (ap->frames != NULL) {
			nframes += ap->nframes;
			ap->frames = NULL;
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
			Atom *ap = ATOM_AT_INDEX(mp->atoms, i);
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
			Atom *ap = ATOM_AT_INDEX(mp->atoms, i);
			if (ap->frames != NULL) {
				memmove(p, ap->frames, sizeof(Vector) * ap->nframes);
				p += sizeof(Vector) * ap->nframes;
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
	
	/*  Expanded atoms  */
	if (mp->nexatoms > 0) {
		len = 8 + sizeof(Int) + sizeof(ExAtom) * mp->nexatoms;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "EXATOM\0\0", 8);
		*((Int *)(p + 8)) = sizeof(ExAtom) * mp->nexatoms;
		p += 8 + sizeof(Int);
		memmove(p, mp->exatoms, sizeof(ExAtom) * mp->nexatoms);
		for (i = 0; i < mp->nexatoms; i++) {
			/*  Clear label id  */
			((ExAtom *)p)[i].labelid = 0;
		}
		len_all += len;
	}
	
	/*  Expanded bonds  */
	if (mp->nexbonds > 0) {
		len = 8 + sizeof(Int) + sizeof(Int) * 2 * mp->nexbonds;
		ptr = (char *)realloc(ptr, len_all + len);
		if (ptr == NULL)
			goto out_of_memory;
		p = ptr + len_all;
		memmove(p, "EXBOND\0\0", 8);
		*((Int *)(p + 8)) = sizeof(Int) * 2 * mp->nexbonds;
		p += 8 + sizeof(Int);
		memmove(p, mp->exbonds, sizeof(Int) * 2 * mp->nexbonds);
		len_all += len;
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

/*  Guess the bonds from the coordinates  */
int
MoleculeGuessBonds(Molecule *mp, Double limit, Int *outNbonds, Int **outBonds)
{
	int i, j, n1, n2;
	Int nbonds, *bonds;
	Atom *ap, *bp;
	Vector r1, r2, dr;
	Double a1, a2, alim;
	Int newbond[2];
	ElementPar *p = gElementParameters;
	nbonds = 0;
	bonds = NULL;
	for (i = 0; i < mp->natoms; i++) {
		ap = ATOM_AT_INDEX(mp->atoms, i);
		n1 = ap->atomicNumber;
		if (n1 >= 0 && n1 < gCountElementParameters)
			a1 = p[n1].radius;
		else a1 = p[6].radius;
		r1 = ap->r;
	/*	if (mp->is_xtal_coord)
			TransformVec(&r1, mp->cell->tr, &r1); */
		for (j = 0; j < i; j++) {
			bp = ATOM_AT_INDEX(mp->atoms, j);
			n2 = bp->atomicNumber;
			if (n2 >= 0 && n2 < gCountElementParameters)
				a2 = p[n2].radius;
			else a2 = p[6].radius;
			r2 = bp->r;
		/*	if (mp->is_xtal_coord)
				TransformVec(&r2, mp->cell->tr, &r2); */
			VecSub(dr, r1, r2);
			alim = limit * (a1 + a2);
			if (VecLength2(dr) < alim * alim) {
				newbond[0] = i;
				newbond[1] = j;
			/*	MoleculeAddBonds(mp, 1, newbonds); */
				AssignArray(&bonds, &nbonds, sizeof(Int) * 2, nbonds, newbond);
			}
		}
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
	Int ibuf[6];
	
	__MoleculeLock(mp);

	/*  Find bonds   */
	if (mp->nbonds == 0) {
		for (i = 0; i < mp->natoms; i++) {
			ap = ATOM_AT_INDEX(mp->atoms, i);
			for (j = 0; j < ap->nconnects; j++) {
				k = ap->connects[j];
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
			for (j = 0; j < ap->nconnects; j++) {
				for (k = j + 1; k < ap->nconnects; k++) {
					ibuf[0] = ap->connects[j];
					ibuf[1] = i;
					ibuf[2] = ap->connects[k];
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
			for (j = 0; j < ap->nconnects; j++) {
				int jj, kk, mm, m;
				Atom *apjj;
				jj = ap->connects[j];
				if (i >= jj)
					continue;
				apjj = ATOM_AT_INDEX(mp->atoms, jj);
				for (k = 0; k < ap->nconnects; k++) {
					if (k == j)
						continue;
					kk = ap->connects[k];
					for (m = 0; m < apjj->nconnects; m++) {
						mm = apjj->connects[m];
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
			for (i1 = 0; i1 < ap->nconnects; i1++) {
				n1 = ap->connects[i1];
				for (i2 = i1 + 1; i2 < ap->nconnects; i2++) {
					n2 = ap->connects[i2];
					for (i4 = i2 + 1; i4 < ap->nconnects; i4++) {
						n4 = ap->connects[i4];
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

int
MoleculeAreAtomsConnected(Molecule *mp, int n1, int n2)
{
	Atom *ap;
	int i;
	if (mp == NULL || n1 < 0 || n1 >= mp->natoms || n2 < 0 || n2 >= mp->natoms)
		return 0;
	ap = ATOM_AT_INDEX(mp->atoms, n1);
	for (i = 0; i < ap->nconnects; i++)
		if (ap->connects[i] == n2)
			return 1;
	return 0;
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
	int ii, jj, ni, nj;
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
	for (ni = 0; ni < api->nconnects; ni++) {
		ii = api->connects[ni];
		if (ig != NULL && IntGroupLookupPoint(ig, ii) < 0)
			continue;
		if (sExistInEqList(ii, 0, list1))
			continue;
		list2 = NULL;
		for (nj = 0; nj < apj->nconnects; nj++) {
			jj = apj->connects[nj];
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
	int i, j, k, ii, jj, kk;
	if (mol == NULL || mol->natoms == 0)
		return NULL;
	db = (Int **)calloc(sizeof(Int *), mol->natoms);

	/*  Find the equivalent univalent atoms  */
	for (i = 0, api = mol->atoms; i < mol->natoms; i++, api = ATOM_NEXT(api)) {
		if (api->nconnects < 2)
			continue;
		for (j = 0; j < api->nconnects; j++) {
			Int ibuf[ATOMS_MAX_CONNECTS], n;
			memset(ibuf, 0, sizeof(ibuf));
			n = 0;
			jj = api->connects[j];
			if (ig != NULL && IntGroupLookupPoint(ig, jj) < 0)
				continue;
			ibuf[n++] = jj;
			apj = ATOM_AT_INDEX(mol->atoms, jj);
			if (apj->nconnects != 1 || db[jj] != NULL)
				continue;
			for (k = j + 1; k < api->nconnects; k++) {
				kk = api->connects[k];
				if (ig != NULL && IntGroupLookupPoint(ig, kk) < 0)
					continue;
				apk = ATOM_AT_INDEX(mol->atoms, kk);
				if (apk->nconnects != 1 || db[kk] != NULL)
					continue;
				if (apj->atomicNumber == apk->atomicNumber)
					ibuf[n++] = kk;
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
MoleculeTransformBySymop(Molecule *mp, const Vector *vpin, Vector *vpout, Symop symop)
{
	if (mp == NULL)
		return 1;
	if (symop.sym >= mp->nsyms)
		return 2;
	if (mp->cell != NULL /* && !mp->is_xtal_coord */) {
		TransformVec(vpout, mp->cell->rtr, vpin);
		TransformVec(vpout, mp->syms[symop.sym], vpout);
		vpout->x += symop.dx;
		vpout->y += symop.dy;
		vpout->z += symop.dz;
		TransformVec(vpout, mp->cell->tr, vpout);
	} else {
		TransformVec(vpout, mp->syms[symop.sym], vpin);
		vpout->x += symop.dx;
		vpout->y += symop.dy;
		vpout->z += symop.dz;
	}
	return 0;
}

/*  TODO: Decide whether the symmetry operations are expressed in internal coordinates or cartesian coordinates. */
int
MoleculeAddExpandedAtoms(Molecule *mp, Symop symop, IntGroup *group)
{
	int i, n, n0, n1, count, *table;
	Atom *ap;
	IntGroupIterator iter;
	int debug = 0;
	if (mp == NULL || mp->natoms == 0 || group == NULL || (count = IntGroupGetCount(group)) == 0)
		return -1;
	if (symop.sym >= mp->nsyms)
		return -2;
	fprintf(stderr, "symop = {%d %d %d %d}\n", symop.sym, symop.dx, symop.dy, symop.dz);

	/*  Create atoms, with avoiding duplicates  */
	n0 = n1 = mp->natoms;
	table = (int *)malloc(sizeof(int) * n0);
	if (table == NULL)
		return -3;
	for (i = 0; i < n0; i++)
		table[i] = -1;
	IntGroupIteratorInit(group, &iter);
	__MoleculeLock(mp);
	for (i = 0; i < count; i++) {
		int n2;
		Atom *ap2;
		Vector nr, dr;
		n = IntGroupIteratorNext(&iter);
		ap = ATOM_AT_INDEX(mp->atoms, n);
		if (SYMOP_ALIVE(ap->symop)) {
			/*  Skip if the atom is expanded  */
			continue;
		}
		/*  Is this expansion already present?  */
		for (n2 = 0, ap2 = mp->atoms; n2 < n0; n2++, ap2 = ATOM_NEXT(ap2)) {
			if (ap2->symbase == n && SYMOP_EQUAL(symop, ap2->symop))
				break;
		}
		if (n2 < n0) {
			/*  If yes, then skip it  */
			continue;
		}
		/*  Is the expanded position coincides with itself?  */
		MoleculeTransformBySymop(mp, &(ap->r), &nr, symop);
		VecSub(dr, ap->r, nr);
		if (VecLength2(dr) < 1e-6) {
			/*  If yes, then this atom is included but no new atom is created  */
			table[n] = n;
		} else {
			/*  Create a new atom  */
			Atom newAtom;
			AtomDuplicate(&newAtom, ap);
			MoleculeCreateAnAtom(mp, &newAtom, -1);
			AtomClean(&newAtom);
			ap2 = ATOM_AT_INDEX(mp->atoms, mp->natoms - 1);
			ap2->r = nr;
			ap2->symbase = n;
			ap2->symop = symop;
			ap2->symop.alive = (symop.dx != 0 || symop.dy != 0 || symop.dz != 0 || symop.sym != 0);
			table[n] = n1;  /*  The index of the new atom  */
			n1++;
		}
	}
	IntGroupIteratorRelease(&iter);

	/*  Create bonds  */
	for (i = 0; i < n0; i++) {
		int b[2];
		b[0] = table[i];
		if (b[0] < 0 || b[0] == i)
			continue;
		ap = ATOM_AT_INDEX(mp->atoms, i);
		for (n = 0; n < ap->nconnects; n++) {
			b[1] = table[ap->connects[n]];
			if (b[1] < 0)
				continue;
			if (b[1] > n0 && b[0] > b[1])
				continue;
			MoleculeAddBonds(mp, 1, b);
		}
	}
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);
	free(table);
	return n1 - n0;  /*  The number of added atoms  */
}

/*  Recalculate the coordinates of symmetry expanded atoms.
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
	if (mp == NULL || mp->natoms == 0 || mp->nsyms == 0)
		return 0;
	if (groupout != NULL && vpout != NULL) {
		*groupout = IntGroupNew();
		if (*groupout == NULL)
			return -1;
		*vpout = (Vector *)malloc(sizeof(Vector) * mp->natoms);
		if (*vpout == NULL) {
			IntGroupRelease(*groupout);
			return -1;
		}
	} else groupout = NULL; /* To simplify test for validity of groupout/vpout */
	
	__MoleculeLock(mp);
	count = 0;
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		Vector nr, dr;
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
			(*vpout)[count] = ap->r;
			IntGroupAdd(*groupout, i, 1);
		}
		ap->r = nr;
		count++;
	}
	mp->needsMDCopyCoordinates = 1;
	__MoleculeUnlock(mp);
	if (groupout != NULL) {
		if (count == 0) {
			free(*vpout);
			*vpout = NULL;
			IntGroupRelease(*groupout);
			*groupout = NULL;
		} else {
			*vpout = (Vector *)realloc(*vpout, sizeof(Vector) * count);
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
	if (objs == NULL || where == NULL)
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
	ap1->nconnects = 0;
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
			for (j = 0; j < api->nconnects; j++) {
				if (api->connects[j] >= pos)
					api->connects[j]++;
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

/*  Merge two molecules. We use this procedure for all add-atom operations.  */
/*  resSeqOffset is an offset to add to the (non-zero) residue numbers in src. */
int
MoleculeMerge(Molecule *dst, Molecule *src, IntGroup *where, int resSeqOffset)
{
	int nsrc, ndst;
	int i, j, n1, n2, n3, n4;
	Int *new2old, *old2new;
	Atom *ap;
	if (dst == NULL || src == NULL || src->natoms == 0 || (where != NULL && IntGroupGetIntervalCount(where) == 0))
		return 0;  /*  Do nothing  */

	if (dst->noModifyTopology)
		return 1;  /*  Prohibited operation  */

	if (where != NULL && IntGroupGetCount(where) != src->natoms)
		return 1;  /*  Bad parameter  */

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
			if (AtomDuplicate(ATOM_AT_INDEX(dst->atoms, ndst + i), ATOM_AT_INDEX(src->atoms, i)) == NULL)
				goto panic;
		}
	//	memmove(ATOM_AT_INDEX(dst->atoms, ndst), src->atoms, gSizeOfAtomRecord * nsrc);
	} else {
		/*  Duplicate to a temporary storage and then insert  */
		Atom *tempatoms = (Atom *)malloc(gSizeOfAtomRecord * nsrc);
		if (tempatoms == NULL)
			goto panic;
		for (i = 0; i < nsrc; i++) {
			if (AtomDuplicate(ATOM_AT_INDEX(tempatoms, i), ATOM_AT_INDEX(src->atoms, i)) == NULL)
				goto panic;
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
		for (j = 0; j < ap->nconnects; j++)
			ap->connects[j] = old2new[ap->connects[j] + n1];
		if (SYMOP_ALIVE(ap->symop))
			ap->symbase = old2new[ap->symbase + n1];
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
		/*  Keep the old number of entries in dst, because it is updated by AssignArray()  */
		n1 = *nitems;
		/*  Also keep the old number of entries in src, in case src and dst point the same molecule  */
		n2 = *nitems_src;
		/*  Expand the array  */
		if (AssignArray(items, nitems, sizeof(Int) * nsize, *nitems + *nitems_src - 1, NULL) == NULL)
			goto panic;
		/*  Copy the items  */
		memmove(*items + n1 * nsize, *items_src, sizeof(Int) * nsize * n2);
		/*  Renumber  */
		for (j = 0; j < n1 * nsize; j++)
			(*items)[j] = old2new[(*items)[j]];
		for (j = n1 * nsize; j < (n1 + n2) * nsize; j++)
			(*items)[j] = old2new[(*items)[j] + ndst];
	}
	
	/*  Copy the residues if necessary  */
	/*  src[1..src->nresidues-1] should become dst[1+resSeqOffset..src->nresidues+resSeqOffset-1];
	    However, 1+resSeqOffset should not overwrite the existing residue in dst;
		i.e. if 1+resSeqOffset is less than dst->nresidues, copy should start from src[dst->nresidues-resSeqOffset] instead of src[1].  */
	n1 = dst->nresidues;
	if (1 + resSeqOffset < n1) {
		n2 = n1;
	} else n2 = 1 + resSeqOffset; /* n2 is the start index of residues from src[] */
	if (src->nresidues > 1 && n1 < src->nresidues + resSeqOffset) {
		if (AssignArray(&dst->residues, &dst->nresidues, sizeof(dst->residues[0]), src->nresidues + resSeqOffset - 1, NULL) == NULL)
			goto panic;
		memmove(dst->residues + n2, src->residues + n2 - resSeqOffset, sizeof(dst->residues[0]) * (src->nresidues - (n2 - resSeqOffset)));
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

static int
sMoleculeUnmergeSub(Molecule *src, Molecule **dstp, IntGroup *where, int resSeqOffset, int moveFlag)
{
	int nsrc, ndst, nsrcnew;
	int i, j, n1, n2, n3, n4;
	Int *new2old, *old2new;
	IntGroup *move_g, *del_g;
	Molecule *dst;
	Atom *ap, *dst_ap;
	
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

	/*  Make sure that no bond between the two fragments exists  */
/*	for (i = 0; i < src->nbonds; i++) {
		n1 = old2new[src->bonds[i * 2]];
		n2 = old2new[src->bonds[i * 2 + 1]];
		if ((n1 < nsrcnew && n2 >= nsrcnew) || (n1 >= nsrcnew && n2 < nsrcnew)) {
			free(new2old);
			return 2;
		}
	} */
	
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
#if 0
		if (sCopyElementsFromArrayAtPositions(src->atoms, src->natoms, dst_ap, gSizeOfAtomRecord, where) != 0)
			goto panic;
		if (dst != NULL) {
			/*  The atom record must be deep-copied correctly  */
			for (i = 0; i < ndst; i++) {
				if (AtomDuplicate(ATOM_AT_INDEX(dst_ap, i), ATOM_AT_INDEX(src->atoms, i)) == NULL)
					goto panic;
			}
		}
#endif
	}
	
	if (dst == NULL) {
		/*  The dummy destination array is no longer needed  */
		free(dst_ap);
		dst_ap = NULL;
	}
	
	/*  Renumber the atom indices in connect[]  */
	if (moveFlag) {
		for (i = 0, ap = src->atoms; i < src->natoms; i++, ap = ATOM_NEXT(ap)) {
			for (j = n1 = 0; j < ap->nconnects; j++) {
				n2 = old2new[ap->connects[j]];
				if (n2 < nsrcnew)
					ap->connects[n1++] = n2;
			}
			ap->nconnects = n1;
		}
	}
	
	/*  Renumber the atom indices in connect[] and the residue indices  */
	if (dst != NULL) {
		for (i = 0, ap = dst->atoms; i < dst->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ap->resSeq != 0 && ap->resSeq - resSeqOffset >= 0)
				ap->resSeq -= resSeqOffset;
			else ap->resSeq = 0;
			for (j = n1 = 0; j < ap->nconnects; j++) {
				n2 = old2new[ap->connects[j]] - nsrcnew;
				if (n2 >= 0)
					ap->connects[n1++] = n2;
			}
			ap->nconnects = n1;
		}
	}

	/*  Separate the bonds, angles, dihedrals, impropers  */
	move_g = IntGroupNew();
	del_g = IntGroupNew();
	if (move_g == NULL || del_g == NULL)
		goto panic;
	for (i = 0; i < 4; i++) {
		Int *nitems, *nitems_dst;
		Int **items, **items_dst;
		Int nsize;  /*  Number of Ints in one element  */
		unsigned char *counts;
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
		/*	if (n1 >= nsrcnew) {
				n1 -= nsrcnew;
				if (j % nsize == 0) {
					if (IntGroupAdd(sep, j / nsize, 1) != 0)
						goto panic;
					n2++;
				}
			}
			(*items)[j] = n1; */
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
			}
			/*  Remove from src  */
			if (moveFlag) {
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
		IntGroupClear(del_g);
	}
	IntGroupRelease(move_g);
	IntGroupRelease(del_g);
	
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

#if 0
	/*  Copy the parameters  */
	if (dst != NULL && src->par != NULL) {
		UnionPar *up;
		for (i = 0; i < nsrc; i++) {
			old2new[i] -= nsrcnew;  /*  new indices for atoms in dst; otherwise negative numbers  */
		}
		move_g = IntGroupNew();
		for (n1 = kFirstParType; n1 <= kLastParType; n1++) {
			n2 = ParameterGetCountForType(src->par, n1);
			if (n2 == 0)
				continue;
			/*  Find parameters to be copied to dst  */
			for (i = 0; i < n2; i++) {
				up = ParameterGetUnionParFromTypeAndIndex(src->par, n1, i);
				for (j = 0, ap = dst->atoms; j < dst->natoms; j++, ap = ATOM_NEXT(ap)) {
					if (ParameterDoesContainAtom(n1, up, new2old[j + nsrcnew], kParameterLookupNoWildcard) || ParameterDoesContainAtom(n1, up, ap->type, kParameterLookupNoWildcard)) {
						IntGroupAdd(move_g, i, 1);
						break;
					}
				}
			}
			n2 = IntGroupGetCount(move_g);
			if (n2 == 0)
				continue;
			up = (UnionPar *)calloc(sizeof(UnionPar), n2);
			if (up == NULL)
				goto panic;
			/*  Renumber indices if necessary  */
			for (i = 0; i < n2; i++)
				ParameterRenumberAtoms(n1, up + i, nsrc, old2new);
			IntGroupClear(move_g);
		}
		for (i = 0; i < nsrc; i++) {
			old2new[i] += nsrcnew;  /*  Restore indices  */
		}
		IntGroupRelease(move_g);
	}
#endif
	
	/*  Clean up  */
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
MoleculeUnmerge(Molecule *src, Molecule **dstp, IntGroup *where, int resSeqOffset)
{
	return sMoleculeUnmergeSub(src, dstp, where, resSeqOffset, 1);
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
	retval = sMoleculeUnmergeSub(src, dstp, where, 0, 0);
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
				MoleculeAddBonds(*dstp, 1, nn);
			}
			IntGroupIteratorRelease(&iter);
			IntGroupRelease(ig);
		}
	}
	
	return 0;
}

int
MoleculeAddBonds(Molecule *mp, Int nbonds, const Int *bonds)
{
	int i, j, n1, n2, n;
	Atom *ap;
	Int *bonds_tmp;

	if (mp == NULL || bonds == NULL || nbonds <= 0)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */

	/*  Check the bonds  */
	bonds_tmp = (Int *)malloc(sizeof(Int) * nbonds * 2);
	if (bonds_tmp == NULL)
		return -4;  /*  Out of memory  */
	n = 0;
	for (i = 0; i < nbonds; i++) {
		n1 = bonds[i * 2];
		n2 = bonds[i * 2 + 1];
		if (n1 == n2 || n1 < 0 || n1 >= mp->natoms || n2 < 0 || n2 >= mp->natoms)
			return -1;  /*  Bad bond specification  */
		ap = ATOM_AT_INDEX(mp->atoms, n1);
		if (ap->nconnects >= ATOMS_MAX_CONNECTS - 1 || ATOM_AT_INDEX(mp->atoms, n2)->nconnects >= ATOMS_MAX_CONNECTS - 1)
			return -2;  /*  Too many bonds  */
		/*  Check duplicates  */
		for (j = 0; j < ap->nconnects; j++) {
			if (ap->connects[j] == n2)
				break;
		}
		if (j == ap->nconnects) {
			bonds_tmp[n * 2] = n1;
			bonds_tmp[n * 2 + 1] = n2;
			n++;
		}
	}
	if (n == 0) {
		/*  No bonds to add  */
		free(bonds_tmp);
		return 0;
	}
	
	__MoleculeLock(mp);

	/*  Add connects[]  */
	for (i = 0; i < n; i++) {
		n1 = bonds_tmp[i * 2];
		n2 = bonds_tmp[i * 2 + 1];
		ap = ATOM_AT_INDEX(mp->atoms, n1);
		ap->connects[ap->nconnects++] = n2;
		ap = ATOM_AT_INDEX(mp->atoms, n2);
		ap->connects[ap->nconnects++] = n1;
	}
	
	/*  Expand the array and insert  */
	n1 = mp->nbonds;
/*	if (AssignArray(&(mp->bonds), &(mp->nbonds), sizeof(Int) * 2, mp->nbonds + nb - 1, NULL) == NULL
	|| sInsertElementsToArrayAtPositions(mp->bonds, n1, bonds, nb, sizeof(Int) * 2, where) != 0) */
	if (AssignArray(&(mp->bonds), &(mp->nbonds), sizeof(Int) * 2, mp->nbonds + n - 1, NULL) == NULL)
		goto panic;
	memmove(mp->bonds + n1 * 2, bonds_tmp, sizeof(Int) * 2 * n);

	/*  Add angles, dihedrals, impropers  */
	{
		Int nangles, ndihedrals, nimpropers;
		Int *angles, *dihedrals, *impropers;
		Int k, n3, n4;
		Int *ip;
		Int temp[4];
		Atom *ap1, *ap2;

		angles = dihedrals = impropers = NULL;
		nangles = ndihedrals = nimpropers = 0;

		for (i = 0; i < n; i++) {
			n1 = bonds_tmp[i * 2];
			n2 = bonds_tmp[i * 2 + 1];
			ap1 = ATOM_AT_INDEX(mp->atoms, n1);
			ap2 = ATOM_AT_INDEX(mp->atoms, n2);
			/*  Angles X-n1-n2  */
			for (j = 0; j < ap1->nconnects; j++) {
				n3 = ap1->connects[j];
				if (n3 == n2)
					continue;
				temp[0] = n3;
				temp[1] = n1;
				temp[2] = n2;
				for (k = 0; k < nangles; k++) {
					ip = angles + k * 3;
					if (ip[1] == n1 && ((ip[0] == n3 && ip[2] == n2) || (ip[0] == n2 && ip[2] == n3)))
						break;
				}
				if (k == nangles) {
					if (AssignArray(&angles, &nangles, sizeof(Int) * 3, nangles, temp) == NULL)
						goto panic;
				}
				/*  Dihedrals X-n1-n2-X  */
				for (k = 0; k < ap2->nconnects; k++) {
					n4 = ap2->connects[k];
					if (n4 == n1 || n4 == n3)
						continue;
					temp[3] = n4;
					if (AssignArray(&dihedrals, &ndihedrals, sizeof(Int) * 4, ndihedrals, temp) == NULL)
						goto panic;
				}
				/*  Impropers X-n2-n1-X  */
			/*	temp[1] = n2;
				temp[2] = n1;
				for (k = 0; k < ap1->nconnects; k++) {
					n4 = ap1->connects[k];
					if (n4 == n2 || n4 <= n3)
						continue;
					temp[3] = n4;
					if (AssignArray(&impropers, &nimpropers, sizeof(Int) * 4, nimpropers, temp) == NULL)
						goto panic;
				} */
			}
			/*  Angles X-n2-n1  */
			for (j = 0; j < ap2->nconnects; j++) {
				n3 = ap2->connects[j];
				if (n3 == n1)
					continue;
				temp[0] = n1;
				temp[1] = n2;
				temp[2] = n3;
				for (k = 0; k < nangles; k++) {
					ip = angles + k * 3;
					if (ip[1] == n2 && ((ip[0] == n3 && ip[2] == n1) || (ip[0] == n1 && ip[2] == n3)))
						break;
				}
				if (k == nangles) {
					if (AssignArray(&angles, &nangles, sizeof(Int) * 3, nangles, temp) == NULL)
						goto panic;
				}
			}
		}
		temp[0] = kInvalidIndex;
		if (AssignArray(&angles, &nangles, sizeof(Int) * 3, nangles, temp) == NULL)
			goto panic;
		if (AssignArray(&dihedrals, &ndihedrals, sizeof(Int) * 4, ndihedrals, temp) == NULL)
			goto panic;
		if (AssignArray(&impropers, &nimpropers, sizeof(Int) * 4, nimpropers, temp) == NULL)
			goto panic;
		MoleculeAddAngles(mp, angles, NULL);
		MoleculeAddDihedrals(mp, dihedrals, NULL);
		MoleculeAddImpropers(mp, impropers, NULL);
		if (angles != NULL)
			free(angles);
		if (dihedrals != NULL)
			free(dihedrals);
		if (impropers != NULL)
			free(impropers);
	}
	
	MoleculeIncrementModifyCount(mp);
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);

	free(bonds_tmp);
	return n;	

  panic:
	__MoleculeUnlock(mp);
	Panic("Low memory while adding bonds");
	return -1;  /*  Not reached  */
}

int
MoleculeDeleteBonds(Molecule *mp, Int nbonds, const Int *bonds)
{
	int i, j, n1, n2;
	Atom *ap;

	if (mp == NULL || nbonds <= 0)
		return 0;
	if (mp->noModifyTopology)
		return -4;  /*  Prohibited operation  */

	__MoleculeLock(mp);

	/*  Update connects[]  */
	for (i = 0; i < nbonds; i++) {
		n1 = bonds[i * 2];
		n2 = bonds[i * 2 + 1];
		ap = ATOM_AT_INDEX(mp->atoms, n1);
		for (j = 0; j < ap->nconnects; j++) {
			if (ap->connects[j] == n2) {
				memmove(&ap->connects[j], &ap->connects[j + 1], sizeof(Int) * (ap->nconnects - j - 1));
				ap->nconnects--;
				break;
			}
		}
		ap = ATOM_AT_INDEX(mp->atoms, n2);
		for (j = 0; j < ap->nconnects; j++) {
			if (ap->connects[j] == n1) {
				memmove(&ap->connects[j], &ap->connects[j + 1], sizeof(Int) * (ap->nconnects - j - 1));
				ap->nconnects--;
				break;
			}
		}
	}

	/*  Remove bonds, angles, dihedrals, impropers  */	
	{
		IntGroup *bg, *ag, *dg, *ig;
		Int *ip;

		bg = IntGroupNew();
		ag = IntGroupNew();
		dg = IntGroupNew();
		ig = IntGroupNew();
		if (bg == NULL || ag == NULL || dg == NULL || ig == NULL)
			goto panic;
		for (i = 0; i < nbonds; i++) {
			n1 = bonds[i * 2];
			n2 = bonds[i * 2 + 1];
			for (j = 0; j < mp->nbonds; j++) {
				ip = mp->bonds + j * 2;
				if ((ip[0] == n1 && ip[1] == n2)
				 || (ip[1] == n1 && ip[0] == n2)) {
					if (IntGroupAdd(bg, j, 1) != 0)
						goto panic;
				}
			}
			for (j = 0; j < mp->nangles; j++) {
				ip = mp->angles + j * 3;
				if ((ip[0] == n1 && ip[1] == n2)
				 || (ip[1] == n1 && ip[0] == n2)
				 || (ip[1] == n1 && ip[2] == n2)
				 || (ip[2] == n1 && ip[1] == n2)) {
					if (IntGroupAdd(ag, j, 1) != 0)
						goto panic;
				}
			}
			for (j = 0; j < mp->ndihedrals; j++) {
				ip = mp->dihedrals + j * 4;
				if ((ip[1] == n1 && ip[2] == n2)
				 || (ip[2] == n1 && ip[1] == n2)) {
					if (IntGroupAdd(dg, j, 1) != 0)
						goto panic;
				}
			}
			for (j = 0; j < mp->nimpropers; j++) {
				ip = mp->impropers + j * 4;
				if ((ip[0] == n1 && ip[2] == n2)
				 || (ip[1] == n1 && ip[2] == n2)
				 || (ip[3] == n1 && ip[2] == n2)
				 || (ip[0] == n2 && ip[2] == n1)
				 || (ip[1] == n2 && ip[2] == n1)
				 || (ip[3] == n2 && ip[2] == n1)) {
					if (IntGroupAdd(ig, j, 1) != 0)
						goto panic;
				}
			}
		}
		if (sRemoveElementsFromArrayAtPositions(mp->bonds, mp->nbonds, NULL, sizeof(Int) * 2, bg) != 0)
			goto panic;
		mp->nbonds -= IntGroupGetCount(bg);
		
		if (IntGroupGetCount(ag) > 0)
			MoleculeDeleteAngles(mp, NULL, ag);
		if (IntGroupGetCount(dg) > 0)
			MoleculeDeleteDihedrals(mp, NULL, dg);
		if (IntGroupGetCount(ig) > 0)
			MoleculeDeleteImpropers(mp, NULL, ig);
		IntGroupRelease(bg);
		IntGroupRelease(ag);
		IntGroupRelease(dg);
		IntGroupRelease(ig);
	}

	MoleculeIncrementModifyCount(mp);
	mp->needsMDRebuild = 1;
	__MoleculeUnlock(mp);

	return nbonds;

  panic:
	__MoleculeUnlock(mp);
	Panic("Low memory while removing bonds");
	return -1;  /*  Not reached  */
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
		Panic("Low memory while adding angles");
	}
	mp->nangles -= (nc = IntGroupGetCount(where));
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
		Panic("Low memory while adding dihedrals");
	}
	mp->ndihedrals -= (nc = IntGroupGetCount(where));
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
		Panic("Low memory while adding impropers");
	}
	mp->nimpropers -= (nc = IntGroupGetCount(where));
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
		nap->nconnects = 0;
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
/*	ig = IntGroupNewWithPoints(bondIndex, 1, -1);
	if (ig == NULL)
		goto panic;
	MoleculeDeleteBonds(mp, NULL, ig);
	IntGroupRelease(ig); */
	MoleculeDeleteBonds(mp, 1, roots);
	newBonds[0] = roots[0];
	newBonds[1] = dummyIndices[0];
	newBonds[2] = roots[1];
	newBonds[3] = dummyIndices[1];
	newBonds[4] = kInvalidIndex;
	
	i = (MoleculeAddBonds(mp, 2, newBonds) < 0 ? -1 : 0);
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
		Int *cp = ap->connects;
		for (j = 0; j < ap->nconnects; j++) {
			Int j0 = cp[j];
			for (k = j + 1; k < ap->nconnects; k++) {
				Int k0 = cp[k];
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
	Int n1, n2, n3, n4, *ip;
	Atom *ap2, *ap3;
	Int ndihedrals;
	Int *dihedrals;
	
	if (mol == NULL || mol->natoms == 0 || mol->atoms == NULL)
		return 0;  /*  molecule is empty  */
	ndihedrals = 0;
	dihedrals = NULL;
	for (n2 = 0, ap2 = mol->atoms; n2 < mol->natoms; n2++, ap2 = ATOM_NEXT(ap2)) {
		Int i1, i3, i4, *ip;
		for (i3 = 0; i3 < ap2->nconnects; i3++) {
			n3 = ap2->connects[i3];
			if (n2 > n3)
				continue;
			ap3 = ATOM_AT_INDEX(mol->atoms, n3);
			for (i1 = 0; i1 < ap2->nconnects; i1++) {
				n1 = ap2->connects[i1];
				if (n1 == n3)
					continue;
				for (i4 = 0; i4 < ap3->nconnects; i4++) {
					n4 = ap3->connects[i4];
					if (n2 == n4 || n1 == n4)
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
	Int n1, n2, n3, n4, t1, t2, t3, t4, *ip;
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
		for (i1 = 0; i1 < ap3->nconnects; i1++) {
			n1 = ap3->connects[i1];
			t1 = ATOM_AT_INDEX(ap, n1)->type;
			for (i2 = i1 + 1; i2 < ap3->nconnects; i2++) {
				n2 = ap3->connects[i2];
				t2 = ATOM_AT_INDEX(ap, n2)->type;
				for (i4 = i2 + 1; i4 < ap3->nconnects; i4++) {
					n4 = ap3->connects[i4];
					t4 = ATOM_AT_INDEX(ap, n4)->type;
					found = 0;
					if (ParameterLookupImproperPar(par, n1, n2, n3, n4, 1) != NULL)
						found = 1;
					else if (ParameterLookupImproperPar(par, t1, t2, t3, t4, 0) != NULL)
						found = 1;
					else if (ParameterLookupImproperPar(gBuiltinParameters, t1, t2, t3, t4, 0) != NULL)
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
		for (j = 0, ip = apArray[i]->connects; j < apArray[i]->nconnects; j++, ip++)
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
	for (i = 0; i < mp->natoms; i++) {
		Atom *ap = ATOM_AT_INDEX(saveAtoms, i);
		Int *ip;
		for (j = 0, ip = ap->connects; j < ap->nconnects; j++, ip++)
			*ip = old2new[*ip];
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

void
MoleculeTranslate(Molecule *mp, const Vector *vp, IntGroup *group)
{
	int i;
	Atom *ap;
	if (mp == NULL || vp == NULL)
		return;
	__MoleculeLock(mp);
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (group != NULL && IntGroupLookup(group, i, NULL) == 0)
			continue;
		VecInc(ap->r, *vp);
	}
	mp->needsMDCopyCoordinates = 1;
	__MoleculeUnlock(mp);
	sMoleculeNotifyChangeAppearance(mp);
}

void
MoleculeRotate(Molecule *mp, const Vector *axis, Double angle, const Vector *center, IntGroup *group)
{
	int i;
	Double w;
	Transform tr;
	Vector cv;
	Atom *ap;
	if (mp == NULL || axis == NULL)
		return;
	w = VecLength(*axis);
	if (w < 1e-7)
		return;
	__MoleculeLock(mp);
	/*  Construct a rotation transform: p' = c + A * (p - c)  */
	if (center == NULL)
		cv.x = cv.y = cv.z = 0.0;
	else
		cv = *center;
	TransformForRotation(tr, axis, angle, &cv);

	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (group != NULL && IntGroupLookup(group, i, NULL) == 0)
			continue;
		TransformVec(&ap->r, tr, &ap->r);
	}
	mp->needsMDCopyCoordinates = 1;
	__MoleculeUnlock(mp);
	sMoleculeNotifyChangeAppearance(mp);
}

void
MoleculeReaxis(Molecule *mp, const Vector *xaxis, const Vector *yaxis, const Vector *zaxis, IntGroup *group)
{
	int i;
	Atom *ap;
	Vector v;
	if (mp == NULL || xaxis == NULL || yaxis == NULL || zaxis == NULL)
		return;
	__MoleculeLock(mp);
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (group != NULL && IntGroupLookup(group, i, NULL) == 0)
			continue;
		v.x = VecDot(ap->r, *xaxis);
		v.y = VecDot(ap->r, *yaxis);
		v.z = VecDot(ap->r, *zaxis);
		ap->r = v;
	}
	mp->needsMDCopyCoordinates = 1;
	__MoleculeUnlock(mp);
	sMoleculeNotifyChangeAppearance(mp);
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
		cp->tr[n1] = vp1->x;
		cp->tr[n1 + 3] = vp1->y;
		cp->tr[n1 + 6] = vp1->z;
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
			cp->tr[n2] = v1.x;
			cp->tr[n2 + 3] = v1.y;
			cp->tr[n2 + 6] = v1.z;
			cp->tr[n3] = v2.x;
			cp->tr[n3 + 3] = v2.y;
			cp->tr[n3 + 6] = v2.z;
		} else {
			VecCross(v1, *vp1, *vp2);
			if (fabs(VecDot(v1, *vp3)) < 1e-7)
				return -1;  /*  Non-regular transform  */
			cp->tr[n2] = vp2->x;
			cp->tr[n2 + 3] = vp2->y;
			cp->tr[n2 + 6] = vp2->z;
			cp->tr[n3] = vp3->x;
			cp->tr[n3 + 3] = vp3->y;
			cp->tr[n3 + 6] = vp3->z;
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
		cp->cell[0] = sqrt(cp->tr[0] * cp->tr[0] + cp->tr[3] * cp->tr[3] + cp->tr[6] * cp->tr[6]);
		cp->cell[1] = sqrt(cp->tr[1] * cp->tr[1] + cp->tr[4] * cp->tr[4] + cp->tr[7] * cp->tr[7]);
		cp->cell[2] = sqrt(cp->tr[2] * cp->tr[2] + cp->tr[5] * cp->tr[5] + cp->tr[8] * cp->tr[8]);
		cp->cell[3] = acos((cp->tr[1] * cp->tr[2] + cp->tr[4] * cp->tr[5] + cp->tr[7] * cp->tr[8]) / (cp->cell[1] * cp->cell[2])) * kRad2Deg;
		cp->cell[4] = acos((cp->tr[2] * cp->tr[0] + cp->tr[5] * cp->tr[3] + cp->tr[8] * cp->tr[6]) / (cp->cell[2] * cp->cell[0])) * kRad2Deg;
		cp->cell[5] = acos((cp->tr[0] * cp->tr[1] + cp->tr[3] * cp->tr[4] + cp->tr[6] * cp->tr[7]) / (cp->cell[0] * cp->cell[1])) * kRad2Deg;
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
	if (a == 0.0) {
		if (mp->cell != NULL) {
			memmove(&cmat, &(mp->cell->tr), sizeof(Transform));
			free(mp->cell);
		} else {
			cmat[0] = cmat[4] = cmat[8] = 1.0;
		}
		mp->cell = NULL;
	/*	mp->is_xtal_coord = 0; */
	} else {
		cp = mp->cell;
		if (cp == NULL) {
			cp = (XtalCell *)malloc(sizeof(XtalCell));
			if (cp == NULL)
				Panic("Low memory during setting cell parameters");
			mp->cell = cp;
			cmat[0] = cmat[4] = cmat[8] = 1.0;
		} else {
		/*	if (mp->is_xtal_coord)
				memmove(&cmat, &(cp->tr), sizeof(Transform)); */
		}
	/*	mp->is_xtal_coord = 1; */
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
		TransformMul(cmat, cp->rtr, cmat);
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
			MoleculeSetAniso(mp, i, 0, anp->bij[0], anp->bij[1], anp->bij[2], anp->bij[3], anp->bij[4], anp->bij[5]);
		}
	}
	__MoleculeUnlock(mp);
	sMoleculeNotifyChangeAppearance(mp);
}

void
MoleculeSetAniso(Molecule *mp, int n1, int type, Double x11, Double x22, Double x33, Double x12, Double x23, Double x31)
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
		anp = (Aniso *)malloc(sizeof(Aniso));
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
	anp->bij[4] = x23 * dx;
	anp->bij[5] = x31 * dx;
	cp = mp->cell;
	if (cp != NULL && u == 1) {
		anp->bij[0] *= cp->rcell[0] * cp->rcell[0];
		anp->bij[1] *= cp->rcell[1] * cp->rcell[1];
		anp->bij[2] *= cp->rcell[2] * cp->rcell[2];
		anp->bij[3] *= cp->rcell[0] * cp->rcell[1]; /* * cos(cp->rcell[5] * kDeg2Rad); */
		anp->bij[4] *= cp->rcell[1] * cp->rcell[2]; /* * cos(cp->rcell[3] * kDeg2Rad); */
		anp->bij[5] *= cp->rcell[2] * cp->rcell[0]; /* * cos(cp->rcell[4] * kDeg2Rad); */
	}
	
	/*  Calculate the principal axes (in Cartesian coordinates)  */
	/*  The principal axes are the eigenvectors of matrix At(B^-1)A, where
		B is (bij) and A is the reciprocal conversion matrix, i.e. x = Az
		in which x and z are the crystal-space and cartesian coordinates. */
	m1[0] = anp->bij[0] / pi22;
	m1[4] = anp->bij[1] / pi22;
	m1[8] = anp->bij[2] / pi22;
	m1[1] = m1[3] = anp->bij[3] / pi22;
	m1[5] = m1[7] = anp->bij[4] / pi22;
	m1[2] = m1[6] = anp->bij[5] / pi22;
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
		anp->pmat[u] = axis[u].x * val[u];
		anp->pmat[u+3] = axis[u].y * val[u];
		anp->pmat[u+6] = axis[u].z * val[u];
	}
	__MoleculeUnlock(mp);
}

int
MoleculeSetPeriodicBox(Molecule *mp, const Vector *ax, const Vector *ay, const Vector *az, const Vector *ao, const char *periodic)
{
	static Vector zeroVec = {0, 0, 0};
/*	static Transform identityTransform = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0}; */
	XtalCell b;
	int n;
	if (mp == NULL)
		return 0;
	if (ax == NULL) {
		if (mp->cell != NULL)
			free(mp->cell);
		mp->cell = NULL;
	/*	mp->is_xtal_coord = 0; */
		return 0;
	}
	b.axes[0] = (ax != NULL ? *ax : zeroVec);
	b.axes[1] = (ay != NULL ? *ay : zeroVec);
	b.axes[2] = (az != NULL ? *az : zeroVec);
	b.origin = *ao;
	memmove(b.flags, periodic, 3);
	if (MoleculeCalculateCellFromAxes(&b, 1) < 0)
		return -1;
	__MoleculeLock(mp);
	if (mp->cell != NULL)
		free(mp->cell);
	mp->cell = (XtalCell *)malloc(sizeof(XtalCell));
	if (mp->cell != NULL) {
		memmove(mp->cell, &b, sizeof(XtalCell));
		n = 0;
	} else n = -2;  /*  Out of memory  */
	__MoleculeUnlock(mp);
	return n;
}

#pragma mark ====== Fragment manipulation ======

static void
sMoleculeFragmentSub(Molecule *mp, int idx, IntGroup *result, IntGroup *exatoms)
{
	Atom *ap;
	int i;
	if (exatoms != NULL && IntGroupLookup(exatoms, idx, NULL))
		return;
	IntGroupAdd(result, idx, 1);
	ap = ATOM_AT_INDEX(mp->atoms, idx);
	for (i = 0; i < ap->nconnects; i++) {
		int idx2 = ap->connects[i];
		if (IntGroupLookup(result, idx2, NULL))
			continue;
		sMoleculeFragmentSub(mp, idx2, result, exatoms);
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
	int i, i1, i2, j, k, bond_count, nval1, nval2;
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
			for (k = 0; k < ap->nconnects; k++) {
				if (IntGroupLookup(group, ap->connects[k], NULL) == 0) {
					bond_count++;
					nval1 = j;
					nval2 = ap->connects[k];
					if (bond_count > 1)
						return 0;  /*  Too many bonds  */
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
	if (mp->nframes < 0) {
		/*  Recalculate  */
		int i, n;
		Atom *ap;
		for (i = n = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ap->nframes > n)
				n = ap->nframes;
		}
		mp->nframes = n;
	}
	return mp->nframes;
}

int
MoleculeInsertFrames(Molecule *mp, IntGroup *group, const Vector *inFrame, const Vector *inFrameCell)
{
	int i, j, count, n_new, n_old, natoms, exframes, last_inserted, old_count;
	Vector *tempv, *vp;
	Atom *ap;
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
	/*  Expand ap->frames for all atoms  */
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (ap->frames == NULL)
			vp = (Vector *)calloc(sizeof(Vector), n_new);
		else
			vp = (Vector *)realloc(ap->frames, sizeof(Vector) * n_new);
		if (vp == NULL) {
			__MoleculeUnlock(mp);
			return -1;
		}
		for (j = ap->nframes; j < n_new; j++)
			vp[j] = ap->r;
		ap->frames = vp;
	}
	if (mp->cell != NULL && (mp->frame_cells != NULL || inFrameCell != NULL)) {
		if (mp->frame_cells == NULL)
			vp = (Vector *)calloc(sizeof(Vector), n_new * 4);
		else
			vp = (Vector *)realloc(mp->frame_cells, sizeof(Vector) * 4 * n_new);
		if (vp == NULL) {
			__MoleculeUnlock(mp);
			return -1;
		}
		mp->frame_cells = vp;
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
	for (i = 0, ap = mp->atoms; i <= mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		int s, t, ns, ne, mult;
		Vector cr;
		ne = s = t = 0;
		if (i == mp->natoms) {
			if (mp->cell == NULL || mp->frame_cells == NULL)
				break;
			vp = mp->frame_cells;
			mult = 4;
		} else {
			cr = ap->r;
			vp = ap->frames;
			mult = 1;
		}
		for (j = 0; (ns = IntGroupGetStartPoint(group, j)) >= 0; j++) {
			if (ns > ne) {
				memmove(tempv + ne * mult, vp + s * mult, sizeof(Vector) * mult * (ns - ne));
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
				} else {
					if (inFrame != NULL)
						tempv[ns] = inFrame[natoms * t + i];
					else
						tempv[ns] = cr;
				}
				t++;
				ns++;
			}
		}
		if (n_new > ne) {
			memmove(tempv + ne * mult, vp + s * mult, sizeof(Vector) * mult * (n_new - ne));
			s += n_new - ne;
		}
		if (i < mp->natoms)
			ap->nframes = n_new;
		memmove(vp, tempv, sizeof(Vector) * mult * n_new);
	}
	free(tempv);
	mp->nframes = n_new;
	MoleculeSelectFrame(mp, last_inserted, 0);
	MoleculeIncrementModifyCount(mp);
	__MoleculeUnlock(mp);
	return count;
}

int
MoleculeRemoveFrames(Molecule *mp, IntGroup *group, Vector *outFrame, Vector *outFrameCell)
{
	int i, count, n_new, n_old, natoms, nframes, old_count;
	Vector *tempv, *vp;
	Atom *ap;
	if (mp == NULL || (natoms = mp->natoms) == 0 || (count = IntGroupGetCount(group)) <= 0)
		return -1;

	/*  outFrame[] should have enough size for Vector * natoms * group.count  */
	memset(outFrame, 0, sizeof(Vector) * natoms * count);
	if (mp->cell != NULL && mp->frame_cells != NULL)
		memset(outFrameCell, 0, sizeof(Vector) * 4 * count);

	n_old = MoleculeGetNumberOfFrames(mp);
	n_new = n_old - count;
	if (n_new < 0)
		n_new = 0;
	tempv = (Vector *)malloc(sizeof(Vector) * n_old * 4);  /*  "*4" for handling cells  */
	if (tempv == NULL)
		return -1;

	/*  group = [n0..n1-1, n2..n3-1, ...]  */
	/*  s = t = 0, */
	/*  tempv[0..n0-1] -> ap[0..n0-1], s += n0,
	    tempv[n0..n1-1] -> outFrame[0..(n1-n0-1)], t += n1-n0,
		tempv[n1..n2-1] -> ap[s..s+(n2-n1-1)], s += n2-n1,
		tempv[n2..n3-1] -> outFrame[t..t+(n3-n2-1)], t += n3-n2,
		...
		tempv[nl..n_old-1] -> ap[s..s+(n_old-nl-1)], s += n_old-nl
		At last, s will become n_new and t will become count.  */
	__MoleculeLock(mp);
	nframes = 0;
	for (i = 0, ap = mp->atoms; i <= mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		int s, t, j, ns, ne;
		int mult;
		/*  if i == mp->natoms, mp->frame_cells is handled  */
		if (i == mp->natoms) {
			if (mp->cell == NULL || mp->frame_cells == NULL)
				break;
			mult = 4;
			vp = mp->frame_cells;
			old_count = n_old;
		} else {
			mult = 1;
			vp = ap->frames;
			if (vp == NULL) {
				ap->frames = vp = (Vector *)calloc(sizeof(Vector), n_old);
				if (vp == NULL) {
					__MoleculeUnlock(mp);
					return -1;
				}
			}
			old_count = ap->nframes;
		}

		/*  Copy vp to tempv  */
		memset(tempv, 0, sizeof(Vector) * mult * n_old);
		memmove(tempv, vp, sizeof(Vector) * mult * (old_count > n_old ? n_old : old_count));
		ne = ns = s = t = 0;
		for (j = 0; ns < n_old && (ns = IntGroupGetStartPoint(group, j)) >= 0; j++) {
			if (ns > n_old)
				ns = n_old;
			if (ns > ne) {
				memmove(vp + s * mult, tempv + ne * mult, sizeof(Vector) * mult * (ns - ne));
				s += ns - ne;
			}
			ne = IntGroupGetEndPoint(group, j);
			if (ne > n_old)
				ne = n_old;
			while (ns < ne) {
				if (i < mp->natoms)
					outFrame[natoms * t + i] = tempv[ns];
				else if (outFrameCell != NULL) {
					outFrameCell[i * 4] = tempv[ns * 4];
					outFrameCell[i * 4 + 1] = tempv[ns * 4 + 1];
					outFrameCell[i * 4 + 2] = tempv[ns * 4 + 2];
					outFrameCell[i * 4 + 3] = tempv[ns * 4 + 3];
				}
				t++;
				ns++;
			}
		}
		if (n_old > ne) {
			memmove(vp + s * mult, tempv + ne * mult, sizeof(Vector) * mult * (n_old - ne));
			s += n_old - ne;
		}
		if (i < mp->natoms)
			ap->nframes = s;
		if (nframes < s)
			nframes = s;
		if (s == 0) {
			if (i < mp->natoms) {
				free(ap->frames);
				ap->frames = NULL;
			} else {
				free(mp->frame_cells);
				mp->frame_cells = NULL;
			}
		} else {
			if (i < mp->natoms)
				ap->frames = (Vector *)realloc(ap->frames, sizeof(Vector) * s);
			else
				mp->frame_cells = (Vector *)realloc(mp->frame_cells, sizeof(Vector) * 4 * s);
		}
	}
	free(tempv);
	mp->nframes = nframes;
	
	/*  Select the "last" frame; do not "copy back" the coordinates to the frame table  */
	i = (mp->cframe >= nframes ? nframes - 1 : mp->cframe);
	MoleculeSelectFrame(mp, i, 0);
	
	MoleculeIncrementModifyCount(mp);
	__MoleculeUnlock(mp);
	return count;
}

#if 0
int
MoleculeInsertFrame(Molecule *mp, int index, const Vector *inFrame)
{
	int i;
	Atom *ap;
	if (mp == NULL || mp->natoms == 0)
		return -1;
	i = MoleculeGetNumberOfFrames(mp);
	if (index < 0 || index > i)
		index = i;

	/*  Insert a new frame at index for each atom  */
	__MoleculeLock(mp);
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		Vector *vp;
		int nframes;
		if (ap->nframes > index)
			nframes = ap->nframes + 1;
		else
			nframes = index + 1;
		if (ap->frames == NULL)
			vp = (Vector *)malloc(sizeof(Vector) * nframes);
		else
			vp = (Vector *)realloc(ap->frames, sizeof(Vector) * nframes);
		if (vp == NULL) {
			__MoleculeUnlock(mp);
			return -1;
		}
		ap->frames = vp;
		/*  Clear newly allocated memory  */
		memset(ap->frames + ap->nframes, 0, sizeof(Vector) * (nframes - ap->nframes));
		/*  Move backward frames[index..nframes-1] by one  */
		memmove(ap->frames + index + 1, ap->frames + index, sizeof(Vector) * (nframes - index - 1));
		ap->nframes = nframes;
		ap->frames[index] = (inFrame != NULL ? inFrame[i] : ap->r);
	}
	mp->nframes = -1;
	MoleculeIncrementModifyCount(mp);
	__MoleculeUnlock(mp);
	return index;
}

int
MoleculeRemoveFrame(Molecule *mp, int frame, Vector *outFrame)
{
	int i, ok, n;
	Atom *ap;
	static Vector zero = {0, 0, 0};
	if (mp == NULL || mp->natoms == 0 || frame < 0)
		return -1;
	ok = 0;
	n = 0;
	__MoleculeLock(mp);
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (frame < ap->nframes) {
			/*  If outFrame != NULL, then copy the coordinates to outFrame[i]  */
			if (outFrame != NULL)
				outFrame[i] = ap->frames[frame];
			/*  Remove this frame  */
			memmove(ap->frames + frame, ap->frames + frame + 1, sizeof(Vector) * (ap->nframes - frame - 1));
			ap->nframes--;
			if (ap->nframes == 0) {
				free(ap->frames);
				ap->frames = NULL;
			}
			ok = 1;
		} else {
			if (outFrame != NULL)
				outFrame[i] = zero;
		}
		if (n < ap->nframes)
			n = ap->nframes;
	}
	mp->nframes = n;
	if (mp->cframe >= n)
		mp->cframe = (n > 0 ? n - 1 : 0);
	else if (mp->cframe > frame)
		mp->cframe--;
	MoleculeIncrementModifyCount(mp);
	__MoleculeUnlock(mp);
	return frame;
}
#endif

int
MoleculeSelectFrame(Molecule *mp, int frame, int copyback)
{
	int i, cframe, ok;
	Atom *ap;
	if (mp == NULL || mp->natoms == 0)
		return -1;
	cframe = mp->cframe;
	ok = 0;
	__MoleculeLock(mp);
	for (i = 0, ap = mp->atoms; i < mp->natoms; i++, ap = ATOM_NEXT(ap)) {
		if (copyback && cframe >= 0 && cframe < ap->nframes) {
			/*  Write the current coordinate back to the frame array  */
			/*  TODO: This behavior may interfere with undo. How to work around this?  */
			ap->frames[cframe] = ap->r;
		}
		if (frame >= 0 && frame < ap->nframes) {
			/*  Read the coordinate from the frame array  */
			ap->r = ap->frames[frame];
			ok = 1;
		}
	}

	if (mp->cell != NULL && mp->frame_cells != NULL) {
		MoleculeSetPeriodicBox(mp, &mp->frame_cells[frame * 4], &mp->frame_cells[frame * 4 + 1], &mp->frame_cells[frame * 4 + 2], &mp->frame_cells[frame * 4 + 3], mp->cell->flags);
	}
	mp->needsMDCopyCoordinates = 1;
	__MoleculeUnlock(mp);
	if (ok) {
		mp->cframe = frame;
		sMoleculeNotifyChangeAppearance(mp);
		return frame;
	} else return -1;
}

#pragma mark ====== MO calculation ======

/*  Calculate an MO value for a single point.  */
/*  Index is the MO number (1-based)  */
/*  tmp is an array of (natoms * 4) atoms, and used to store dr and |dr|^2 for each atom.  */
static Double
sCalcMOPoint(const BasisSet *bset, Int index, const Vector *vp, Double *tmp)
{
	ShellInfo *sp;
	PrimInfo *pp;
	Double val, tval, *cnp, *tmpp, *mobasep, *mop;
	Int i, j;
	/*  Cache dr and |dr|^2  */
	for (i = 0; i < bset->natoms; i++) {
		Vector r = bset->pos[i];
		tmp[i * 4] = r.x = vp->x - r.x;
		tmp[i * 4 + 1] = r.y = vp->y - r.y;
		tmp[i * 4 + 2] = r.z = vp->z - r.z;
		tmp[i * 4 + 3] = r.x * r.x + r.y * r.y + r.z * r.z;
	}
	/*  Iterate over all shells  */
	val = 0.0;
	mobasep = bset->mo + (index - 1) * bset->nmos;
	for (i = 0, sp = bset->shells; i < bset->nshells; i++, sp++) {
		pp = bset->priminfos + sp->p_idx;
		cnp = bset->cns + sp->cn_idx;
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
		}
	}
	return val;
}

/*  Calculate one MO. The input vectors should be in bohr unit (angstrom * 1.889725989 = kAngstrom2Bohr).  */
/*  mono is the MO number (1-based)  */
int
MoleculeCalcMO(Molecule *mp, Int mono, const Vector *op, const Vector *dxp, const Vector *dyp, const Vector *dzp, Int nx, Int ny, Int nz, int (*callback)(double progress, void *ref), void *ref)
{
	int ix, iy, iz, n, nn;
	Cube *cp;
	Double *tmp;
	if (mp == NULL || mp->bset == NULL)
		return -1;
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
	tmp = (Double *)calloc(sizeof(Double), mp->bset->natoms * 4);
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
				cp->dp[n++] = sCalcMOPoint(mp->bset, mono, &p, tmp);
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

int
MoleculeGetDefaultMOGrid(Molecule *mp, Int npoints, Vector *op, Vector *xp, Vector *yp, Vector *zp, Int *nx, Int *ny, Int *nz)
{
	int i;
	Vector rmin, rmax, *vp;
	Double dr, dx, dy, dz;
	if (mp == NULL || mp->bset == NULL || mp->bset->natoms == 0)
		return -1;
	if (npoints <= 0)
		npoints = 1000000;
	rmin.x = rmin.y = rmin.z = 1e10;
	rmax.x = rmax.y = rmax.z = -1e10;
	for (i = 0, vp = mp->bset->pos; i < mp->bset->natoms; i++, vp++) {
		dr = RadiusForAtomicNumber(ATOM_AT_INDEX(mp->atoms, i)->atomicNumber);
		if (dr == 0.0)
			dr = 1.0;
		dr = dr * kAngstrom2Bohr * 3.0 + 2.0;
		if (rmin.x > vp->x - dr)
			rmin.x = vp->x - dr;
		if (rmin.y > vp->y - dr)
			rmin.y = vp->y - dr;
		if (rmin.z > vp->z - dr)
			rmin.z = vp->z - dr;
		if (rmax.x < vp->x + dr)
			rmax.x = vp->x + dr;
		if (rmax.y < vp->y + dr)
			rmax.y = vp->y + dr;
		if (rmax.z < vp->z + dr)
			rmax.z = vp->z + dr;
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
	
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", -(mp->bset->natoms), cp->origin.x, cp->origin.y, cp->origin.z);
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", cp->nx, cp->dx.x, cp->dx.y, cp->dx.z);
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", cp->ny, cp->dy.x, cp->dy.y, cp->dy.z);
	fprintf(fp, "%5d %11.6f %11.6f %11.6f\n", cp->nz, cp->dz.x, cp->dz.y, cp->dz.z);
	
	/*  Atomic information  */
	for (i = 0; i < mp->bset->natoms; i++) {
		Vector *vp = mp->bset->pos + i;
		Atom *ap = ATOM_AT_INDEX(mp->atoms, i);
		/*  The second number should actually be the effective charge  */
		fprintf(fp, "%5d %11.6f %11.6f %11.6f %11.6f\n", ap->atomicNumber, (double)ap->atomicNumber, vp->x, vp->y, vp->z);
	}
	fprintf(fp, "%5d%5d\n", 1, 1);
	
	/*  3D data  */
	for (i = n = 0; i < cp->nx; i++) {
		for (j = 0; j < cp->ny; j++) {
			for (k = 0; k < cp->nz; k++) {
				fprintf(fp, " %12.5e", cp->dp[n++]);
				if (k == cp->nz - 1 || k % 6 == 5)
					fprintf(fp, "\n");
			}
		}
	}
	fclose(fp);
	return 0;
}
