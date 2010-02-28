/*
 *  Parameter.c
 *
 *  Created by Toshi Nagata on 06/03/11.
 *  Copyright 2006-2008 Toshi Nagata. All rights reserved.
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
#include <ctype.h>
#include <stdarg.h>
#include <math.h>

/*  Global parameter: it is initialized by the first call to ParameterReadFromFile()  */
Parameter *gBuiltinParameters = NULL;

/*  Global parameter  */
ElementPar *gElementParameters = NULL;
Int gCountElementParameters = 0;

static Parameter *sParameterRoot = NULL;
static int sParameterUntitledCount = 0;

/*  Global struct for storing information for global parameter files  */
GlobalParInfoRecord gGlobalParInfo;

#pragma mark ====== Parameter Alloc/Release ======

Parameter *
ParameterNew(void)
{
	char name[40];
	Parameter *par = (Parameter *)calloc(sizeof(Parameter), 1);
	if (par == NULL)
		Panic("Cannot allocate new parameter record");
	snprintf(name, sizeof name, "Untitled %d", sParameterUntitledCount++);
	ObjectInit((Object *)par, (Object **)&sParameterRoot, name);
	return par;
}

Parameter *
ParameterDuplicate(const Parameter *par)
{
	Parameter *npar = ParameterNew();
	if (par->bondPars != NULL) {
		npar->bondPars = (BondPar *)malloc(sizeof(BondPar) * par->nbondPars);
		if (npar->bondPars == NULL)
			goto release;
		memmove(npar->bondPars, par->bondPars, sizeof(BondPar) * par->nbondPars);
		npar->nbondPars = par->nbondPars;
	}
	if (par->anglePars != NULL) {
		npar->anglePars = (AnglePar *)malloc(sizeof(AnglePar) * par->nanglePars);
		if (npar->anglePars == NULL)
			goto release;
		memmove(npar->anglePars, par->anglePars, sizeof(AnglePar) * par->nanglePars);
		npar->nanglePars = par->nanglePars;
	}
	if (par->dihedralPars != NULL) {
		npar->dihedralPars = (TorsionPar *)malloc(sizeof(TorsionPar) * par->ndihedralPars);
		if (npar->dihedralPars == NULL)
			goto release;
		memmove(npar->dihedralPars, par->dihedralPars, sizeof(TorsionPar) * par->ndihedralPars);
		npar->ndihedralPars = par->ndihedralPars;
	}
	if (par->improperPars != NULL) {
		npar->improperPars = (TorsionPar *)malloc(sizeof(TorsionPar) * par->nimproperPars);
		if (npar->improperPars == NULL)
			goto release;
		memmove(npar->improperPars, par->improperPars, sizeof(TorsionPar) * par->nimproperPars);
		npar->nimproperPars = par->nimproperPars;
	}
	if (par->vdwPars != NULL) {
		npar->vdwPars = (VdwPar *)malloc(sizeof(VdwPar) * par->nvdwPars);
		if (npar->vdwPars == NULL)
			goto release;
		memmove(npar->vdwPars, par->vdwPars, sizeof(VdwPar) * par->nvdwPars);
		npar->nvdwPars = par->nvdwPars;
	}
	if (par->vdwpPars != NULL) {
		npar->vdwpPars = (VdwPairPar *)malloc(sizeof(VdwPairPar) * par->nvdwpPars);
		if (npar->vdwpPars == NULL)
			goto release;
		memmove(npar->vdwpPars, par->vdwpPars, sizeof(VdwPairPar) * par->nvdwpPars);
		npar->nvdwpPars = par->nvdwpPars;
	}
	if (par->vdwCutoffPars != NULL) {
		npar->vdwCutoffPars = (VdwCutoffPar *)malloc(sizeof(VdwCutoffPar) * par->nvdwCutoffPars);
		if (npar->vdwCutoffPars == NULL)
			goto release;
		memmove(npar->vdwCutoffPars, par->vdwCutoffPars, sizeof(VdwCutoffPar) * par->nvdwCutoffPars);
		npar->nvdwCutoffPars = par->nvdwCutoffPars;
	}
/*	if (par->atomPars != NULL) {
		npar->atomPars = (ElementPar *)malloc(sizeof(ElementPar) * par->natomPars);
		if (npar->atomPars == NULL)
			goto release;
		memmove(npar->atomPars, par->atomPars, sizeof(ElementPar) * par->natomPars);
		npar->natomPars = par->natomPars;
	} */
	return npar;
release:
	ParameterRelease(npar);
	return NULL;
}

Parameter *
ParameterWithName(const char *name)
{
	return (Parameter *)ObjectWithName(name, (Object *)sParameterRoot);
}

/*  Assign a unique name to this parameter record  */
void
ParameterSetName(Parameter *par, const char *name)
{
	ObjectSetName((Object *)par, name, (Object *)sParameterRoot);
}

const char *
ParameterGetName(Parameter *par)
{
	return ObjectGetName((Object *)par);
}

void
ParameterRetain(Parameter *par)
{
	ObjectIncrRefCount((Object *)par);
}

void
ParameterRelease(Parameter *par)
{
	if (ObjectDecrRefCount((Object *)par) == 0) {
		if (par->bondPars != NULL)
			free(par->bondPars);
		if (par->anglePars != NULL)
			free(par->anglePars);
		if (par->dihedralPars != NULL)
			free(par->dihedralPars);
		if (par->improperPars != NULL)
			free(par->improperPars);
		if (par->vdwPars != NULL)
			free(par->vdwPars);
		if (par->vdwpPars != NULL)
			free(par->vdwpPars);
		if (par->vdwCutoffPars != NULL)
			free(par->vdwCutoffPars);
	/*	if (par->atomPars != NULL)
			free(par->atomPars); */
		ObjectDealloc((Object *)par, (Object **)&sParameterRoot);
	}
}

#pragma mark ====== ParameterRef Definitions ======

ParameterRef *
ParameterRefNew(struct Molecule *mol, int type, int idx)
{
	ParameterRef *pref = (ParameterRef *)calloc(sizeof(ParameterRef), 1);
	if (pref != NULL) {
		pref->mol = mol;
		if (mol != NULL)
			MoleculeRetain(mol);
		pref->parType = type;
		pref->idx = idx;
	}
	return pref;
}

void
ParameterRefRelease(ParameterRef *pref)
{
	if (pref != NULL) {
		if (pref->mol != NULL)
			MoleculeRelease(pref->mol);
		free(pref);
	}
}

UnionPar *
ParameterGetUnionParFromTypeAndIndex(Parameter *par, int type, int index)
{
	if (par == NULL)
		return NULL;
	switch (type) {
		case kBondParType:
			if (index >= 0 && index < par->nbondPars)
				return (UnionPar *)(par->bondPars + index);
			else return NULL;
		case kAngleParType:
			if (index >= 0 && index < par->nanglePars)
				return (UnionPar *)(par->anglePars + index);
			else return NULL;
		case kDihedralParType:
			if (index >= 0 && index < par->ndihedralPars)
				return (UnionPar *)(par->dihedralPars + index);
			else return NULL;
		case kImproperParType:
			if (index >= 0 && index < par->nimproperPars)
				return (UnionPar *)(par->improperPars + index);
			else return NULL;
		case kVdwParType:
			if (index >= 0 && index < par->nvdwPars)
				return (UnionPar *)(par->vdwPars + index);
			else return NULL;
		case kVdwPairParType:
			if (index >= 0 && index < par->nvdwpPars)
				return (UnionPar *)(par->vdwpPars + index);
			else return NULL;
		case kVdwCutoffParType:
			if (index >= 0 && index < par->nvdwCutoffPars)
				return (UnionPar *)(par->vdwCutoffPars + index);
			else return NULL;
	/*	case kElementParType:
			if (index >= 0 && index < par->natomPars)
				return (UnionPar *)(par->atomPars + index);
			else return NULL; */
		default:
			return NULL;
	}
}

Int
ParameterGetCountForType(Parameter *par, int type)
{
	if (par == NULL)
		return 0;
	switch (type) {
		case kBondParType:
			return par->nbondPars;
		case kAngleParType:
			return par->nanglePars;
		case kDihedralParType:
			return par->ndihedralPars;
		case kImproperParType:
			return par->nimproperPars;
		case kVdwParType:
			return par->nvdwPars;
		case kVdwPairParType:
			return par->nvdwpPars;
		case kVdwCutoffParType:
			return par->nvdwCutoffPars;
	/*	case kElementParType:
			return par->natomPars; */
		default:
			return 0;
	}
}

Int
ParameterGetSizeForType(int type)
{
	switch (type) {
		case kBondParType:
			return sizeof(BondPar);
		case kAngleParType:
			return sizeof(AnglePar);
		case kDihedralParType:
			return sizeof(TorsionPar);
		case kImproperParType:
			return sizeof(TorsionPar);
		case kVdwParType:
			return sizeof(VdwPar);
		case kVdwPairParType:
			return sizeof(VdwPairPar);
		case kVdwCutoffParType:
			return sizeof(VdwCutoffPar);
			/*	case kElementParType:
			 return par->natomPars; */
		default:
			return 0;
	}
}

UnionPar *
ParameterRefGetPar(ParameterRef *pref)
{
	Parameter *par;
	if (pref == NULL)
		return NULL;
	if (pref->mol != NULL)
		par = pref->mol->par;
	else par = gBuiltinParameters;
	if (par == NULL)
		return NULL;
	return ParameterGetUnionParFromTypeAndIndex(par, pref->parType, pref->idx);
}

#pragma mark ====== Insert/Delete (for MolAction) ======

int
ParameterInsert(Parameter *par, Int type, const UnionPar *up, struct IntGroup *where)
{
	Int i, n1, n2, size, *ip;
	void *p;
	if (par == NULL)
		return -1;
	switch (type) {
		case kBondParType:
			p = &par->bondPars;
			ip = &par->nbondPars;
			size = sizeof(BondPar);
			break;
		case kAngleParType:
			p = &par->anglePars;
			ip = &par->nanglePars;
			size = sizeof(AnglePar);
			break;
		case kDihedralParType:
			p = &par->dihedralPars;
			ip = &par->ndihedralPars;
			size = sizeof(TorsionPar);
			break;
		case kImproperParType:
			p = &par->improperPars;
			ip = &par->nimproperPars;
			size = sizeof(TorsionPar);
			break;
		case kVdwParType:
			p = &par->vdwPars;
			ip = &par->nvdwPars;
			size = sizeof(VdwPar);
			break;
		case kVdwPairParType:
			p = &par->vdwpPars;
			ip = &par->nvdwpPars;
			size = sizeof(VdwPairPar);
			break;
		case kVdwCutoffParType:
			p = &par->vdwCutoffPars;
			ip = &par->nvdwCutoffPars;
			size = sizeof(VdwCutoffPar);
			break;
	/*	case kElementParType:
			p = &par->atomPars;
			ip = &par->natomPars;
			size = sizeof(ElementPar);
			break; */
		default:
			return -1;
	}
	n2 = 0;
	for (i = 0; (n1 = IntGroupGetNthPoint(where, i)) >= 0; i++) {
		if (InsertArray(p, ip, size, n1, 1, (up ? up + i : NULL)) == NULL)
			return n2;
		n2++;
	}
	return n2;
}

static int
sParameterDeleteOrCopy(Parameter *par, Int type, UnionPar *up, struct IntGroup *where, Int copyflag)
{
	Int i, n1, n2, size, *ip;
	char **p;
	if (par == NULL)
		return -1;
	switch (type) {
		case kBondParType:
			p = (char **)&par->bondPars;
			ip = &par->nbondPars;
			size = sizeof(BondPar);
			break;
		case kAngleParType:
			p = (char **)&par->anglePars;
			ip = &par->nanglePars;
			size = sizeof(AnglePar);
			break;
		case kDihedralParType:
			p = (char **)&par->dihedralPars;
			ip = &par->ndihedralPars;
			size = sizeof(TorsionPar);
			break;
		case kImproperParType:
			p = (char **)&par->improperPars;
			ip = &par->nimproperPars;
			size = sizeof(TorsionPar);
			break;
		case kVdwParType:
			p = (char **)&par->vdwPars;
			ip = &par->nvdwPars;
			size = sizeof(VdwPar);
			break;
		case kVdwPairParType:
			p = (char **)&par->vdwpPars;
			ip = &par->nvdwpPars;
			size = sizeof(VdwPairPar);
			break;
		case kVdwCutoffParType:
			p = (char **)&par->vdwCutoffPars;
			ip = &par->nvdwCutoffPars;
			size = sizeof(VdwCutoffPar);
			break;
	/*	case kElementParType:
			p = (char **)&par->atomPars;
			ip = &par->natomPars;
			size = sizeof(ElementPar);
			break; */
		default:
			return -1;
	}
	n2 = 0;
	for (i = IntGroupGetCount(where) - 1; i >= 0 && (n1 = IntGroupGetNthPoint(where, i)) >= 0; i--) {
		if (n1 >= *ip)
			return -1; /*  Internal error  */
		if (copyflag) {
			if (up != NULL)
				memmove(up + i, *p + size * n1, size);
		} else {
			if (DeleteArray(p, ip, size, n1, 1, (up ? up + i : NULL)) == NULL)
				return n2;
		}
		n2++;
	}
	return n2;
}

int
ParameterDelete(Parameter *par, Int type, UnionPar *up, struct IntGroup *where)
{
	return sParameterDeleteOrCopy(par, type, up, where, 0);
}

int
ParameterCopy(Parameter *par, Int type, UnionPar *up, struct IntGroup *where)
{
	return sParameterDeleteOrCopy(par, type, up, where, 1);
}

/*  Renumber the atom index field according to the conversion table.  If the atom index
    points to a non-existing atom, returns non-zero.  */
int
ParameterRenumberAtoms(Int type, UnionPar *up, Int oldnatoms, const Int *old2new)
{
#define RENUMBER_FIELD(_tp, _field) \
  (((_tp *)up)->_field >= 0 && ((_tp *)up)->_field < oldnatoms && \
	(old2new[((_tp *)up)->_field] >= 0 ? \
      ((((_tp *)up)->_field = old2new[((_tp *)up)->_field]), 0) : \
      ((((_tp *)up)->_field = kAtomTypeWildcard), 1)))
	switch (type) {
		case kBondParType:
			return RENUMBER_FIELD(BondPar, type1) + RENUMBER_FIELD(BondPar, type2);
		case kAngleParType:
			return RENUMBER_FIELD(AnglePar, type1) + RENUMBER_FIELD(AnglePar, type2) + RENUMBER_FIELD(AnglePar, type3);
		case kDihedralParType:
		case kImproperParType:
			return RENUMBER_FIELD(TorsionPar, type1) + RENUMBER_FIELD(TorsionPar, type2) + RENUMBER_FIELD(TorsionPar, type3) + RENUMBER_FIELD(TorsionPar, type4);
		case kVdwParType:
			return RENUMBER_FIELD(VdwPar, type1);
		case kVdwPairParType:
			return RENUMBER_FIELD(VdwPairPar, type1) + RENUMBER_FIELD(VdwPairPar, type2);
		case kVdwCutoffParType:
			if (((VdwCutoffPar *)up)->type == 1)
				return RENUMBER_FIELD(VdwCutoffPar, n1) + RENUMBER_FIELD(VdwCutoffPar, n2) + RENUMBER_FIELD(VdwCutoffPar, n3) || RENUMBER_FIELD(VdwCutoffPar, n4);
			else return 0;
		default:
			return -1;
	}
}

#pragma mark ====== Load from Files ======

static int
s_AppendWarning(char **ptr, const char *fmt, ...)
{
	char *buf;
	int len, n, nn;
	va_list ap;
	if (ptr == NULL)
		return 0;
	if (*ptr == NULL) {
		*ptr = (char *)malloc(128);
		if (*ptr == NULL)
			return -1;
		(*ptr)[0] = 0;
		len = 128;
		n = 0;
	} else {
		n = strlen(*ptr);
		for (len = 128; len <= n; len *= 2);
	}
	va_start(ap, fmt);
	nn = vasprintf(&buf, fmt, ap);
	if (nn < 0)
		return nn;
	if (len <= n + nn) {
		while (len <= n + nn)
			len *= 2;
		*ptr = (char *)realloc(*ptr, len);
		if (*ptr == NULL)
			return -1;
	}
	strncpy(*ptr + n, buf, nn);
	*(*ptr + n + nn) = 0;
	free(buf);
	return nn;
}

int
ParameterGlobalParIndexForSrcIndex(int src)
{
	int i;
	if (src == 0 || src == -1)
		return src - 1;
	for (i = 0; i < gGlobalParInfo.count; i++) {
		if (gGlobalParInfo.files[i].src == src)
			return i;
	}
	return -2;
}

int
ParameterCommentIndexForGlobalFileName(const char *p)
{
	const char *pp = p + strlen(p), *p1 = NULL, *p2 = NULL;
	char buf[128];
	while (--pp >= p) {
		if (p2 == NULL && *pp == '.')
			p2 = pp;
		if (p1 == NULL && (*pp == '\'' || *pp == '/'))
			p1 = pp + 1;
	}
	if (p1 == NULL)
		p1 = p;
	if (p2 == NULL)
		p2 = p + strlen(p);
	snprintf(buf, sizeof(buf), "%.*s", (int)(p2 - p1), p1);
	return ParameterCommentIndex(buf);
}

int
ParameterCompare(const UnionPar *up1, const UnionPar *up2, int type)
{
	switch (type) {
		case kBondParType: {
			const BondPar *bp1 = &up1->bond;
			const BondPar *bp2 = &up2->bond;
			return (((bp1->type1 == bp2->type1 && bp1->type2 == bp2->type2)
					 || (bp1->type1 == bp2->type2 && bp1->type2 == bp2->type1))
					&& fabs(bp1->k - bp2->k) < 1e-8 && fabs(bp1->r0 - bp2->r0) < 1e-8);
		}
		case kAngleParType: {
			const AnglePar *ap1 = &up1->angle;
			const AnglePar *ap2 = &up2->angle;
			return (ap1->type2 == ap2->type2
					&& ((ap1->type1 == ap2->type1 && ap1->type3 == ap2->type3)
						|| (ap1->type1 == ap2->type3 && ap1->type3 == ap2->type1))
					&& fabs(ap1->k - ap2->k) < 1e-8 && fabs(ap1->a0 - ap2->a0) < 1e-8);
		}
		case kDihedralParType:
		case kImproperParType: {
			const TorsionPar *tp1 = &up1->torsion;
			const TorsionPar *tp2 = &up2->torsion;
			int i;
			if (tp1->mult != tp2->mult)
				return 0;
			if (type == kDihedralParType) {
				if ((tp1->type1 != tp2->type1 || tp1->type2 != tp2->type2 || tp1->type3 != tp2->type3 || tp1->type4 != tp2->type4) 
					&& (tp1->type1 != tp2->type4 || tp1->type2 != tp2->type3 || tp1->type3 != tp2->type2 || tp1->type4 != tp2->type1))
					return 0;
			} else {
				if (tp1->type3 != tp2->type3)
					return 0;
				if ((tp1->type1 != tp2->type1 || tp1->type2 != tp2->type2 || tp1->type4 != tp2->type4)
					&& (tp1->type1 != tp2->type1 || tp1->type2 != tp2->type4 || tp1->type4 != tp2->type2)
					&& (tp1->type1 != tp2->type2 || tp1->type2 != tp2->type1 || tp1->type4 != tp2->type4)
					&& (tp1->type1 != tp2->type2 || tp1->type2 != tp2->type4 || tp1->type4 != tp2->type1)
					&& (tp1->type1 != tp2->type4 || tp1->type2 != tp2->type1 || tp1->type4 != tp2->type2)
					&& (tp1->type1 != tp2->type4 || tp1->type2 != tp2->type2 || tp1->type4 != tp2->type1))
					return 0;
			}
			for (i = 0; i < tp1->mult; i++) {
				if (tp1->period[i] != tp2->period[i] || fabs(tp1->k[i] - tp2->k[i]) > 1e-8 || fabs(tp1->phi0[i] - tp2->phi0[i]) > 1e-8)
					return 0;
			}
			return 1;
		}
		case kVdwParType: {
			const VdwPar *vp1 = &up1->vdw;
			const VdwPar *vp2 = &up2->vdw;
			return (vp1->type1 == vp2->type1 && fabs(vp1->eps - vp2->eps) < 1e-6 && fabs(vp1->r_eq - vp2->r_eq) < 1e-6 && fabs(vp1->eps14 - vp2->eps14) < 1e-4 && fabs(vp1->r_eq14 - vp2->r_eq14) < 1e-4);
		}
		case kVdwPairParType: {
			const VdwPairPar *vp1 = &up1->vdwp;
			const VdwPairPar *vp2 = &up2->vdwp;
			return (vp1->type1 == vp2->type1 && fabs(vp1->eps - vp2->eps) < 1e-6 && fabs(vp1->r_eq - vp2->r_eq) < 1e-6 && fabs(vp1->eps14 - vp2->eps14) < 1e-4 && fabs(vp1->r_eq14 - vp2->r_eq14) < 1e-4);
		}
		case kVdwCutoffParType: {
			const VdwCutoffPar *vp1 = &up1->vdwcutoff;
			const VdwCutoffPar *vp2 = &up2->vdwcutoff;
			if (vp1->type != vp2->type)
				return 0;
			if (vp1->type == 0) {
				if ((vp1->n1 != vp2->n1 || vp1->n2 != vp2->n2) && (vp1->n1 != vp2->n2 || vp1->n2 != vp2->n1))
					return 0;
			} else {
				if (vp1->n1 != vp2->n1 || vp1->n2 != vp2->n2 || vp1->n3 != vp2->n3 || vp1->n4 != vp2->n4)
					return 0;
			}
			return fabs(vp1->cutoff - vp2->cutoff) < 1e-8;
		}
	}
	return 0;
}

/*  bp can also be AnglePar *, TorsionPar *, etc.  */
static void
s_StoreComment(int parType, BondPar *bp, char *p, const char *src)
{
	char *s, *pp, *pcom, buf[40];
	int embedded = 0, src_idx;
	if (p == NULL)
		return;
	while (isspace(*p) || *p == ';' || *p == '!')
		p++;
	pcom = p;
	if (src == NULL && *p == '[') {
		/*  Source is embedded  */
		int len;
		s = p + 1;
		p += 2;
		while (*p != ']' && *p != 0)
			p++;
		len = p - s;
		if (len >= sizeof(buf))
			len = sizeof(buf) - 1;
		strncpy(buf, s, len);
		buf[len] = 0;
		src = buf;
		if (*p != 0) {
			p++;
			while (isspace(*p))
				p++;
		}
		embedded = 1;
	}
	pp = p;
	while (*pp != 0 && *pp != '\r' && *pp != '\n')
		pp++;
	*pp = 0;
	if (src != NULL && *src != 0) {
		src_idx = ParameterCommentIndex(src);
		if (embedded) {
			/*  Compare with already known global parameters, and src is set if the same one is found  */
			int i;
			UnionPar *up2;
			for (i = 0; (up2 = ParameterGetUnionParFromTypeAndIndex(gBuiltinParameters, parType, i)) != NULL; i++) {
				if (up2->bond.src == src_idx && ParameterCompare((UnionPar *)bp, up2, parType)) {
					bp->src = src_idx;
					break;
				}
			}
			if (up2 == NULL) {
				/*  Not found: embedded source is retained, and this entry is regarded as "no source"  */
				bp->src = 0;
				p = pcom;
			}
		} else bp->src = src_idx;
	}
	if (*p != 0)
		bp->com = ParameterCommentIndex(p);
}

/*  bp can also be BondPar*, AnglePar *, TorsionPar *, etc.  */
static void
s_CommentToString(char *buf, int bufsize, void *bp)
{
	const char *src, *com;
	src = ParameterGetComment(((BondPar *)bp)->src);
	com = ParameterGetComment(((BondPar *)bp)->com);
	if (src == NULL && com == NULL)
		buf[0] = 0;
	else if (src != NULL)
		snprintf(buf, bufsize, " ![%s] %s", src, (com == NULL ? "" : com));
	else
		snprintf(buf, bufsize, " ! %s", com);
}

int
ParameterReadFromString(Parameter *par, char *buf, char **wbufp, const char *fname, int lineNumber, int src_idx)
{
	int i, len, options;
	char com[12], com2[12], type[4][8];
	Int itype[4];
	float val[6];  /*  not Double  */
	Int ival[12];
	int n;
	int retval = 0;
	const char *src;
	if (sscanf(buf, " %11s", com) <= 0 || !isalpha(com[0]))
		return 0;
	len = strlen(com);
	for (i = 0; i < len; i++)
		com[i] = tolower(com[i]);
	if (strncmp(com, "include", len) == 0) {
		char c, *p, *pp;
		int wc1;
		char *wbuf1;
		i = 0;
		for ( ; (c = buf[len]) != 0; len++) {
			if (c == ' ' || c == '\t')
				continue;
			if (c == '\"') {
				if (i == 0)
					continue;
				else break;
			}
			if (c == '\\') {
				len++;
				buf[i] = buf[len];
				if (buf[len] == 0)
					break;
				i++;
				continue;
			}
			buf[i++] = c;
		}
		buf[i] = 0;
		p = (char *)malloc(strlen(fname) + i + 1);
		if (p == NULL) {
			return -1;
		}
		strcpy(p, fname);
		pp = strrchr(p, '/');
		if (pp == NULL)
			strcpy(p, buf);
		else
			strcpy(p + 1, buf);
		i = ParameterReadFromFile(par, p, &wbuf1, &wc1);
		if (wbuf1 != NULL) {
			s_AppendWarning(wbufp, "In included file %s:\n%s", p, wbuf1);
			free(wbuf1);
			retval = wc1;
		}
		free(p);
		if (i < 0) {
			retval = i;
		}
		return retval;
	}
	if (par == gBuiltinParameters && src_idx == 0) {
		/*  Comes here only when the reading file is "default.par" at the initialization of the built-in parameters. */
		/*  In this case, only "include" statements are allowed.  */
		return -1;
	}
	if (src_idx == 0)
		src = NULL;
	else src = ParameterGetComment(src_idx);
	options = kParameterLookupNoWildcard | kParameterLookupNoBaseAtomType;
	if (strncmp(com, "bond", len) == 0) {
		BondPar *bp;
		if (sscanf(buf, " %11s %4s %4s %f %f %n", com2, type[0], type[1], &val[0], &val[1], &n) < 5) {
			s_AppendWarning(wbufp, "%s:%d: missing parameter in BOND record\n", fname, lineNumber);
			return 1;
		}
		itype[0] = AtomTypeEncodeToUInt(type[0]);
		itype[1] = AtomTypeEncodeToUInt(type[1]);
		if (itype[0] > itype[1]) {
			i = itype[0];
			itype[0] = itype[1];
			itype[1] = i;
		}
		val[0] *= KCAL2INTERNAL;
		bp = ParameterLookupBondPar(par, itype[0], itype[1], options);
		if (bp != NULL) {
			if (bp->k != val[0] || bp->r0 != val[1]) {
				s_AppendWarning(wbufp, "%s:%d: The BOND %s-%s parameter appeared twice; the values (%f, %f) are used\n", fname, lineNumber, AtomTypeDecodeToString(itype[0], type[0]), AtomTypeDecodeToString(itype[1], type[1]), val[0], val[1]);
				retval = 1;
			}
		}
		bp = AssignArray(&par->bondPars, &par->nbondPars, sizeof(*bp), par->nbondPars, NULL);
		bp->type1 = itype[0];
		bp->type2 = itype[1];
		bp->k = val[0];
		bp->r0 = val[1];
		s_StoreComment(kBondParType, bp, buf + n, src);
	} else if (strncmp(com, "angle", len) == 0) {
		AnglePar *ap;
		if (sscanf(buf, " %11s %4s %4s %4s %f %f %n", com2, type[0], type[1], type[2], &val[0], &val[1], &n) < 6) {
			s_AppendWarning(wbufp, "%s:%d: missing parameter in ANGLE record\n", fname, lineNumber);
			return 1;
		}
		itype[0] = AtomTypeEncodeToUInt(type[0]);
		itype[1] = AtomTypeEncodeToUInt(type[1]);
		itype[2] = AtomTypeEncodeToUInt(type[2]);
		if (itype[0] > itype[2]) {
			i = itype[0];
			itype[0] = itype[2];
			itype[2] = i;
		}
		val[0] *= KCAL2INTERNAL;
		val[1] *= (3.14159265358979 / 180.0);
		ap = ParameterLookupAnglePar(par, itype[0], itype[1], itype[2], options);
		if (ap != NULL) {
			if (ap->k != val[0] || ap->a0 != val[1]) {
				s_AppendWarning(wbufp, "%s:%d: The ANGLE %s-%s-%s parameter appeared twice; the values (%f, %f) are used\n", fname, lineNumber, AtomTypeDecodeToString(itype[0], type[0]), AtomTypeDecodeToString(itype[1], type[1]), AtomTypeDecodeToString(itype[2], type[2]), val[0], val[1]);
				retval = 1;
			}
		}
		ap = AssignArray(&par->anglePars, &par->nanglePars, sizeof(*ap), par->nanglePars, NULL);
		ap->type1 = itype[0];
		ap->type2 = itype[1];
		ap->type3 = itype[2];
		ap->k = val[0];
		ap->a0 = val[1];
		s_StoreComment(kAngleParType, (BondPar *)ap, buf + n, src);
	} else if (strncmp(com, "dihedral", len) == 0) {
		TorsionPar *dp;
		if (sscanf(buf, " %11s %4s %4s %4s %4s %f %d %f %n", com2, type[0], type[1], type[2], type[3], &val[0], &ival[0], &val[1], &n) < 8) {
			s_AppendWarning(wbufp, "%s:%d: missing parameter in DIHEDRAL record\n", fname, lineNumber);
			return 1;
		}
		itype[0] = AtomTypeEncodeToUInt(type[0]);
		itype[1] = AtomTypeEncodeToUInt(type[1]);
		itype[2] = AtomTypeEncodeToUInt(type[2]);
		itype[3] = AtomTypeEncodeToUInt(type[3]);
		if (itype[0] > itype[3]) {
			i = itype[0];
			itype[0] = itype[3];
			itype[3] = i;
			i = itype[1];
			itype[1] = itype[2];
			itype[2] = i;
		}
		dp = ParameterLookupDihedralPar(par, itype[0], itype[1], itype[2], itype[3], options);
		val[0] *= KCAL2INTERNAL;
		val[1] *= 3.14159265358979 / 180.0;
		if (dp != NULL) {
			if (dp->mult != 1 || dp->k[0] != val[0] || dp->period[0] != ival[0] || dp->phi0[0] != val[1]) {
				s_AppendWarning(wbufp, "%s:%d: The DIHEDRAL %s-%s-%s-%s parameter appeared twice; the values (%f, %d, %f) are used\n", fname, lineNumber, AtomTypeDecodeToString(itype[0], type[0]), AtomTypeDecodeToString(itype[1], type[1]), AtomTypeDecodeToString(itype[2], type[2]), AtomTypeDecodeToString(itype[3], type[3]), val[0], ival[0], val[1]);
				retval = 1;
			}
		}
		dp = AssignArray(&par->dihedralPars, &par->ndihedralPars, sizeof(*dp), par->ndihedralPars, NULL);
		dp->type1 = itype[0];
		dp->type2 = itype[1];
		dp->type3 = itype[2];
		dp->type4 = itype[3];
		dp->k[0] = val[0];
		dp->period[0] = ival[0];
		dp->phi0[0] = val[1];
		dp->mult = 1;
		s_StoreComment(kDihedralParType, (BondPar *)dp, buf + n, src);
	} else if (strncmp(com, "improper", len) == 0) {
		TorsionPar *ip;
		if (sscanf(buf, " %11s %4s %4s %4s %4s %f %d %f %n", com2, type[0], type[1], type[2], type[3], &val[0], &ival[0], &val[1], &n) < 8) {
			s_AppendWarning(wbufp, "%s:%d: missing parameter in IMPROPER record\n", fname, lineNumber);
			return 1;
		}
		itype[0] = AtomTypeEncodeToUInt(type[0]);
		itype[1] = AtomTypeEncodeToUInt(type[1]);
		itype[2] = AtomTypeEncodeToUInt(type[2]);
		itype[3] = AtomTypeEncodeToUInt(type[3]);
		if (itype[0] > itype[1]) {
			i = itype[0];
			itype[0] = itype[1];
			itype[1] = i;
		}
		if (itype[0] > itype[3]) {
			i = itype[0];
			itype[0] = itype[3];
			itype[3] = i;
		}
		if (itype[1] > itype[3]) {
			i = itype[1];
			itype[1] = itype[3];
			itype[3] = i;
		}
		ip = ParameterLookupImproperPar(par, itype[0], itype[1], itype[2], itype[3], options);
		val[0] *= KCAL2INTERNAL;
		val[1] *= 3.14159265358979 / 180.0;
		if (ip != NULL) {
			if (ip->mult != 1 || ip->k[0] != val[0] || ip->period[0] != ival[0] || ip->phi0[0] != val[1]) {
				s_AppendWarning(wbufp, "%s:%d: The IMPROPER %s-%s-%s-%s parameter appeared twice; the values (%f, %d, %f) are used\n", fname, lineNumber, AtomTypeDecodeToString(itype[0], type[0]), AtomTypeDecodeToString(itype[1], type[1]), AtomTypeDecodeToString(itype[2], type[2]), AtomTypeDecodeToString(itype[3], type[3]), val[0], ival[0], val[1]);
				retval = 1;
			}
		}
		ip = AssignArray(&par->improperPars, &par->nimproperPars, sizeof(*ip), par->nimproperPars, NULL);
		ip->type1 = itype[0];
		ip->type2 = itype[1];
		ip->type3 = itype[2];
		ip->type4 = itype[3];
		ip->k[0] = val[0];
		ip->period[0] = ival[0];
		ip->phi0[0] = val[1];
		ip->mult = 1;	
		s_StoreComment(kImproperParType, (BondPar *)ip, buf + n, src);
	} else if (strncmp(com, "nonbonded", len) == 0 || strncmp(com, "vdw", len) == 0) {
		VdwPar *vp, vtemp;
		char *p;
		/*  NOTE: the nonbonded record lists "2*sigma", not "sigma"!  */
		int flag = (com[0] == 'v');
		if (sscanf(buf, " %11s %4s %f %f %f %f %n", com2, type[0], &val[0], &val[1], &val[2], &val[3], &n) < 6) {
			s_AppendWarning(wbufp, "%s:%d: missing parameter in %s record\n", fname, lineNumber, (flag ? "VDW" : "NONBONDED"));
			return 1;
		}
		itype[0] = AtomTypeEncodeToUInt(type[0]);
		memset(&vtemp, 0, sizeof(vtemp));
		vtemp.type1 = itype[0];
		vtemp.atomicNumber = 0;  /*  No definition given  */
		vtemp.eps = val[0] * KCAL2INTERNAL;
		vtemp.r_eq = val[1] * (flag ? 1.0 : 0.561231024154687); /* 1/2 * 2**(1/6)  */
		vtemp.A = pow(vtemp.r_eq * 2, 12) * vtemp.eps;
		vtemp.B = 2 * pow(vtemp.r_eq * 2, 6) * vtemp.eps;
		vtemp.eps14 = val[2] * KCAL2INTERNAL;
		vtemp.r_eq14 = val[3] * (flag ? 1.0 : 0.561231024154687); /* 1/2 * 2**(1/6)  */
		vtemp.A14 = pow(vtemp.r_eq14 * 2, 12) * vtemp.eps14;
		vtemp.B14 = 2 * pow(vtemp.r_eq14 * 2, 6) * vtemp.eps14;
		p = buf + n;
		ival[0] = 0;
		val[4] = val[5] = 0.0;
		if (sscanf(p, "%d %n", &ival[0], &n) == 1) {
			vtemp.atomicNumber = ival[0];
			p += n;
		}
		if (sscanf(p, "%f %n", &val[4], &n) == 1) {
			vtemp.weight = val[4];
			p += n;
		}
		if (val[4] == 0.0 && ival[0] != 0)
			vtemp.weight = WeightForAtomicNumber(ival[0]);
		if (sscanf(p, "%f %n", &val[5], &n) == 1) {
			vtemp.polarizability = val[5];
			p += n;
		}
		n = p - buf;
		if (ival[0] == 0 && val[4] != 0.0) {
			for (i = 1; (val[5] = WeightForAtomicNumber(i)) != 0.0; i++) {
				if (fabs(val[4] - val[5]) < 0.1) {
					vtemp.atomicNumber = i;
					break;
				}
			}
		}
		vp = NULL;
		for (i = 0; i < par->nvdwPars; i++) {
			if (itype[0] == par->vdwPars[i].type1) {
				vp = par->vdwPars + i;
				if (vp->A != vtemp.A || vp->B != vtemp.B || vp->A14 != vtemp.A14 || vp->B14 != vtemp.B14) {
					s_AppendWarning(wbufp, "%s:%d: The %s %s parameter appeared twice; the values (%f, %f, %f, %f, %d, %f, %f) are used\n", fname, lineNumber, (flag ? "VDW" : "NONBONDED"), AtomTypeDecodeToString(itype[0], type[0]), val[0], val[1], val[2], val[3], ival[0], val[4], val[5]);
					retval = 1;
				}
				break;
			}
		}
		vp = AssignArray(&par->vdwPars, &par->nvdwPars, sizeof(*vp), i, NULL);
		vtemp.com = vp->com;  /*  Keep comment field  */
		*vp = vtemp;
		s_StoreComment(kVdwParType, (BondPar *)vp, buf + n, src);
	} else if (strncmp(com, "nbfi", len) == 0 || strncmp(com, "vdwpair", len) == 0) {
		VdwPairPar *vp, vtemp;
		int flag = (com[0] == 'v');
		if (sscanf(buf, " %11s %4s %4s %f %f %f %f %n", com2, type[0], type[1], &val[0], &val[1], &val[2], &val[3], &n) < 6) {
			s_AppendWarning(wbufp, "%s:%d: missing parameter in %s record\n", fname, lineNumber, (flag ? "VDWP" : "NBFI"));
			return 1;
		}
		itype[0] = AtomTypeEncodeToUInt(type[0]);
		itype[1] = AtomTypeEncodeToUInt(type[1]);
		if (itype[0] > itype[1]) {
			i = itype[0];
			itype[0] = itype[1];
			itype[1] = i;
		}
		vtemp.type1 = itype[0];
		vtemp.type2 = itype[1];
		if (flag) {  /*  eps/r_eq representation  */
			vtemp.eps = val[0] * KCAL2INTERNAL;
			vtemp.r_eq = val[1];
			vtemp.A = pow(val[1] * 2, 12) * vtemp.eps;
			vtemp.B = 2 * pow(val[1] * 2, 6) * vtemp.eps;
			vtemp.eps14 = val[2] * KCAL2INTERNAL;
			vtemp.r_eq14 = val[3];
			vtemp.A14 = pow(val[3] * 2, 12) * vtemp.eps14;
			vtemp.B14 = 2 * pow(val[3] * 2, 6) * vtemp.eps14;
		} else {    /*  A/B representation  */
			vtemp.A = val[0] * KCAL2INTERNAL;
			vtemp.B = val[1] * KCAL2INTERNAL;
			vtemp.eps = pow(0.25 * vtemp.B * vtemp.B / vtemp.A, 0.16666666667);
			vtemp.r_eq = pow(vtemp.A / vtemp.B * 2.0, 0.16666666667) * 0.5;
			vtemp.A14 = val[2] * KCAL2INTERNAL;
			vtemp.B14 = val[3] * KCAL2INTERNAL;
			vtemp.eps14 = pow(0.25 * vtemp.B14 * vtemp.B14 / vtemp.A14, 0.16666666667);
			vtemp.r_eq14 = pow(vtemp.A14 / vtemp.B14 * 2.0, 0.16666666667) * 0.5;
		}
		vp = NULL;
		for (i = 0; i < par->nvdwpPars; i++) {
			if (itype[0] == par->vdwpPars[i].type1 && itype[1] == par->vdwpPars[i].type2) {
				vp = par->vdwpPars + i;
				if (vp->A != vtemp.A || vp->B != vtemp.B || vp->A14 != vtemp.A14 || vp->B14 != vtemp.B14) {
					s_AppendWarning(wbufp, "%s:%d: The %s %s-%s parameter appeared twice; the values (%f, %f, %f, %f) are used\n", fname, lineNumber, (flag ? "VDWP" : "NBFI"), AtomTypeDecodeToString(itype[0], type[0]), AtomTypeDecodeToString(itype[1], type[1]), val[0], val[1], val[2], val[3]);
					retval = 1;
				}
				break;
			}
		}
		vp = AssignArray(&par->vdwpPars, &par->nvdwpPars, sizeof(*vp), i, NULL);
		vtemp.com = vp->com;  /*  Keep comment field  */
		*vp = vtemp;
		s_StoreComment(kVdwPairParType, (BondPar *)vp, buf + n, src);
	} else {
		s_AppendWarning(wbufp, "%s:%d: unknown keyword %s\n", fname, lineNumber, com);
		return 1;
	}
	return retval;
}

int
ParameterReadFromFile(Parameter *par, const char *fname, char **outWarningMessage, int *outWarningCount)
{
	char *wbuf;
	char buf[1024];
	FILE *fp = NULL;
	int first = 0;
	int wcount;
	int lineNumber;
	int src_idx;
	
	if (par == NULL) {
		par = gBuiltinParameters;
		if (par == NULL) {
			par = ParameterNew();
			gBuiltinParameters = par;
			first = 1;
		}
	}
	
	wcount = 0;
	wbuf = NULL;

	fp = fopen(fname, "rb");
	if (fp == NULL) {
		s_AppendWarning(&wbuf, "Cannot open parameter file %s\n", fname);
		wcount = -1;
		goto exit;
	}

	if (par != gBuiltinParameters || first)
		src_idx = 0;
	else
		src_idx = ParameterCommentIndexForGlobalFileName(fname);

	if (par == gBuiltinParameters && !first) {
		/*  Ensure the "source" field is unique  */
		int i;
		ParFileInfo *ip;
		const char *p;
		for (i = 0; i < gGlobalParInfo.count; i++) {
			if (gGlobalParInfo.files[i].src == src_idx)
				break;
		}
		if (i < gGlobalParInfo.count) {
			/*  Delete the existing Parameters from the same source  */
			ParameterDeleteAllEntriesForSource(par, src_idx);
		}
		
		/*  Register the global file info  */
		ip = AssignArray(&(gGlobalParInfo.files), &(gGlobalParInfo.count), sizeof(ParFileInfo), gGlobalParInfo.count, NULL);
		for (p = fname + strlen(fname) - 1; p >= fname; p--) {
			if (*p == '/' || *p == '\\')
				break;
		}
		if (p < fname)
			ip->dir = NULL;
		else {
			i = p - fname;
			ip->dir = (char *)malloc(i + 1);
			if (ip->dir != NULL) {
				strncpy(ip->dir, fname, i);
				ip->dir[i] = 0;
			}
		}
		p++;
		ip->name = strdup(p);
		ip->src = src_idx;
	}
	
	lineNumber = 0;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		int wc1 = ParameterReadFromString(par, buf, &wbuf, fname, lineNumber, src_idx);
		if (wc1 >= 0)
			wcount += wc1;
		else {
			wcount = wc1;
			break;
		}
	}
	if (first)
		gGlobalParInfo.builtinCount = gGlobalParInfo.count;

exit:
	if (fp != NULL)
		fclose(fp);
	if (outWarningMessage != NULL)
		*outWarningMessage = wbuf;
	else if (wbuf != NULL)
		free(wbuf);
	if (outWarningCount != NULL)
		*outWarningCount = (wcount >= 0 ? wcount : 0);

	return (wcount >= 0 ? 0 : wcount);
}

int
ParameterAppendToFile(Parameter *par, FILE *fp)
{
	int i, n;
	char cbuf[6][8];
	char buf[120];
	int bufsize = sizeof(buf);

	if (par == NULL)
		return 0;
	
	n = 0;
	if (par->nbondPars > 0)
		fprintf(fp, "! type1 type2 k r0\n");
	for (i = 0; i < par->nbondPars; i++) {
		BondPar *bp = par->bondPars + i;
		AtomTypeDecodeToString(bp->type1, cbuf[0]);
		AtomTypeDecodeToString(bp->type2, cbuf[1]);
		s_CommentToString(buf, bufsize, bp);
		fprintf(fp, "bond %s %s %.6f %f%s\n", cbuf[0], cbuf[1], bp->k * INTERNAL2KCAL, bp->r0, buf);
		n++;
	}
	if (par->nanglePars > 0)
		fprintf(fp, "! type1 type2 type3 k a0\n");
	for (i = 0; i < par->nanglePars; i++) {
		AnglePar *ap = par->anglePars + i;
		AtomTypeDecodeToString(ap->type1, cbuf[0]);
		AtomTypeDecodeToString(ap->type2, cbuf[1]);
		AtomTypeDecodeToString(ap->type3, cbuf[2]);
		s_CommentToString(buf, bufsize, ap);
		fprintf(fp, "angle %s %s %s %.6f %f%s\n", cbuf[0], cbuf[1], cbuf[2], ap->k * INTERNAL2KCAL, ap->a0 * kRad2Deg, buf);
		n++;
	}
	if (par->ndihedralPars > 0)
		fprintf(fp, "! type1 type2 type3 type4 k periodicity phi0\n");
	for (i = 0; i < par->ndihedralPars; i++) {
		TorsionPar *tp = par->dihedralPars + i;
		AtomTypeDecodeToString(tp->type1, cbuf[0]);
		AtomTypeDecodeToString(tp->type2, cbuf[1]);
		AtomTypeDecodeToString(tp->type3, cbuf[2]);
		AtomTypeDecodeToString(tp->type4, cbuf[3]);
		s_CommentToString(buf, bufsize, tp);
		fprintf(fp, "dihe %s %s %s %s %.6f %d %f%s\n", cbuf[0], cbuf[1], cbuf[2], cbuf[3], tp->k[0] * INTERNAL2KCAL, tp->period[0], tp->phi0[0] * kRad2Deg, buf);
		n++;
	}
	if (par->nimproperPars > 0)
		fprintf(fp, "! type1 type2 type3 type4 k periodicity phi0\n");
	for (i = 0; i < par->nimproperPars; i++) {
		TorsionPar *tp = par->improperPars + i;
		AtomTypeDecodeToString(tp->type1, cbuf[0]);
		AtomTypeDecodeToString(tp->type2, cbuf[1]);
		AtomTypeDecodeToString(tp->type3, cbuf[2]);
		AtomTypeDecodeToString(tp->type4, cbuf[3]);
		s_CommentToString(buf, bufsize, tp);
		fprintf(fp, "impr %s %s %s %s %.6f %d %f%s\n", cbuf[0], cbuf[1], cbuf[2], cbuf[3], tp->k[0] * INTERNAL2KCAL, tp->period[0], tp->phi0[0] * kRad2Deg, buf);
		n++;
	}
	if (par->nvdwPars > 0)
		fprintf(fp, "! type eps r_eq eps14 r_eq14 atomic_number weight\n");
	for (i = 0; i < par->nvdwPars; i++) {
		VdwPar *vp = par->vdwPars + i;
	/*	Double eps, eps14;  */
		AtomTypeDecodeToString(vp->type1, cbuf[0]);
	/*	eps = (vp->A == 0.0 ? 0.0 : vp->B * vp->B / vp->A * 0.25 * INTERNAL2KCAL);
		eps14 = (vp->A14 == 0.0 ? 0.0 : vp->B14 * vp->B14 / vp->A14 * 0.25 * INTERNAL2KCAL);  */
		s_CommentToString(buf, bufsize, vp);
		fprintf(fp, "vdw %s %.6f %.6f %.6f %.6f %d %f%s\n", cbuf[0], vp->eps * INTERNAL2KCAL, vp->r_eq, vp->eps14 * INTERNAL2KCAL, vp->r_eq14, vp->atomicNumber, vp->weight, buf);  /*  polarizability is not written because it is not used now  */
		n++;
	}
	if (par->nvdwpPars > 0)
		fprintf(fp, "! type1 type2 eps r_eq eps14 r_eq14\n");
	for (i = 0; i < par->nvdwpPars; i++) {
		VdwPairPar *vpp = par->vdwpPars + i;
	/*	Double eps, eps14;  */
		AtomTypeDecodeToString(vpp->type1, cbuf[0]);
		AtomTypeDecodeToString(vpp->type2, cbuf[1]);
	/*	eps = (vpp->A == 0.0 ? 0.0 : vpp->B * vpp->B / vpp->A * 0.25 * INTERNAL2KCAL);
		eps14 = (vpp->A14 == 0.0 ? 0.0 : vpp->B14 * vpp->B14 / vpp->A14 * 0.25 * INTERNAL2KCAL); */
		s_CommentToString(buf, bufsize, vpp);
		fprintf(fp, "vdwp %s %s %.6f %.6f %.6f %.6f%s\n", cbuf[0], cbuf[1], vpp->eps * INTERNAL2KCAL, vpp->r_eq, vpp->eps14 * INTERNAL2KCAL, vpp->r_eq14, buf);
		n++;
	}
/*	if (par->natomPars > 0)
		fprintf(fp, "! name atomic_number radius red green blue weight\n");
	for (i = 0; i < par->natomPars; i++) {
		ElementPar *app = par->atomPars + i;
		s_CommentToString(buf, bufsize, app);
		fprintf(fp, "element %.4s %d %f %f %f %f %f%s\n", app->name, app->number, app->radius, app->r, app->g, app->b, app->weight, buf);
		n++;
	} */
	return n;
}

int
ParameterDeleteAllEntriesForSource(Parameter *par, int src_idx)
{
	int i, n, type;
	UnionPar *up;
	IntGroup *ig;
/*	if (fname == NULL)
		return -1;
	src_idx = ParameterCommentIndexForGlobalFileName(fname); */
	if (src_idx == 0)
		return -1;
	n = 0;
	for (type = kFirstParType; type <= kLastParType; type++) {
		ig = IntGroupNew();
		for (i = ParameterGetCountForType(par, type) - 1; i >= 0; i--) {
			up = ParameterGetUnionParFromTypeAndIndex(par, type, i);
			if (up != NULL && up->bond.src == src_idx)
				IntGroupAdd(ig, i, 1);
		}
		i = IntGroupGetCount(ig);
		if (i > 0) {
			ParameterDelete(par, type, NULL, ig);
			n += i;
		}
		IntGroupRelease(ig);
	}
	if (par == gBuiltinParameters) {
		/*  Unregister from the global info  */
		for (i = gGlobalParInfo.builtinCount; i < gGlobalParInfo.count; i++) {
			if (gGlobalParInfo.files[i].src == src_idx) {
				DeleteArray(&gGlobalParInfo.files, &gGlobalParInfo.count, sizeof(ParFileInfo), i, 1, NULL);
				break;
			}
		}
	}
	return n;
}

#pragma mark ====== Parameter Comments ======

static const char **sParameterComments;
static Int sNumberOfParameterComments;

int
ParameterCommentIndex(const char *comment)
{
	int i;
	char *p, *pp;
	if (comment == NULL || comment[0] == 0)
		return 0;
	/*  Duplicate the comment, convert whitespaces to ' ', and chop trailing whitespaces  */
	p = strdup(comment);
	i = 0;
	for (pp = p + strlen(p) - 1; pp >= p; pp--) {
		if (isspace(*pp)) {
			if (i == 0)
				*pp = 0;
			else *pp = ' ';
		} else i++;
	}
	for (i = 1; i < sNumberOfParameterComments; i++) {
		if (strcmp(sParameterComments[i], p) == 0) {
			free(p);
			return i;
		}
	}
	if (sNumberOfParameterComments == 0) {
		/*  Index 0 is skipped  */
		AssignArray(&sParameterComments, &sNumberOfParameterComments, sizeof(char *), 1, &p);
	} else {
		AssignArray(&sParameterComments, &sNumberOfParameterComments, sizeof(char *), sNumberOfParameterComments, &p);
	}
	return sNumberOfParameterComments - 1;
}

const char *
ParameterGetComment(int idx)
{
	if (idx <= 0 || idx >= sNumberOfParameterComments)
		return NULL;  /*  No such number  */
	return sParameterComments[idx];
}

#pragma mark ====== Parameter Lookup ======

#define s_ParMatch(_t1, _t2, _nowildcard) (_t1 == _t2 || (!_nowildcard && _t1 == kAtomTypeWildcard))

/*  Returns non-zero if the parameter record contains designated atom_type.
 The atom_type can also be an atom index.  */
int
ParameterDoesContainAtom(Int type, UnionPar *up, UInt atom_type, Int options)
{
#define CHECK_FIELD(_tp, _field) s_ParMatch((((_tp *)up)->_field), atom_type, nowildcard)
	Int nowildcard = (options & kParameterLookupNoWildcard);
	switch (type) {
		case kBondParType:
			return CHECK_FIELD(BondPar, type1) || CHECK_FIELD(BondPar, type2);
		case kAngleParType:
			return CHECK_FIELD(AnglePar, type1) || CHECK_FIELD(AnglePar, type2) || CHECK_FIELD(AnglePar, type3);
		case kDihedralParType:
		case kImproperParType:
			return CHECK_FIELD(TorsionPar, type1) || CHECK_FIELD(TorsionPar, type2) || CHECK_FIELD(TorsionPar, type3) || CHECK_FIELD(TorsionPar, type4);
		case kVdwParType:
			return CHECK_FIELD(VdwPar, type1);
		case kVdwPairParType:
			return CHECK_FIELD(VdwPairPar, type1) || CHECK_FIELD(VdwPairPar, type2);
		case kVdwCutoffParType:
			if (((VdwCutoffPar *)up)->type == 1)
				return CHECK_FIELD(VdwCutoffPar, n1) || CHECK_FIELD(VdwCutoffPar, n2) || CHECK_FIELD(VdwCutoffPar, n3) || CHECK_FIELD(VdwCutoffPar, n4);
			else return 0;
		default: return 0;
	}
}

BondPar *
ParameterLookupBondPar(Parameter *par, UInt t1, UInt t2, Int options)
{
	int i;
	BondPar *bp;
	Int nowildcard = (options & kParameterLookupNoWildcard);
	if (par == NULL)
		return NULL;
	if ((options & (kParameterLookupGlobal | kParameterLookupLocal | kParameterLookupMissing)) == 0)
		options |= kParameterLookupGlobal | kParameterLookupLocal;
	for (i = par->nbondPars - 1, bp = par->bondPars + i; i >= 0; i--, bp--) {
		if ((bp->src > 0 && !(options & kParameterLookupGlobal))
			|| (bp->src == 0 && !(options & kParameterLookupLocal))
			|| (bp->src < 0 && !(options & kParameterLookupMissing)))
			continue;
		if (s_ParMatch(bp->type1, t1, nowildcard) && s_ParMatch(bp->type2, t2, nowildcard))
			return bp;
		if (s_ParMatch(bp->type1, t2, nowildcard) && s_ParMatch(bp->type2, t1, nowildcard))
			return bp;		
	}
	if (options & kParameterLookupNoBaseAtomType)
		return NULL;
	return ParameterLookupBondPar(par, t1 % kAtomTypeVariantBase, t2 % kAtomTypeVariantBase, options | kParameterLookupNoBaseAtomType);
}

AnglePar *
ParameterLookupAnglePar(Parameter *par, UInt t1, UInt t2, UInt t3, Int options)
{
	int i;
	AnglePar *ap;
	Int nowildcard = (options & kParameterLookupNoWildcard);
	if (par == NULL)
		return NULL;
	if ((options & (kParameterLookupGlobal | kParameterLookupLocal | kParameterLookupMissing)) == 0)
		options |= kParameterLookupGlobal | kParameterLookupLocal;
	for (i = par->nanglePars - 1, ap = par->anglePars + i; i >= 0; i--, ap--) {
		if ((ap->src > 0 && !(options & kParameterLookupGlobal))
			|| (ap->src == 0 && !(options & kParameterLookupLocal))
			|| (ap->src < 0 && !(options & kParameterLookupMissing)))
			continue;
		if (s_ParMatch(ap->type1, t1, nowildcard) && s_ParMatch(ap->type2, t2, nowildcard) && s_ParMatch(ap->type3, t3, nowildcard))
			return ap;
		if (s_ParMatch(ap->type1, t3, nowildcard) && s_ParMatch(ap->type2, t2, nowildcard) && s_ParMatch(ap->type3, t1, nowildcard))
			return ap;
	}
	if (options & kParameterLookupNoBaseAtomType)
		return NULL;
	return ParameterLookupAnglePar(par, t1 % kAtomTypeVariantBase, t2 % kAtomTypeVariantBase, t3 % kAtomTypeVariantBase, options | kParameterLookupNoBaseAtomType);
}

TorsionPar *
ParameterLookupDihedralPar(Parameter *par, UInt t1, UInt t2, UInt t3, UInt t4, Int options)
{
	int i;
	TorsionPar *tp;
	Int nowildcard = (options & kParameterLookupNoWildcard);
	if (par == NULL)
		return NULL;
	if ((options & (kParameterLookupGlobal | kParameterLookupLocal | kParameterLookupMissing)) == 0)
		options |= kParameterLookupGlobal | kParameterLookupLocal;
	for (i = par->ndihedralPars - 1, tp = par->dihedralPars + i; i >= 0; i--, tp--) {
		if ((tp->src > 0 && !(options & kParameterLookupGlobal))
			|| (tp->src == 0 && !(options & kParameterLookupLocal))
			|| (tp->src < 0 && !(options & kParameterLookupMissing)))
			continue;
		if (s_ParMatch(tp->type1, t1, nowildcard) && s_ParMatch(tp->type2, t2, nowildcard) && s_ParMatch(tp->type3, t3, nowildcard) && s_ParMatch(tp->type4, t4, nowildcard))
			return tp;
		if (s_ParMatch(tp->type1, t4, nowildcard) && s_ParMatch(tp->type2, t3, nowildcard) && s_ParMatch(tp->type3, t2, nowildcard) && s_ParMatch(tp->type4, t1, nowildcard))
			return tp;
	}
	if (options & kParameterLookupNoBaseAtomType)
		return NULL;
	return ParameterLookupDihedralPar(par, t1 % kAtomTypeVariantBase, t2 % kAtomTypeVariantBase, t3 % kAtomTypeVariantBase, t4 % kAtomTypeVariantBase, options | kParameterLookupNoBaseAtomType);
}

TorsionPar *
ParameterLookupImproperPar(Parameter *par, UInt t1, UInt t2, UInt t3, UInt t4, Int options)
{
	int i;
	TorsionPar *tp;
	Int nowildcard = (options & kParameterLookupNoWildcard);
	if (par == NULL)
		return NULL;
	if ((options & (kParameterLookupGlobal | kParameterLookupLocal | kParameterLookupMissing)) == 0)
		options |= kParameterLookupGlobal | kParameterLookupLocal;
	for (i = par->nimproperPars - 1, tp = par->improperPars + i; i >= 0; i--, tp--) {
		if ((tp->src > 0 && !(options & kParameterLookupGlobal))
			|| (tp->src == 0 && !(options & kParameterLookupLocal))
			|| (tp->src < 0 && !(options & kParameterLookupMissing)))
			continue;
		if (!s_ParMatch(tp->type3, t3, nowildcard))
			continue;
		if ((s_ParMatch(tp->type1, t1, nowildcard) && s_ParMatch(tp->type2, t2, nowildcard) && s_ParMatch(tp->type4, t4, nowildcard))
			|| (s_ParMatch(tp->type1, t1, nowildcard) && s_ParMatch(tp->type2, t4, nowildcard) && s_ParMatch(tp->type4, t2, nowildcard))
			|| (s_ParMatch(tp->type1, t2, nowildcard) && s_ParMatch(tp->type2, t1, nowildcard) && s_ParMatch(tp->type4, t4, nowildcard))
			|| (s_ParMatch(tp->type1, t2, nowildcard) && s_ParMatch(tp->type2, t4, nowildcard) && s_ParMatch(tp->type4, t1, nowildcard))
			|| (s_ParMatch(tp->type1, t4, nowildcard) && s_ParMatch(tp->type2, t1, nowildcard) && s_ParMatch(tp->type4, t2, nowildcard))
			|| (s_ParMatch(tp->type1, t4, nowildcard) && s_ParMatch(tp->type2, t2, nowildcard) && s_ParMatch(tp->type4, t1, nowildcard)))
			return tp;
	}
	if (options & kParameterLookupNoBaseAtomType)
		return NULL;
	return ParameterLookupImproperPar(par, t1 % kAtomTypeVariantBase, t2 % kAtomTypeVariantBase, t3 % kAtomTypeVariantBase, t4 % kAtomTypeVariantBase, options | kParameterLookupNoBaseAtomType);
}

VdwPar *
ParameterLookupVdwPar(Parameter *par, UInt t1, Int options)
{
	int i;
	VdwPar *vp;
	Int nowildcard = (options & kParameterLookupNoWildcard);
	if (par == NULL)
		return NULL;
	if ((options & (kParameterLookupGlobal | kParameterLookupLocal | kParameterLookupMissing)) == 0)
		options |= kParameterLookupGlobal | kParameterLookupLocal;
	for (i = par->nvdwPars - 1, vp = par->vdwPars + i; i >= 0; i--, vp--) {
		if ((vp->src > 0 && !(options & kParameterLookupGlobal))
			|| (vp->src == 0 && !(options & kParameterLookupLocal))
			|| (vp->src < 0 && !(options & kParameterLookupMissing)))
			continue;
		if (s_ParMatch(vp->type1, t1, nowildcard))
			return vp;
	}
	if (options & kParameterLookupNoBaseAtomType)
		return NULL;
	return ParameterLookupVdwPar(par, t1 % kAtomTypeVariantBase, options | kParameterLookupNoBaseAtomType);
}

VdwPairPar *
ParameterLookupVdwPairPar(Parameter *par, UInt t1, UInt t2, Int options)
{
	int i;
	VdwPairPar *vp;
	Int nowildcard = (options & kParameterLookupNoWildcard);
	if (par == NULL)
		return NULL;
	if ((options & (kParameterLookupGlobal | kParameterLookupLocal | kParameterLookupMissing)) == 0)
		options |= kParameterLookupGlobal | kParameterLookupLocal;
	for (i = par->nvdwpPars - 1, vp = par->vdwpPars + i; i >= 0; i--, vp--) {
		if ((vp->src > 0 && !(options & kParameterLookupGlobal))
			|| (vp->src == 0 && !(options & kParameterLookupLocal))
			|| (vp->src < 0 && !(options & kParameterLookupMissing)))
			continue;
		if ((s_ParMatch(vp->type1, t1, nowildcard) && s_ParMatch(vp->type2, t2, nowildcard))
			|| (s_ParMatch(vp->type1, t2, nowildcard) && s_ParMatch(vp->type2, t1, nowildcard)))
			return vp;
	}
	if (options & kParameterLookupNoBaseAtomType)
		return NULL;
	return ParameterLookupVdwPairPar(par, t1 % kAtomTypeVariantBase, t2 % kAtomTypeVariantBase, options | kParameterLookupNoBaseAtomType);
}

VdwCutoffPar *
ParameterLookupVdwCutoffPar(Parameter *par, UInt t1, UInt t2, Int options)
{
	int i;
	VdwCutoffPar *vp;
	Int nowildcard = (options & kParameterLookupNoWildcard);
	if (par == NULL)
		return NULL;
	if ((options & (kParameterLookupGlobal | kParameterLookupLocal | kParameterLookupMissing)) == 0)
		options |= kParameterLookupGlobal | kParameterLookupLocal;
	for (i = par->nvdwCutoffPars - 1, vp = par->vdwCutoffPars + i; i >= 0; i--, vp--) {
		if ((vp->src > 0 && !(options & kParameterLookupGlobal))
			|| (vp->src == 0 && !(options & kParameterLookupLocal))
			|| (vp->src < 0 && !(options & kParameterLookupMissing)))
			continue;
		if (vp->type == 0) {
			if (s_ParMatch(vp->n1, t1, nowildcard) && s_ParMatch(vp->n2, t2, nowildcard))
				return vp;
			if (s_ParMatch(vp->n1, t2, nowildcard) && s_ParMatch(vp->n2, t1, nowildcard))
				return vp;
		} else {
			if (vp->n1 <= t1 && vp->n2 >= t1 && vp->n3 <= t2 && vp->n4 <= t2)
				return vp;
			if (vp->n1 <= t2 && vp->n2 >= t2 && vp->n3 <= t1 && vp->n4 >= t1)
				return vp;
		}
	}
	if (options & kParameterLookupNoBaseAtomType)
		return NULL;
	return ParameterLookupVdwCutoffPar(par, t1 % kAtomTypeVariantBase, t2 % kAtomTypeVariantBase, options | kParameterLookupNoBaseAtomType);
}

#pragma mark ====== Table View Support ======

int
ParameterTableNumberOfRows(Parameter *par)
{
	if (par == NULL)
		return 0;
	return par->nvdwPars + par->nbondPars + par->nanglePars + par->ndihedralPars + par->nimproperPars + par->nvdwpPars + 6;
}

int
ParameterTableGetItemIndex(Parameter *par, int row, int *type)
{
	if (par == NULL || row < 0) {
		*type = kInvalidParType;
		return -1;
	}
	if (--row < par->nvdwPars) {
		*type = kVdwParType;
	} else if ((row -= par->nvdwPars + 1) < par->nbondPars) {
		*type = kBondParType;
	} else if ((row -= par->nbondPars + 1) < par->nanglePars) {
		*type = kAngleParType;
	} else if ((row -= par->nanglePars + 1) < par->ndihedralPars) {
		*type = kDihedralParType;
	} else if ((row -= par->ndihedralPars + 1) < par->nimproperPars) {
		*type = kImproperParType;
	} else if ((row -= par->nimproperPars + 1) < par->nvdwpPars) {
		*type = kVdwPairParType;
	} else {
		*type = kInvalidParType;
		return -1;
	}
	return row;
}

UnionPar *
ParameterTableGetItemPtr(Parameter *par, int row, int *type)
{
	if (par == NULL || row < 0) {
		*type = kInvalidParType;
		return NULL;
	}
	if (--row < par->nvdwPars) {
		*type = kVdwParType;
		return (UnionPar *)(row >= 0 ? par->vdwPars + row : NULL);
	} else if ((row -= par->nvdwPars + 1) < par->nbondPars) {
		*type = kBondParType;
		return (UnionPar *)(row >= 0 ? par->bondPars + row : NULL);
	} else if ((row -= par->nbondPars + 1) < par->nanglePars) {
		*type = kAngleParType;
		return (UnionPar *)(row >= 0 ? par->anglePars + row : NULL);
	} else if ((row -= par->nanglePars + 1) < par->ndihedralPars) {
		*type = kDihedralParType;
		return (UnionPar *)(row >= 0 ? par->dihedralPars + row : NULL);
	} else if ((row -= par->ndihedralPars + 1) < par->nimproperPars) {
		*type = kImproperParType;
		return (UnionPar *)(row >= 0 ? par->improperPars + row : NULL);
	} else if ((row -= par->nimproperPars + 1) < par->nvdwpPars) {
		*type = kVdwPairParType;
		return (UnionPar *)(row >= 0 ? par->vdwpPars + row : NULL);
/*	} else if ((row -= par->nvdwpPars + 1) < par->natomPars) {
		*type = kElementParType;
		return (UnionPar *)(row >= 0 ? par->atomPars + row : NULL); */
	} else {
		*type = kInvalidParType;
		return NULL;
	}
}

void
ParameterTableGetItemText(Parameter *par, int column, int row, char *buf, int bufsize)
{
	static char *sBondParTitles[] = {"", "Bonds", "k", "r0"};
	static char *sAngleParTitles[] = {"", "Angles", "k", "a0"};
	static char *sDihedralParTitles[] = {"", "Dihedrals", "k", "period", "phi0"};
	static char *sImproperParTitles[] = {"", "Impropers", "k", "period", "phi0"};
	static char *sVdwParTitles[] = {"", "VDWs", "eps", "r", "eps14", "r14", "atomNo", "weight"};
	static char *sVdwPairParTitles[] = {"", "VDW Pairs", "eps", "r", "eps14", "r14"};
/*	static char *sAtomParTitles[] = {"", "Atom Display", "atomNo", "radius", "red", "green", "blue", "weight"}; */
	const char *p;
	int type;
	UnionPar *up = NULL;
	char types[4][8];
	buf[0] = 0;
	if (par == NULL || row < 0)
		return;
	up = ParameterTableGetItemPtr(par, row, &type);
	switch (type) {
		case kVdwParType: {
			VdwPar *vp = (VdwPar *)up;
			if (vp == NULL) {
				if (column >= 0 && column < 8)
					snprintf(buf, bufsize, "%s", sVdwParTitles[column]);
				return;
			}
			switch (column) {
				case 0: snprintf(buf, bufsize, "vdw"); break;
				case 1:
					AtomTypeDecodeToString(vp->type1, types[0]);
					snprintf(buf, bufsize, "%s", types[0]);
					break;
				case 2:
					snprintf(buf, bufsize, "%.5f", vp->eps * INTERNAL2KCAL);
					break;
				case 3:
					snprintf(buf, bufsize, "%.5f", vp->r_eq);
					break;
				case 4:
					snprintf(buf, bufsize, "%.5f", vp->eps14 * INTERNAL2KCAL);
					break;
				case 5:
					snprintf(buf, bufsize, "%.5f", vp->r_eq14);
					break;
				case 6:
					snprintf(buf, bufsize, "%d", vp->atomicNumber);
					break;
				case 7:
					snprintf(buf, bufsize, "%.3f", vp->weight);
					break;
			}
			break;
		}
		case kBondParType: {
			BondPar *bp = (BondPar *)up;
			if (bp == NULL) {
				if (column >= 0 && column < 4)
					snprintf(buf, bufsize, "%s", sBondParTitles[column]);
				return;
			}
			switch (column) {
				case 0: snprintf(buf, bufsize, "bond"); break;
				case 1:
					AtomTypeDecodeToString(bp->type1, types[0]);
					AtomTypeDecodeToString(bp->type2, types[1]);
					snprintf(buf, bufsize, "%s-%s", types[0], types[1]);
					break;
				case 2:
				case 3:
					snprintf(buf, bufsize, "%.3f", (column == 2 ? bp->k * INTERNAL2KCAL : bp->r0));
					break;
			}
			break;
		}
		case kAngleParType: {
			AnglePar *ap = (AnglePar *)up;
			if (ap == NULL) {
				if (column >= 0 && column < 4)
					snprintf(buf, bufsize, "%s", sAngleParTitles[column]);
				return;
			}
			switch (column) {
				case 0: snprintf(buf, bufsize, "angle"); break;
				case 1:
					AtomTypeDecodeToString(ap->type1, types[0]);
					AtomTypeDecodeToString(ap->type2, types[1]);
					AtomTypeDecodeToString(ap->type3, types[2]);
					snprintf(buf, bufsize, "%s-%s-%s", types[0], types[1], types[2]);
					break;
				case 2:
				case 3:
					snprintf(buf, bufsize, "%.3f", (column == 2 ? ap->k * INTERNAL2KCAL : ap->a0 * kRad2Deg));
					break;
			}
			break;
		}
		case kDihedralParType: {
			TorsionPar *tp = (TorsionPar *)up;
			if (tp == NULL) {
				if (column >= 0 && column < 5)
					snprintf(buf, bufsize, "%s", sDihedralParTitles[column]);
				return;
			}
			switch (column) {
				case 0: snprintf(buf, bufsize, "dihe"); break;
				case 1:
					AtomTypeDecodeToString(tp->type1, types[0]);
					AtomTypeDecodeToString(tp->type2, types[1]);
					AtomTypeDecodeToString(tp->type3, types[2]);
					AtomTypeDecodeToString(tp->type4, types[3]);
					snprintf(buf, bufsize, "%s-%s-%s-%s", types[0], types[1], types[2], types[3]);
					break;
				case 3:
					snprintf(buf, bufsize, "%d", tp->period[0]);
					break;
				case 2:
				case 4:
					snprintf(buf, bufsize, "%.3f", (column == 2 ? tp->k[0] * INTERNAL2KCAL : tp->phi0[0] * kRad2Deg));
					break;
			}
			break;
		}
		case kImproperParType: {
			TorsionPar *tp = (TorsionPar *)up;
			if (tp == NULL) {
				if (column >= 0 && column < 5)
					snprintf(buf, bufsize, "%s", sImproperParTitles[column]);
				return;
			}		
			switch (column) {
				case 0: snprintf(buf, bufsize, "impr"); break;
				case 1:
					AtomTypeDecodeToString(tp->type1, types[0]);
					AtomTypeDecodeToString(tp->type2, types[1]);
					AtomTypeDecodeToString(tp->type3, types[2]);
					AtomTypeDecodeToString(tp->type4, types[3]);
					snprintf(buf, bufsize, "%s-%s-%s-%s", types[0], types[1], types[2], types[3]);
					break;
				case 3:
					snprintf(buf, bufsize, "%d", tp->period[0]);
					break;
				case 2:
				case 4:
					snprintf(buf, bufsize, "%.3f", (column == 2 ? tp->k[0] * INTERNAL2KCAL : tp->phi0[0] * kRad2Deg));
					break;
			}
			break;
		}
		case kVdwPairParType: {
			VdwPairPar *vp = (VdwPairPar *)up;
			if (vp == NULL) {
				if (column >= 0 && column < 6)
					snprintf(buf, bufsize, "%s", sVdwPairParTitles[column]);
				return;
			}
			switch (column) {
				case 0: snprintf(buf, bufsize, "pvdw"); break;
				case 1:
					AtomTypeDecodeToString(vp->type1, types[0]);
					AtomTypeDecodeToString(vp->type2, types[1]);
					snprintf(buf, bufsize, "%s-%s", types[0], types[1]);
					break;
				case 2:
					snprintf(buf, bufsize, "%.6f", vp->eps * INTERNAL2KCAL);
					break;
				case 3:
					snprintf(buf, bufsize, "%.6f", vp->r_eq);
					break;
				case 4:
					snprintf(buf, bufsize, "%.6f", (vp->A14 == 0.0 ? 0.0 : vp->B14 * vp->B14 / vp->A14 * 0.25 * INTERNAL2KCAL));
					break;
				case 5:
					snprintf(buf, bufsize, "%.6f", vp->eps14 * INTERNAL2KCAL);
					break;
			}
			break;
		}
/*		case kElementParType: {
			ElementPar *ap = (ElementPar *)up;
			if (ap == NULL) {
				if (column >= 0 && column < 8)
					snprintf(buf, bufsize, "%s", sAtomParTitles[column]);
				return;
			}
			switch (column) {
				case 0: snprintf(buf, bufsize, "disp"); break;
				case 1: snprintf(buf, bufsize, "%.4s", ap->name); break;
				case 2: snprintf(buf, bufsize, "%d", ap->number); break;
				case 3: snprintf(buf, bufsize, "%.2f", ap->radius); break;
				case 4: snprintf(buf, bufsize, "%.3f", ap->r); break;
				case 5: snprintf(buf, bufsize, "%.3f", ap->g); break;
				case 6: snprintf(buf, bufsize, "%.3f", ap->b); break;
				case 7: snprintf(buf, bufsize, "%.3f", ap->weight); break;
			}
			break;
		} */
		default: return;
	}
	if (up != NULL && (column == 8 || column == 9)) {
		if (column == 8 && ((BondPar *)up)->src == -1)
			snprintf(buf, bufsize, "!NONE!");
		else if ((p = ParameterGetComment(column == 8 ? ((BondPar *)up)->src : ((BondPar *)up)->com)) != NULL)
			snprintf(buf, bufsize, "%s", p);
	}
}

/*  Return values:
    -3: invalid
    -2: separator item
    -1: missing parameter
     0: molecule-local parameter
     1 and larger: global parameter values (gGlobalParInfo index + 1)  */
int
ParameterTableGetItemSource(Parameter *par, int row)
{
	int src, type;
	UnionPar *up;
	if (par == NULL || row < 0)
		return -3;
	up = ParameterTableGetItemPtr(par, row, &type);
	src = (type == kInvalidParType ? -3 : (up == NULL ? -2 : ((BondPar *)up)->src));
	if (type == kInvalidParType)
		src = -3;
	else if (up == NULL)
		src = -2;
	else src = ((BondPar *)up)->src;
	if (src > 0) {
		/*  Search src in gGlobalParInfo  */
		int i;
		for (i = 0; i < gGlobalParInfo.count; i++) {
			if (gGlobalParInfo.files[i].src == src)
				return i + 1;
		}
		return -3;  /*  Must not happen  */
	}
	return src;
}

int
ParameterTableIsItemEditable(Parameter *par, int column, int row)
{
	UnionPar *up;
	int type, f;
	if (par == NULL || row < 0)
		return 0;
	up = ParameterTableGetItemPtr(par, row, &type);
	if (type != kInvalidParType && up != NULL) {
		/*  Valid type, not separator row, molecule-local value  */
		int src = ((BondPar *)up)->src;
		f = (src == 0 || src == -1);
	} else f = 0;
	switch (type) {
		case kVdwParType: return (f && column > 0 && column < 8);
		case kBondParType: return (f && column > 0 && column < 4);
		case kAngleParType: return (f && column > 0 && column < 4);
		case kDihedralParType: return (f && column > 0 && column < 5);
		case kImproperParType: return (f && column > 0 && column < 5);
		case kVdwPairParType: return (f && column > 0 && column < 5);
	/*	case kElementParType: return (f && column > 0 && column < 7); */
		default: return 0;
	}
}

#pragma mark ====== Utility Functions ======

int
ElementParameterInitialize(const char *fname, char **outWarningMessage)
{
	char buf[1024], name[6], fullname[16];
	float val[6];
	FILE *fp = NULL;
	int i, lineNumber, retval = 0;
	char *wbuf = NULL;

	fp = fopen(fname, "rb");
	if (fp == NULL) {
		retval = -1;
		goto exit;
	}
	lineNumber = 0;
	while (ReadLine(buf, sizeof buf, fp, &lineNumber) > 0) {
		ElementPar *ep;
		if (strncmp(buf, "element ", 8) != 0)
			continue;  /*  Skip non-relevant lines  */
		fullname[0] = 0;
		if (sscanf(buf + 8, " %4s %f %f %f %f %f %f %15s", name, &val[0], &val[1], &val[2], &val[3], &val[4], &val[5], fullname) < 7) {
			asprintf(&wbuf, "%s:%d: missing parameter in ELEMENT record", fname, lineNumber);
			retval = 1;
			goto exit;
		}
		i = (int)val[0];
		if (i < 0 || i >= 200) {
			asprintf(&wbuf, "%s:%d: The atomic number (%d) in ELEMENT record is out of range", fname, lineNumber, i);
			retval = 2;
			goto exit;
		}
		ep = AssignArray(&gElementParameters, &gCountElementParameters, sizeof(ElementPar), i, NULL);
		memmove(ep->name, name, 4);
		ep->number = i;
		ep->radius = val[1];
		ep->r = val[2];
		ep->g = val[3];
		ep->b = val[4];
		ep->weight = val[5];
		fullname[15] = 0;
		memmove(ep->fullname, fullname, 16);
	}
exit:
	if (fp != NULL)
		fclose(fp);
	if (outWarningMessage != NULL)
		*outWarningMessage = wbuf;
	return retval;
}

UInt
AtomTypeEncodeToUInt(const char *s)
{
	Int i;
	UInt n, t;
	if ((s[0] == 'x' || s[0] == 'X') && s[1] == 0)
		return kAtomTypeWildcard;
	if (s[0] >= '0' && s[0] <= '9')
		return atoi(s);
	for (i = t = 0; i < 4; i++, s++) {
		/*  Encode: variant*96*96*96*96 + a[0]*96*96*96 + a[1]*96*96 + a[2] * 96 + a[3]  */
		static const UInt s_coeff[4] = {96*96*96, 96*96, 96, 1};
		while (*s == ' ')
			s++;
		if (*s == 0)
			break;
		if (*s == ',' || *s == '.' || *s == ';' || *s == ':') {
			/*  Variant (only [0-9a-z] are allowed) */
			s++;
			if (*s >= '0' && *s <= '9')
				n = *s - '0' + 1;
			else if (*s >= 'A' && *s <= 'Z')
				n = *s - 'A' + 11;
			else if (*s >= 'a' && *s <= 'z')
				n = *s - 'a' + 11;
			else n = (*s % 36) + 1;  /*  Map to something non-zero  */
			t += n * (96*96*96*96);
			break;
		}
		n = (*s - 0x20) % 96;  /*  Map to 1..95  */
		if (i == 0 && n < 32)
			n = 32;
		t += n * s_coeff[i];
	}
	return t;
}

char *
AtomTypeDecodeToString(UInt type, char *s)
{
	static const UInt s_coeff[4] = {96*96*96, 96*96, 96, 1};
	static char buf[8];
	int i;
	UInt variant, n;
	if (s == NULL) {
		s = buf;
		buf[6] = 0;
	}
	if (type == kAtomTypeWildcard) {
		s[0] = 'X';
		s[1] = 0;
		return s;
	}
	if (type < kAtomTypeMinimum) {
		snprintf(s, 6, "%d", type);
		return s;
	}
	for (i = 0; i < 4; i++) {
		s[i] = (type / s_coeff[i]) % 96;
		if (s[i] != 0)
			s[i] += 0x20;
	}
	s[4] = 0;
	n = strlen(s);
	if ((variant = (type / (96*96*96*96)) % 96) != 0) {
		s[n] = '.';
		s[n + 1] = (variant <= 10 ? '0' + variant - 1 : 'a' + variant - 11);
		s[n + 2] = 0;
	}
	return s;
}

Int
ElementToInt(const char *s)
{
	int i;
	ElementPar *p;
	for (i = 0, p = gElementParameters; i < gCountElementParameters; i++, p++) {
		if (p->name[0] == toupper(s[0]) && p->name[1] == tolower(s[1]))
			return p->number;
	}
	return 0;
}

char *
ElementToString(Int elem, char *s)
{
	if (elem >= 0 && elem < gCountElementParameters) {
		const char *cs = gElementParameters[elem].name;
		s[0] = cs[0];
		s[1] = cs[1];
		s[2] = 0;
		return s;
	} else return NULL;
}

Int
AtomNameToElement(const char *s)
{
	char element[4];
	const char *p;
	/*  $element = ($name =~ /([A-Za-z]{1,2})/); # in Perl  */
	element[0] = 0;
	for (p = s; *p != 0; p++) {
		if (isalpha(*p) && *p != '_') {
			element[0] = toupper(*p);
			if (isalpha(p[1]) && p[1] != '_') {
				element[1] = toupper(p[1]);
				element[2] = 0;
			} else {
				element[1] = 0;
			}
			break;
		}
	}
	return ElementToInt(element);
}

Int
GuessAtomicNumber(const char *name, Double weight)
{
	int i;
	ElementPar *p;
	char buf[4];
	const char *cp;
	for (i = 0, p = gElementParameters; i < gCountElementParameters; i++, p++) {
		if (p->weight > 0.0 && fabs(weight - p->weight) < 0.1)
			return p->number;
	}
	buf[0] = 0;
	for (cp = name; *cp != 0 && cp < name + 4; cp++) {
		if (isalpha(*cp) && *cp != '_') {
			buf[0] = toupper(*cp);
			if (isalpha(cp[1]) && cp[1] != '_') {
				buf[1] = toupper(cp[1]);
				buf[2] = 0;
			} else {
				buf[1] = 0;
			}
			break;
		}
	}
	return ElementToInt(buf);
}

Double
WeightForAtomicNumber(Int elem)
{
	if (elem >= 1 && elem < gCountElementParameters)
		return gElementParameters[elem].weight;
	else return 0.0;
}

Double
RadiusForAtomicNumber(Int elem)
{
	if (elem >= 1 && elem < gCountElementParameters)
		return gElementParameters[elem].radius;
	else return 0.0;
}

