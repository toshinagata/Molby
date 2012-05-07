/*
 *  Parameter.h
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

#ifndef __Parameter_h__
#define __Parameter_h__

#include "Types.h"
#include "Object.h"
#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif

/*  Information record for global parameter files  */
typedef struct ParFileInfo {
	char *name;   /* basename */
	char *dir;    /* directory */
	int src;      /* comment index */
} ParFileInfo;

typedef struct GlobalParInfoRecord {
	Int count;
	Int builtinCount;  /*  Number of 'built-in' global parameter files  */
	ParFileInfo *files;
} GlobalParInfoRecord;

extern struct GlobalParInfoRecord gGlobalParInfo;
	
/*  Parameter types  */
enum {
	kInvalidParType,
	kBondParType, kAngleParType, kDihedralParType, kImproperParType,
	kVdwParType, kVdwPairParType, kVdwCutoffParType,
	kFirstParType = kBondParType,
	kLastParType = kVdwCutoffParType,
	kElementParType = 100  /*  Only used in rb_cParameterRef  */
};
	
/*  Parameter Lookup options  */
enum {
	kParameterLookupNoWildcard = 1,
	kParameterLookupGlobal = 2,
	kParameterLookupLocal = 4,
	kParameterLookupMissing = 8,
	kParameterLookupNoBaseAtomType = 16
};

/*  The atom type is a 32-bit signed integer; a "variant" number and 4 ascii characters are
    encoded as: variant*96*96*96*96 + a[0]*96*96*96 + a[1]*96*96 + a[2]*96 + a[3].
    variant is [0-9a-z] mapped to [1..36].
    a[n] (n = 0..3) is a printable ASCII character mapped to [1..95].
    The first character must be no smaller than '@'.  */

#define kAtomTypeMinimum 28311552  /* 32*96*96*96 */
#define kAtomTypeWildcard 0x7fffffff
#define kAtomTypeVariantBase 84934656  /* 96*96*96*96 */

	
/*  Parameters are stored in separate arrays for each type. However, it is sometimes more
    convenient when a single number can represent a specific parameter regardless of type.
    kParameterIndexOffset is used for that purpose. Thus, the i-th bond parameter corresponds
    to i, the i-th angle parameter to i + kParameterIndexOffset, the i-th dihedral parameter
    to i + 2 * kParameterIndexOffset, etc.  */
#define kParameterIndexOffset 10000000

/*  Parameters  */
typedef struct BondPar {
	Int    com, src;   /*  Index to the comment array  */
	UInt    type1, type2;
	Double k, r0;
} BondPar;

typedef struct AnglePar {
	Int    com, src;   /*  Index to the comment array  */
	UInt    type1, type2, type3;
	Double k, a0;
} AnglePar;

typedef struct TorsionPar {
	Int    com, src;   /*  Index to the comment array  */
	UInt    type1, type2, type3, type4;
	Int    mult;  /*  The multiple term for CHARMm parameter sets  */
	Double k[3];
	Int    period[3];
	Double phi0[3];
} TorsionPar;

typedef struct VdwPar {
	Int    com, src;   /*  Index to the comment array  */
	UInt    type1;
	Int    atomicNumber;
	Double A, B;  /* A = 4*(sigma**12)*eps = ((2*r_eq)**12)*eps, B = 4*(sigma**6)*eps = 2*((2*r_eq)**6)*eps, sigma = (2*r_eq)*((1/2)**(1/6)), where r_eq is a van der Waals radius (not interatomic distance at the lowest energy!)  */
	Double eps, r_eq;   /*  These are redundant, but useful from the user's point of view  */
	Double A14, B14;
	Double eps14, r_eq14;
	Double weight; /*  Atomic weight  */
	Double polarizability;  /*  Not used at present  */
} VdwPar;

typedef struct VdwPairPar {
	Int    com, src;   /*  Index to the comment array  */
	UInt    type1, type2;
	Double A, B;  /* A = 4*(sigma**12)*eps = ((2*r_eq)**12)*eps, B = 4*(sigma**6)*eps = 2*((2*r_eq)**6)*eps, sigma = (2*r_eq)*((1/2)**(1/6))  */
	Double eps, r_eq;   /*  These are redundant, but useful from the user's point of view  */
	Double A14, B14;
	Double eps14, r_eq14;
} VdwPairPar;

typedef struct VdwCutoffPar {
	Int    com, src;   /*  Index to the comment array  */
	UInt   type1, type2;
/*  The fields type and n1-n4 are now obsolete; atom index range is no longer supported  */
/*	Byte   type;          *//* 0: atom type specific, 1: atom number specific */
/*	UInt   n1, n2, n3, n4;*//*  atom type specific: n1 and n2 are atom types  */
	                        /*  atom number specific: (n1,n2) and (n3,n4) are atom ranges */
							/*  (n1<=n2, n3<=n4, including the bounds)  */
	Double cutoff;
} VdwCutoffPar;

/*  Display parameters (defined for elements)  */
typedef struct ElementPar {
	Int    com, src;   /*  Index to the comment array  */
	Int    number;   /*  Atomic number  */
	char   name[4];
	Double radius;
	Double r, g, b;  /*  Color: [0.0, 1.0] for each component  */
	Double weight;
	char   fullname[16];
} ElementPar;
	
typedef struct Parameter {
	Object base;
	Int    nbondPars;
	BondPar *bondPars;
	Int    nanglePars;
	AnglePar *anglePars;
	Int    ndihedralPars;
	TorsionPar *dihedralPars;
	Int    nimproperPars;
	TorsionPar *improperPars;
	Int    nvdwPars;
	VdwPar *vdwPars;
	Int    nvdwpPars;
	VdwPairPar *vdwpPars;
	Int    nvdwCutoffPars;
	VdwCutoffPar *vdwCutoffPars;
/*	Int    natomPars;
	ElementPar *atomPars;
	Int    ncomments;
	char **comments;  */
} Parameter;

/*  Parameter reference (for MolAction and Ruby types)  */
typedef struct ParameterRef {
	struct Molecule *mol;  /*  If non-NULL, points to mol->par; otherwise, points to gBuiltinParameters  */
	Int parType;
	Int idx;
} ParameterRef;

/*  For return value for parameter pointer  */
typedef union UnionPar {
	BondPar bond;
	AnglePar angle;
	TorsionPar torsion;
	VdwPar vdw;
	VdwPairPar vdwp;
	VdwCutoffPar vdwcutoff;
	ElementPar atom;
} UnionPar;

Parameter *ParameterNew(void);
Parameter *ParameterDuplicate(const Parameter *par);
Parameter *ParameterWithName(const char *name);
void ParameterSetName(Parameter *par, const char *name);
const char *ParameterGetName(Parameter *par);
void ParameterRetain(Parameter *par);
void ParameterRelease(Parameter *par);
int ParameterReadFromFile(Parameter *par, const char *fname, char **outWarningMessage, int *outWarningCount);
int ParameterReadFromString(Parameter *par, char *buf, char **wbufp, const char *fname, int lineNumber, int src_idx);
int ParameterAppendToFile(Parameter *par, FILE *fp);
int ParameterDeleteAllEntriesForSource(Parameter *par, int src_idx);

int ParameterGlobalParIndexForSrcIndex(int src);
int ParameterCommentIndexForGlobalFileName(const char *p);
int ParameterCommentIndex(const char *comment);
const char *ParameterGetComment(int idx);
	
UnionPar *ParameterGetUnionParFromTypeAndIndex(Parameter *par, int type, int index);
Int ParameterGetCountForType(Parameter *par, int type);
Int ParameterGetSizeForType(int type);
	
ParameterRef *ParameterRefNew(struct Molecule *mol, int type, int idx);
void ParameterRefRelease(ParameterRef *pref);
UnionPar *ParameterRefGetPar(ParameterRef *pref);
	
struct IntGroup;  /*  forward declaration  */
struct Atom;      /*  forward declaration  */
int ParameterInsert(Parameter *par, Int type, const UnionPar *up, struct IntGroup *where);
int ParameterDelete(Parameter *par, Int type, UnionPar *up, struct IntGroup *where);
int ParameterCopy(Parameter *par, Int type, UnionPar *up, struct IntGroup *where);

/*  Caution!  When up is a UnionPar pointer given by ParameterGetUnionParFromTypeAndIndex() etc and
    u is a UnionPar variable, u = *up will cause Bad Address exception, if up is at the last
    of the allocated array and sizeof(UnionPar) is larger than the size of actual parameter record.
    This copy function does take care of such case.  */
void ParameterCopyOneWithType(UnionPar *dst, const UnionPar *src, int type);

int ParameterGetAtomTypes(Int type, const UnionPar *up, UInt *outTypes);
int ParameterRenumberAtoms(Int type, UnionPar *up, Int oldnatoms, const Int *old2new);
int ParameterDoesContainAtom(Int type, const UnionPar *up, UInt atom_type, Int options);
int ParameterIsRelevantToAtomGroup(Int type, const UnionPar *up, const struct Atom *ap, struct IntGroup *ig);

BondPar *ParameterLookupBondPar(Parameter *par, UInt t1, UInt t2, Int options);
AnglePar *ParameterLookupAnglePar(Parameter *par, UInt t1, UInt t2, UInt t3, Int options);
TorsionPar *ParameterLookupDihedralPar(Parameter *par, UInt t1, UInt t2, UInt t3, UInt t4, Int options);
TorsionPar *ParameterLookupImproperPar(Parameter *par, UInt t1, UInt t2, UInt t3, UInt t4, Int options);
VdwPar *ParameterLookupVdwPar(Parameter *par, UInt t1, Int options);
VdwPairPar *ParameterLookupVdwPairPar(Parameter *par, UInt t1, UInt t2, Int options);
VdwCutoffPar *ParameterLookupVdwCutoffPar(Parameter *par, UInt t1, UInt t2, Int options);

int ElementParameterInitialize(const char *fname, char **outWarningMessage);
UInt AtomTypeEncodeToUInt(const char *s);
char *AtomTypeDecodeToString(UInt type, char *s);
Int ElementToInt(const char *s);
char *ElementToString(Int type, char *s);
Int AtomNameToElement(const char *s);

Int GuessAtomicNumber(const char *name, Double weight);
Double WeightForAtomicNumber(Int elem);
Double RadiusForAtomicNumber(Int elem);
	
int ParameterTableNumberOfRows(Parameter *par);
int ParameterTableGetItemIndex(Parameter *par, int row, int *type);
int ParameterTableGetRowFromTypeAndIndex(Parameter *par, int type, int idx);
UnionPar *ParameterTableGetItemPtr(Parameter *par, int row, int *type);
void ParameterTableGetItemText(Parameter *par, int column, int row, char *buf, int bufsize);
int ParameterTableGetItemSource(Parameter *par, int row);
int ParameterTableIsItemEditable(Parameter *par, int column, int row);

extern Parameter *gBuiltinParameters;
extern ElementPar *gElementParameters;
extern Int gCountElementParameters;
	
#ifdef __cplusplus
}
#endif
		
#endif /* __Parameter_h__ */
