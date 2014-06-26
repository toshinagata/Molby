/*
 *  Molecule.h
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

#ifndef __Molecule_h__
#define __Molecule_h__

#include "Types.h"
#include "Object.h"
#include "IntGroup.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif
	
#define ATOMS_MAX_SYMMETRY 12
#define ATOMS_MAX_NUMBER 100000000  /*  Sufficiently large value  */

/*  Conversion between kcal/mol and internal energy unit (am*ang^2/fs^2, am = atomic mass) */
#define KCAL2INTERNAL (4.184e-4)
#define INTERNAL2KCAL (1.0/KCAL2INTERNAL)
#define J2INTERNAL (1e-4)
#define INTERNAL2J (1.0/J2INTERNAL)
	
#define BOLTZMANN (8.31441e-3*J2INTERNAL)
#define PI 3.14159265358979
#define PI2R 0.564189583547756    /*  1.0/sqrt(PI)  */
	
/*  Anisotropic thermal parameter  */
typedef struct Aniso {
	Double  bij[6];    /*  b11, b22, b33, b12, b13, b23 (ORTEP type 0) */
	char has_bsig;     /*  Has sigma values?  */
	Double  bsig[6];   /*  sigma values  */
	Mat33  pmat;      /*  A 3x3 matrix whose three column vectors are the principal axes of the ellipsoid. Note: If the B matrix is not positive definite, the axis length corresponding to the negative eigenvalue is replaced with 0.001.  */
} Aniso;

/*  Symmetry operation  */
/*  If periodic box is defined, dx/dy/dz denote multiples of the axes of the periodic box.
    Otherwise, dx/dy/dz denote offset to the x/y/z coordinates of the atoms.  */
typedef struct Symop {
	signed int dx : 4;
	signed int dy : 4;
	signed int dz : 4;
	unsigned int sym : 8;
	unsigned int alive: 1;
} Symop;

/*  Exflags  */
enum {
	kAtomHiddenFlag = 1
};

/*  Atom connection record  */
/*  If nconnects <= ATOM_CONNECT_LIMIT, data[] field is used. Otherwise,
    memory is allocated by malloc().  */
#define ATOM_CONNECT_LIMIT 6
typedef struct AtomConnect {
	Int    count;  /*  Number of connections  */
	union {
		Int *ptr;
		Int data[ATOM_CONNECT_LIMIT];
	} u;
} AtomConnect;

typedef struct PiAnchor {
	AtomConnect connect;
	Int ncoeffs;
	Double *coeffs;
} PiAnchor;

/*  Atom record  */
typedef struct Atom {
	Int    segSeq;
	char   segName[4];
	Int    resSeq;
	char   resName[4];
	char   aname[4];
	UInt   type;
	Double  charge;
	Double  weight;
	char   element[4];
	Int    atomicNumber;
	AtomConnect connect;
	Vector r;  /*  position  */
	Vector v;  /*  velocity  */
	Vector f;  /*  force  */
	Vector sigma;   /*  For crystallographic data only; sigma for each crystallographic coordinates  */
					/*  (Unlike r, these are not converted to the cartesian system)  */
	Double  occupancy;
	Double  tempFactor;
	Aniso  *aniso;
	Int    intCharge;
	Int    exflags;
	Int    nframes;  /*  Multiple frames  */
	Vector *frames;
	Symop  symop;    /*  For symmetry-expanded atom  */
	Int    symbase;  /*  The index of original atom for symmetry-expansion  */
	PiAnchor *anchor;  /*  Non-NULL if this atom is a pi-anchor  */
	Int    labelid;  /*  The label ID; 0 for no label  */
	short  wrap_dx, wrap_dy, wrap_dz; /*  Calculated by md_wrap_coordinates; used only in wrapped output.  */
	Double fix_force; /*  0: no fix, >0: fix at fix_pos with harmonic potential, <0: fix at fix_pos without force  */
	Vector fix_pos;
	Byte   mm_exclude;        /*  If nonzero, then this atom is excluded from MM/MD calculations  */
	Byte   periodic_exclude;  /*  If nonzero, then this atom is excluded from periodic calculations  */
	char   uff_type[6]; /*  UFF type string  */
} Atom;

extern Int gSizeOfAtomRecord;

#define ATOM_AT_INDEX(p, i)  ((Atom *)((char *)(p) + (i) * gSizeOfAtomRecord))
#define ATOM_NEXT(p)         ((Atom *)((char *)(p) + gSizeOfAtomRecord))
#define ATOM_PREV(p)         ((Atom *)((char *)(p) - gSizeOfAtomRecord))
#define SYMOP_ALIVE(s) ((s.dx || s.dy || s.dz || s.sym) != 0)
#define SYMOP_EQUAL(s1, s2) (s1.dx == s2.dx && s1.dy == s2.dy && s1.dz == s2.dz && s1.sym == s2.sym)
#define SYMMETRY_AT_INDEX(p, i) (*((i) == 0 ? &gIdentityTransform : &p[i]))

/*  atom.connects is a union entry, including direct data for nconnects <= ATOM_CONNECT_LIMIT
    and malloc()'ed entry for nconnects > ATOM_CONNECT_LIMIT. The following functions
	automatically take care of the memory allocation/deallocation.  */
Int *AtomConnectData(AtomConnect *ac);
void AtomConnectResize(AtomConnect *ac, Int nconnects);
void AtomConnectInsertEntry(AtomConnect *ac, Int idx, Int connect);
void AtomConnectDeleteEntry(AtomConnect *ac, Int idx);

#define ATOM_CONNECT_PTR(ac) ((ac)->count > ATOM_CONNECT_LIMIT ? (ac)->u.ptr : (ac)->u.data)

/*  Duplicate an atom. If dst is non-NULL, *src is copied to *dst and dst is returned. If dst is NULL, a new atom is allocated by malloc() and that atom is returned. It is the called's responsibility to release the returned memory. */
extern Atom *AtomDuplicate(Atom *dst, const Atom *src);
	
/*  Duplicate an atom, except for the frame entry  */
extern Atom *AtomDuplicateNoFrame(Atom *dst, const Atom *src);
	
/*  Clean the content of an atom record  */
extern void AtomClean(Atom *ap);

/*  MolEnumerable type code  */
enum {
	kAtomKind = 0,
	kBondKind = 1,
	kAngleKind = 2,
	kDihedralKind = 3,
	kImproperKind = 4,
	kResidueKind = 5,
	kEndKind
};

/*  Enumerable class to access to atoms, bonds, etc.  */
typedef struct MolEnumerable {
	struct Molecule *mol;
	int    kind;
} MolEnumerable;

/*  Atom reference  */
typedef struct AtomRef {
	struct Molecule *mol;
	int idx;
} AtomRef;

/*  Crystallographic cell parameter (also used as periodic box in MD) */
typedef struct XtalCell {
	Double  cell[6];     /*  a, b, c, alpha, beta, gamma (in degree)  */
	Double  rcell[6];    /*  Reciprocal cell  */
	Vector  axes[3];     /*  Cartesian unit vectors along the three axis  */
	Vector  origin;      /*  Cartesian origin of the periodic box  */
	char    flags[3];    /*  1 for periodic, 0 for non-periodic  */
	char    has_sigma;   /*  Has sigma?  */
	Transform tr;        /*  Crystal coord -> cartesian  */
	Transform rtr;       /*  Cartesian -> crystal coord  */
	Double  cellsigma[6];  /*  For crystallographic data; sigma for the cell parameters  */
} XtalCell;

/*  3-Dimensional distribution  */
typedef struct Cube {
	Int idn;             /*  Integer identifier (such as MO number)  */
	Vector origin;
	Vector dx, dy, dz;
	Int nx, ny, nz;
	Double *dp;          /*  Value for point (ix, iy, iz) is in dp[(ix*ny+iy)*nz+iz]  */
} Cube;

/*  Gaussian orbital symmetry types  */
enum {
	kGTOType_S,
	kGTOType_SP,
	kGTOType_P,
	kGTOType_D,
	kGTOType_D5,
	kGTOType_F,
	kGTOType_F7,
	kGTOType_G,
	kGTOType_G9,
	kGTOType_UU
};

/*  Exponent/coefficient info for a single gaussian primitive  */
typedef struct PrimInfo {
	Double A;            /*  Exponent  */
	Double C;            /*  Contraction coefficient  */
	Double Csp;          /*  P(S=P) contraction coefficient  */
} PrimInfo;

/*  Gaussian orbital shell information  */
typedef struct ShellInfo {
	signed char sym;     /*  Symmetry of the basis; S, P, ... */
	signed char ncomp;   /*  Number of components (S: 1, P: 3, SP: 4, etc.)  */
	signed char nprim;   /*  Number of primitives for this shell  */
	Int p_idx;           /*  Index to the PrimInfo (exponent/coefficient) table  */
	Int cn_idx;          /*  Index to the normalized (cached) contraction coefficient table  */
	Int a_idx;           /*  Index to the atom which this primitive belongs to */
	Int m_idx;           /*  Index to the MO matrix  */
} ShellInfo;

/*  Basis set and MO information  */
typedef struct BasisSet {
	Int nshells;         /*  Number of gaussian orbital shells  */
	ShellInfo *shells;   /*  Gaussian orbital shells  */
	Int npriminfos;      /*  Size of primitive information table  */
	PrimInfo *priminfos; /*  Primitive information table  */
	Int ncns;            /*  Number of normalized (cached) contraction coefficient values  */
	Double *cns;         /*  Normalized (cached) contraction coefficients; (up to 10 values for each primitive)  */
	Int natoms_bs;       /*  Number of atoms; separately cached here because MO info should be invariant during editing */
	Double *nuccharges;  /*  Nuclear charges (for ECP atoms)  */
	Int ne_alpha, ne_beta;  /*  Number of alpha/beta electrons  */
	Int rflag;           /*  0: UHF, 1: RHF, 2:ROHF  */
	Int ncomps;          /*  Number of AO components; equal to sum of shells[i].ncomp  */
	Int nmos;            /*  Number of MOs; equal to ncomps if close shell, ncomps*2 if open shell */
	Double *mo;          /*  MO matrix (mo[i][j] represents the j-th AO coefficient for the i-th MO)  */
						 /*  Memory are allocated for (2*nmos+1) entries; the last entry is for displaying arbitrary vector  */
	Double *moenergies;  /*  MO energies  */
	Double *scfdensities; /*  SCF densities; lower triangle of a symmetric matrix (size nmos*(nmos+1)/2)  */
	Int ncubes;          /*  Number of calculated MOs  */
	Cube **cubes;        /*  Calculated MOs (an array of pointers to Cubes)  */
} BasisSet;

/*  Marching Cube (for drawing isosurface)  */
typedef struct MCubePoint {
	Int key;       /*  key = ((ix*ny+iy)*nz+iz)*3+ii, ii=0/1/2 for x/y/z direction, respectively */
	float d;       /*  offset toward the direction; 0 <= d < 1  */
	float pos[3];  /*  cartesian coordinate of the point  */
	float grad[3]; /*  gradient vector  */
} MCubePoint;
	
typedef struct MCube {
	Int idn;             /*  MO number  */
	Vector origin;       /*  Cube origin */
	Double dx, dy, dz;   /*  Cube steps */
	Int nx, ny, nz;      /*  Cube dimension (must be multiples of 8)  */
	Double thres;        /*  Threshold value  */
	Double *dp;          /*  Value for point (ix, iy, iz) is in dp[(ix*ny+iy)*nz+iz]  */
	Int nradii;
	Double *radii;       /*  Estimated radius (with margin) for each atom  */
	Double expand;       /*  Expand the estimated radius by this value (default: 1.0)  */
	struct {
		/*  Flags for cube (ix, iy, iz)-(ix+1, iy+1, iz+1). It is an 8-bit */
		/*  integer representing whether the values at the 8 corners are */
		/*  larger than the threshold value or not. As special cases,  */
		/*  the values 0 and 255 (all corners are below or above the threshold) */
		/*  are represented as 255, and the value 0 is used to indicate "yet undefined". */
		unsigned char *fp;
		/*  Cube points and triangles: for positive and negative surfaces  */
		Int ncubepoints;
		MCubePoint *cubepoints;
		Int ntriangles;
		Int *triangles;  /*  Triangles; indices to cubepoints[]  */
		float rgba[4];   /*  Surface color  */
	} c[2];
} MCube;

/*  Electrostatic potential  */
typedef struct Elpot {
	Vector pos;
	Double esp;
} Elpot;

/*  Properties (total energy etc.; specific for each frame)  */
typedef struct MolProp {
	char *propname;
	Double *propvals;
} MolProp;
	
/*  Molecule record  */
typedef struct Molecule {
	Object base;
	Int    natoms;
	Atom   *atoms;
	Int    nbonds;
	Int    *bonds;       /*  The size of array is 2*nbonds  */
	Int    nangles;
	Int    *angles;      /*  The size of array is 3*nangles  */
	Int    ndihedrals;
	Int    *dihedrals;   /*  The size of array is 4*ndihedrals  */
	Int    nimpropers;
	Int    *impropers;   /*  The size of array is 4*nimpropers  */
	Int    nresidues;    /*  Number of residues; maximum residue number + 1 (because residue 0 is 'undefined residue')  */
	char   (*residues)[4];
	XtalCell   *cell;
	Int    nsyms;        /*  Symmetry operations; syms are always described in crystallographic units (even when the unit cell is not defined)  */
	Transform *syms;

	IntGroup *selection;
	Int    nframes;      /*  The number of frames (>= 1). This is a cached value, and should be
							 recalculated from the atoms if it is -1  */
	Int    cframe;       /*  The current frame number  */

	Int    nframe_cells;
	Vector *frame_cells; /*  The cell vectors for frames; (nframe_cells*4) array of Vectors  */

	struct MainView *mview;  /*  Reference to the MainView object if present (no retain)  */
	Int    modifyCount;  /*  Internal counter for modification. This value is not to be modified
	                         manually; instead, call MoleculeIncrementModifyCount() whenever
						     modification is done, which also takes care necessary notification
							 to the other part of the application (system dependent)  */

	struct MDArena *arena;  /*  Reference to the MDArena record during MM/MD run (no retain)  */

	const char *path;     /*  The full path of the molecule, when this molecule is last saved to/loaded from file. Only used in the command-line version. (In GUI version, the path is obtained by the Document mechanism) */

	/*  Information from the dcd files  */
	Int    startStep;     /*  the timestep for frame 0  */
	Int    stepsPerFrame; /*  the number of timesteps between neighboring frames  */
	Double psPerStep;     /*  picosecond per step  */

	/*  Information for basis sets and MOs  */
	BasisSet *bset;
	
	/*  Marching cube  */
	MCube *mcube;

	/*  Electrostatic potential  */
	Int    nelpots;
	Elpot  *elpots;

	/*  Properties  */
	Int    nmolprops;
	MolProp *molprops;

	/*  Parameters specific for this molecule  */
	struct Parameter *par;
	
	/*  Bond order (not to be used in MM/MD, but may be necessary to hold this info)  */
	Int    nbondOrders;
	Double *bondOrders;

	/*  Flag to request rebuilding MD internal information  */
	Byte   needsMDRebuild;
	
	/*  Flag to clear selection of the parameter table  */
	Byte   parameterTableSelectionNeedsClear;
	
	/*  Flag to request copying coordinates to MD arena  */
	Byte   needsMDCopyCoordinates;

	/*  Prohibit modification of the topology (to avoid interfering MD) */
	Byte   noModifyTopology;
	
	/*  Flag to request aborting a subthread  */
	Byte   requestAbortThread;

	/*  Flag to signal that a subthread is terminated  */
	Byte   threadTerminated;

	/*  Mutex object. If non-NULL, it should be locked before modifying molecule  */
	void *mutex;
	
	/*  Flag to prohibit modification from user interface  */
	Byte   dontModifyFromGUI;

	/*  Ruby pointer (see ruby_bind.c)  */
	void *exmolobj;
	Byte exmolobjProtected;

} Molecule;

int strlen_limit(const char *s, int limit);

void BasisSetRelease(BasisSet *bset);

Molecule *MoleculeNew(void);
int MoleculeLoadFile(Molecule *mp, const char *fname, const char *ftype, char **errbuf);
int MoleculeLoadPsfFile(Molecule *mp, const char *fname, char **errbuf);
int MoleculeLoadTepFile(Molecule *mp, const char *fname, char **errbuf);
int MoleculeLoadShelxFile(Molecule *mp, const char *fname, char **errbuf);
int MoleculeLoadGaussianFchkFile(Molecule *mp, const char *fname, char **errbuf);
int MoleculeLoadMbsfFile(Molecule *mp, const char *fname, char **errbuf);
Molecule *MoleculeNewWithName(const char *name);
Molecule *MoleculeInitWithAtoms(Molecule *mp, const Atom *atoms, int natoms);
Molecule *MoleculeInitWithMolecule(Molecule *mp2, Molecule *mp);
void MoleculeSetName(Molecule *par, const char *name);
const char *MoleculeGetName(Molecule *mp);
void MoleculeSetPath(Molecule *mol, const char *fname);
const char *MoleculeGetPath(Molecule *mol);
Molecule *MoleculeWithName(const char *name);
Molecule *MoleculeRetain(Molecule *mp);
void MoleculeRelease(Molecule *mp);
void MoleculeExchange(Molecule *mp1, Molecule *mp2);

int MoleculeAddGaussianOrbitalShell(Molecule *mol, Int a_idx, Int sym, Int nprims);
int MoleculeAddGaussianPrimitiveCoefficients(Molecule *mol, Double exponent, Double contraction, Double contraction_sp);
int MoleculeGetGaussianComponentInfo(Molecule *mol, Int comp_idx, Int *outAtomIdx, char *outLabel, Int *outShellIdx);
int MoleculeSetMOCoefficients(Molecule *mol, Int idx, Double energy, Int ncomps, Double *coeffs);
int MoleculeGetMOCoefficients(Molecule *mol, Int idx, Double *energy, Int *ncoeffs, Double **coeffs);
int MoleculeSetMOInfo(Molecule *mol, Int rflag, Int ne_alpha, Int ne_beta);

void MoleculeIncrementModifyCount(Molecule *mp);
void MoleculeClearModifyCount(Molecule *mp);

MolEnumerable *MolEnumerableNew(Molecule *mol, int kind);
void MolEnumerableRelease(MolEnumerable *mseq);
AtomRef *AtomRefNew(Molecule *mol, int idx);
void AtomRefRelease(AtomRef *aref);

void MoleculeSetCell(Molecule *mp, Double a, Double b, Double c, Double alpha, Double beta, Double gamma, int convertCoordinates);
void MoleculeSetAniso(Molecule *mp, int n1, int type, Double x11, Double x22, Double x33, Double x12, Double x13, Double x23, const Double *sigmap);
void MoleculeSetAnisoBySymop(Molecule *mp, int idx);
int MoleculeSetPeriodicBox(Molecule *mp, const Vector *ax, const Vector *ay, const Vector *az, const Vector *ao, const char *periodic, int convertCoordinates);
int MoleculeCalculateCellFromAxes(XtalCell *cp, int calc_abc);

int MoleculeReadCoordinatesFromFile(Molecule *mp, const char *fname, const char *ftype, char **errbuf);
int MoleculeReadCoordinatesFromPdbFile(Molecule *mp, const char *fname, char **errbuf);
int MoleculeReadCoordinatesFromDcdFile(Molecule *mp, const char *fname, char **errbuf);

int MoleculeLoadGamessDatFile(Molecule *mol, const char *fname, char **errbuf);

int MoleculeReadExtendedInfo(Molecule *mp, const char *fname, char **errbuf);
int MoleculeWriteExtendedInfo(Molecule *mp, const char *fname, char **errbuf);

int MoleculeWriteToFile(Molecule *mp, const char *fname, const char *ftype, char **errbuf);
int MoleculeWriteToPsfFile(Molecule *mp, const char *fname, char **errbuf);
int MoleculeWriteToPdbFile(Molecule *mp, const char *fname, char **errbuf);
int MoleculeWriteToDcdFile(Molecule *mp, const char *fname, char **errbuf);
int MoleculeWriteToTepFile(Molecule *mp, const char *fname, char **errbuf);
int MoleculeWriteToMbsfFile(Molecule *mp, const char *fname, char **errbuf);
void MoleculeDump(Molecule *mol);

int MoleculePrepareMDArena(Molecule *mol, int check_only, char **retmsg);

char *MoleculeSerialize(Molecule *mp, Int *outLength, Int *timep);
Molecule *MoleculeDeserialize(const char *data, Int length, Int *timep);

void MoleculeCleanUpResidueTable(Molecule *mp);
int MoleculeChangeNumberOfResidues(Molecule *mp, int nresidues);
int MoleculeChangeResidueNumberWithArray(Molecule *mp, IntGroup *group, Int *resSeqs);
int MoleculeChangeResidueNumber(Molecule *mp, IntGroup *group, int resSeq);
int MoleculeOffsetResidueNumbers(Molecule *mp, IntGroup *group, int offset, int nresidues);
int MoleculeChangeResidueNames(Molecule *mp, int argc, Int *resSeqs, char *names);
int MoleculeMaximumResidueNumber(Molecule *mp, IntGroup *group);
int MoleculeMinimumResidueNumber(Molecule *mp, IntGroup *group);

struct MolAction;
#if defined(DEBUG)
	int MoleculeCheckSanity(Molecule *mp);
#else
#define MoleculeCheckSanity(mp)
#endif

int MoleculeCreateAnAtom(Molecule *mp, const Atom *ap, int pos);
int MoleculeMerge(Molecule *dst, Molecule *src, IntGroup *where, int resSeqOffset, Int *nactions, struct MolAction ***actions, Int forUndo);
int MoleculeUnmerge(Molecule *src, Molecule **dstp, IntGroup *where, int resSeqOffset, Int *nactions, struct MolAction ***actions, Int forUndo);
int MoleculeExtract(Molecule *src, Molecule **dstp, IntGroup *where, int dummyFlag);
int MoleculeAddBonds(Molecule *mp, Int nbonds, const Int *bonds, IntGroup *where, Int autoGenerate);
int MoleculeDeleteBonds(Molecule *mp, Int *bonds, IntGroup *where, Int **outRemoved, IntGroup **outRemovedPos);
int MoleculeAssignBondOrders(Molecule *mp, const Double *orders, IntGroup *where);
int MoleculeGetBondOrders(Molecule *mp, Double *outOrders, IntGroup *where);
int MoleculeAddAngles(Molecule *mp, const Int *angles, IntGroup *where);
int MoleculeDeleteAngles(Molecule *mp, Int *angles, IntGroup *where);
int MoleculeAddDihedrals(Molecule *mp, const Int *dihedrals, IntGroup *where);
int MoleculeDeleteDihedrals(Molecule *mp, Int *dihedrals, IntGroup *where);
int MoleculeAddImpropers(Molecule *mp, const Int *impropers, IntGroup *where);
int MoleculeDeleteImpropers(Molecule *mp, Int *impropers, IntGroup *where);
int MoleculeLookupBond(Molecule *mp, Int n1, Int n2);
int MoleculeLookupAngle(Molecule *mp, Int n1, Int n2, Int n3);
int MoleculeLookupDihedral(Molecule *mp, Int n1, Int n2, Int n3, Int n4);
int MoleculeLookupImproper(Molecule *mp, Int n1, Int n2, Int n3, Int n4);

Int MoleculeFindMissingAngles(Molecule *mol, Int **outAngles);
Int MoleculeFindMissingDihedrals(Molecule *mol, Int **outDihedrals);
Int MoleculeFindMissingImpropers(Molecule *mol, Int **outImpropers);
	
IntGroup *MoleculeSearchBondsIncludingAtoms(Molecule *mp, IntGroup *atomgroup);
IntGroup *MoleculeSearchAnglesIncludingAtoms(Molecule *mp, IntGroup *atomgroup);
IntGroup *MoleculeSearchDihedralsIncludingAtoms(Molecule *mp, IntGroup *atomgroup);
IntGroup *MoleculeSearchImpropersIncludingAtoms(Molecule *mp, IntGroup *atomgroup);

IntGroup *MoleculeSearchBondsAcrossAtomGroup(Molecule *mp, IntGroup *atomgroup);

IntGroup *MoleculeSearchAnglesIncludingBond(Molecule *mp, int n1, int n2);
IntGroup *MoleculeSearchDihedralsIncludingBond(Molecule *mp, int n1, int n2);
IntGroup *MoleculeSearchImpropersIncludingBond(Molecule *mp, int n1, int n2);

int MoleculeLookupAtomInResidue(Molecule *mp, int n1, int resno);
int MoleculeAnalyzeAtomName(const char *s, char *resName, int *resSeq, char *atomName);
int MoleculeAtomIndexFromString(Molecule *mp, const char *s);

int MoleculeFindCloseAtoms(Molecule *mp, const Vector *vp, Double radius, Double limit, Int *outNbonds, Int **outBonds, Int triangle);
int MoleculeGuessBonds(Molecule *mp, Double limit, Int *outNbonds, Int **outBonds);
int MoleculeRebuildTablesFromConnects(Molecule *mp);
int MoleculeAreAtomsConnected(Molecule *mol, int idx1, int idx2);
	
void MoleculeGetAtomName(Molecule *mp, int index, char *buf, int bufsize);

void MoleculeSetSelection(Molecule *mp, IntGroup *select);
IntGroup *MoleculeGetSelection(Molecule *mp);
void MoleculeSelectAtom(Molecule *mp, int n1, int extending);
void MoleculeUnselectAtom(Molecule *mp, int n1);
void MoleculeToggleSelectionOfAtom(Molecule *mp, int n1);
int MoleculeIsAtomSelected(Molecule *mp, int n1);
int MoleculeIsBondSelected(Molecule *mp, int n1, int n2);
IntGroup *MoleculeModifySelectionByRemovingAtoms(Molecule *mp, IntGroup *selection, IntGroup *remove);

int MoleculeGetTransformForSymop(Molecule *mp, Symop symop, Transform *tf, int is_cartesian);
int MoleculeGetSymopForTransform(Molecule *mp, const Transform tf, Symop *symop, int is_cartesian);

int MoleculeTransformBySymop(Molecule *mp, const Vector *vpin, Vector *vpout, Symop symop);
int MoleculeAddExpandedAtoms(Molecule *mp, Symop symop, IntGroup *group, Int *indices, Int allowOverlap);
int MoleculeAmendBySymmetry(Molecule *mp, IntGroup *group, IntGroup **groupout, Vector **vpout);

int MoleculeShowAllAtoms(Molecule *mp);
int MoleculeShowReverse(Molecule *mp);
int MoleculeHideAtoms(Molecule *mp, IntGroup *ig);

int MoleculeRenumberAtoms(Molecule *mp, const Int *new2old, Int *old2new_out, Int isize);

void MoleculeTransform(Molecule *mp, Transform tr, IntGroup *group);
void MoleculeTranslate(Molecule *mp, const Vector *vp, IntGroup *group);
void MoleculeRotate(Molecule *mp, const Vector *axis, Double angle, const Vector *center, IntGroup *group);
int MoleculeCenterOfMass(Molecule *mp, Vector *center, IntGroup *group);
int MoleculeBounds(Molecule *mp, Vector *min, Vector *max, IntGroup *group);
	
Int *MoleculeSearchEquivalentAtoms(Molecule *mol, IntGroup *ig);
	
void MoleculeAddExpansion(Molecule *mp, Vector dr, Int symop, IntGroup *group, Double limit);
void MoleculeClearExpansion(Molecule *mp, IntGroup *group);
void MoleculeRemoveExpansion(Molecule *mp, Vector dr, Int symop, IntGroup *group);
void MoleculeAutoExpansion(Molecule *mp, const float *boxstart, const float *boxend, IntGroup *group, Double limit);

void MoleculeXtalToCartesian(Molecule *mp, Vector *dst, const Vector *src);
void MoleculeCartesianToXtal(Molecule *mp, Vector *dst, const Vector *src);
Double MoleculeMeasureBond(Molecule *mp, const Vector *vp1, const Vector *vp2);
Double MoleculeMeasureAngle(Molecule *mp, const Vector *vp1, const Vector *vp2, const Vector *vp3);
Double MoleculeMeasureDihedral(Molecule *mp, const Vector *vp1, const Vector *vp2, const Vector *vp3, const Vector *vp4);

IntGroup *MoleculeFragmentExcludingAtomGroup(Molecule *mp, int n1, IntGroup *exatoms);
IntGroup *MoleculeFragmentExcludingAtoms(Molecule *mp, int n1, int argc, int *argv);
IntGroup *MoleculeFragmentWithAtomGroups(Molecule *mp, IntGroup *inatoms, IntGroup *exatoms);
int MoleculeIsFragmentDetachable(Molecule *mp, IntGroup *group, int *n1, int *n2);
int MoleculeIsFragmentRotatable(Molecule *mp, IntGroup *group, int *n1, int *n2, IntGroup **rotGroup);

int MoleculeGetNumberOfFrames(Molecule *mp);
int MoleculeInsertFrames(Molecule *mp, IntGroup *group, const Vector *inFrame, const Vector *inFrameCell);
int MoleculeRemoveFrames(Molecule *mp, IntGroup *group, Vector *outFrame, Vector *outFrameCell);
int MoleculeSelectFrame(Molecule *mp, int frame, int copyback);
int MoleculeFlushFrames(Molecule *mp);
int MoleculeReorderFrames(Molecule *mp, const Int *old_idx);

int MoleculeCreateProperty(Molecule *mp, const char *name);
int MoleculeLookUpProperty(Molecule *mp, const char *name);
int MoleculeDeletePropertyAtIndex(Molecule *mp, int idx);
int MoleculeSetProperty(Molecule *mp, int idx, IntGroup *ig, const Double *values);
int MoleculeGetProperty(Molecule *mp, int idx, IntGroup *ig, Double *outValues);

void MoleculeUpdatePiAnchorPositions(Molecule *mol);
void MoleculeCalculatePiAnchorPosition(Molecule *mol, int idx);
int MoleculeSetPiAnchorList(Molecule *mol, Int idx, Int nentries, Int *entries, Double *weights, Int *nUndoActions, struct MolAction ***undoActions);
	
int MoleculeCalcMO(Molecule *mp, Int mono, const Vector *op, const Vector *dxp, const Vector *dyp, const Vector *dzp, Int nx, Int ny, Int nz, int (*callback)(double progress, void *ref), void *ref);
int MoleculeGetDefaultMOGrid(Molecule *mp, Int npoints, Vector *op, Vector *xp, Vector *yp, Vector *zp, Int *nx, Int *ny, Int *nz);
const Cube *MoleculeGetCubeAtIndex(Molecule *mp, Int index);
int MoleculeLookUpCubeWithMONumber(Molecule *mp, Int mono);
int MoleculeClearCubeAtIndex(Molecule *mp, Int index);
int MoleculeOutputCube(Molecule *mp, Int index, const char *fname, const char *comment);

MCube *MoleculeClearMCube(Molecule *mol, Int nx, Int ny, Int nz, const Vector *origin, Double dx, Double dy, Double dz);
int MoleculeUpdateMCube(Molecule *mol, int idn);
void MoleculeDeallocateMCube(MCube *mcube);

extern char *gMoleculePasteboardType;
extern char *gParameterPasteboardType;
extern char *gLoadSaveErrorMessage;
	
STUB void MoleculeRetainExternalObj(Molecule *mol);
STUB void MoleculeReleaseExternalObj(Molecule *mol);

STUB int MoleculeCallback_writeToPasteboard(const char *type, const void *data, int length);
STUB int MoleculeCallback_readFromPasteboard(const char *type, void **dptr, int *length);
STUB int MoleculeCallback_isDataInPasteboard(const char *type);

STUB Molecule *MoleculeCallback_openNewMolecule(const char *fname);
STUB void MoleculeCallback_notifyModification(Molecule *mp, int now_flag);
STUB Molecule *MoleculeCallback_currentMolecule(void);
STUB Molecule *MoleculeCallback_moleculeAtIndex(int idx);
STUB Molecule *MoleculeCallback_moleculeAtOrderedIndex(int idx);
STUB void MoleculeCallback_displayName(Molecule *mol, char *buf, int bufsize);
STUB void MoleculeCallback_pathName(Molecule *mol, char *buf, int bufsize);
STUB int MoleculeCallback_setDisplayName(Molecule *mol, const char *name);

STUB void MoleculeCallback_lockMutex(void *mutex);
STUB void MoleculeCallback_unlockMutex(void *mutex);
STUB void MoleculeCallback_disableModificationFromGUI(Molecule *mol);
STUB void MoleculeCallback_enableModificationFromGUI(Molecule *mol);
	
STUB void MoleculeCallback_cannotModifyMoleculeDuringMDError(Molecule *mol);

STUB int MoleculeCallback_callSubProcessAsync(Molecule *mol, const char *cmd, int (*callback)(Molecule *, int), int (*timerCallback)(Molecule *, int), FILE *output, FILE *errout);

/*  This is also defined in Molby_extern.h, but it may be called from functions in Molecule.c  */
STUB int MyAppCallback_checkInterrupt(void);
	
void MoleculeLock(Molecule *mol);
void MoleculeUnlock(Molecule *mol);

#if 0
#define __MoleculeLock(mol) MoleculeLock(mol)
#define __MoleculeUnlock(mol) MoleculeUnlock(mol)
#else
#define __MoleculeLock(mol)
#define __MoleculeUnlock(mol)
#endif
	
#ifdef __cplusplus
}
#endif
		
#endif /* __Molecule_h__ */
