/*
 *  MolAction.h
 *
 *  Created by Toshi Nagata on 07/06/23.
 *  Copyright 2007-2008 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#ifndef __MolAction_h__
#define __MolAction_h__

#include "Types.h"

#ifdef __cplusplus
extern "C" {
#endif
	
/*  Action names with signatures; definitions are in MolAction.c  */

extern const char *gMolActionNone;
extern const char *gMolActionAddAnAtom;
extern const char *gMolActionDeleteAnAtom;
extern const char *gMolActionMergeMolecule;
extern const char *gMolActionUnmergeMolecule;
extern const char *gMolActionAddBonds;
extern const char *gMolActionDeleteBonds;
extern const char *gMolActionAddAngles;
extern const char *gMolActionDeleteAngles;
extern const char *gMolActionAddDihedrals;
extern const char *gMolActionDeleteDihedrals;
extern const char *gMolActionAddImpropers;
extern const char *gMolActionDeleteImpropers;
extern const char *gMolActionTranslateAtoms;
extern const char *gMolActionRotateAtoms;
extern const char *gMolActionTransformAtoms;
extern const char *gMolActionSetAtomPositions;
extern const char *gMolActionSetAtomVelocities;
extern const char *gMolActionSetAtomForces;
extern const char *gMolActionInsertFrames;
extern const char *gMolActionRemoveFrames;
extern const char *gMolActionSetSelection;
extern const char *gMolActionChangeResidueNumber;
extern const char *gMolActionChangeResidueNumberForUndo;
extern const char *gMolActionChangeResidueNames;
extern const char *gMolActionOffsetResidueNumbers;
extern const char *gMolActionChangeNumberOfResidues;
extern const char *gMolActionRenumberAtoms;
extern const char *gMolActionExpandBySymmetry;
extern const char *gMolActionDeleteSymmetryOperation;
extern const char *gMolActionAddSymmetryOperation;
extern const char *gMolActionSetCell;
extern const char *gMolActionSetCellPeriodicity;
extern const char *gMolActionSetBox;
extern const char *gMolActionClearBox;
/*extern const char *gMolActionSetBoxForFrames; */
/*extern const char *gMolActionSetCellFlexibility; */
extern const char *gMolActionAddParameters;
extern const char *gMolActionDeleteParameters;
extern const char *gMolActionAmendBySymmetry;
extern const char *gMolActionInsertOnePiAtom;
extern const char *gMolActionReplaceOnePiAtom;
extern const char *gMolActionRemoveOnePiAtom;
extern const char *gMolActionInsertPiBonds;
extern const char *gMolActionRemovePiBonds;
	
/*  Special action signatures to invoke the Ruby script. Used as follows:
 *  MolActionCreateAndPerform(mol, SCRIPT_ACTION("vd"), "rotate", vec, angle);
 *    (Will perform 'mol.rotate(vec, angle)')
 *  or:
 *  MolActionCreateAndPerform(mol, SCRIPT_ACTION("vd"), "proc {|v,d| rotate(v,d)}", vec, deg)
 *    (Will perform '(mol.instance_eval "proc {...}").call(vec, deg)')
 */
#define kMolActionPerformScript "script:s"
#define SCRIPT_ACTION(sig) (kMolActionPerformScript sig)

/*  Action record for reversible editing  */
typedef struct MolAction {
	int refCount;
	const char *name;
	int frame;  /*  The frame number which the action should be performed on.
	                Usually -1 = no specific frame, and set by MolActionSetFrame() if necessary. */
	int nargs;
	struct MolActionArg *args;
} MolAction;

MolAction *MolActionNew(const char *name, ...);
MolAction *MolActionRetain(MolAction *action);
void MolActionRelease(MolAction *action);
void MolActionSetFrame(MolAction *action, int frame);

/*  Perform a MolAction, and register undo action through MolActionCallback_registerUndo(). */
int MolActionPerform(Molecule *mol, MolAction *action);

/*  A convenient function, which creates a MolAction, perform it, and release it.  */
int MolActionCreateAndPerform(Molecule *mol, const char *name, ...);

/*  Show an error dialog saying Ruby is already running  */
void MolActionAlertRubyIsRunning(void);
	
STUB void MolActionCallback_registerUndo(Molecule *mol, MolAction *action);
STUB int MolActionCallback_setUndoRegistrationEnabled(Molecule *mol, int flag);
STUB int MolActionCallback_isUndoRegistrationEnabled(Molecule *mol);

#ifdef __cplusplus
}
#endif
		
#endif /* __MolAction_h__ */
