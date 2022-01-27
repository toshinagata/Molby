/*
 *  Molby.h
 *
 *  Created by Toshi Nagata on 2005/06/03.
 *  Copyright 2005-2008 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
*/

#ifndef __Molby_h__
#define __Molby_h__

#include <ruby.h>
#include "../MolLib.h"

#ifndef RSTRING_PTR
#define RSTRING_PTR(_s) (RSTRING(_s)->ptr)
#endif
#ifndef RSTRING_LEN
#define RSTRING_LEN(_s) (RSTRING(_s)->len)
#endif

#ifndef RARRAY_PTR
#define RARRAY_PTR(_a) (RARRAY(_a)->ptr)
#endif
#ifndef RARRAY_LEN
#define RARRAY_LEN(_a) (RARRAY(_a)->len)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*extern VALUE Ruby_IncrementInterruptLevel(VALUE dummy);
extern VALUE Ruby_DecrementInterruptLevel(VALUE dummy); */

/*extern VALUE Ruby_EnableInterrupt(VALUE dummy);
extern VALUE Ruby_DisableInterrupt(VALUE dummy); */

/*  Flag to avoid multiple running ruby interpreter  */
extern int gMolbyRunLevel;
extern int gMolbyIsCheckingInterrupt;
extern VALUE gMolbyBacktrace;

extern VALUE Ruby_CallMethodWithInterrupt(VALUE receiver, ID method_id, VALUE args, int *status);
extern VALUE Ruby_GetInterruptFlag(void);
extern VALUE Ruby_SetInterruptFlag(VALUE val);

extern VALUE Ruby_funcall2_protect(VALUE recv, ID mid, int argc, VALUE *argv, int *status);
extern VALUE Ruby_ObjectAtIndex(VALUE ary, int idx);
	
extern VALUE rb_eMolbyError;
extern VALUE rb_mMolby;
extern VALUE rb_cMolecule, rb_cMolEnumerable, rb_cAtomRef, rb_cIntGroup;
extern VALUE rb_cVector3D, rb_cTransform, rb_cLAMatrix;
extern VALUE rb_cMDArena;
	
extern void Init_MolbyTypes(void);
extern void Init_MolbyMDTypes(void);
	
extern VALUE ValueFromMolecule(Molecule *mol);
extern VALUE ValueFromIntGroup(IntGroup *ig);
extern VALUE ValueFromVector(const Vector *vp);
extern VALUE ValueFromTransform(Transform *tp);
extern VALUE ValueFromMoleculeWithParameterTypeAndIndex(Molecule *mol, int type, int idx1);

extern void VectorFromValue(VALUE val, Vector *vp);
extern void TransformFromValue(VALUE val, Transform *tp);
extern LAMatrix *LAMatrixFromValue(VALUE val, int *needsRelease, int rowHint, int columnHint);
extern IntGroup *IntGroupFromValue(VALUE val);
extern Molecule *MoleculeFromValue(VALUE val);

extern struct MDArena *MDArenaFromValue(VALUE val);
extern VALUE ValueFromMDArena(struct MDArena *arena);

extern void IntGroup_RaiseIfError(int err);
extern VALUE IntGroup_Alloc(VALUE klass);
extern VALUE ValueFromVector(const Vector *vp);

#define FileStringValuePtr(val) Ruby_FileStringValuePtr(&val)
extern char *Ruby_FileStringValuePtr(VALUE *valp);
extern VALUE Ruby_NewFileStringValue(const char *fstr);

#define EncodedStringValuePtr(val) Ruby_EncodedStringValuePtr(&val)
extern char *Ruby_EncodedStringValuePtr(VALUE *valp);
extern VALUE Ruby_NewEncodedStringValue(const char *str, int len);
extern VALUE Ruby_NewEncodedStringValue2(const char *str);

/*
STUB VALUE MyAppCallback_getGlobalSettings(const char *key);
STUB void MyAppCallback_setGlobalSettings(const char *key, VALUE value);
STUB int MyAppCallback_showScriptMessage(const char *fmt, ...);
STUB void MyAppCallback_showScriptError(int status);
STUB int MyAppCallback_checkInterrupt(void);
STUB void MyAppCallback_showProgressPanel(const char *msg);
STUB void MyAppCallback_hideProgressPanel(void);
STUB void MyAppCallback_setProgressValue(double dval);
STUB void MyAppCallback_setProgressMessage(const char *msg);
STUB int MyAppCallback_processUIWithTimeout(double seconds);
STUB int MyAppCallback_getTextWithPrompt(const char *prompt, char *buf, int bufsize);
STUB int MyAppCallback_registerScriptMenu(const char *title);
*/

#ifdef __cplusplus
}
#endif
		
#include "Molby_extern.h"

#endif /* __Molby_h__ */
