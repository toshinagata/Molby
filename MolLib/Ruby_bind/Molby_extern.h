/*
 *  Molby_extern.h
 *
 *  Created by Toshi Nagata on 2008/11/05.
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

#ifndef __Molby_extern_h__
#define __Molby_extern_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "../MolLib.h"

/*  This definition is to work around 'VALUE' type in sources without "ruby.h"  */
typedef void *RubyValue;
	
extern char *gRubyVersion;
extern char *gRubyCopyright;

extern void Molby_startup(const char *script_path, const char *dir);
extern void Molby_showError(int status);
extern RubyValue Molby_evalRubyScriptOnMolecule(const char *script, Molecule *mol, const char *fname, int *status);
/* extern RubyValue Molby_evalRubyScript(const char *script, int *status);
extern RubyValue Molby_evalRubyScriptOnActiveMoleculeWithInterrupt(const char *script, int *status); */
extern void Molby_showRubyValue(RubyValue value);
extern int Ruby_methodType(const char *className, const char *methodName);
extern void Molby_buildARGV(int argc, const char **argv);
	
STUB char *MyAppCallback_getGlobalSettings(const char *key);
STUB void MyAppCallback_setGlobalSettings(const char *key, const char *value);
STUB int MyAppCallback_getGlobalSettingsWithType(const char *key, int type, void *ptr);
STUB int MyAppCallback_setGlobalSettingsWithType(const char *key, int type, const void *ptr);
STUB int MyAppCallback_showScriptMessage(const char *fmt, ...);
STUB void MyAppCallback_setConsoleColor(int color);
STUB void MyAppCallback_showRubyPrompt(void);
STUB int MyAppCallback_checkInterrupt(void);
STUB void MyAppCallback_showProgressPanel(const char *msg);
STUB void MyAppCallback_hideProgressPanel(void);
STUB void MyAppCallback_setProgressValue(double dval);
STUB void MyAppCallback_setProgressMessage(const char *msg);
//STUB int MyAppCallback_processUIWithTimeout(double seconds);
STUB int MyAppCallback_getTextWithPrompt(const char *prompt, char *buf, int bufsize);
STUB int MyAppCallback_messageBox(const char *message, const char *title, int flags, int icon);
STUB void MyAppCallback_errorMessageBox(const char *fmt, ...);
STUB char *MyAppCallback_getHomeDir(void);
STUB char *MyAppCallback_getDocumentHomeDir(void);
STUB void MyAppCallback_registerScriptMenu(const char *cmd, const char *title);
STUB int MyAppCallback_lookupScriptMenu(const char *title);
STUB RubyValue MyAppCallback_executeScriptFromFile(const char *path, int *status);
STUB int MyAppCallback_callSubProcess(const char *cmdline, const char *procname);
STUB void MyAppCallback_beginUndoGrouping(void);
STUB void MyAppCallback_endUndoGrouping(void);

#ifdef __cplusplus
}
#endif

#endif /* __Molby_h__ */
