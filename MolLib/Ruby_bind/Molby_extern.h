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
#ifndef RubyValue_is_Defined
#define RubyValue_is_Defined
typedef void *RubyValue;
#define RubyNil ((RubyValue)4)
#endif

extern int gMolbyRunLevel;
extern int gMolbyIsCheckingInterrupt;

extern int gSuppressConsole;
extern int gUseGUI;
    
extern char *gRubyVersion;
extern char *gRubyCopyright;

extern int Ruby_WasInterruptRaised(void);
extern int Ruby_getErrorMessage(int status, char **title, char **msg1, char **msg2);
extern void Ruby_showError(int status);
extern void Ruby_clearError(void);
extern int Ruby_showValue(RubyValue value, char **outValueString);
extern int Ruby_UpdateUI(int index, Molecule *mol, int *outChecked, char **outTitle);

extern int Molby_loadScript(const char *script, int from_file);
extern void Molby_startup(const char *script_path, const char *dir);
extern void Molby_getDescription(char **versionString, char **auxString);
extern RubyValue Molby_evalRubyScriptOnMolecule(const char *script, Molecule *mol, const char *fname, int *status);
/* extern RubyValue Molby_evalRubyScript(const char *script, int *status);
extern RubyValue Molby_evalRubyScriptOnActiveMoleculeWithInterrupt(const char *script, int *status); */
/*extern int Ruby_methodType(const char *className, const char *methodName);*/
extern void Molby_buildARGV(int argc, const char **argv);

extern int Molby_updateNamedFragments(int *count, char ***ary);
extern int Molby_evalStringAsType(const char *str, int type, void *ptr);
extern char *Molby_inspectedValueOfType(int type, const void *ptr, int *status);

/*  RubyValue version of Ruby_funcall2_protect()  */
extern RubyValue Ruby_funcall2_protect_extern(RubyValue recv, int mid, int argc, RubyValue *argv, int *status);

extern int g_RubyID_call;  /*  rb_intern("call") for external use  */
	
STUB char *MyAppCallback_getGUIDescriptionString(void);
STUB char *MyAppCallback_getGlobalSettings(const char *key);
STUB void MyAppCallback_setGlobalSettings(const char *key, const char *value);
STUB int MyAppCallback_getGlobalSettingsWithType(const char *key, int type, void *ptr);
STUB int MyAppCallback_setGlobalSettingsWithType(const char *key, int type, const void *ptr);
STUB int MyAppCallback_showScriptMessage(const char *fmt, ...);
STUB void MyAppCallback_setConsoleColor(int color);
STUB void MyAppCallback_showRubyPrompt(void);
STUB int MyAppCallback_checkInterrupt(int id);
STUB int MyAppCallback_showProgressPanel(const char *msg);
STUB void MyAppCallback_hideProgressPanel(int id);
STUB void MyAppCallback_setProgressValue(double dval, int id);
STUB void MyAppCallback_setProgressMessage(const char *msg, int id);
//STUB int MyAppCallback_processUIWithTimeout(double seconds);
STUB int MyAppCallback_getTextWithPrompt(const char *prompt, char *buf, int bufsize);
STUB int MyAppCallback_messageBox(const char *message, const char *title, int flags, int icon);
STUB void MyAppCallback_errorMessageBox(const char *fmt, ...);
STUB char *MyAppCallback_getHomeDir(void);
STUB char *MyAppCallback_getDocumentHomeDir(void);
STUB int MyAppCallback_registerScriptMenu(const char *title);
STUB int MyAppCallback_lookupScriptMenu(const char *title);
STUB RubyValue MyAppCallback_executeScriptFromFile(const char *path, int *status);
STUB int MyAppCallback_callSubProcess(const char **argv, const char *procname, int (*callback)(void *), void *callback_data, FILE *output, FILE *errout, int *exitstatus_p, int *pid_p);
STUB void MyAppCallback_beginUndoGrouping(void);
STUB void MyAppCallback_endUndoGrouping(void);
STUB int MyAppCallback_switchToFilterMode(void);
STUB void MyAppCallback_showConsoleWindow(void);
STUB void MyAppCallback_hideConsoleWindow(void);
STUB void MyAppCallback_bell(void);
STUB int MyAppCallback_playSound(const char *filename, int flag);
STUB void MyAppCallback_stopSound(void);
STUB void MyAppCallback_initImageHandlers(void);
	
#define DUMMY_CALLBACK ((int (*)(void *))1)

#ifdef __cplusplus
}
#endif

#endif /* __Molby_h__ */
