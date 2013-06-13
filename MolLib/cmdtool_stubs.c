/*
 *  cmdtool_stubs.c
 *
 *  Created by Toshi Nagata on Sun Jun 17 2001.
 
 Copyright (c) 2010 Toshi Nagata. All rights reserved.
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

/*  Empty stub functions for command-line tool  */

#include "MolLib.h"
#include "Molby_extern.h"
#include "Missing.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <fcntl.h>

Molecule *
MoleculeCallback_moleculeAtOrderedIndex(int idx)
{
	return NULL;
}

void
MainViewCallback_setNeedsDisplay(MainView *mview, int flag)
{
}

/*
void
MainViewCallback_moleculeReplaced(MainView *mview, struct Molecule *mol)
{
}
*/

void
MyAppCallback_beginUndoGrouping(void)
{
}

Molecule *
MoleculeCallback_currentMolecule(void)
{
	return NULL;
}

void
MoleculeCallback_lockMutex(void *mutex)
{
}

void
MoleculeCallback_unlockMutex(void *mutex)
{
}

void
MyAppCallback_setGlobalSettings(const char *key, const char *value)
{
}

void
MyAppCallback_hideProgressPanel(void)
{
}

int
MyAppCallback_showScriptMessage(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	return vprintf(fmt, ap);
}

char *
MyAppCallback_getDocumentHomeDir(void)
{
	char *s;
	s = getenv("HOME");
	return (s == NULL ? NULL : strdup(s));
}

RubyValue
MyAppCallback_executeScriptFromFile(const char *path, int *status)
{
	return 0;
}

void
MyAppCallback_setProgressMessage(const char *msg)
{
}

int
MainView_isAtomHidden(MainView *mview, int index)
{
	return 0;
}

int
MolActionCallback_isUndoRegistrationEnabled(Molecule *mol)
{
	return 0;
}

int
MyAppCallback_messageBox(const char *message, const char *title, int flags, int icon)
{
	return printf("%s\n%s\n", title, message);
}

void
MainViewCallback_drawInfoText(MainView *mview, const char *label)
{
}

int
MyAppCallback_checkInterrupt(void)
{
	return 0;
}

char *
MyAppCallback_getGlobalSettings(const char *key)
{
	return NULL;
}

void
MyAppCallback_setProgressValue(double dval)
{
}

void
MyAppCallback_showProgressPanel(const char *msg)
{
}

int
MyAppCallback_getTextWithPrompt(const char *prompt, char *buf, int bufsize)
{
	buf[0] = 0;
	return 0;
}

int
MyAppCallback_callSubProcess(const char *cmdline, const char *procname, int (*callback)(void *), void *callback_data)
{
	return system(cmdline);
}

void
MainViewCallback_display(MainView *mview)
{
}

void
MainViewCallback_makeFront(MainView *mview)
{
}

int
MolActionCallback_setUndoRegistrationEnabled(Molecule *mol, int flag)
{
	return 0;
}

void
MoleculeCallback_notifyModification(Molecule *mp, int now_flag)
{
}

void
MyAppCallback_setConsoleColor(int color)
{
}

void
MyAppCallback_errorMessageBox(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
}

void
MyAppCallback_registerScriptMenu(const char *cmd, const char *title)
{
}

int
MyAppCallback_lookupScriptMenu(const char *title)
{
	return 0;
}

void
MyAppCallback_endUndoGrouping(void)
{
}

void
MyAppCallback_showConsoleWindow(void)
{
}

void
MyAppCallback_hideConsoleWindow(void)
{
}

Molecule *
MoleculeCallback_openNewMolecule(const char *fname)
{
	return NULL;
}

void
MainView_getCamera(MainView *mview, Vector *outCamera, Vector *outLookAt, Vector *outUp)
{
}

void
RubyDialogInitClass(void)
{
}

Molecule *
MoleculeCallback_moleculeAtIndex(int idx)
{
	return NULL;
}

void
MoleculeCallback_displayName(Molecule *mol, char *buf, int bufsize)
{
	buf[0] = 0;
}

int
MoleculeCallback_setDisplayName(Molecule *mol, const char *name)
{
	return 0;
}

void
MoleculeCallback_pathName(Molecule *mol, char *buf, int bufsize)
{
	if (mol != NULL && mol->path != NULL) {
		strncpy(buf, mol->path, bufsize - 1);
		buf[bufsize - 1] = 0;
	} else buf[0] = 0;
}

IntGroup *
MainView_selectedMO(MainView *mview)
{
	return NULL;
}

void
MainView_resizeToFit(MainView *mview)
{
}

void
MolActionCallback_registerUndo(Molecule *mol, MolAction *action)
{
}

int
main(int argc, const char **argv)
{
	int fd;
	char *scriptdir;
	static const char fname[] = "startup.rb";
	char *molbydir = getenv("MOLBYDIR");
	if (molbydir == NULL) {
		fprintf(stderr, "Please define the environmental variable MOLBYDIR to specify the location in which the 'Scripts' directory is present.\n");
		exit(1);
	}
	asprintf(&scriptdir, "%s%cScripts", molbydir, PATH_SEPARATOR);
	fd = open(".", O_RDONLY);
	chdir(scriptdir);
	
	Molby_startup(fname, scriptdir);

	fchdir(fd);
	close(fd);
	
	free(scriptdir);
	
	ruby_options(argc, argv);
	ruby_run();
	return 0;
}
