/*
 *  wxKillAddition.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 10/04/29.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#if defined(__WXMSW__)

#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include <wx/msw/missing.h>     // CHARSET_HANGUL
#include <wx/utils.h>
#include <wx/app.h>
#include <wx/intl.h>
#include <wx/log.h>
#include <wx/timer.h>
#endif  //WX_PRECOMP

#include <wx/msw/registry.h>
#include <wx/apptrait.h>
#include <wx/dynlib.h>
#include <wx/dynload.h>
#include <wx/scopeguard.h>

#include <wx/confbase.h>        // for wxExpandEnvVars()

#include <wx/msw/private.h>     // includes <windows.h>

#if defined(__CYGWIN__)
//CYGWIN gives annoying warning about runtime stuff if we don't do this
#   define USE_SYS_TYPES_FD_SET
#   include <sys/types.h>
#endif

// Doesn't work with Cygwin at present
#if wxUSE_SOCKETS && (defined(__GNUWIN32_OLD__) || defined(__WXWINCE__) || defined(__CYGWIN32__))
// apparently we need to include winsock.h to get WSADATA and other stuff
// used in wxGetFullHostName() with the old mingw32 versions
#include <winsock.h>
#endif

#if !defined(__GNUWIN32__) && !defined(__SALFORDC__) && !defined(__WXMICROWIN__) && !defined(__WXWINCE__)
#include <direct.h>

#ifndef __MWERKS__
#include <dos.h>
#endif
#endif  //GNUWIN32

#if defined(__CYGWIN__)
#include <sys/unistd.h>
#include <sys/stat.h>
#include <sys/cygwin.h> // for cygwin_conv_to_full_win32_path()
#endif  //GNUWIN32

#ifdef __BORLANDC__ // Please someone tell me which version of Borland needs
// this (3.1 I believe) and how to test for it.
// If this works for Borland 4.0 as well, then no worries.
#include <dir.h>
#endif

// VZ: there is some code using NetXXX() functions to get the full user name:
//     I don't think it's a good idea because they don't work under Win95 and
//     seem to return the same as wxGetUserId() under NT. If you really want
//     to use them, just #define USE_NET_API
#undef USE_NET_API

#ifdef USE_NET_API
#include <lm.h>
#endif // USE_NET_API

#if defined(__WIN32__) && !defined(__WXMICROWIN__) && !defined(__WXWINCE__)
#ifndef __UNIX__
#include <io.h>
#endif

#ifndef __GNUWIN32__
#include <shellapi.h>
#endif
#endif

#ifndef __WATCOMC__
#if !(defined(_MSC_VER) && (_MSC_VER > 800))
#include <errno.h>
#endif
#endif

// For wxKillAllChildren
#include <tlhelp32.h>

//  Taken from wxWidgets source, and modified so that recursive killing of
//  child processes can be done

typedef HANDLE (WINAPI *CreateToolhelp32Snapshot_t)(DWORD,DWORD);
typedef BOOL (WINAPI *Process32_t)(HANDLE,LPPROCESSENTRY32);

static CreateToolhelp32Snapshot_t lpfCreateToolhelp32Snapshot;
static Process32_t lpfProcess32First, lpfProcess32Next;

static void InitToolHelp32()
{
    static bool s_initToolHelpDone = false;
	
    if (s_initToolHelpDone)
        return;
	
    s_initToolHelpDone = true;
	
    lpfCreateToolhelp32Snapshot = NULL;
    lpfProcess32First = NULL;
    lpfProcess32Next = NULL;
	
#if wxUSE_DYNLIB_CLASS
	
    wxDynamicLibrary dllKernel(_T("kernel32.dll"), wxDL_VERBATIM);
	
    // Get procedure addresses.
    // We are linking to these functions of Kernel32
    // explicitly, because otherwise a module using
    // this code would fail to load under Windows NT,
    // which does not have the Toolhelp32
    // functions in the Kernel 32.
    lpfCreateToolhelp32Snapshot =
	(CreateToolhelp32Snapshot_t)dllKernel.RawGetSymbol(_T("CreateToolhelp32Snapshot"));
	
    lpfProcess32First =
	(Process32_t)dllKernel.RawGetSymbol(_T("Process32First"));
	
    lpfProcess32Next =
	(Process32_t)dllKernel.RawGetSymbol(_T("Process32Next"));
	
#endif // wxUSE_DYNLIB_CLASS
}

typedef struct {
	long pid;
	long p_pid;
	char flag;
} sPidListEntry;

/*  Mark the process with pid and all its children  */
static void
sMarkProcessAndChildren(long pid, sPidListEntry *myList, int count, long *toKill, int *countToKill)
{
	int i;
	sPidListEntry *ep;
	//  Add the process with this pid to the list
	for (i = 0, ep = myList; i < count; i++, ep++) {
		if (ep->flag)
			continue;
		if (ep->pid == pid && *countToKill < count) {
			toKill[*countToKill] = pid;
			(*countToKill)++;
			ep->flag = 1;
			break;
		}
	}
	//  Look for the child process
	for (i = 0, ep = myList; i < count; i++, ep++) {
		if (ep->flag)
			continue;
		if (ep->p_pid == pid)
			sMarkProcessAndChildren(ep->pid, myList, count, toKill, countToKill);
	}
}

int 
myKillAllChildren(long pid, wxSignal sig, wxKillError *krc)
{
    InitToolHelp32();
    if (krc)
        *krc = wxKILL_OK;
	
    // If not implemented for this platform (e.g. NT 4.0), silently ignore
    if (!lpfCreateToolhelp32Snapshot || !lpfProcess32First || !lpfProcess32Next)
        return 0;
	
    // Take a snapshot of all processes in the system.
    HANDLE hProcessSnap = lpfCreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);
    if (hProcessSnap == INVALID_HANDLE_VALUE) {
        if (krc)
            *krc = wxKILL_ERROR;
        return -1;
    }
	
    //Fill in the size of the structure before using it.
    PROCESSENTRY32 pe;
    wxZeroMemory(pe);
    pe.dwSize = sizeof(PROCESSENTRY32);
	
    if (!lpfProcess32First(hProcessSnap, &pe)) {
        // Can't get first process.
        if (krc)
            *krc = wxKILL_ERROR;
        CloseHandle (hProcessSnap);
        return -1;
    }
	
	//Register the pids and parent pids
	int count = 4, i = 0;
	sPidListEntry *myList = (sPidListEntry *)malloc(sizeof(sPidListEntry) * count);
	if (myList == NULL)
		return -1;
    do {
		if (i >= count) {
			count += 4;
			myList = (sPidListEntry *)realloc(myList, sizeof(sPidListEntry) * count);
			if (myList == NULL)
				return -1;
		}
		myList[i].pid = pe.th32ProcessID;
		myList[i].p_pid = pe.th32ParentProcessID;
		myList[i].flag = 0;
		i++;
    } while (lpfProcess32Next (hProcessSnap, &pe));
	CloseHandle (hProcessSnap);
	count = i;
	
	//  Mark the processes to kill
	long *toKill = (long *)malloc(sizeof(long) * count);
	i = 0;
	sMarkProcessAndChildren(pid, myList, count, toKill, &i);
	count = i;
	free(myList);
	
	//  Kill the processes
	for (i = 0; i < count; i++) {
		if (wxKill(toKill[i], sig, krc))
			return -1;
	}
	free(toKill);
	
    return 0;
}
#endif
