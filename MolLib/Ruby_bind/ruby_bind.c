/*
 *  ruby_bind.c
 *  Ruby binding
 *
 *  Created by Toshi Nagata on 07/11/09.
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

#include "Molby.h"
#include "ruby_dialog.h"
#include "../MD/MDCore.h"

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include <version.h>  /*  for Ruby version  */
#include <node.h>     /*  for rb_add_event_hook()  */

#if defined(__WXMAC__) || defined(__CMDMAC__)
#define USE_PTHREAD_FOR_TIMER 1
#endif

#if !__WXMSW__
#if USE_PTHREAD_FOR_TIMER
#include <unistd.h>   /*  for usleep()  */
#include <pthread.h>  /*  for pthread  */
#else
#include <signal.h>   /*  for sigaction()  */
#endif
#endif

#include "../Missing.h"

#pragma mark ====== Global Values ======

VALUE rb_eMolbyError;
VALUE rb_mMolby;
VALUE rb_cMolecule, rb_cMolEnumerable, rb_cAtomRef;
VALUE rb_cParameter, rb_cParEnumerable, rb_cParameterRef;

VALUE gMolbyBacktrace;

int gMolbyRunLevel = 0;
int gMolbyIsCheckingInterrupt = 0;

char *gRubyVersion, *gRubyCopyright;

/*  For convenience  */
static ID s_ID_equal;  /*  rb_intern("==")  */

int g_RubyID_call;

/*  Symbols for atom attributes  */
static VALUE
	s_IndexSym, s_SegSeqSym, s_SegNameSym, s_ResSeqSym,
	s_ResNameSym, s_NameSym, s_AtomTypeSym, s_ChargeSym,
	s_WeightSym, s_ElementSym, s_AtomicNumberSym, s_ConnectsSym,
	s_RSym, s_XSym, s_YSym, s_ZSym,
	s_FractRSym, s_FractXSym, s_FractYSym, s_FractZSym,
	s_SigmaSym, s_SigmaXSym, s_SigmaYSym, s_SigmaZSym,
	s_VSym, s_FSym, s_OccupancySym, s_TempFactorSym,
	s_AnisoSym, s_SymopSym, s_IntChargeSym, s_FixForceSym,
	s_FixPosSym, s_ExclusionSym, s_MMExcludeSym, s_PeriodicExcludeSym,
	s_HiddenSym, s_AnchorListSym, s_UFFTypeSym;

/*  Symbols for parameter attributes  */
static VALUE
	s_ParTypeSym, s_AtomTypesSym, s_KSym, s_R0Sym,
	s_A0Sym, s_MultSym, s_PeriodSym, s_Phi0Sym,
	/* s_ASym, s_BSym, */
	s_ReqSym, s_EpsSym,
	/* s_A14Sym, s_B14Sym, */
	s_Req14Sym, s_Eps14Sym,
	s_CutoffSym, s_RadiusSym, s_ColorSym, s_FullNameSym,
	s_CommentSym, s_SourceSym;

/*
 *  Utility function
 *  Get ary[i] by calling "[]" method
 */
VALUE
Ruby_ObjectAtIndex(VALUE ary, int idx)
{
	static ID index_method = 0;
	if (TYPE(ary) == T_ARRAY) {
		int len = RARRAY_LEN(ary);
		if (idx >= 0 && idx < len)
			return (RARRAY_PTR(ary))[idx];
		else return Qnil;
	}
	if (index_method == 0)
		index_method = rb_intern("[]");
	return rb_funcall(ary, index_method, 1, INT2NUM(idx));
}

char *
Ruby_FileStringValuePtr(VALUE *valp)
{
#if __WXMSW__
	char *p = strdup(StringValuePtr(*valp));
	translate_char(p, '/', '\\');
	*valp = rb_str_new2(p);
	free(p);
	return StringValuePtr(*valp);
#else
	return StringValuePtr(*valp);
#endif
}

VALUE
Ruby_NewFileStringValue(const char *fstr)
{
#if __WXMSW__
	VALUE retval;
	char *p = strdup(fstr);
	translate_char(p, '\\', '/');
	retval = rb_str_new2(p);
	free(p);
	return retval;
#else
	return rb_str_new2(fstr);
#endif
}

VALUE
Ruby_ObjToStringObj(VALUE val)
{
	switch (TYPE(val)) {
		case T_STRING:
			return val;
		case T_SYMBOL:
			return rb_str_new2(rb_id2name(SYM2ID(val)));
		default:
			return rb_str_to_str(val);
	}
}

#pragma mark ====== Message input/output ======

/*
 *  call-seq:
 *     message_box(str, title, button = nil, icon = :info)
 *
 *  Show a message box.
 *  Buttons: nil (ok and cancel), :ok (ok only), :cancel (cancel only)
 *  Icon: :info, :warning, :error
 */
static VALUE
s_Kernel_MessageBox(int argc, VALUE *argv, VALUE self)
{
	char *str, *title, *s;
	int buttons, icon;
	VALUE sval, tval, bval, ival;
	rb_scan_args(argc, argv, "22", &sval, &tval, &bval, &ival);
	str = StringValuePtr(sval);
	title = StringValuePtr(tval);
	if (bval != Qnil) {
		bval = Ruby_ObjToStringObj(bval);
		s = RSTRING_PTR(bval);
		if (strncmp(s, "ok", 2) == 0)
			buttons = 1;
		else if (strncmp(s, "cancel", 6) == 0)
			buttons = 2;
		else
			rb_raise(rb_eMolbyError, "the button specification should be either nil, :ok or :cancel");
	} else buttons = 3;
	if (ival != Qnil) {
		ival = Ruby_ObjToStringObj(ival);
		s = RSTRING_PTR(ival);
		if (strncmp(s, "info", 4) == 0)
			icon = 1;
		else if (strncmp(s, "warn", 4) == 0)
			icon = 2;
		else if (strncmp(s, "err", 3) == 0)
			icon = 3;
		else
			rb_raise(rb_eMolbyError, "the icon specification should be either :info, :warning or :error");
	} else icon = 1;
	MyAppCallback_messageBox(str, title, buttons, icon);
	return Qnil;
}

/*
 *  call-seq:
 *     error_message_box(str)
 *
 *  Show an error message box.
 */
static VALUE
s_Kernel_ErrorMessageBox(VALUE self, VALUE sval)
{
	char *str = StringValuePtr(sval);
	MyAppCallback_errorMessageBox("%s", str);
	return Qnil;
}

/*
 *  call-seq:
 *     ask(prompt, default = nil) -> string
 *
 *  Open a modal dialog and get a line of text.
 */
static VALUE
s_Kernel_Ask(int argc, VALUE *argv, VALUE self)
{
	volatile VALUE prompt, message;
	char buf[1024];
	int retval;
	rb_scan_args(argc, argv, "11", &prompt, &message);
	if (message != Qnil) {
		strncpy(buf, StringValuePtr(message), sizeof buf - 1);
		buf[sizeof buf - 1] = 0;
	} else buf[0] = 0;
	retval = MyAppCallback_getTextWithPrompt(StringValuePtr(prompt), buf, sizeof buf);
	if (retval)
		return rb_str_new2(buf);
	else
		return Qnil;	
}

/*
 *  call-seq:
 *     show_console_window
 *
 *  Show the console window and bring to the front.
 */
static VALUE
s_Kernel_ShowConsoleWindow(VALUE self)
{
	MyAppCallback_showConsoleWindow();
	return Qnil;
}

/*
 *  call-seq:
 *     hide_console_window
 *
 *  Hide the console window.
 */
static VALUE
s_Kernel_HideConsoleWindow(VALUE self)
{
	MyAppCallback_hideConsoleWindow();
	return Qnil;
}

/*
 *  call-seq:
 *     bell
 *
 *  Ring the system bell.
 */
static VALUE
s_Kernel_Bell(VALUE self)
{
	MyAppCallback_bell();
	return Qnil;
}

/*
 *  call-seq:
 *     play_sound(filename, flag = 0)
 *
 *  Play the sound (a WAV file) in the file. Flag: 0, pause until sound ends;
 *  1, play the sound asynchronously; 3, play the sound with loop asynchronously
 */
static VALUE
s_Kernel_PlaySound(int argc, VALUE *argv, VALUE self)
{
	VALUE fnval, flval;
	int flag, retval;
	char *fname;
	rb_scan_args(argc, argv, "11", &fnval, &flval);
	if (flval == Qnil)
		flag = 0;
	else flag = NUM2INT(rb_Integer(flval));
	fnval = rb_funcall(rb_cFile, rb_intern("expand_path"), 1, fnval);
	fname = StringValuePtr(fnval);
	retval = MyAppCallback_playSound(fname, flag);
	return (retval ? Qtrue : Qnil);
}

/*
 *  call-seq:
 *     stop_sound
 *
 *  Stop the sound if playing.
 */
static VALUE
s_Kernel_StopSound(VALUE self)
{
	MyAppCallback_stopSound();
	return Qnil;
}

/*
 *  call-seq:
 *     export_to_clipboard(str)
 *
 *  Export the given string to clipboard.
 */
static VALUE
s_Kernel_ExportToClipboard(VALUE self, VALUE sval)
{
#if !defined(__CMDMAC__)
	const char *s = StringValuePtr(sval);
	char *ns;
#if __WXMSW__
	/*  Convert the end-of-line characters  */
	{	const char *p; int nc; char *np;
		nc = 0;
		for (p = s; *p != 0; p++) {
			if (*p == '\n')
				nc++;
		}	
		ns = (char *)malloc(strlen(s) + nc + 1);
		for (np = ns, p = s; *p != 0; p++, np++) {
			if (*p == '\n')
				*np++ = '\r';
			*np = *p;
		}
		*np = 0;
	}
#else
	ns = (char *)malloc(strlen(s) + 1);
	strcpy(ns, s);
#if __WXMAC__
	{	char *np;
		/*  wxMac still has Carbon code. Oops.  */
		for (np = ns; *np != 0; np++) {
			if (*np == '\n')
				*np = '\r';
		}
	}
#endif
#endif
	if (MoleculeCallback_writeToPasteboard("TEXT", ns, strlen(ns) + 1))
		rb_raise(rb_eMolbyError, "Cannot export string to clipboard");
#endif
	return Qnil;
}

/*
 *  call-seq:
 *     stdout.write(str)
 *
 *  Put the message in the main text view in black color.
 */
static VALUE
s_StandardOutput(VALUE self, VALUE str)
{
	int n;
	MyAppCallback_setConsoleColor(0);
	n = MyAppCallback_showScriptMessage("%s", StringValuePtr(str));
	return INT2NUM(n);
}

/*
 *  call-seq:
 *     stderr.write(str)
 *
 *  Put the message in the main text view in red color.
 */
static VALUE
s_StandardErrorOutput(VALUE self, VALUE str)
{
	int n;
	MyAppCallback_setConsoleColor(1);
	n = MyAppCallback_showScriptMessage("%s", StringValuePtr(str));
	MyAppCallback_setConsoleColor(0);
	return INT2NUM(n);
}

/*
 *  call-seq:
 *     stdin.gets(rs = $/)
 *
 *  Read one line message via dialog box.
 */
static VALUE
s_StandardInputGets(int argc, VALUE *argv, VALUE self)
{
	VALUE pval, rval;
	pval = rb_str_new2("Enter a line:");
	rval = s_Kernel_Ask(1, &pval, self);
	if (rval == Qnil)
		rb_interrupt();
	rb_str_cat2(rval, "\n");
	return rval;
}

/*
 *  call-seq:
 *     stdin.method_missing(name, args, ...)
 *
 *  Throw an exception, noting only gets and readline are defined.
 */
static VALUE
s_StandardInputMethodMissing(int argc, VALUE *argv, VALUE self)
{
	VALUE nval;
	rb_scan_args(argc, argv, "10", &nval);
	rb_raise(rb_eMolbyError, "'%s' is undefined. Only 'gets' and 'readline' can be used for stdin within Molby.", rb_id2name(SYM2ID(nval)));
	return Qnil;  /*  Not reached  */
}

#pragma mark ====== Track key events ======

/*  User interrupt handling
 *  User interrupt (command-period on Mac OS) is handled by periodic polling of
 *  key events. This polling should only be enabled during "normal" execution
 *  of scripts and must be disabled when the rest of the application (or Ruby
 *  script itself) is handling GUI. This is ensured by appropriate calls to
 *  enable_interrupt and disable_interrupt.  */

static VALUE s_interrupt_flag = Qfalse;

static VALUE
s_ShowProgressPanel(int argc, VALUE *argv, VALUE self)
{
	volatile VALUE message;
	const char *p;
	if (Ruby_GetInterruptFlag() == Qtrue) {
		rb_scan_args(argc, argv, "01", &message);
		if (message != Qnil)
			p = StringValuePtr(message);
		else
			p = NULL;
		MyAppCallback_showProgressPanel(p);
	}
	return Qnil;
}

static VALUE
s_HideProgressPanel(VALUE self)
{
	MyAppCallback_hideProgressPanel();
	return Qnil;
}

static VALUE
s_SetProgressValue(VALUE self, VALUE val)
{
	double dval = NUM2DBL(rb_Float(val));
	MyAppCallback_setProgressValue(dval);
	return Qnil;
}

static VALUE
s_SetProgressMessage(VALUE self, VALUE msg)
{
	const char *p;
	if (msg == Qnil)
		p = NULL;
	else p = StringValuePtr(msg);
	MyAppCallback_setProgressMessage(p);
	return Qnil;
}

static VALUE
s_SetInterruptFlag(VALUE self, VALUE val)
{
	VALUE oldval;
	if (val != Qundef) {
		if (val == Qfalse || val == Qnil)
			val = Qfalse;
		else val = Qtrue;
	}
	oldval = s_interrupt_flag;
	if (val != Qundef) {
		s_interrupt_flag = val;
		if (val == Qfalse) {
			s_HideProgressPanel(self);
		}
	}
	return oldval;
}

static VALUE
s_GetInterruptFlag(VALUE self)
{
	return s_SetInterruptFlag(self, Qundef);
}

#if 0
static VALUE
s_Ruby_CallMethod(VALUE val)
{
	void **ptr = (void **)val;
	VALUE receiver = (VALUE)ptr[0];
	ID method_id = (ID)ptr[1];
	VALUE args = (VALUE)ptr[2];
	VALUE retval;
	if (method_id == 0) {
		/*  args should be a string, which is evaluated  */
		if (receiver == Qnil) {
			retval = rb_eval_string(StringValuePtr(args));
		} else {
			retval = rb_obj_instance_eval(1, &args, receiver);
		}
	} else {
		/*  args should be an array of arguments  */
		retval = rb_apply(receiver, method_id, args);
	}
	return retval;
}

VALUE
Ruby_CallMethodWithInterrupt(VALUE receiver, ID method_id, VALUE args, int *status)
{
	VALUE retval, save_interrupt_flag;
	void *ptr[3];
	save_interrupt_flag = s_SetInterruptFlag(Qnil, Qtrue);
	ptr[0] = (void *)receiver;
	ptr[1] = (void *)method_id;
	ptr[2] = (void *)args;
	MyAppCallback_beginUndoGrouping();
	retval = rb_protect(s_Ruby_CallMethod, (VALUE)ptr, status);
	MyAppCallback_endUndoGrouping();
	s_SetInterruptFlag(Qnil, save_interrupt_flag);
	MyAppCallback_hideProgressPanel();  /*  In case when the progress panel is still onscreen */
	return retval;
}
#endif

VALUE
Ruby_SetInterruptFlag(VALUE val)
{
	return s_SetInterruptFlag(Qnil, val);
}

VALUE
Ruby_GetInterruptFlag(void)
{
	return s_SetInterruptFlag(Qnil, Qundef);
}

/*
 *  call-seq:
 *     check_interrupt -> integer
 *
 *  Returns 1 if interrupted, 0 if not, -1 if interrupt is disabled.
 */
static VALUE
s_Kernel_CheckInterrupt(VALUE self)
{
	if (Ruby_GetInterruptFlag() == Qfalse)
		return INT2NUM(-1);
	else if (MyAppCallback_checkInterrupt())
		return INT2NUM(1);
	else return INT2NUM(0);
}

static volatile unsigned long sITimerCount = 0;

#if __WXMSW__
static HANDLE sITimerEvent;
static HANDLE sITimerThread;
static int sITimerInterval;

static __stdcall unsigned
s_ITimerThreadFunc(void *p)
{
	while (WaitForSingleObject(sITimerEvent, sITimerInterval) == WAIT_TIMEOUT) {
		sITimerCount++;
	}
	return 0;
}

#elif USE_PTHREAD_FOR_TIMER

/*  Timer thread  */
static pthread_t sTimerThread;

/*  -1: uninitiated; 0: active, 1: inactive, -2: request to terminate  */
static volatile signed char sTimerFlag = -1;
static volatile int sTimerIntervalMicrosec = 0;

static void *
s_TimerThreadEntry(void *param)
{
	while (1) {
		usleep(sTimerIntervalMicrosec);
		if (sTimerFlag == 0)
			sITimerCount++;
		else if (sTimerFlag == -2)
			break;
	}
	return NULL;	
}

#endif

static void
s_SignalAction(int n)
{
	sITimerCount++;
}

static void
s_SetIntervalTimer(int n, int msec)
{
#if __WXMSW__
	if (n == 0) {
		/*  Start interval timer  */
		sITimerEvent = CreateEvent(NULL, FALSE, FALSE, NULL);
		sITimerInterval = msec;
		if (sITimerEvent) {
			sITimerThread = (HANDLE)_beginthreadex(NULL, 0, s_ITimerThreadFunc, NULL, 0, NULL);
		}
	} else {
		/*  Stop interval timer  */
		if (sITimerEvent)
			SetEvent(sITimerEvent);   /*  Tell thread to terminate  */
		if (sITimerThread) {
			WaitForSingleObject(sITimerThread, 1000);
			CloseHandle(sITimerThread);
		}
		if (sITimerEvent)
			CloseHandle(sITimerEvent);
		sITimerEvent = NULL;
		sITimerThread = NULL;
	}
#elif USE_PTHREAD_FOR_TIMER
	if (n == 0) {
		if (sTimerFlag == -1) {
			int status = pthread_create(&sTimerThread, NULL, s_TimerThreadEntry, NULL);
			if (status != 0) {
				fprintf(stderr, "pthread_create failed while setting Ruby interval timer: status = %d\n", status);
			}
		}
		sTimerFlag = 0;  /*  Active  */
		sTimerIntervalMicrosec = msec * 1000;
	} else if (sTimerFlag != -1)
		sTimerFlag = 1;  /*  Inactive  */	
#else
	static struct itimerval sOldValue;
	static struct sigaction sOldAction;
	struct itimerval val;
	struct sigaction act;
	if (n == 0) {
		sITimerCount = 0;
		act.sa_handler = s_SignalAction;
		act.sa_mask = 0;
		act.sa_flags = 0;
		sigaction(SIGALRM, &act, &sOldAction);
		val.it_value.tv_sec = 0;
		val.it_value.tv_usec = msec * 1000;
		val.it_interval.tv_sec = 0;
		val.it_interval.tv_usec = msec * 1000;
		setitimer(ITIMER_REAL, &val, &sOldValue);
	} else {
		setitimer(ITIMER_REAL, &sOldValue, &val);
		sigaction(SIGALRM, &sOldAction, &act);
	}
#endif
}

static unsigned long
s_GetTimerCount(void)
{
	return sITimerCount;
}

static void
s_Event_Callback(rb_event_t event, NODE *node, VALUE self, ID rid, VALUE klass)
{
	if (s_interrupt_flag != Qfalse) {
		static unsigned long sLastTime = 0;
		unsigned long currentTime;
		int flag;
		currentTime = s_GetTimerCount();
		if (currentTime != sLastTime) {
			sLastTime = currentTime;
			gMolbyIsCheckingInterrupt = 1;
			flag = MyAppCallback_checkInterrupt();
			gMolbyIsCheckingInterrupt = 0;
			if (flag) {
				s_SetInterruptFlag(Qnil, Qfalse);
				rb_interrupt();
			}
		}
	}
}

#pragma mark ====== Menu handling ======

/*
 *  call-seq:
 *     register_menu(title, method)
 *
 *  Register the method (specified as a symbol) in the script menu.
 *  The method must be either an instance method of Molecule with no argument,
 *  or a class method of Molecule with one argument (the current molecule).
 *  The menu associated with the class method can be invoked even when no document
 *  is open (the argument is set to Qnil in this case). On the other hand, the
 *  menu associated with the instance method can only be invoked when at least one 
 *  document is active.
 */
static VALUE
s_Kernel_RegisterMenu(VALUE self, VALUE title, VALUE method)
{
	if (TYPE(method) == T_SYMBOL) {
		method = rb_funcall(method, rb_intern("to_s"), 0);
	}
	MyAppCallback_registerScriptMenu(StringValuePtr(method), StringValuePtr(title));
	return self;
}

static VALUE
s_Kernel_LookupMenu(VALUE self, VALUE title)
{
	int n = MyAppCallback_lookupScriptMenu(StringValuePtr(title));
	return INT2NUM(n);
}

static VALUE
s_Ruby_methodType_sub(VALUE data)
{
	const char **p = (const char **)data;
	VALUE klass = rb_const_get(rb_cObject, rb_intern(p[0]));
	ID mid = rb_intern(p[1]);
	int ival;
	if (rb_funcall(klass, rb_intern("method_defined?"), 1, ID2SYM(mid)) != Qfalse)
		ival = 1;
	else if (rb_respond_to(klass, mid))
		ival = 2;
	else ival = 0;
	return INT2FIX(ival);
}
	
/*  Returns 1 if the class defines the instance method with the given name, 2 if the class
    has the singleton method (class method) with the given name, 0 otherwise.  */
int
Ruby_methodType(const char *className, const char *methodName)
{
	int status;
	VALUE retval;
	const char *p[2];
	p[0] = className;
	p[1] = methodName;
	retval = rb_protect(s_Ruby_methodType_sub, (VALUE)p, &status);
	if (status == 0)
		return FIX2INT(retval);
	else return 0;
}

/*
 *  call-seq:
 *     execute_script_file(fname)
 *
 *  Execute the script in the given file. If a molecule is active, then
 *  the script is evaluated as Molecule.current.instance_eval(script).
 *  Before entering the script, the current directory is set to the parent
 *  directory of the script.
 */
static VALUE
s_Kernel_ExecuteScript(VALUE self, VALUE fname)
{
	int status;
	VALUE retval = (VALUE)MyAppCallback_executeScriptFromFile(StringValuePtr(fname), &status);
	if (retval == (VALUE)6 && status == -1)
		rb_raise(rb_eMolbyError, "Cannot open script file: %s", StringValuePtr(fname));
	if (status != 0)
		rb_jump_tag(status);
	return retval;
}

/*
 *  call-seq:
 *     document_home
 *
 *  Get the directory suitable for storing user documents. On Windows
 *  it is the home directory + "My Documents". On other platforms
 *  it is the home directory.
 */
static VALUE
s_Kernel_DocumentHome(VALUE self)
{
	char *s = MyAppCallback_getDocumentHomeDir();
	VALUE retval = Ruby_NewFileStringValue(s);
	free(s);
	return retval;
}

/*  The callback function for call_subprocess  */
static int
s_Kernel_CallSubProcess_Callback(void *data)
{
	int status;
	VALUE retval = Ruby_funcall2_protect((VALUE)data, rb_intern("call"), 0, NULL, &status);
	if (status != 0 || retval == Qnil || retval == Qfalse)
		return 1;
	else return 0;
}

/*
 *  call-seq:
 *     call_subprocess(cmd, process_name [, callback_proc [, stdout_file [, stderr_file]]])
 *
 *  Call subprocess. A progress dialog window is displayed, with a message
 *  "Running #{process_name}...".
 *  A callback proc can be given, which is called periodically during execution. If the proc returns
 *  nil or false, then the execution will be interrupted.
 *  If stdout_file or stderr_file is a filename, then the message will be sent to the file; if the
 *  filename begins with ">>", then the message will be appended to the file.
 *  If the filename is "/dev/null" or "NUL", then the message will be lost.
 *  If the argument is nil, then the message will be sent to the Ruby console.
 */
static VALUE
s_Kernel_CallSubProcess(int argc, VALUE *argv, VALUE self)
{
	VALUE cmd, procname, cproc, stdout_val, stderr_val;
	int n;
	char *sout, *serr;
	FILE *fpout, *fperr;

	rb_scan_args(argc, argv, "23", &cmd, &procname, &cproc, &stdout_val, &stderr_val);

	if (stdout_val == Qnil) {
		fpout = (FILE *)1;
	} else {
		sout = StringValuePtr(stdout_val);
		if (strcmp(sout, "/dev/null") == 0 || strcmp(sout, "NUL") == 0)
			fpout = NULL;
		else {
			if (strncmp(sout, ">>", 2) == 0) {
				sout += 2;
				fpout = fopen(sout, "a");
			} else {
				if (*sout == '>')
					sout++;
				fpout = fopen(sout, "w");
			}
			if (fpout == NULL)
				rb_raise(rb_eMolbyError, "Cannot open file for standard output (%s)", sout);
		}
	}
	if (stderr_val == Qnil) {
		fperr = (FILE *)1;
	} else {
		serr = StringValuePtr(stderr_val);
		if (strcmp(serr, "/dev/null") == 0 || strcmp(serr, "NUL") == 0)
			fperr = NULL;
		else {
			if (strncmp(serr, ">>", 2) == 0) {
				serr += 2;
				fpout = fopen(serr, "a");
			} else {
				if (*serr == '>')
					serr++;
				fperr = fopen(serr, "w");
			}
			if (fperr == NULL)
				rb_raise(rb_eMolbyError, "Cannot open file for standard output (%s)", serr);
		}
	}

	n = MyAppCallback_callSubProcess(StringValuePtr(cmd), StringValuePtr(procname), (cproc == Qnil ? NULL : s_Kernel_CallSubProcess_Callback), (cproc == Qnil ? NULL : (void *)cproc), fpout, fperr);
	
	if (fpout != NULL && fpout != (FILE *)1)
		fclose(fpout);
	if (fperr != NULL && fperr != (FILE *)1)
		fclose(fperr);

	return INT2NUM(n);

	
}

/*
 *  call-seq:
 *     backquote(cmd)
 *
 *  Same as the builtin backquote, except that, under Windows, no console window gets opened.
 */
static VALUE
s_Kernel_Backquote(VALUE self, VALUE cmd)
{
	char *buf;
	int n;
	VALUE val;
	n = MyAppCallback_callSubProcess(StringValuePtr(cmd), NULL, DUMMY_CALLBACK, &buf, NULL, NULL);
	if (n != 0)
		rb_raise(rb_eMolbyError, "Cannot invoke command '%s'", StringValuePtr(cmd));
	if (buf != NULL) {
		val = rb_str_new2(buf);
		free(buf);
	} else {
		val = rb_str_new2("");
	}
	return val;
}

#pragma mark ====== User defaults ======

/*
 *  call-seq:
 *     get_global_settings(key)
 *
 *  Get a setting data for key from the application preferences.
 */
static VALUE
s_Kernel_GetGlobalSettings(VALUE self, VALUE key)
{
	char *p = MyAppCallback_getGlobalSettings(StringValuePtr(key));
	if (p != NULL) {
		VALUE retval = rb_eval_string(p);
		free(p);
		return retval;
	} else return Qnil;
}

/*
 *  call-seq:
 *     set_global_settings(key, value)
 *
 *  Set a setting data for key to the application preferences.
 */
static VALUE
s_Kernel_SetGlobalSettings(VALUE self, VALUE key, VALUE value)
{
	VALUE sval = rb_inspect(value);
	MyAppCallback_setGlobalSettings(StringValuePtr(key), StringValuePtr(sval));
	return value;
}

#pragma mark ====== Utility functions (protected funcall) ======

struct Ruby_funcall2_record {
	VALUE recv;
	ID mid;
	int argc;
	VALUE *argv;
};

static VALUE
s_Ruby_funcall2_sub(VALUE data)
{
	struct Ruby_funcall2_record *rp = (struct Ruby_funcall2_record *)data;
	return rb_funcall2(rp->recv, rp->mid, rp->argc, rp->argv);
}

VALUE
Ruby_funcall2_protect(VALUE recv, ID mid, int argc, VALUE *argv, int *status)
{
	struct Ruby_funcall2_record rec;
	rec.recv = recv;
	rec.mid = mid;
	rec.argc = argc;
	rec.argv = argv;
	return rb_protect(s_Ruby_funcall2_sub, (VALUE)&rec, status);
}

RubyValue
Ruby_funcall2_protect_extern(RubyValue recv, int mid, int argc, RubyValue *argv, int *status)
{
	return (RubyValue)Ruby_funcall2_protect((VALUE)recv, mid, argc, (VALUE *)argv, status);
}

#pragma mark ====== ParameterRef Class ======

static UnionPar *
s_UnionParFromValue(VALUE self, Int *typep, Int checkEditable)
{
	UnionPar *up;
	ParameterRef *pref;
	Data_Get_Struct(self, ParameterRef, pref);
	if (typep != NULL)
		*typep = pref->parType;
	if (pref->parType == kElementParType) {
		up = (UnionPar *)&gElementParameters[pref->idx];
	} else {
		up = ParameterRefGetPar(pref);
		if (checkEditable) {
			if (pref->idx < 0)
				rb_raise(rb_eMolbyError, "Cannot modify parameter because it is internally cached in the MDArena");
			if (up->bond.src != 0 && up->bond.src != -1)
				rb_raise(rb_eMolbyError, "Cannot modify parameter because it is not molecule-local");
		}
	}
	return up;
}

static void
s_RegisterUndoForParameterAttrChange(VALUE self, VALUE key, VALUE val, VALUE oldval, int oldsrc)
{
	UnionPar *up;
	ParameterRef *pref;
	Data_Get_Struct(self, ParameterRef, pref);
	if (pref->mol == NULL)
		return;
	up = ParameterRefGetPar(pref);
	if (key != s_SourceSym)
		up->bond.src = 0;  /*  Becomes automatically molecule-local  */
	if (MolActionCallback_isUndoRegistrationEnabled(pref->mol)) {
		/*  Register undo  */
		MolAction *act;
		act = MolActionNew(SCRIPT_ACTION("iirri"), "set_parameter_attr", pref->parType, pref->idx, key, oldval, oldsrc);
		MolActionCallback_registerUndo(pref->mol, act);
		MoleculeCallback_notifyModification(pref->mol, 0);
		pref->mol->needsMDRebuild = 1;
	}
}

VALUE
ValueFromMoleculeWithParameterTypeAndIndex(Molecule *mol, int type, int idx1)
{
	ParameterRef *pref = ParameterRefNew(mol, type, idx1);
	if (pref != NULL)
		return Data_Wrap_Struct(rb_cParameterRef, 0, (void (*)(void *))ParameterRefRelease, pref);
	else
		rb_raise(rb_eMolbyError, "Cannot create parameter reference");
}

static int
s_AtomTypeIndexFromValue(VALUE val)
{
	if (rb_obj_is_kind_of(val, rb_cNumeric))
		return NUM2INT(val);
	else
		return AtomTypeEncodeToUInt(StringValuePtr(val));
}

static const char *s_ParameterTypeNames[] = {
	"bond", "angle", "dihedral", "improper", "vdw", "vdw_pair", "vdw_cutoff", "element"
};
static ID s_ParameterTypeIDs[8] = {0, 0, 0, 0, 0, 0, 0, 0};

static int
s_ParTypeFromValue(VALUE val)
{
	int i, n;
	ID valid;
	n = sizeof(s_ParameterTypeNames) / sizeof(s_ParameterTypeNames[0]);
	if (s_ParameterTypeIDs[0] == 0) {
		for (i = 0; i < n; i++)
			s_ParameterTypeIDs[i] = rb_intern(s_ParameterTypeNames[i]);
	}
	valid = rb_to_id(val);
	for (i = 0; i < n; i++) {
		if (valid == s_ParameterTypeIDs[i]) {
			if (i == 7)
				return kElementParType;
			else return kFirstParType + i;
		}
	}
	return kInvalidParType;
}

/*
 *  call-seq:
 *     index -> Integer
 *
 *  Get the index in the parameter list.
 */
static VALUE s_ParameterRef_GetIndex(VALUE self) {
	ParameterRef *pref;
	Data_Get_Struct(self, ParameterRef, pref);
	return INT2NUM(pref->idx);
}

/*
 *  call-seq:
 *     par_type -> String
 *
 *  Get the parameter type, like "bond", "angle", etc.
 */
static VALUE s_ParameterRef_GetParType(VALUE self) {
	Int tp;
	s_UnionParFromValue(self, &tp, 0);
	if (tp == kElementParType)
		return rb_str_new2("element");
	tp -= kFirstParType;
	if (tp >= 0 && tp < sizeof(s_ParameterTypeNames) / sizeof(s_ParameterTypeNames[0]))
		return rb_str_new2(s_ParameterTypeNames[tp]);
	else rb_raise(rb_eMolbyError, "Internal error: parameter type tag is out of range (%d)", tp);
}

/*
 *  call-seq:
 *     atom_type -> String or Array of String
 *     atom_types -> String or Array of String
 *
 *  Get the atom types. For a bond parameter, an array of two strings (like ["ca", "ha"])
 *  is returned. For an angle parameter, an array of three strings (like ["ha", "ca", "ha"])
 *  is returned. For a dihedral or improper parameter, an array of four strings is returned.
 *  The atom type may be "X", which is a wildcard that matches any atom type.
 */
static VALUE s_ParameterRef_GetAtomTypes(VALUE self) {
	UnionPar *up;
	Int tp, i, n;
	UInt types[4];
	VALUE vals[4];
	up = s_UnionParFromValue(self, &tp, 0);
	n = ParameterGetAtomTypes(tp, up, types);
	if (n == 0)
		rb_raise(rb_eMolbyError, "invalid member atom_types");
	for (i = 0; i < n; i++) {
		if (types[i] >= 0 && types[i] < kAtomTypeMinimum)
			vals[i] = INT2NUM(types[i]);
		else
			vals[i] = rb_str_new2(AtomTypeDecodeToString(types[i], NULL));
	}
	if (n == 1)
		return vals[0];
	else
		return rb_ary_new4(n, vals);
}

/*
 *  call-seq:
 *     k -> Float
 *
 *  Get the force constant. Available for bond, angle, dihedral, and improper parameters.
 */
static VALUE s_ParameterRef_GetK(VALUE self) {
	UnionPar *up;
	Int tp, i, n;
	VALUE vals[3];
	up = s_UnionParFromValue(self, &tp, 0);
	switch (tp) {
		case kBondParType:
			return rb_float_new(up->bond.k * INTERNAL2KCAL);
		case kAngleParType:
			return rb_float_new(up->angle.k * INTERNAL2KCAL);
		case kDihedralParType:
		case kImproperParType:
			if (up->torsion.mult == 1)
				return rb_float_new(up->torsion.k[0] * INTERNAL2KCAL);
			n = up->torsion.mult;
			if (n > 3)
				n = 3;
			for (i = 0; i < n; i++)
				vals[i] = rb_float_new(up->torsion.k[i] * INTERNAL2KCAL);
			return rb_ary_new4(n, vals);
		default:
			rb_raise(rb_eMolbyError, "invalid member k");
	}
}

/*
 *  call-seq:
 *     r0 -> Float
 *
 *  Get the equilibrium bond length. Only available for bond parameters.
 */
static VALUE s_ParameterRef_GetR0(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kBondParType)
		return rb_float_new(up->bond.r0);
	else rb_raise(rb_eMolbyError, "invalid member r0");
}

/*
 *  call-seq:
 *     a0 -> Float
 *
 *  Get the equilibrium angle (in degree). Only available for angle parameters.
 */
static VALUE s_ParameterRef_GetA0(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kAngleParType)
		return rb_float_new(up->angle.a0 * kRad2Deg);
	else rb_raise(rb_eMolbyError, "invalid member a0");
}

/*
 *  call-seq:
 *     mult -> Float
 *
 *  Get the multiplicity. Only available for dihedral and improper parameters.
 *  (Note: Implementation of multiple dihedral/improper parameters is not well tested)
 */
static VALUE s_ParameterRef_GetMult(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kDihedralParType || tp == kImproperParType)
		return rb_float_new(up->torsion.mult);
	else rb_raise(rb_eMolbyError, "invalid member mult");
}

/*
 *  call-seq:
 *     period -> Integer or Array of Integers
 *
 *  Get the periodicity. Only available for dihedral and improper parameters.
 *  If the multiplicity is larger than 1, then an array of integers is returned. 
 *  (Note: Implementation of multiple dihedral/improper parameters is not well tested)
 */
static VALUE s_ParameterRef_GetPeriod(VALUE self) {
	UnionPar *up;
	Int tp, i, n;
	VALUE vals[3];
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kDihedralParType || tp == kImproperParType) {
		if (up->torsion.mult == 1)
			return INT2NUM(up->torsion.period[0]);
		n = up->torsion.mult;
		if (n > 3)
			n = 3;
		for (i = 0; i < n; i++)
			vals[i] = INT2NUM(up->torsion.period[i]);
		return rb_ary_new4(n, vals);
	} else rb_raise(rb_eMolbyError, "invalid member period");
}

/*
 *  call-seq:
 *     phi0 -> Float or Array of Floats
 *
 *  Get the equilibrium dihedral angle. Only available for dihedral and improper parameters.
 *  If the multiplicity is larger than 1, then an array of floats is returned. 
 *  (Note: Implementation of multiple dihedral/improper parameters is not well tested)
 */
static VALUE s_ParameterRef_GetPhi0(VALUE self) {
	UnionPar *up;
	Int tp, i, n;
	VALUE vals[3];
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kDihedralParType || tp == kImproperParType) {
		if (up->torsion.mult == 1)
			return rb_float_new(up->torsion.phi0[0] * kRad2Deg);
		n = up->torsion.mult;
		if (n > 3)
			n = 3;
		for (i = 0; i < n; i++)
			vals[i] = rb_float_new(up->torsion.phi0[i] * kRad2Deg);
		return rb_ary_new4(n, vals);
	} else rb_raise(rb_eMolbyError, "invalid member phi0");
}

/*
 *  call-seq:
 *     A -> Float
 *
 *  Get the "A" value for the van der Waals parameter.
 */
/*
 static VALUE s_ParameterRef_GetA(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kVdwParType)
		return rb_float_new(up->vdw.A);
	else if (tp == kVdwPairParType)
		return rb_float_new(up->vdwp.A);
	else rb_raise(rb_eMolbyError, "invalid member A");
}
*/

/*
 *  call-seq:
 *     B -> Float
 *
 *  Get the "B" value for the van der Waals parameter.
 */
/*
static VALUE s_ParameterRef_GetB(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kVdwParType)
		return rb_float_new(up->vdw.B);
	else if (tp == kVdwPairParType)
		return rb_float_new(up->vdwp.B);
	else rb_raise(rb_eMolbyError, "invalid member B");
}
*/

/*
 *  call-seq:
 *     r_eq -> Float
 *
 *  Get the equilibrium radius (half of the minimum energy distance) for the van der Waals parameter.
 */
static VALUE s_ParameterRef_GetReq(VALUE self) {
	UnionPar *up;
	Int tp;
/*	Double a, b, r; */
	Double r;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kVdwParType) {
	/*	a = up->vdw.A;
		b = up->vdw.B;  */
		r = up->vdw.r_eq;
	} else if (tp == kVdwPairParType) {
	/*	a = up->vdwp.A;
		b = up->vdwp.B;  */
		r = up->vdwp.r_eq;
	} else rb_raise(rb_eMolbyError, "invalid member r_eq");
/*	if (a == 0.0 || b == 0.0) */
	return rb_float_new(r);
/*	else return rb_float_new(pow(2*a/b, 1.0/6.0)); */
}

/*
 *  call-seq:
 *     eps -> Float
 *
 *  Get the minimum energy for the van der Waals parameter.
 */
static VALUE s_ParameterRef_GetEps(VALUE self) {
	UnionPar *up;
	Int tp;
/*	Double a, b; */
	Double eps;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kVdwParType) {
	/*	a = up->vdw.A;
		b = up->vdw.B;  */
		eps = up->vdw.eps;
	} else if (tp == kVdwPairParType) {
	/*	a = up->vdwp.A;
		b = up->vdwp.B; */
		eps = up->vdwp.eps;
	} else rb_raise(rb_eMolbyError, "invalid member eps");
/*	if (a == 0.0 || b == 0.0)  */
		return rb_float_new(eps * INTERNAL2KCAL);
/*	else return rb_float_new(b*b/a/4.0 * INTERNAL2KCAL);  */
}

/*
 *  call-seq:
 *     A14 -> Float
 *
 *  Get the "A" value for the 1-4 van der Waals parameter.
 */
/*
static VALUE s_ParameterRef_GetA14(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kVdwParType)
		return rb_float_new(up->vdw.A14);
	else if (tp == kVdwPairParType)
		return rb_float_new(up->vdwp.A14);
	else rb_raise(rb_eMolbyError, "invalid member A14");
}
*/

/*
 *  call-seq:
 *     B14 -> Float
 *
 *  Get the "B" value for the 1-4 van der Waals parameter.
 */
/*
static VALUE s_ParameterRef_GetB14(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kVdwParType)
		return rb_float_new(up->vdw.B14);
	else if (tp == kVdwPairParType)
		return rb_float_new(up->vdwp.B14);
	else rb_raise(rb_eMolbyError, "invalid member B14");
}
*/

/*
 *  call-seq:
 *     r_eq14 -> Float
 *
 *  Get the equilibrium radius (half of the minimum energy distance) for the 1-4 van der Waals parameter.
 */
static VALUE s_ParameterRef_GetReq14(VALUE self) {
	UnionPar *up;
	Int tp;
/*	Double a, b, r; */
	Double r;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kVdwParType) {
	/*	a = up->vdw.A14;
		b = up->vdw.B14; */
		r = up->vdw.r_eq14;
	} else if (tp == kVdwPairParType) {
	/*	a = up->vdwp.A14;
		b = up->vdwp.B14;  */
		r = up->vdwp.r_eq14;
	} else rb_raise(rb_eMolbyError, "invalid member r_eq14");
/*	if (a == 0.0 || b == 0.0)  */
	return rb_float_new(r);
/*	else return rb_float_new(pow(2*a/b, 1.0/6.0));  */
}

/*
 *  call-seq:
 *     eps14 -> Float
 *
 *  Get the minimum energy for the 1-4 van der Waals parameter.
 */
static VALUE s_ParameterRef_GetEps14(VALUE self) {
	UnionPar *up;
	Int tp;
/*	Double a, b;  */
	Double eps;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kVdwParType) {
	/*	a = up->vdw.A14;
		b = up->vdw.B14;  */
		eps = up->vdw.eps14;
	} else if (tp == kVdwPairParType) {
	/*	a = up->vdwp.A14;
		b = up->vdwp.B14; */
		eps = up->vdwp.eps14;
	} else rb_raise(rb_eMolbyError, "invalid member eps14");
/*	if (a == 0.0 || b == 0.0) */
	return rb_float_new(eps * INTERNAL2KCAL);
/*	else return rb_float_new(b*b/a/4.0 * INTERNAL2KCAL);  */
}

/*
 *  call-seq:
 *     cutoff -> Float
 *
 *  Get the cutoff distance for the van der Waals pair-specific cutoff parameter.
 */
static VALUE s_ParameterRef_GetCutoff(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kVdwCutoffParType)
		return rb_float_new(up->vdwcutoff.cutoff);
	else rb_raise(rb_eMolbyError, "invalid member cutoff");
}

/*
 *  call-seq:
 *     radius -> Float
 *
 *  Get the atomic radius for the atom display parameter.
 */
static VALUE s_ParameterRef_GetRadius(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kElementParType)
		return rb_float_new(up->atom.radius);
	else rb_raise(rb_eMolbyError, "invalid member radius");
}

/*
 *  call-seq:
 *     color -> [Float, Float, Float]
 *
 *  Get the rgb color for the atom display parameter.
 */
static VALUE s_ParameterRef_GetColor(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kElementParType)
		return rb_ary_new3(3, rb_float_new(up->atom.r), rb_float_new(up->atom.g), rb_float_new(up->atom.b));
	else rb_raise(rb_eMolbyError, "invalid member color");
}

/*
 *  call-seq:
 *     atomic_number -> Integer
 *
 *  Get the atomic number for the vdw or atom parameter.
 */
static VALUE s_ParameterRef_GetAtomicNumber(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kElementParType)
		return INT2NUM(up->atom.number);
	else if (tp == kVdwParType)
		return INT2NUM(up->vdw.atomicNumber);
	else rb_raise(rb_eMolbyError, "invalid member atomic_number");
}

/*
 *  call-seq:
 *     name -> String
 *
 *  Get the name for the atom display parameter.
 */
static VALUE s_ParameterRef_GetName(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kElementParType) {
		char name[5];
		strncpy(name, up->atom.name, 4);
		name[4] = 0;
		return rb_str_new2(name);
	} else rb_raise(rb_eMolbyError, "invalid member name");
}

/*
 *  call-seq:
 *     weight -> Float
 *
 *  Get the atomic weight for the atom display parameter.
 */
static VALUE s_ParameterRef_GetWeight(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kElementParType)
		return rb_float_new(up->atom.weight);
	else if (tp == kVdwParType)
		return rb_float_new(up->vdw.weight);
	else rb_raise(rb_eMolbyError, "invalid member weight");
}

/*
 *  call-seq:
 *     fullname -> String
 *
 *  Get the full name for the atom display parameter.
 */
static VALUE s_ParameterRef_GetFullName(VALUE self) {
	UnionPar *up;
	Int tp;
	up = s_UnionParFromValue(self, &tp, 0);
	if (tp == kElementParType) {
		char fullname[16];
		strncpy(fullname, up->atom.fullname, 15);
		fullname[15] = 0;
		return rb_str_new2(fullname);
	} else rb_raise(rb_eMolbyError, "invalid member fullname");
}

/*
 *  call-seq:
 *     comment -> String
 *
 *  Get the comment for the parameter.
 */
static VALUE s_ParameterRef_GetComment(VALUE self) {
	UnionPar *up;
	Int tp, com;
	up = s_UnionParFromValue(self, &tp, 0);
	com = up->bond.com;
	if (com == 0)
		return Qnil;
	else return rb_str_new2(ParameterGetComment(com));
}

/*
 *  call-seq:
 *     source -> String
 *
 *  Get the source string for the parameter. Returns false for undefined parameter,
 *  and nil for "local" parameter that is specific for the molecule.
 */
static VALUE s_ParameterRef_GetSource(VALUE self) {
	UnionPar *up;
	Int tp, src;
	up = s_UnionParFromValue(self, &tp, 0);
	src = up->bond.src;
	if (src < 0)
		return Qfalse;  /* undefined */
	else if (src == 0)
		return Qnil;  /*  local  */
	else return rb_str_new2(ParameterGetComment(src));
}

static void
s_ScanAtomTypes(VALUE val, Int n, UInt *types)
{
	VALUE *valp;
	int i;
	if (n == 1)
		valp = &val;
	else {
		if (rb_obj_is_kind_of(val, rb_cString)) {
			char *s = StringValuePtr(val);
			char *p;
			for (i = 0; i < n; i++) {
				char buf[40];
				int len;
				/*  Skip leading separaters  */
				while (*s == '-' || *s == ' ' || *s == '\t')
					s++;
				for (p = s; *p != 0; p++) {
					if (*p == '-' || *p == ' ' || *p == '\t')
						break;
				}
				len = p - s;
				if (len >= sizeof(buf))
					len = sizeof(buf) - 1;
				strncpy(buf, s, len);
				buf[len] = 0;
				/*  Skip trailing blanks  */
				while (--len >= 0 && (buf[len] == ' ' || buf[len] == '\t'))
					buf[len] = 0;
				if (buf[0] == 0)
					rb_raise(rb_eMolbyError, "Bad atom type specification: %s", StringValuePtr(val));
				if (buf[0] >= '0' && buf[0] <= '9')
					types[i] = atoi(buf);
				else
					types[i] = AtomTypeEncodeToUInt(buf);
				if (p == NULL || *p == 0) {
					i++;
					break;
				} else s = p + 1;
			}
			if (i < n)
				rb_raise(rb_eMolbyError, "%d atom types are required but only %d are given; %s", n, i, StringValuePtr(val));
			return;
		}
		val = rb_ary_to_ary(val);
		if (RARRAY_LEN(val) != n)
			rb_raise(rb_eMolbyError, "an array of %d atom types is required", n);
		valp = RARRAY_PTR(val);
	}
	for (i = 0; i < n; i++) {
		if (rb_obj_is_kind_of(valp[i], rb_cNumeric))
			types[i] = NUM2INT(rb_Integer(valp[i]));
		else {
			VALUE sval = valp[i];
			types[i] = AtomTypeEncodeToUInt(StringValuePtr(sval));
		}
	}
}

static VALUE s_ParameterRef_SetAtomTypes(VALUE self, VALUE val) {
	UnionPar *up;
	VALUE oldval;
	Int oldsrc, tp;
	UInt types[4];
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetAtomTypes(self);
	oldsrc = up->bond.src;
	switch (tp) {
		case kBondParType:
			s_ScanAtomTypes(val, 2, types);
			up->bond.type1 = types[0];
			up->bond.type2 = types[1];
			break;
		case kAngleParType:
			s_ScanAtomTypes(val, 3, types);
			up->angle.type1 = types[0];
			up->angle.type2 = types[1];
			up->angle.type3 = types[2];
			break;
		case kDihedralParType:
		case kImproperParType:
			s_ScanAtomTypes(val, 4, types);
			up->torsion.type1 = types[0];
			up->torsion.type2 = types[1];
			up->torsion.type3 = types[2];
			up->torsion.type4 = types[3];
			break;
		case kVdwParType:
			s_ScanAtomTypes(val, 1, types);
			up->vdw.type1 = types[0];
			break;
		case kVdwPairParType:
			s_ScanAtomTypes(val, 2, types);
			up->vdwp.type1 = types[0];
			up->vdwp.type2 = types[1];
			break;
		case kVdwCutoffParType:
			s_ScanAtomTypes(val, 2, types);
			up->vdwcutoff.type1 = types[0];
			up->vdwcutoff.type2 = types[1];
			break;
		default:
			return Qnil;
	}
	s_RegisterUndoForParameterAttrChange(self, s_AtomTypesSym, val, oldval, oldsrc);
	return val;
}

static VALUE s_ParameterRef_SetK(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, i, n, oldsrc;
	VALUE *valp, oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetK(self);
	oldsrc = up->bond.src;
	switch (tp) {
		case kBondParType:
			val = rb_Float(val);
			up->bond.k = NUM2DBL(val) * KCAL2INTERNAL;
			break;
		case kAngleParType:
			val = rb_Float(val);
			up->angle.k = NUM2DBL(val) * KCAL2INTERNAL;
			break;
		case kDihedralParType:
		case kImproperParType:
			if (up->torsion.mult == 1 || up->torsion.mult == 0) {
				up->torsion.mult = 1;
				val = rb_Float(val);
				up->torsion.k[0] = NUM2DBL(val) * KCAL2INTERNAL;
				break;
			}
			n = up->torsion.mult;
			if (n > 3)
				n = 3;
			val = rb_ary_to_ary(val);
			if (RARRAY_LEN(val) != n)
				rb_raise(rb_eMolbyError, "the value should be an array of %d floats", n);
			valp = RARRAY_PTR(val);
			for (i = 0; i < n; i++) {
				up->torsion.k[i] = NUM2DBL(rb_Float(valp[i])) * KCAL2INTERNAL;
			}
			break;
		default:
			rb_raise(rb_eMolbyError, "invalid member k");
	}
	s_RegisterUndoForParameterAttrChange(self, s_KSym, val, oldval, oldsrc);
	return val;
}

static VALUE s_ParameterRef_SetR0(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetR0(self);
	oldsrc = up->bond.src;
	if (tp == kBondParType) {
		val = rb_Float(val);
		up->bond.r0 = NUM2DBL(val);
	} else rb_raise(rb_eMolbyError, "invalid member r0");
	s_RegisterUndoForParameterAttrChange(self, s_R0Sym, val, oldval, oldsrc);
	return val;
}

static VALUE s_ParameterRef_SetA0(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetA0(self);
	oldsrc = up->bond.src;
	if (tp == kAngleParType) {
		val = rb_Float(val);
		up->angle.a0 = NUM2DBL(val) * kDeg2Rad;
	} else rb_raise(rb_eMolbyError, "invalid member a0");
	s_RegisterUndoForParameterAttrChange(self, s_A0Sym, val, oldval, oldsrc);
	return val;
}

static VALUE s_ParameterRef_SetMult(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetMult(self);
	oldsrc = up->bond.src;
	if (tp == kDihedralParType || tp == kImproperParType) {
		int i;
		val = rb_Integer(val);
		i = NUM2INT(val);
		if (i < 0 || i > 3)
			rb_raise(rb_eMolbyError, "torsion multiplicity should be 0..3");
		up->torsion.mult = i;
	} else rb_raise(rb_eMolbyError, "invalid member mult");
	s_RegisterUndoForParameterAttrChange(self, s_MultSym, val, oldval, oldsrc);
	return val;
}

static VALUE s_ParameterRef_SetPeriod(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, i, n, oldsrc;
	VALUE *valp, oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetPeriod(self);
	oldsrc = up->bond.src;
	if (tp == kDihedralParType || tp == kImproperParType) {
		if (up->torsion.mult == 1 || up->torsion.mult == 0) {
			up->torsion.mult = 1;
			val = rb_Integer(val);
			up->torsion.period[0] = NUM2INT(val);
		} else {
			n = up->torsion.mult;
			if (n > 3)
				n = 3;
			val = rb_ary_to_ary(val);
			if (RARRAY_LEN(val) != n)
				rb_raise(rb_eMolbyError, "the value should be an array of %d integers", n);
			valp = RARRAY_PTR(val);
			for (i = 0; i < n; i++) {
				up->torsion.period[i] = NUM2INT(rb_Integer(valp[i]));
			}
		}
	} else rb_raise(rb_eMolbyError, "invalid member period");
	s_RegisterUndoForParameterAttrChange(self, s_PeriodSym, val, oldval, oldsrc);
	return val;
}

static VALUE s_ParameterRef_SetPhi0(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, i, n, oldsrc;
	VALUE *valp, oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetPhi0(self);
	oldsrc = up->bond.src;
	if (tp == kDihedralParType || tp == kImproperParType) {
		if (up->torsion.mult == 1 || up->torsion.mult == 0) {
			up->torsion.mult = 1;
			val = rb_Float(val);
			up->torsion.phi0[0] = NUM2DBL(val) * kDeg2Rad;
		} else {
			n = up->torsion.mult;
			if (n > 3)
				n = 3;
			val = rb_ary_to_ary(val);
			if (RARRAY_LEN(val) != n)
				rb_raise(rb_eMolbyError, "the value should be an array of %d floats", n);
			valp = RARRAY_PTR(val);
			for (i = 0; i < n; i++)
				up->torsion.phi0[i] = NUM2DBL(rb_Float(valp[i])) * kDeg2Rad;
		}
	} else rb_raise(rb_eMolbyError, "invalid member phi0");
	s_RegisterUndoForParameterAttrChange(self, s_Phi0Sym, val, oldval, oldsrc);
	return val;
}

/*
static VALUE s_ParameterRef_SetA(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	double d;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetA(self);
	oldsrc = up->bond.src;
	val = rb_Float(val);
	d = NUM2DBL(val);
	if (tp == kVdwParType)
		up->vdw.A = d;
	else if (tp == kVdwPairParType)
		up->vdwp.A = d;
	else rb_raise(rb_eMolbyError, "invalid member A");
	s_RegisterUndoForParameterAttrChange(self, s_ASym, val, oldval, oldsrc);
	return val;
}

static VALUE s_ParameterRef_SetB(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	double d;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetB(self);
	oldsrc = up->bond.src;
	val = rb_Float(val);
	d = NUM2DBL(val);
	if (tp == kVdwParType)
		up->vdw.B = d;
	else if (tp == kVdwPairParType)
		up->vdwp.B = d;
	else rb_raise(rb_eMolbyError, "invalid member B");
	s_RegisterUndoForParameterAttrChange(self, s_BSym, val, oldval, oldsrc);
	return val;
}
*/

static VALUE s_ParameterRef_SetReq(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	Double r;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetReq(self);
	oldsrc = up->bond.src;
	val = rb_Float(val);
	r = NUM2DBL(val);
	if (tp == kVdwParType) {
		up->vdw.r_eq = r;
		up->vdw.A = pow(r * 2, 12.0) * up->vdw.eps;
		up->vdw.B = 2 * pow(r * 2, 6.0) * up->vdw.eps;
	} else if (tp == kVdwPairParType) {
		up->vdwp.r_eq = r;
		up->vdwp.A = pow(r * 2, 12.0) * up->vdwp.eps;
		up->vdwp.B = 2 * pow(r * 2, 6.0) * up->vdwp.eps;
	} else rb_raise(rb_eMolbyError, "invalid member r_eq");
	s_RegisterUndoForParameterAttrChange(self, s_ReqSym, val, oldval, oldsrc);
	return val;
}

static VALUE s_ParameterRef_SetEps(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	Double e;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetEps(self);
	oldsrc = up->bond.src;
	val = rb_Float(val);
	e = NUM2DBL(val) * KCAL2INTERNAL;
	if (tp == kVdwParType) {
		up->vdw.eps = e;
		up->vdw.A = pow(up->vdw.r_eq * 2, 12.0) * up->vdw.eps;
		up->vdw.B = 2 * pow(up->vdw.r_eq * 2, 6.0) * up->vdw.eps;
	} else if (tp == kVdwPairParType) {
		up->vdwp.eps = e;
		up->vdwp.A = pow(up->vdwp.r_eq * 2, 12.0) * up->vdwp.eps;
		up->vdwp.B = 2 * pow(up->vdwp.r_eq * 2, 6.0) * up->vdwp.eps;
	} else rb_raise(rb_eMolbyError, "invalid member eps");
	s_RegisterUndoForParameterAttrChange(self, s_EpsSym, val, oldval, oldsrc);
	return val;
}

/*
static VALUE s_ParameterRef_SetA14(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	double d;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetA14(self);
	oldsrc = up->bond.src;
	val = rb_Float(val);
	d = NUM2DBL(val);
	if (tp == kVdwParType)
		up->vdw.A14 = d;
	else if (tp == kVdwPairParType)
		up->vdwp.A14 = d;
	else rb_raise(rb_eMolbyError, "invalid member A14");
	s_RegisterUndoForParameterAttrChange(self, s_A14Sym, val, oldval, oldsrc);	
	return val;
}

static VALUE s_ParameterRef_SetB14(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	double d;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetB14(self);
	oldsrc = up->bond.src;
	val = rb_Float(val);
	d = NUM2DBL(val);
	if (tp == kVdwParType)
		up->vdw.B14 = d;
	else if (tp == kVdwPairParType)
		up->vdwp.B14 = d;
	else rb_raise(rb_eMolbyError, "invalid member B14");
	s_RegisterUndoForParameterAttrChange(self, s_B14Sym, val, oldval, oldsrc);	
	return val;
}
*/

static VALUE s_ParameterRef_SetReq14(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	Double r;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetReq14(self);
	oldsrc = up->bond.src;
	val = rb_Float(val);
	r = NUM2DBL(val);
	if (tp == kVdwParType) {
		up->vdw.r_eq14 = r;
		up->vdw.A14 = pow(up->vdw.r_eq14 * 2, 12.0) * up->vdw.eps14;
		up->vdw.B14 = 2 * pow(up->vdw.r_eq14 * 2, 6.0) * up->vdw.eps14;
	} else if (tp == kVdwPairParType) {
		up->vdwp.r_eq14 = r;
		up->vdwp.A14 = pow(up->vdwp.r_eq14 * 2, 12.0) * up->vdwp.eps14;
		up->vdwp.B14 = 2 * pow(up->vdwp.r_eq14 * 2, 6.0) * up->vdwp.eps14;
	} else rb_raise(rb_eMolbyError, "invalid member r_eq14");
	s_RegisterUndoForParameterAttrChange(self, s_Req14Sym, val, oldval, oldsrc);	
	return val;
}

static VALUE s_ParameterRef_SetEps14(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	Double e;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetEps14(self);
	oldsrc = up->bond.src;
	val = rb_Float(val);
	e = NUM2DBL(val) * KCAL2INTERNAL;
	if (tp == kVdwParType) {
		up->vdw.eps14 = e;
		up->vdw.A14 = pow(up->vdw.r_eq14 * 2, 12.0) * up->vdw.eps14;
		up->vdw.B14 = 2 * pow(up->vdw.r_eq14 * 2, 6.0) * up->vdw.eps14;
	} else if (tp == kVdwPairParType) {
		up->vdwp.eps14 = e;
		up->vdwp.A14 = pow(up->vdwp.r_eq14 * 2, 12.0) * up->vdwp.eps14;
		up->vdwp.B14 = 2 * pow(up->vdwp.r_eq14 * 2, 6.0) * up->vdwp.eps14;
	} else rb_raise(rb_eMolbyError, "invalid member eps14");
	s_RegisterUndoForParameterAttrChange(self, s_Eps14Sym, val, oldval, oldsrc);	
	return val;
}

static VALUE s_ParameterRef_SetCutoff(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE oldval;
	oldval = s_ParameterRef_GetCutoff(self);
	oldsrc = up->bond.src;
	up = s_UnionParFromValue(self, &tp, 1);
	val = rb_Float(val);
	if (tp == kVdwCutoffParType) {
		up->vdwcutoff.cutoff = NUM2DBL(val);
	} else rb_raise(rb_eMolbyError, "invalid member cutoff");
	s_RegisterUndoForParameterAttrChange(self, s_CutoffSym, val, oldval, oldsrc);	
	return val;
}

static VALUE s_ParameterRef_SetRadius(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetRadius(self);
	oldsrc = up->bond.src;
	val = rb_Float(val);
	if (tp == kElementParType) {
		up->atom.radius = NUM2DBL(val);
	} else rb_raise(rb_eMolbyError, "invalid member radius");
	s_RegisterUndoForParameterAttrChange(self, s_RadiusSym, val, oldval, oldsrc);	
	return val;
}

static VALUE s_ParameterRef_SetColor(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE *valp, oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetColor(self);
	oldsrc = up->bond.src;
	val = rb_ary_to_ary(val);
	if (RARRAY_LEN(val) != 3)
		rb_raise(rb_eMolbyError, "value should be an array of three floats (r, g, b)");
	valp = RARRAY_PTR(val);
	if (tp == kElementParType) {
		up->atom.r = NUM2DBL(rb_Float(valp[0]));
		up->atom.g = NUM2DBL(rb_Float(valp[1]));
		up->atom.b = NUM2DBL(rb_Float(valp[2]));
	} else rb_raise(rb_eMolbyError, "invalid member color");
	s_RegisterUndoForParameterAttrChange(self, s_ColorSym, val, oldval, oldsrc);	
	return val;
}

static VALUE s_ParameterRef_SetAtomicNumber(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetAtomicNumber(self);
	oldsrc = up->bond.src;
	val = rb_Integer(val);
	if (tp == kElementParType)
		up->atom.number = NUM2INT(val);
	else if (tp == kVdwParType) {
		up->vdw.atomicNumber = NUM2INT(val);
		up->vdw.weight = WeightForAtomicNumber(up->vdw.atomicNumber);
	} else rb_raise(rb_eMolbyError, "invalid member atomic_number");
	s_RegisterUndoForParameterAttrChange(self, s_AtomicNumberSym, val, oldval, oldsrc);	
	return val;
}

static VALUE s_ParameterRef_SetName(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetName(self);
	oldsrc = up->bond.src;
	if (tp == kElementParType) {
		strncpy(up->atom.name, StringValuePtr(val), 4);
	} else rb_raise(rb_eMolbyError, "invalid member name");
	s_RegisterUndoForParameterAttrChange(self, s_NameSym, val, oldval, oldsrc);	
	return val;
}

static VALUE s_ParameterRef_SetWeight(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE oldval;
	val = rb_Float(val);
	oldval = s_ParameterRef_GetWeight(self);
	up = s_UnionParFromValue(self, &tp, 1);
	oldsrc = up->bond.src;
	if (tp == kElementParType)
		up->atom.weight = NUM2DBL(val);
	else if (tp == kVdwParType)
		up->vdw.weight = NUM2DBL(val);
	else rb_raise(rb_eMolbyError, "invalid member weight");
	s_RegisterUndoForParameterAttrChange(self, s_WeightSym, val, oldval, oldsrc);	
	return val;
}

static VALUE s_ParameterRef_SetFullName(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetFullName(self);
	oldsrc = up->bond.src;
	if (tp == kElementParType) {
		strncpy(up->atom.fullname, StringValuePtr(val), 15);
		up->atom.fullname[15] = 0;
	} else rb_raise(rb_eMolbyError, "invalid member fullname");
	s_RegisterUndoForParameterAttrChange(self, s_FullNameSym, val, oldval, oldsrc);	
	return val;
}

static VALUE s_ParameterRef_SetComment(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, com, oldsrc;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	oldval = s_ParameterRef_GetComment(self);
	oldsrc = up->bond.src;
	if (val == Qnil)
		up->bond.com = 0;
	else {
		com = ParameterCommentIndex(StringValuePtr(val));
		up->bond.com = com;
	}
	s_RegisterUndoForParameterAttrChange(self, s_CommentSym, val, oldval, oldsrc);	
	return val;	
}

/*  Only false (undefined) and nil (local) can be set  */
static VALUE s_ParameterRef_SetSource(VALUE self, VALUE val) {
	UnionPar *up;
	Int tp, oldsrc;
	VALUE oldval;
	up = s_UnionParFromValue(self, &tp, 1);
	if (val != Qfalse && val != Qnil)
		rb_raise(rb_eMolbyError, "set source: only false (undefined parameter) or nil (local parameter) is allowed");
	oldval = s_ParameterRef_GetSource(self);
	oldsrc = up->bond.src;
	if (oldsrc != 0 && oldsrc != -1)
		rb_raise(rb_eMolbyError, "source information of global parameter cannot be modified");
	up->bond.src = (val == Qfalse ? -1 : 0);
	s_RegisterUndoForParameterAttrChange(self, s_SourceSym, val, oldval, oldsrc);	
	return val;	
}

static struct s_ParameterAttrDef {
	char *name;
	VALUE *symref;  /*  Address of s_IndexSymbol etc. */
	ID id;			/*  Will be set within InitMolby()  */
	VALUE (*getter)(VALUE);
	VALUE (*setter)(VALUE, VALUE);
} s_ParameterAttrDefTable[] = {
	{"index",        &s_IndexSym,        0, s_ParameterRef_GetIndex,        NULL},
	{"par_type",     &s_ParTypeSym,      0, s_ParameterRef_GetParType,      NULL},
	{"atom_types",   &s_AtomTypesSym,    0, s_ParameterRef_GetAtomTypes,    s_ParameterRef_SetAtomTypes},
	{"atom_type",    &s_AtomTypeSym,     0, s_ParameterRef_GetAtomTypes,    s_ParameterRef_SetAtomTypes},
	{"k",            &s_KSym,            0, s_ParameterRef_GetK,            s_ParameterRef_SetK},
	{"r0",           &s_R0Sym,           0, s_ParameterRef_GetR0,           s_ParameterRef_SetR0},
	{"a0",           &s_A0Sym,           0, s_ParameterRef_GetA0,           s_ParameterRef_SetA0},
	{"mult",         &s_MultSym,         0, s_ParameterRef_GetMult,         s_ParameterRef_SetMult},
	{"period",       &s_PeriodSym,       0, s_ParameterRef_GetPeriod,       s_ParameterRef_SetPeriod},
	{"phi0",         &s_Phi0Sym,         0, s_ParameterRef_GetPhi0,         s_ParameterRef_SetPhi0},
/*	{"A",            &s_ASym,            0, s_ParameterRef_GetA,            NULL},
	{"B",            &s_BSym,            0, s_ParameterRef_GetB,            NULL}, */
	{"r_eq",         &s_ReqSym,          0, s_ParameterRef_GetReq,          s_ParameterRef_SetReq},
	{"eps",          &s_EpsSym,          0, s_ParameterRef_GetEps,          s_ParameterRef_SetEps},
/*	{"A14",          &s_A14Sym,          0, s_ParameterRef_GetA14,          NULL},
	{"B14",          &s_B14Sym,          0, s_ParameterRef_GetB14,          NULL}, */
	{"r_eq14",       &s_Req14Sym,        0, s_ParameterRef_GetReq14,        s_ParameterRef_SetReq14},
	{"eps14",        &s_Eps14Sym,        0, s_ParameterRef_GetEps14,        s_ParameterRef_SetEps14},
	{"cutoff",       &s_CutoffSym,       0, s_ParameterRef_GetCutoff,       s_ParameterRef_SetCutoff},
	{"radius",       &s_RadiusSym,       0, s_ParameterRef_GetRadius,       s_ParameterRef_SetRadius},
	{"color",        &s_ColorSym,        0, s_ParameterRef_GetColor,        s_ParameterRef_SetColor},
	{"atomic_number",&s_AtomicNumberSym, 0, s_ParameterRef_GetAtomicNumber, s_ParameterRef_SetAtomicNumber},
	{"name",         &s_NameSym,         0, s_ParameterRef_GetName,         s_ParameterRef_SetName},
	{"weight",       &s_WeightSym,       0, s_ParameterRef_GetWeight,       s_ParameterRef_SetWeight},
	{"fullname",     &s_FullNameSym,     0, s_ParameterRef_GetFullName,     s_ParameterRef_SetFullName},
	{"comment",      &s_CommentSym,      0, s_ParameterRef_GetComment,      s_ParameterRef_SetComment},
	{"source",       &s_SourceSym,       0, s_ParameterRef_GetSource,       s_ParameterRef_SetSource},
	{NULL} /* Sentinel */
};

static VALUE
s_ParameterRef_SetAttr(VALUE self, VALUE key, VALUE value)
{
	int i;
	ID kid;
	if (TYPE(key) != T_SYMBOL) {
		kid = rb_intern(StringValuePtr(key));
		key = ID2SYM(kid);
	} else kid = SYM2ID(key);
	for (i = 0; s_ParameterAttrDefTable[i].name != NULL; i++) {
		if (s_ParameterAttrDefTable[i].id == kid) {
			if (value == Qundef)
				return (*(s_ParameterAttrDefTable[i].getter))(self);
			else if (s_ParameterAttrDefTable[i].setter == NULL)
				rb_raise(rb_eMolbyError, "the attribute \"%s\" is read-only", rb_id2name(kid));
			else
				return (*(s_ParameterAttrDefTable[i].setter))(self, value);
		}
	}
	rb_raise(rb_eMolbyError, "unknown parameter attribute \"%s\"", rb_id2name(kid));
	return Qnil; /* not reached */
}

static VALUE
s_ParameterRef_GetAttr(VALUE self, VALUE key)
{
	return s_ParameterRef_SetAttr(self, key, Qundef);
}

/*
 *  call-seq:
 *     keys(idx)          -> array of valid parameter attributes
 *  
 *  Returns an array of valid parameter attributes (as Symbols).
 */
static VALUE
s_ParameterRef_Keys(VALUE self)
{
	ParameterRef *pref;
	Data_Get_Struct(self, ParameterRef, pref);
	switch (pref->parType) {
		case kBondParType:
			return rb_ary_new3(7, s_IndexSym, s_ParTypeSym, s_AtomTypesSym, s_KSym, s_R0Sym, s_CommentSym, s_SourceSym);
		case kAngleParType:
			return rb_ary_new3(7, s_IndexSym, s_ParTypeSym, s_AtomTypesSym, s_KSym, s_A0Sym, s_CommentSym, s_SourceSym);
		case kDihedralParType:
		case kImproperParType:
			return rb_ary_new3(9, s_IndexSym, s_ParTypeSym, s_AtomTypesSym, s_MultSym, s_KSym, s_PeriodSym, s_Phi0Sym, s_CommentSym, s_SourceSym);
		case kVdwParType:
			return rb_ary_new3(11, s_IndexSym, s_ParTypeSym, s_AtomTypesSym, s_AtomicNumberSym, s_ReqSym, s_EpsSym, s_Req14Sym, s_Eps14Sym, s_WeightSym, s_CommentSym, s_SourceSym);
		case kVdwPairParType:
			return rb_ary_new3(9, s_IndexSym, s_ParTypeSym, s_AtomTypesSym, s_ReqSym, s_EpsSym, s_Req14Sym, s_Eps14Sym, s_CommentSym, s_SourceSym);
		case kVdwCutoffParType:
			return rb_ary_new3(6, s_IndexSym, s_ParTypeSym, s_AtomTypesSym, s_CutoffSym, s_CommentSym, s_SourceSym);
		case kElementParType:
			return rb_ary_new3(10, s_IndexSym, s_ParTypeSym, s_AtomicNumberSym, s_NameSym, s_FullNameSym, s_RadiusSym, s_ColorSym, s_WeightSym, s_CommentSym, s_SourceSym);
		default:
			rb_raise(rb_eMolbyError, "internal error: invalid parameter type");
	}
	return Qnil;  /*  Not reached  */
}

/*
 *  call-seq:
 *     to_hash(idx)          -> Hash
 *  
 *  Returns a hash containing valid parameter names and values
 */
static VALUE
s_ParameterRef_ToHash(VALUE self)
{
	VALUE keys = s_ParameterRef_Keys(self);
	VALUE retval;
	int i;
	if (keys == Qnil)
		return Qnil;
	retval = rb_hash_new();
	for (i = 0; i < RARRAY_LEN(keys); i++) {
		VALUE key = RARRAY_PTR(keys)[i];
		VALUE val = s_ParameterRef_GetAttr(self, key);
		rb_hash_aset(retval, key, val);
	}
	return retval;
}

/*
 *  call-seq:
 *     parameter.to_s(idx)          -> String
 *  
 *  Returns a string representation of the given parameter
 */
static VALUE
s_ParameterRef_ToString(VALUE self)
{
	Int tp, i, n;
	char buf[1024], types[4][8];
	UnionPar *up = s_UnionParFromValue(self, &tp, 0);
	switch (tp) {
		case kBondParType:
			snprintf(buf, sizeof buf, "bond %4.6s %4.6s %8.2f %8.3f", AtomTypeDecodeToString(up->bond.type1, types[0]), AtomTypeDecodeToString(up->bond.type2, types[1]), up->bond.k * INTERNAL2KCAL, up->bond.r0);
			break;
		case kAngleParType:
			snprintf(buf, sizeof buf, "angle %4.6s %4.6s %4.6s %8.2f %8.3f", AtomTypeDecodeToString(up->angle.type1, types[0]), AtomTypeDecodeToString(up->angle.type2, types[1]), AtomTypeDecodeToString(up->angle.type3, types[2]), up->angle.k * INTERNAL2KCAL, up->angle.a0 * kRad2Deg);
			break;
		case kDihedralParType:
		case kImproperParType:
			snprintf(buf, sizeof buf, "%s %4.6s %4.6s %4.6s %4.6s", (tp == kDihedralParType ? "dihe" : "impr"), AtomTypeDecodeToString(up->torsion.type1, types[0]), AtomTypeDecodeToString(up->torsion.type2, types[1]), AtomTypeDecodeToString(up->torsion.type3, types[2]), AtomTypeDecodeToString(up->torsion.type4, types[3]));
			n = strlen(buf);
			for (i = 0; i < up->torsion.mult; i++) {
				snprintf(buf + n, sizeof buf - n, " %8.2f %2d %8.3f", up->torsion.k[i] * INTERNAL2KCAL, up->torsion.period[i], up->torsion.phi0[i] * kRad2Deg);
				n = strlen(buf);
			}
			break;
		case kVdwParType:
			snprintf(buf, sizeof buf, "nonbonded %4.6s %8.4f %8.4f %8.4f %8.4f", AtomTypeDecodeToString(up->vdw.type1, types[0]), up->vdw.A * INTERNAL2KCAL / pow(up->vdw.r_eq, 12.0), up->vdw.r_eq / 1.12246204830937, up->vdw.A14 * INTERNAL2KCAL / pow(up->vdw.r_eq14, 12.0), up->vdw.r_eq14 / 1.12246204830937);
			break;
		case kVdwPairParType:
			snprintf(buf, sizeof buf, "nbfi %4.6s %4.6s %12.8e %12.8e %12.8e %12.8e", AtomTypeDecodeToString(up->vdwp.type1, types[0]), AtomTypeDecodeToString(up->vdwp.type2, types[1]), up->vdwp.A * INTERNAL2KCAL, up->vdwp.B * INTERNAL2KCAL, up->vdwp.A14 * INTERNAL2KCAL, up->vdwp.B14 * INTERNAL2KCAL);
			break;
		case kVdwCutoffParType:
			snprintf(buf, sizeof buf, "vdwcutoff %4.6s %4.6s %8.4f", AtomTypeDecodeToString(up->vdwcutoff.type1, types[0]), AtomTypeDecodeToString(up->vdwcutoff.type2, types[1]), up->vdwcutoff.cutoff);
			break;
		case kElementParType:
			snprintf(buf, sizeof buf, "element %2.2s %3d %6.3f %6.3f %6.3f %6.3f %8.4f %s", up->atom.name, up->atom.number, up->atom.radius, up->atom.r, up->atom.g, up->atom.b, up->atom.weight, up->atom.fullname);
			break;
	}
	return rb_str_new2(buf);
}

/*
 *  call-seq:
 *     self == parameterRef -> boolean
 *  
 *  True if the parameters point to the same parameter record.
 */
static VALUE
s_ParameterRef_Equal(VALUE self, VALUE val)
{
	Int tp1, tp2;
	if (rb_obj_is_kind_of(val, rb_cParameterRef)) {
		return (s_UnionParFromValue(self, &tp1, 0) == s_UnionParFromValue(val, &tp2, 0) ? Qtrue : Qfalse);
	} else return Qfalse;
}
	
#pragma mark ====== Parameter Class ======

/*  The Parameter class actually encapsulate Molecule record. If the record pointer
 *  is NULL, then the global parameters are looked for.  */

/*  Rebuild the MD parameter record if necessary: may throw an exception  */
/*  The second parameter is passed to md_arena.prepare; if true, then check only  */
static void
s_RebuildMDParameterIfNecessary(VALUE val, VALUE cval)
{
	Molecule *mol;
	Data_Get_Struct(val, Molecule, mol);
	if (mol == NULL)
		rb_raise(rb_eMolbyError, "the molecule is empty");
	if (mol->par == NULL || mol->arena == NULL || mol->arena->is_initialized == 0 || mol->needsMDRebuild) {
		/*  Do self.md_arena.prepare  */
		VALUE val2 = rb_funcall(val, rb_intern("md_arena"), 0);
		if (val2 != Qnil)
			val2 = rb_funcall(val2, rb_intern("prepare"), 1, cval);
	}
}

static VALUE
s_NewParameterValueFromValue(VALUE val)
{
	Molecule *mol;
	if (rb_obj_is_kind_of(val, rb_cMolecule)) {
		Data_Get_Struct(val, Molecule, mol);
		s_RebuildMDParameterIfNecessary(val, Qtrue);
		MoleculeRetain(mol);
		return Data_Wrap_Struct(rb_cParameter, 0, (void (*)(void *))MoleculeRelease, mol);
	} else {
		mol = NULL;
		return Data_Wrap_Struct(rb_cParameter, 0, NULL, mol);
	}
}

static Molecule *
s_MoleculeFromParameterValue(VALUE val)
{
	Molecule *mol;
	Data_Get_Struct(val, Molecule, mol);
	return mol;
}

static Parameter *
s_ParameterFromParameterValue(VALUE val)
{
	Molecule *mol;
	Data_Get_Struct(val, Molecule, mol);
	if (mol != NULL)
		return mol->par;
	return gBuiltinParameters;
}

/*  Forward declarations  */
static VALUE s_NewParEnumerableValueFromMoleculeAndType(Molecule *mol, Int parType);
static Molecule *s_MoleculeFromParEnumerableValue(VALUE val);

static Molecule *
s_MoleculeFromParameterOrParEnumerableValue(VALUE val)
{
	if (val == rb_cParameter) {
		return NULL;  /*  Parameter class method: builtin parameters  */
	} else if (rb_obj_is_kind_of(val, rb_cParameter)) {
		return s_MoleculeFromParameterValue(val);
	} else if (rb_obj_is_kind_of(val, rb_cParEnumerable)) {
		return s_MoleculeFromParEnumerableValue(val);
	} else return NULL;
}

/*
 *  call-seq:
 *     builtin    -> Parameter
 *  
 *  Returns a parameter value that points to the global (builtin) parameters.
 *  Equivalent to Parameter::Builtin (constant).
 */
static VALUE
s_Parameter_Builtin(VALUE self)
{
	static ID s_builtin_id = 0;
	if (s_builtin_id == 0)
		s_builtin_id = rb_intern("Builtin");
	return rb_const_get(rb_cParameter, s_builtin_id);
}

/*
 *  call-seq:
 *     bond(idx)          -> ParameterRef
 *  
 *  The index-th bond parameter record is returned.
 */
static VALUE
s_Parameter_Bond(VALUE self, VALUE ival)
{
	Molecule *mol;
	int idx, n;
	mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	idx = NUM2INT(rb_Integer(ival));
	if (mol == NULL)
		n = gBuiltinParameters->nbondPars;
	else if (mol->par != NULL)
		n = mol->par->nbondPars;
	else n = 0;
	if (idx < -n || idx >= n)
		rb_raise(rb_eMolbyError, "Bond parameter index (%d) out of range", idx);
	if (idx < 0)
		idx += n;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kBondParType, idx);
}

/*
 *  call-seq:
 *     angle(idx)          -> ParameterRef
 *  
 *  The index-th angle parameter record is returned.
 */
static VALUE
s_Parameter_Angle(VALUE self, VALUE ival)
{
	Molecule *mol;
	int idx, n;
	mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	idx = NUM2INT(rb_Integer(ival));
	if (mol == NULL)
		n = gBuiltinParameters->nanglePars;
	else if (mol->par != NULL)
		n = mol->par->nanglePars;
	else n = 0;
	if (idx < -n || idx >= n)
		rb_raise(rb_eMolbyError, "Angle parameter index (%d) out of range", idx);
	if (idx < 0)
		idx += n;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kAngleParType, idx);
}

/*
 *  call-seq:
 *     dihedral(idx)          -> ParameterRef
 *  
 *  The index-th dihedral parameter record is returned.
 */
static VALUE
s_Parameter_Dihedral(VALUE self, VALUE ival)
{
	Molecule *mol;
	int idx, n;
	mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	idx = NUM2INT(rb_Integer(ival));
	if (mol == NULL)
		n = gBuiltinParameters->ndihedralPars;
	else if (mol->par != NULL)
		n = mol->par->ndihedralPars;
	else n = 0;
	if (idx < -n || idx >= n)
		rb_raise(rb_eMolbyError, "Dihedral parameter index (%d) out of range", idx);
	if (idx < 0)
		idx += n;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kDihedralParType, idx);
}

/*
 *  call-seq:
 *     improper(idx)          -> ParameterRef
 *  
 *  The index-th improper parameter record is returned.
 */
static VALUE
s_Parameter_Improper(VALUE self, VALUE ival)
{
	Molecule *mol;
	int idx, n;
	mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	idx = NUM2INT(rb_Integer(ival));
	if (mol == NULL)
		n = gBuiltinParameters->nimproperPars;
	else if (mol->par != NULL)
		n = mol->par->nimproperPars;
	else n = 0;
	if (idx < -n || idx >= n)
		rb_raise(rb_eMolbyError, "Improper parameter index (%d) out of range", idx);
	if (idx < 0)
		idx += n;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kImproperParType, idx);
}

/*
 *  call-seq:
 *     vdw(idx)          -> ParameterRef
 *  
 *  The index-th vdw parameter record is returned.
 */
static VALUE
s_Parameter_Vdw(VALUE self, VALUE ival)
{
	Molecule *mol;
	int idx, n;
	mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	idx = NUM2INT(rb_Integer(ival));
	if (mol == NULL)
		n = gBuiltinParameters->nvdwPars;
	else if (mol->par != NULL)
		n = mol->par->nvdwPars;
	else n = 0;
	if (idx < -n || idx >= n)
		rb_raise(rb_eMolbyError, "Vdw parameter index (%d) out of range", idx);
	if (idx < 0)
		idx += n;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kVdwParType, idx);
}

/*
 *  call-seq:
 *     vdw_pair(idx)          -> ParameterRef
 *  
 *  The index-th vdw pair parameter record is returned.
 */
static VALUE
s_Parameter_VdwPair(VALUE self, VALUE ival)
{
	Molecule *mol;
	int idx, n;
	mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	idx = NUM2INT(rb_Integer(ival));
	if (mol == NULL)
		n = gBuiltinParameters->nvdwpPars;
	else if (mol->par != NULL)
		n = mol->par->nvdwpPars;
	else n = 0;
	if (idx < -n || idx >= n)
		rb_raise(rb_eMolbyError, "Vdw pair parameter index (%d) out of range", idx);
	if (idx < 0)
		idx += n;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kVdwPairParType, idx);
}

/*
 *  call-seq:
 *     vdw_cutoff(idx)          -> ParameterRef
 *  
 *  The index-th vdw cutoff parameter record is returned.
 */
static VALUE
s_Parameter_VdwCutoff(VALUE self, VALUE ival)
{
	Molecule *mol;
	int idx, n;
	mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	idx = NUM2INT(rb_Integer(ival));
	if (mol == NULL)
		n = gBuiltinParameters->nvdwCutoffPars;
	else if (mol->par != NULL)
		n = mol->par->nvdwCutoffPars;
	else n = 0;
	if (idx < -n || idx >= n)
		rb_raise(rb_eMolbyError, "Dihedral parameter index (%d) out of range", idx);
	if (idx < 0)
		idx += n;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kVdwCutoffParType, idx);
}

/*
 *  call-seq:
 *     element(idx)            -> ParameterRef
 *     element(t1)             -> ParameterRef
 *  
 *  In the first form, the index-th element parameter record is returned. In the second
 *  form, the element parameter for t1 is looked up (the last index first). t1
 *  is the element name string (up to 4 characters).
 *  Unlike other Parameter methods, this is used only for the global parameter.
 */
static VALUE
s_Parameter_Element(VALUE self, VALUE ival)
{
	Int idx1;
	if (rb_obj_is_kind_of(ival, rb_cNumeric)) {
		int n = gCountElementParameters;
		idx1 = NUM2INT(rb_Integer(ival));
		if (idx1 < -n || idx1 >= n)
			rb_raise(rb_eMolbyError, "Element parameter index (%d) out of range", idx1);
		if (idx1 < 0)
			idx1 += n;
		return ValueFromMoleculeWithParameterTypeAndIndex(NULL, kElementParType, idx1);
	} else {
		ElementPar *ep;
		char name[6];
		int i;
		strncpy(name, StringValuePtr(ival), 4);
		name[4] = 0;
		for (i = gCountElementParameters - 1, ep = gElementParameters + i; i >= 0; i--, ep--) {
			if (strncmp(ep->name, name, 4) == 0)
				return ValueFromMoleculeWithParameterTypeAndIndex(NULL, kElementParType, i);
		}
		return Qnil;
	}
}

/*
 *  call-seq:
 *     nbonds          -> Integer
 *  
 *  Returns the number of bond parameters.
 */
static VALUE
s_Parameter_Nbonds(VALUE self)
{
	Int n;
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	if (mol == NULL)
		n = gBuiltinParameters->nbondPars;
	else if (mol->par != NULL)
		n = mol->par->nbondPars;
	else n = 0;
	return INT2NUM(n);
}

/*
 *  call-seq:
 *     nangles          -> Integer
 *  
 *  Returns the number of angle parameters.
 */
static VALUE
s_Parameter_Nangles(VALUE self)
{
	Int n;
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	if (mol == NULL)
		n = gBuiltinParameters->nanglePars;
	else if (mol->par != NULL)
		n = mol->par->nanglePars;
	else n = 0;
	return INT2NUM(n);
}

/*
 *  call-seq:
 *     ndihedrals          -> Integer
 *  
 *  Returns the number of dihedral parameters.
 */
static VALUE
s_Parameter_Ndihedrals(VALUE self)
{
	Int n;
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	if (mol == NULL)
		n = gBuiltinParameters->ndihedralPars;
	else if (mol->par != NULL)
		n = mol->par->ndihedralPars;
	else n = 0;
	return INT2NUM(n);
}

/*
 *  call-seq:
 *     nimpropers          -> Integer
 *  
 *  Returns the number of improper parameters.
 */
static VALUE
s_Parameter_Nimpropers(VALUE self)
{
	Int n;
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	if (mol == NULL)
		n = gBuiltinParameters->nimproperPars;
	else if (mol->par != NULL)
		n = mol->par->nimproperPars;
	else n = 0;
	return INT2NUM(n);
}

/*
 *  call-seq:
 *     nvdws          -> Integer
 *  
 *  Returns the number of vdw parameters.
 */
static VALUE
s_Parameter_Nvdws(VALUE self)
{
	Int n;
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	if (mol == NULL)
		n = gBuiltinParameters->nvdwPars;
	else if (mol->par != NULL)
		n = mol->par->nvdwPars;
	else n = 0;
	return INT2NUM(n);
}

/*
 *  call-seq:
 *     nvdw_pairs          -> Integer
 *  
 *  Returns the number of vdw pair parameters.
 */
static VALUE
s_Parameter_NvdwPairs(VALUE self)
{
	Int n;
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	if (mol == NULL)
		n = gBuiltinParameters->nvdwpPars;
	else if (mol->par != NULL)
		n = mol->par->nvdwpPars;
	else n = 0;
	return INT2NUM(n);
}

/*
 *  call-seq:
 *     nvdw_cutoffs          -> Integer
 *  
 *  Returns the number of vdw cutoff parameters.
 */
static VALUE
s_Parameter_NvdwCutoffs(VALUE self)
{
	Int n;
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	if (mol == NULL)
		n = gBuiltinParameters->nvdwCutoffPars;
	else if (mol->par != NULL)
		n = mol->par->nvdwCutoffPars;
	else n = 0;
	return INT2NUM(n);
}

/*
 *  call-seq:
 *     nelements          -> Integer
 *  
 *  Returns the number of element parameters.
 */
static VALUE
s_Parameter_Nelements(VALUE self)
{
	return INT2NUM(gCountElementParameters);
}

/*
 *  call-seq:
 *     bonds          -> ParEnumerable
 *  
 *  Returns a ParEnumerable value that (formally) points to the collection of bond parameters.
 *  Parameter.bonds[x] is equivalent to Parameter.bond(x). ParEnumerable class is
 *  useful when all accessible parameters should be examined by use of 'each' method.
 */
static VALUE
s_Parameter_Bonds(VALUE self)
{
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	return s_NewParEnumerableValueFromMoleculeAndType(mol, kBondParType);
}

/*
 *  call-seq:
 *     angles          -> ParEnumerable
 *  
 *  Returns a ParEnumerable value that (formally) points to the collection of angle parameters.
 *  Parameter.angles[x] is equivalent to Parameter.angle(x). ParEnumerable class is
 *  useful when all accessible parameters should be examined by use of 'each' method.
 */
static VALUE
s_Parameter_Angles(VALUE self)
{
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	return s_NewParEnumerableValueFromMoleculeAndType(mol, kAngleParType);
}

/*
 *  call-seq:
 *     dihedrals          -> ParEnumerable
 *  
 *  Returns a ParEnumerable value that (formally) points to the collection of dihedral parameters.
 *  Parameter.dihedrals[x] is equivalent to Parameter.dihedral(x). ParEnumerable class is
 *  useful when all accessible parameters should be examined by use of 'each' method.
 */
static VALUE
s_Parameter_Dihedrals(VALUE self)
{
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	return s_NewParEnumerableValueFromMoleculeAndType(mol, kDihedralParType);
}

/*
 *  call-seq:
 *     impropers          -> ParEnumerable
 *  
 *  Returns a ParEnumerable value that (formally) points to the collection of improper parameters.
 *  Parameter.impropers[x] is equivalent to Parameter.improper(x). ParEnumerable class is
 *  useful when all accessible parameters should be examined by use of 'each' method.
 */
static VALUE
s_Parameter_Impropers(VALUE self)
{
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	return s_NewParEnumerableValueFromMoleculeAndType(mol, kImproperParType);
}

/*
 *  call-seq:
 *     vdws          -> ParEnumerable
 *  
 *  Returns a ParEnumerable value that (formally) points to the collection of vdw parameters.
 *  Parameter.vdws[x] is equivalent to Parameter.vdw(x). ParEnumerable class is
 *  useful when all accessible parameters should be examined by use of 'each' method.
 */
static VALUE
s_Parameter_Vdws(VALUE self)
{
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	return s_NewParEnumerableValueFromMoleculeAndType(mol, kVdwParType);
}

/*
 *  call-seq:
 *     vdw_pairs          -> ParEnumerable
 *  
 *  Returns a ParEnumerable value that (formally) points to the collection of vdw pair parameters.
 *  Parameter.vdw_pairs[x] is equivalent to Parameter.vdw_pair(x). ParEnumerable class is
 *  useful when all accessible parameters should be examined by use of 'each' method.
 */
static VALUE
s_Parameter_VdwPairs(VALUE self)
{
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	return s_NewParEnumerableValueFromMoleculeAndType(mol, kVdwPairParType);
}

/*
 *  call-seq:
 *     vdw_cutoffs          -> ParEnumerable
 *  
 *  Returns a ParEnumerable value that (formally) points to the collection of vdw cutoff parameters.
 *  Parameter.vdw_cutoffs[x] is equivalent to Parameter.vdw_cutoff(x). ParEnumerable class is
 *  useful when all accessible parameters should be examined by use of 'each' method.
 */
static VALUE
s_Parameter_VdwCutoffs(VALUE self)
{
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	return s_NewParEnumerableValueFromMoleculeAndType(mol, kVdwCutoffParType);
}

/*
 *  call-seq:
 *     elements          -> ParEnumerable
 *  
 *  Returns a ParEnumerable value that (formally) points to the collection of element parameters.
 *  Parameter.elements[x] is equivalent to Parameter.element(x). ParEnumerable class is
 *  useful when all accessible parameters should be examined by use of 'each' method.
 */
static VALUE
s_Parameter_Elements(VALUE self)
{
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	return s_NewParEnumerableValueFromMoleculeAndType(mol, kElementParType);
}

static VALUE
s_Parameter_Lookup_sub(int argc, VALUE *argv, int parType, Molecule *mol)
{
	VALUE atval, optval;
	UInt t[4];
	Int ii[4];
	int i, n, idx, flags, is_global;

	rb_scan_args(argc, argv, "1*", &atval, &optval);
	
	/*  Get the atom types  */
	switch (parType) {
		case kBondParType: n = 2; break;
		case kAngleParType: n = 3; break;
		case kDihedralParType: n = 4; break;
		case kImproperParType: n = 4; break;
		case kVdwParType: n = 1; break;
		case kVdwPairParType: n = 2; break;
		default: return Qnil;
	}
	s_ScanAtomTypes(atval, n, t);
	for (i = 0; i < n; i++) {
		if (t[i] < kAtomTypeMinimum) {
			/*  Explicit atom index  */
			if (mol == NULL)
				rb_raise(rb_eMolbyError, "Explicit atom index (%d) is invalid for global parameters", t[i]);
			if (t[i] >= mol->natoms)
				rb_raise(rb_eMolbyError, "Atom index (%d) out of range", t[i]);
			ii[i] = t[i];
			t[i] = ATOM_AT_INDEX(mol->atoms, t[i])->type;
		} else ii[i] = -1;
	}
	
	/*  Analyze options  */
	flags = 0;
	n = RARRAY_LEN(optval);
	for (i = 0; i < n; i++) {
		VALUE oval = RARRAY_PTR(optval)[i];
		if (oval == ID2SYM(rb_intern("global")))
			flags |= kParameterLookupGlobal;
		else if (oval == ID2SYM(rb_intern("local")))
			flags |= kParameterLookupLocal;
		else if (oval == ID2SYM(rb_intern("missing")))
			flags |= kParameterLookupMissing;
		else if (oval == ID2SYM(rb_intern("nowildcard")))
			flags |= kParameterLookupNoWildcard;
		else if (oval == ID2SYM(rb_intern("nobasetype")))
			flags |= kParameterLookupNoBaseAtomType;
		else if (oval == ID2SYM(rb_intern("create")))
			flags |= 256;
	}
	if (flags == 0)
		flags = kParameterLookupGlobal | kParameterLookupLocal;
	
	idx = -1;
	is_global = 0;
	switch (parType) {
		case kBondParType: {
			BondPar *bp;
			if (mol != NULL) {
				bp = ParameterLookupBondPar(mol->par, t[0], t[1], ii[0], ii[1], flags);
				if (bp != NULL) {
					idx = bp - mol->par->bondPars;
					break;
				}
			}
			bp = ParameterLookupBondPar(gBuiltinParameters, t[0], t[1], -1, -1, flags);
			if (bp != NULL) {
				idx = bp - gBuiltinParameters->bondPars;
				is_global = 1;
			}
			break;
		}
		case kAngleParType: {
			AnglePar *ap;
			if (mol != NULL) {
				ap = ParameterLookupAnglePar(mol->par, t[0], t[1], t[2], ii[0], ii[1], ii[2], flags);
				if (ap != NULL) {
					idx = ap - mol->par->anglePars;
					break;
				}
			}
			ap = ParameterLookupAnglePar(gBuiltinParameters, t[0], t[1], t[2], -1, -1, -1, flags);
			if (ap != NULL) {
				idx = ap - gBuiltinParameters->anglePars;
				is_global = 1;
			}
			break;
		}
		case kDihedralParType: {
			TorsionPar *tp;
			if (mol != NULL) {
				tp = ParameterLookupDihedralPar(mol->par, t[0], t[1], t[2], t[3], ii[0], ii[1], ii[2], ii[3], flags);
				if (tp != NULL) {
					idx = tp - mol->par->dihedralPars;
					break;
				}
			}
			tp = ParameterLookupDihedralPar(gBuiltinParameters, t[0], t[1], t[2], t[3], -1, -1, -1, -1, flags);
			if (tp != NULL) {
				idx = tp - gBuiltinParameters->dihedralPars;
				is_global = 1;
			}
			break;
		}
		case kImproperParType: {
			TorsionPar *tp;
			if (mol != NULL) {
				tp = ParameterLookupImproperPar(mol->par, t[0], t[1], t[2], t[3], ii[0], ii[1], ii[2], ii[3], flags);
				if (tp != NULL) {
					idx = tp - mol->par->improperPars;
					break;
				}
			}
			tp = ParameterLookupImproperPar(gBuiltinParameters, t[0], t[1], t[2], t[3], -1, -1, -1, -1, flags);
			if (tp != NULL) {
				idx = tp - gBuiltinParameters->improperPars;
				is_global = 1;
			}
			break;
		}	
		case kVdwParType: {
			VdwPar *vp;
			if (mol != NULL) {
				vp = ParameterLookupVdwPar(mol->par, t[0], ii[0], flags);
				if (vp != NULL) {
					idx = vp - mol->par->vdwPars;
					break;
				}
			}
			vp = ParameterLookupVdwPar(gBuiltinParameters, t[0], -1, flags);
			if (vp != NULL) {
				idx = vp - gBuiltinParameters->vdwPars;
				is_global = 1;
			}
			break;
		}	
		case kVdwPairParType: {
			VdwPairPar *vp;
			if (mol != NULL) {
				vp = ParameterLookupVdwPairPar(mol->par, t[0], t[1], ii[0], ii[1], flags);
				if (vp != NULL) {
					idx = vp - mol->par->vdwpPars;
					break;
				}
			}
			vp = ParameterLookupVdwPairPar(gBuiltinParameters, t[0], t[1], -1, -1, flags);
			if (vp != NULL) {
				idx = vp - gBuiltinParameters->vdwpPars;
				is_global = 1;
			}
			break;
		}
		default:
			return Qnil;
	}
	if (idx < 0) {
		if ((flags & 256) == 0 || mol == NULL || mol->par == NULL)
			return Qnil;		
		else {
			/*  Insert a new parameter record  */
			UnionPar *up;
			Int count = ParameterGetCountForType(mol->par, parType);
			IntGroup *ig = IntGroupNewWithPoints(count, 1, -1);
			MolActionCreateAndPerform(mol, gMolActionAddParameters, parType, ig, 0, NULL);
			IntGroupRelease(ig);
			is_global = 0;
			idx = count;
			/*  Set atom types  */
			up = ParameterGetUnionParFromTypeAndIndex(mol->par, parType, idx);
			if (up == NULL)
				return Qnil;
			switch (parType) {
				case kBondParType:
					up->bond.type1 = t[0];
					up->bond.type2 = t[1];
					break;
				case kAngleParType:
					up->angle.type1 = t[0];
					up->angle.type2 = t[1];
					up->angle.type3 = t[2];
					break;
				case kDihedralParType:
				case kImproperParType:
					up->torsion.type1 = t[0];
					up->torsion.type2 = t[1];
					up->torsion.type3 = t[2];
					up->torsion.type4 = t[3];
					break;
				case kVdwParType:
					up->vdw.type1 = t[0];
					break;
				case kVdwPairParType:
					up->vdwp.type1 = t[0];
					up->vdwp.type2 = t[1];
					break;
				default:
					return Qnil;
			}
		}
	}
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, parType, idx);
}

/*
 *  call-seq:
 *     lookup(par_type, atom_types, options, ...) -> ParameterRef
 *     lookup(par_type, atom_type_string, options, ...) -> ParameterRef
 *
 *  Find the parameter record that matches the given atom types. The atom types are given
 *  either as an array of string, or a single string delimited by whitespaces or hyphens.
 *  Options are given as symbols. Valid values are :global (look for global parameters), :local
 *  (look for local parameters), :missing (look for missing parameters), :nowildcard (do not 
 *  allow wildcard matching), :nobasetype (the base type does not match for the variant types)
 */
static VALUE
s_Parameter_LookUp(int argc, VALUE *argv, VALUE self)
{
	int parType;
	Molecule *mol = s_MoleculeFromParameterOrParEnumerableValue(self);
	if (argc == 0)
		rb_raise(rb_eMolbyError, "parameter type and atom types must be specified");
	parType = s_ParTypeFromValue(argv[0]);
	return s_Parameter_Lookup_sub(argc - 1, argv + 1, parType, mol);
}

/*
 *  call-seq:
 *     self == parameter -> boolean
 *  
 *  True if the parameters point to the same parameter table.
 */
static VALUE
s_Parameter_Equal(VALUE self, VALUE val)
{
	if (rb_obj_is_kind_of(val, rb_cParameter)) {
		return (s_MoleculeFromParameterOrParEnumerableValue(self) == s_MoleculeFromParameterOrParEnumerableValue(val) ? Qtrue : Qfalse);
	} else return Qfalse;
}

#pragma mark ====== ParEnumerable Class ======

/*  The ParEnumerable class encapsulates the Molecule (not Parameter) pointer
 and the parameter type. If the Molecule is NULL, then it refers to the
 global (built-in) parameters. Note that, even when the Molecule is not NULL,
 the global parameters are always accessible. */

typedef struct ParEnumerable {
	Molecule *mol;
	Int parType;   /*  Same as parType in ParameterRef  */
} ParEnumerable;

static ParEnumerable *
s_ParEnumerableNew(Molecule *mol, Int parType)
{
	ParEnumerable *pen = (ParEnumerable *)calloc(sizeof(ParEnumerable), 1);
	if (pen != NULL) {
		pen->mol = mol;
		if (mol != NULL)
			MoleculeRetain(mol);
		pen->parType = parType;
	}
	return pen;
}

static void
s_ParEnumerableRelease(ParEnumerable *pen)
{
	if (pen != NULL) {
		if (pen->mol != NULL)
			MoleculeRelease(pen->mol);
		free(pen);
	}
}

static Molecule *
s_MoleculeFromParEnumerableValue(VALUE val)
{
	ParEnumerable *pen;
    Data_Get_Struct(val, ParEnumerable, pen);
	return pen->mol;
}

static VALUE
s_NewParEnumerableValueFromMoleculeAndType(Molecule *mol, Int parType)
{
	ParEnumerable *pen = s_ParEnumerableNew(mol, parType);
	if (pen == NULL)
		rb_raise(rb_eMolbyError, "cannot allocate ParEnumerable record");
	return Data_Wrap_Struct(rb_cParEnumerable, 0, (void (*)(void *))s_ParEnumerableRelease, pen);
}

/*
 *  call-seq:
 *     par_type -> String
 *
 *  Get the parameter type, like "bond", "angle", etc.
 */
static VALUE
s_ParEnumerable_ParType(VALUE self) {
	ParEnumerable *pen;
	Int tp;
    Data_Get_Struct(self, ParEnumerable, pen);
	tp = pen->parType;
	if (tp == kElementParType)
		return rb_str_new2("element");
	tp -= kFirstParType;
	if (tp >= 0 && tp < sizeof(s_ParameterTypeNames) / sizeof(s_ParameterTypeNames[0]))
		return rb_str_new2(s_ParameterTypeNames[tp]);
	else rb_raise(rb_eMolbyError, "Internal error: parameter type tag is out of range (%d)", tp);
}

/*
 *  call-seq:
 *     self[idx]          -> ParameterRef
 *  
 *  Call the accessor of the Parameter object from which this ParEnumerable object is derived from.
 *  Thus, if self is "bond" type, self[idx] is equivalent to p.bond(idx), where p is the
 *  parent Parameter object of self.
 *
 *  <b>See Also</b>: Parameter#bond, Parameter#angle, Parameter#dihedral, Parameter#improper, 
 *  Parameter#vdw, Parameter#vdw_pair, Parameter#vdw_cutoff, Parameter#element.
 */
static VALUE
s_ParEnumerable_Aref(VALUE self, VALUE ival)
{
	ParEnumerable *pen;
    Data_Get_Struct(self, ParEnumerable, pen);
	switch (pen->parType) {
			/*  s_Parameter_XXXX() also accepts ParEnumerable argument  */
		case kBondParType:      return s_Parameter_Bond(self, ival);
		case kAngleParType:     return s_Parameter_Angle(self, ival);
		case kDihedralParType:  return s_Parameter_Dihedral(self, ival);
		case kImproperParType:  return s_Parameter_Improper(self, ival);
		case kVdwParType:       return s_Parameter_Vdw(self, ival);
		case kVdwPairParType:   return s_Parameter_VdwPair(self, ival);
		case kVdwCutoffParType: return s_Parameter_VdwCutoff(self, ival);
		case kElementParType:   return s_Parameter_Element(self, ival);
		default:
			rb_raise(rb_eMolbyError, "internal error: unknown parameter type (%d)", pen->parType);
	}
	return Qnil;  /*  Not reached  */
}

/*
 *  call-seq:
 *     length          -> Integer
 *  
 *  Returns the number of parameters included in this enumerable.
 */
static VALUE
s_ParEnumerable_Length(VALUE self)
{
	ParEnumerable *pen;
    Data_Get_Struct(self, ParEnumerable, pen);
	switch (pen->parType) {
			/*  s_Parameter_XXXX() also accepts ParEnumerable argument  */
		case kBondParType:      return s_Parameter_Nbonds(self);
		case kAngleParType:     return s_Parameter_Nangles(self); 
		case kDihedralParType:  return s_Parameter_Ndihedrals(self);
		case kImproperParType:  return s_Parameter_Nimpropers(self);
		case kVdwParType:       return s_Parameter_Nvdws(self);
		case kVdwPairParType:   return s_Parameter_NvdwPairs(self);
		case kVdwCutoffParType: return s_Parameter_NvdwCutoffs(self);
		case kElementParType:   return s_Parameter_Nelements(self);
		default:
			rb_raise(rb_eMolbyError, "internal error: unknown parameter type (%d)", pen->parType);
	}
	return Qnil;  /*  Not reached  */
}

/*
 *  call-seq:
 *     each {|pref| ...}
 *  
 *  Call the block for each parameter, passing a ParameterRef object as a block argument.
 */
VALUE
s_ParEnumerable_Each(VALUE self)
{
	VALUE aval;
	ParEnumerable *pen;
	ParameterRef *pref;
	int i, ofs, n;
    Data_Get_Struct(self, ParEnumerable, pen);
	if (pen->parType == kElementParType)
		n = gCountElementParameters;
	else {
		switch (pen->parType) {
			case kBondParType:      ofs = offsetof(Parameter, nbondPars); break;
			case kAngleParType:     ofs = offsetof(Parameter, nanglePars); break;
			case kDihedralParType:  ofs = offsetof(Parameter, ndihedralPars); break;
			case kImproperParType:  ofs = offsetof(Parameter, nimproperPars); break;
			case kVdwParType:       ofs = offsetof(Parameter, nvdwPars); break;
			case kVdwPairParType:   ofs = offsetof(Parameter, nvdwpPars); break;
			case kVdwCutoffParType: ofs = offsetof(Parameter, nvdwCutoffPars); break;
			default:
				rb_raise(rb_eMolbyError, "internal error: unknown parameter type (%d)", pen->parType);
		}
		if (pen->mol == NULL)
			n = *((Int *)((char *)gBuiltinParameters + ofs));
		else if (pen->mol->par != NULL)
			n = *((Int *)((char *)(pen->mol->par) + ofs));
		else return self;
	}		
	aval = ValueFromMoleculeWithParameterTypeAndIndex(pen->mol, pen->parType, 0);
	Data_Get_Struct(aval, ParameterRef, pref);
	for (i = 0; i < n; i++) {
		pref->idx = i;
		rb_yield(aval);
	}
    return self;
}

/*
 *  call-seq:
 *     reverse_each {|pref| ...}
 *  
 *  Call the block for each parameter in the reverse order, passing a ParameterRef object as a block argument.
 */
VALUE
s_ParEnumerable_ReverseEach(VALUE self)
{
	VALUE aval;
	ParEnumerable *pen;
	ParameterRef *pref;
	int i, ofs, n;
    Data_Get_Struct(self, ParEnumerable, pen);
	if (pen->parType == kElementParType)
		n = gCountElementParameters;
	else {
		switch (pen->parType) {
			case kBondParType:      ofs = offsetof(Parameter, nbondPars); break;
			case kAngleParType:     ofs = offsetof(Parameter, nanglePars); break;
			case kDihedralParType:  ofs = offsetof(Parameter, ndihedralPars); break;
			case kImproperParType:  ofs = offsetof(Parameter, nimproperPars); break;
			case kVdwParType:       ofs = offsetof(Parameter, nvdwPars); break;
			case kVdwPairParType:   ofs = offsetof(Parameter, nvdwpPars); break;
			case kVdwCutoffParType: ofs = offsetof(Parameter, nvdwCutoffPars); break;
			default:
				rb_raise(rb_eMolbyError, "internal error: unknown parameter type (%d)", pen->parType);
		}
		if (pen->mol == NULL)
			n = *((Int *)((char *)gBuiltinParameters + ofs));
		else if (pen->mol->par != NULL)
			n = *((Int *)((char *)(pen->mol->par) + ofs));
		else return self;
	}		
	aval = ValueFromMoleculeWithParameterTypeAndIndex(pen->mol, pen->parType, 0);
	Data_Get_Struct(aval, ParameterRef, pref);
	for (i = n - 1; i >= 0; i--) {
		pref->idx = i;
		rb_yield(aval);
	}
    return self;
}

/*
 *  call-seq:
 *     insert(idx = nil, pref = nil)       -> ParameterRef
 *  
 *  Insert a new parameter at the specified position (if idx is nil, then at the end).
 *  If a ParameterRef is given, then the content of the parameter is copied to the new parameter,
 *  and the parameter is set as molecule-local. Otherwise, the new parameter is blank, and the
 *  parameter is left undefined.
 *  Throws an exception if ParEnumerable points to the global parameter.
 */
static VALUE
s_ParEnumerable_Insert(int argc, VALUE *argv, VALUE self)
{
	VALUE ival, pval;
	ParEnumerable *pen;
	int i, n;
	IntGroup *ig;
	UnionPar u;
	MolAction *act;
    Data_Get_Struct(self, ParEnumerable, pen);
	if (pen->mol == NULL)
		rb_raise(rb_eMolbyError, "the global parameters cannot be modified");
	n = ParameterGetCountForType(pen->mol->par, pen->parType);
	rb_scan_args(argc, argv, "02", &ival, &pval);
	if (ival != Qnil) {
		i = NUM2INT(rb_Integer(ival));
		if (i < 0 || i > n)
			rb_raise(rb_eMolbyError, "the parameter index (%d) out of range (should be 0..%d)", i, n);
		n = i;
	}
	if (pval != Qnil) {
		Int type;
		UnionPar *up = s_UnionParFromValue(pval, &type, 0);
		if (up == NULL || type != pen->parType)
			rb_raise(rb_eMolbyError, "the parameter specification is not correct");
		ParameterCopyOneWithType(&u, up, pen->parType);
		u.bond.src = 0;
	} else {
		memset(&u, 0, sizeof(u));
		u.bond.src = 0;
	}
	ig = IntGroupNewWithPoints(n, 1, -1);
	ParameterInsert(pen->mol->par, pen->parType, &u, ig);

	act = MolActionNew(gMolActionDeleteParameters, pen->parType, ig);
	MolActionCallback_registerUndo(pen->mol, act);
	MolActionRelease(act);
	
	IntGroupRelease(ig);
	return ValueFromMoleculeWithParameterTypeAndIndex(pen->mol, pen->parType, n);
}

/*
 *  call-seq:
 *     delete(Integer)
 *     delete(IntGroup)
 *  
 *  Delete the parameter(s) specified by the argument.
 */
static VALUE
s_ParEnumerable_Delete(VALUE self, VALUE ival)
{
	ParEnumerable *pen;
	int i, n;
	IntGroup *ig;
    Data_Get_Struct(self, ParEnumerable, pen);
	if (pen->mol == NULL)
		rb_raise(rb_eMolbyError, "the global parameters cannot be modified");
	n = ParameterGetCountForType(pen->mol->par, pen->parType);
	if (rb_obj_is_kind_of(ival, rb_cInteger)) {
		ig = IntGroupNewWithPoints(NUM2INT(ival), 1, -1);
		i = 1;
	} else {
		ig = IntGroupFromValue(ival);
		if ((i = IntGroupGetCount(ig)) == 0) {
			IntGroupRelease(ig);
			return Qnil;
		}
	}
	if ((i = IntGroupGetNthPoint(ig, i - 1)) >= n)
		rb_raise(rb_eMolbyError, "the parameter index (%d) out of range (should be 0..%d)", i, n);
	n = i;

	MolActionCreateAndPerform(pen->mol, gMolActionDeleteParameters, pen->parType, ig);
	IntGroupRelease(ig);
	return ival;
}

/*
 *  call-seq:
 *     lookup(atom_types, options, ...) -> ParameterRef
 *     lookup(atom_type_string, options, ...) -> ParameterRef
 *
 *  Find the parameter record that matches the given atom types. The arguments are
 *  the same as Parameter#lookup, except for the parameter type which is implicitly
 *  specified.
 */
static VALUE
s_ParEnumerable_LookUp(int argc, VALUE *argv, VALUE self)
{
	ParEnumerable *pen;
    Data_Get_Struct(self, ParEnumerable, pen);
	return s_Parameter_Lookup_sub(argc, argv, pen->parType, pen->mol);
}

/*
 *  call-seq:
 *     self == parEnumerable -> boolean
 *  
 *  True if the arguments point to the same parameter table and type.
 */
static VALUE
s_ParEnumerable_Equal(VALUE self, VALUE val)
{
	if (rb_obj_is_kind_of(val, rb_cParEnumerable)) {
		ParEnumerable *pen1, *pen2;
		Data_Get_Struct(self, ParEnumerable, pen1);
		Data_Get_Struct(val, ParEnumerable, pen2);
		return (pen1->mol == pen2->mol && pen1->parType == pen2->parType) ? Qtrue : Qfalse;
	} else return Qfalse;
}

#pragma mark ====== AtomRef Class ======

/*  Forward declaration for register undo  */
static VALUE s_Molecule_RegisterUndo(int argc, VALUE *argv, VALUE self);

/*  Ruby string "set_atom_attr"  */
static VALUE s_SetAtomAttrString;

static int
s_AtomIndexFromValue(VALUE self, Atom **app, Molecule **mpp)
{
	AtomRef *aref;
	int idx;
	Data_Get_Struct(self, AtomRef, aref);
	idx = (aref->idx >= 0 ? aref->idx : aref->mol->natoms + aref->idx);
	if (idx < 0 || idx >= aref->mol->natoms)
		rb_raise(rb_eMolbyError, "atom index out of range (%d; should be %d..%d)", aref->idx, -aref->mol->natoms, aref->mol->natoms - 1);
	if (app != NULL)
		*app = aref->mol->atoms + idx;
	if (mpp != NULL)
		*mpp = aref->mol;
	return idx;
}

static Atom *
s_AtomFromValue(VALUE self)
{
	Atom *ap;
	s_AtomIndexFromValue(self, &ap, NULL);
	return ap;
}

static Atom *
s_AtomAndMoleculeFromValue(VALUE self, Molecule **mpp)
{
	Atom *ap;
	s_AtomIndexFromValue(self, &ap, mpp);
	return ap;
}

static void
s_NotifyModificationForAtomRef(VALUE self)
{
	AtomRef *aref;
	Data_Get_Struct(self, AtomRef, aref);
	MoleculeIncrementModifyCount(aref->mol);
}

static void
s_RegisterUndoForAtomAttrChange(VALUE self, VALUE key, VALUE val, VALUE oldval)
{
	AtomRef *aref;
	Data_Get_Struct(self, AtomRef, aref);
	if (MolActionCallback_isUndoRegistrationEnabled(aref->mol) && !RTEST(rb_funcall(val, s_ID_equal, 1, oldval))) {
		/*  Register undo  */
		MolAction *act;
		act = MolActionNew(SCRIPT_ACTION("rrr"), "set_atom_attr", INT2NUM(aref->idx), key, oldval);
		MolActionCallback_registerUndo(aref->mol, act);
		MoleculeCallback_notifyModification(aref->mol, 0);
		/*  Request MD rebuilt if necessary  */
		if (key == s_AtomTypeSym || key == s_ChargeSym || key == s_WeightSym || key == s_FixPosSym || key == s_FixForceSym)
			aref->mol->needsMDRebuild = 1;
	}
}

VALUE
ValueFromMoleculeAndIndex(Molecule *mol, int idx)
{
	AtomRef *aref;
	aref = AtomRefNew(mol, idx);
	return Data_Wrap_Struct(rb_cAtomRef, 0, (void (*)(void *))AtomRefRelease, aref);
}

static VALUE
s_AtomRef_GetMolecule(VALUE self)
{
	Molecule *mpp;
	s_AtomIndexFromValue(self, NULL, &mpp);
	return ValueFromMolecule(mpp);
}

static VALUE s_AtomRef_GetIndex(VALUE self) {
	return INT2NUM(s_AtomIndexFromValue(self, NULL, NULL));
}

static VALUE s_AtomRef_GetSegSeq(VALUE self) {
	return INT2NUM(s_AtomFromValue(self)->segSeq);
}

static VALUE s_AtomRef_GetSegName(VALUE self) {
	char *p = s_AtomFromValue(self)->segName;
	return rb_str_new(p, strlen_limit(p, 4));
}

static VALUE s_AtomRef_GetResSeq(VALUE self) {
	return INT2NUM(s_AtomFromValue(self)->resSeq);
}

static VALUE s_AtomRef_GetResName(VALUE self) {
	char *p = s_AtomFromValue(self)->resName;
	return rb_str_new(p, strlen_limit(p, 4));
}

static VALUE s_AtomRef_GetName(VALUE self) {
	char *p = s_AtomFromValue(self)->aname;
	return rb_str_new(p, strlen_limit(p, 4));
}

static VALUE s_AtomRef_GetAtomType(VALUE self) {
	int type = s_AtomFromValue(self)->type;
	char *p = (type == 0 ? "" : AtomTypeDecodeToString(type, NULL));
	return rb_str_new(p, strlen_limit(p, 6));
}

static VALUE s_AtomRef_GetCharge(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->charge);
}

static VALUE s_AtomRef_GetWeight(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->weight);
}

static VALUE s_AtomRef_GetElement(VALUE self) {
	char *p = s_AtomFromValue(self)->element;
	return rb_str_new(p, strlen_limit(p, 4));
}

static VALUE s_AtomRef_GetAtomicNumber(VALUE self) {
	return INT2NUM(s_AtomFromValue(self)->atomicNumber);
}

static VALUE s_AtomRef_GetConnects(VALUE self) {
	VALUE retval;
	Int i, *cp;
	Atom *ap = s_AtomFromValue(self);
	retval = rb_ary_new();
	cp = AtomConnectData(&ap->connect);
	for (i = 0; i < ap->connect.count; i++)
		rb_ary_push(retval, INT2NUM(cp[i]));
	return retval;
}

static VALUE s_AtomRef_GetR(VALUE self) {
	return ValueFromVector(&(s_AtomFromValue(self)->r));
}

static VALUE s_AtomRef_GetX(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->r.x);
}

static VALUE s_AtomRef_GetY(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->r.y);
}

static VALUE s_AtomRef_GetZ(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->r.z);
}

static Vector s_AtomRef_GetFractionalRAsVector(VALUE self) {
	Atom *ap;
	Molecule *mp;
	Vector r1;
	s_AtomIndexFromValue(self, &ap, &mp);
	r1 = ap->r;
	if (mp->cell != NULL)
		TransformVec(&r1, mp->cell->rtr, &r1);
	return r1;
}

static VALUE s_AtomRef_GetFractionalR(VALUE self) {
	Vector r1 = s_AtomRef_GetFractionalRAsVector(self);
	return ValueFromVector(&r1);
}

static VALUE s_AtomRef_GetFractionalX(VALUE self) {
	return rb_float_new(s_AtomRef_GetFractionalRAsVector(self).x);
}

static VALUE s_AtomRef_GetFractionalY(VALUE self) {
	return rb_float_new(s_AtomRef_GetFractionalRAsVector(self).y);
}

static VALUE s_AtomRef_GetFractionalZ(VALUE self) {
	return rb_float_new(s_AtomRef_GetFractionalRAsVector(self).z);
}

static VALUE s_AtomRef_GetSigma(VALUE self) {
	return ValueFromVector(&(s_AtomFromValue(self)->sigma));
}

static VALUE s_AtomRef_GetSigmaX(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->sigma.x);
}

static VALUE s_AtomRef_GetSigmaY(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->sigma.y);
}

static VALUE s_AtomRef_GetSigmaZ(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->sigma.z);
}

static VALUE s_AtomRef_GetV(VALUE self) {
	return ValueFromVector(&(s_AtomFromValue(self)->v));
}

static VALUE s_AtomRef_GetF(VALUE self) {
	Vector v = s_AtomFromValue(self)->f;
	VecScaleSelf(v, INTERNAL2KCAL);
	return ValueFromVector(&v);
}

static VALUE s_AtomRef_GetOccupancy(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->occupancy);
}

static VALUE s_AtomRef_GetTempFactor(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->tempFactor);
}

static VALUE s_AtomRef_GetAniso(VALUE self) {
	VALUE retval;
	int i;
	Atom *ap = s_AtomFromValue(self);
	if (ap->aniso == NULL)
		return Qnil;
	retval = rb_ary_new();
	for (i = 0; i < 6; i++)
		rb_ary_push(retval, rb_float_new(ap->aniso->bij[i]));
	if (ap->aniso->has_bsig) {
		rb_ary_push(retval, INT2NUM(0));
		for (i = 0; i < 6; i++)
			rb_ary_push(retval, rb_float_new(ap->aniso->bsig[i]));
	}
	return retval;
}

static VALUE s_AtomRef_GetSymop(VALUE self) {
	VALUE retval;
	Atom *ap = s_AtomFromValue(self);
	if (!ap->symop.alive)
		return Qnil;
	retval = rb_ary_new();
	rb_ary_push(retval, INT2NUM(ap->symop.sym));
	rb_ary_push(retval, INT2NUM(ap->symop.dx));
	rb_ary_push(retval, INT2NUM(ap->symop.dy));
	rb_ary_push(retval, INT2NUM(ap->symop.dz));
	rb_ary_push(retval, INT2NUM(ap->symbase));
	return retval;
}

static VALUE s_AtomRef_GetIntCharge(VALUE self) {
	return INT2NUM(s_AtomFromValue(self)->intCharge);
}

static VALUE s_AtomRef_GetFixForce(VALUE self) {
	return rb_float_new(s_AtomFromValue(self)->fix_force * INTERNAL2KCAL);
}

static VALUE s_AtomRef_GetFixPos(VALUE self) {
	return ValueFromVector(&(s_AtomFromValue(self)->fix_pos));
}

static VALUE s_AtomRef_GetExclusion(VALUE self) {
	Molecule *mol;
	Atom *ap;
	int idx, i;
	MDExclusion *exinfo;
	Int *exlist;
	VALUE retval, aval;
	idx = s_AtomIndexFromValue(self, &ap, &mol);
	if (mol->par == NULL || mol->arena == NULL || mol->arena->is_initialized == 0 || mol->needsMDRebuild) {
		VALUE mval = ValueFromMolecule(mol);
		s_RebuildMDParameterIfNecessary(mval, Qnil);
	}
	if (mol->arena->exinfo == NULL)
		return Qnil;
	exinfo = mol->arena->exinfo + idx;
	exlist = mol->arena->exlist;
	retval = rb_ary_new();
	aval = rb_ary_new();
	for (i = exinfo->index1; i < exinfo->index2; i++)  /* 1-2 exclusion  */
		rb_ary_push(aval, INT2NUM(exlist[i]));
	rb_ary_push(retval, aval);
	aval = rb_ary_new();
	for (i = exinfo->index2; i < exinfo->index3; i++)  /* 1-3 exclusion  */
		rb_ary_push(aval, INT2NUM(exlist[i]));
	rb_ary_push(retval, aval);
	aval = rb_ary_new();
	for (i = exinfo->index3; i < (exinfo + 1)->index0; i++)  /* 1-4 exclusion  */
		rb_ary_push(aval, INT2NUM(exlist[i]));
	rb_ary_push(retval, aval);
	return retval;
}

static VALUE s_AtomRef_GetMMExclude(VALUE self) {
	return INT2NUM(s_AtomFromValue(self)->mm_exclude);
}

static VALUE s_AtomRef_GetPeriodicExclude(VALUE self) {
	return INT2NUM(s_AtomFromValue(self)->periodic_exclude);
}

static VALUE s_AtomRef_GetHidden(VALUE self) {
	return ((s_AtomFromValue(self)->exflags & kAtomHiddenFlag) != 0 ? Qtrue : Qfalse);
}

static VALUE s_AtomRef_GetAnchorList(VALUE self) {
	VALUE retval;
	Int i, count, *cp;
	Atom *ap = s_AtomFromValue(self);
	if (ap->anchor == NULL)
		return Qnil;
	count = ap->anchor->connect.count;
	retval = rb_ary_new2(count * 2);
	cp = AtomConnectData(&ap->anchor->connect);
	for (i = 0; i < count; i++) {
		rb_ary_store(retval, i, INT2NUM(cp[i]));
		rb_ary_store(retval, i + count, rb_float_new(ap->anchor->coeffs[i]));
	}
	return retval;
}

static VALUE s_AtomRef_GetUFFType(VALUE self) {
	char *p = s_AtomFromValue(self)->uff_type;
	return rb_str_new(p, strlen_limit(p, 5));
}

static VALUE s_AtomRef_SetIndex(VALUE self, VALUE val) {
	rb_raise(rb_eMolbyError, "index cannot be directly set");
	return Qnil;
}

static VALUE s_AtomRef_SetSegSeq(VALUE self, VALUE val) {
	VALUE oval = s_AtomRef_GetSegSeq(self);
	val = rb_Integer(val);
	s_AtomFromValue(self)->segSeq = NUM2INT(val);
	s_RegisterUndoForAtomAttrChange(self, s_SegSeqSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetSegName(VALUE self, VALUE val) {
	char *p = StringValuePtr(val);
	VALUE oval = s_AtomRef_GetSegName(self);
	strncpy(s_AtomFromValue(self)->segName, p, 4);
	s_RegisterUndoForAtomAttrChange(self, s_SegNameSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetResSeqOrResName(VALUE self, VALUE val) {
	rb_raise(rb_eMolbyError, "residue number/names cannot be directly set. Use Molecule#assign_residue instead.");
	return val; /* Not reached */
}

static VALUE s_AtomRef_SetName(VALUE self, VALUE val) {
	Atom *ap = s_AtomFromValue(self);
	char *p = StringValuePtr(val);
	VALUE oval = s_AtomRef_GetName(self);
	if (ap->anchor != NULL && p[0] == '_')
		rb_raise(rb_eMolbyError, "please avoid a name beginning with '_' for a pi anchor, because it causes some internal confusion.");
	strncpy(ap->aname, p, 4);
	s_RegisterUndoForAtomAttrChange(self, s_NameSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetAtomType(VALUE self, VALUE val) {
	Molecule *mp;
	char *p = StringValuePtr(val);
	VALUE oval = s_AtomRef_GetAtomType(self);
	int type = (p[0] == 0 ? 0 : AtomTypeEncodeToUInt(p));
	if (type != 0 && type < kAtomTypeMinimum)
		rb_raise(rb_eMolbyError, "Invalid atom type '%s'", p);
	s_AtomAndMoleculeFromValue(self, &mp)->type = type;
	mp->needsMDRebuild = 1;
	s_RegisterUndoForAtomAttrChange(self, s_AtomTypeSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetCharge(VALUE self, VALUE val) {
	Molecule *mp;
	VALUE oval = s_AtomRef_GetCharge(self);
	val = rb_Float(val);
	s_AtomAndMoleculeFromValue(self, &mp)->charge = NUM2DBL(val);
	mp->needsMDRebuild = 1;
	s_RegisterUndoForAtomAttrChange(self, s_ChargeSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetWeight(VALUE self, VALUE val) {
	Molecule *mp;
	VALUE oval = s_AtomRef_GetWeight(self);
	val = rb_Float(val);
	s_AtomAndMoleculeFromValue(self, &mp)->weight = NUM2DBL(val);
	mp->needsMDRebuild = 1;
	s_RegisterUndoForAtomAttrChange(self, s_WeightSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetElement(VALUE self, VALUE val) {
	Double w;
	Molecule *mp;
	Atom *ap = s_AtomAndMoleculeFromValue(self, &mp);
	char *p = StringValuePtr(val);
	VALUE oval = s_AtomRef_GetElement(self);
	ap->atomicNumber = ElementToInt(p);
	ElementToString(ap->atomicNumber, ap->element);
	if ((w = WeightForAtomicNumber(ap->atomicNumber)) > 0.0)
		ap->weight = w;
	s_RegisterUndoForAtomAttrChange(self, s_ElementSym, val, oval);
	mp->needsMDRebuild = 1;
	return val;
}

static VALUE s_AtomRef_SetAtomicNumber(VALUE self, VALUE val) {
	Double w;
	Molecule *mp;
	Atom *ap = s_AtomAndMoleculeFromValue(self, &mp);
	VALUE oval = s_AtomRef_GetAtomicNumber(self);
	val = rb_Integer(val);
	ap->atomicNumber = NUM2INT(val);
	ElementToString(ap->atomicNumber, ap->element);
	if ((w = WeightForAtomicNumber(ap->atomicNumber)) > 0.0)
		ap->weight = w;
	s_RegisterUndoForAtomAttrChange(self, s_AtomicNumberSym, val, oval);
	mp->needsMDRebuild = 1;
	return val;
}

static VALUE s_AtomRef_SetConnects(VALUE self, VALUE val) {
	rb_raise(rb_eMolbyError, "connection table cannot be directly set. Use Molecule#create_bond instead.");
	return val; /* Not reached */
}

static VALUE s_AtomRef_SetR(VALUE self, VALUE val) {
	Vector v;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetR(self);
	VectorFromValue(val, &v);
	s_AtomAndMoleculeFromValue(self, &mp)->r = v;
	s_RegisterUndoForAtomAttrChange(self, s_RSym, val, oval);
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetX(VALUE self, VALUE val) {
	Double f;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetX(self);
	val = rb_Float(val);
	f = NUM2DBL(val);
 	s_AtomAndMoleculeFromValue(self, &mp)->r.x = f;
	s_RegisterUndoForAtomAttrChange(self, s_XSym, val, oval);
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetY(VALUE self, VALUE val) {
	Double f;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetY(self);
	val = rb_Float(val);
	f = NUM2DBL(val);
	s_AtomAndMoleculeFromValue(self, &mp)->r.y = f;
	s_RegisterUndoForAtomAttrChange(self, s_YSym, val, oval);
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetZ(VALUE self, VALUE val) {
	Double f;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetZ(self);
	val = rb_Float(val);
	f = NUM2DBL(val);
	s_AtomAndMoleculeFromValue(self, &mp)->r.z = f;
	s_RegisterUndoForAtomAttrChange(self, s_ZSym, val, oval);
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetFractionalR(VALUE self, VALUE val) {
	Vector v, ov;
	Atom *ap;
	Molecule *mp;
	s_AtomIndexFromValue(self, &ap, &mp);
	ov = ap->r;
	VectorFromValue(val, &v);
	if (mp->cell != NULL)
		TransformVec(&v, mp->cell->tr, &v);
	ap->r = v;
	s_RegisterUndoForAtomAttrChange(self, s_RSym, Qnil, ValueFromVector(&ov));
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetFractionalX(VALUE self, VALUE val) {
	double f;
	Vector v, ov;
	Atom *ap;
	Molecule *mp;
	s_AtomIndexFromValue(self, &ap, &mp);
	ov = v = ap->r;
	val = rb_Float(val);
	f = NUM2DBL(val);
	if (mp->cell != NULL) {
		TransformVec(&v, mp->cell->rtr, &v);
		v.x = f;
		TransformVec(&v, mp->cell->tr, &v);
	} else v.x = f;
	ap->r = v;
	s_RegisterUndoForAtomAttrChange(self, s_RSym, Qnil, ValueFromVector(&ov));
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetFractionalY(VALUE self, VALUE val) {
	double f;
	Vector v, ov;
	Atom *ap;
	Molecule *mp;
	s_AtomIndexFromValue(self, &ap, &mp);
	ov = v = ap->r;
	val = rb_Float(val);
	f = NUM2DBL(val);
	if (mp->cell != NULL) {
		TransformVec(&v, mp->cell->rtr, &v);
		v.y = f;
		TransformVec(&v, mp->cell->tr, &v);
	} else v.y = f;
	ap->r = v;
	s_RegisterUndoForAtomAttrChange(self, s_RSym, Qnil, ValueFromVector(&ov));
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetFractionalZ(VALUE self, VALUE val) {
	double f;
	Vector v, ov;
	Atom *ap;
	Molecule *mp;
	s_AtomIndexFromValue(self, &ap, &mp);
	ov = v = ap->r;
	val = rb_Float(val);
	f = NUM2DBL(val);
	if (mp->cell != NULL) {
		TransformVec(&v, mp->cell->rtr, &v);
		v.z = f;
		TransformVec(&v, mp->cell->tr, &v);
	} else v.z = f;
	ap->r = v;
	s_RegisterUndoForAtomAttrChange(self, s_RSym, Qnil, ValueFromVector(&ov));
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetSigma(VALUE self, VALUE val) {
	Vector v;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetSigma(self);
	VectorFromValue(val, &v);
	s_AtomAndMoleculeFromValue(self, &mp)->sigma = v;
	s_RegisterUndoForAtomAttrChange(self, s_SigmaSym, val, oval);
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetSigmaX(VALUE self, VALUE val) {
	Double f;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetSigmaX(self);
	val = rb_Float(val);
	f = NUM2DBL(val);
 	s_AtomAndMoleculeFromValue(self, &mp)->sigma.x = f;
	s_RegisterUndoForAtomAttrChange(self, s_SigmaXSym, val, oval);
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetSigmaY(VALUE self, VALUE val) {
	Double f;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetSigmaY(self);
	val = rb_Float(val);
	f = NUM2DBL(val);
	s_AtomAndMoleculeFromValue(self, &mp)->sigma.y = f;
	s_RegisterUndoForAtomAttrChange(self, s_SigmaYSym, val, oval);
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetSigmaZ(VALUE self, VALUE val) {
	Double f;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetSigmaZ(self);
	val = rb_Float(val);
	f = NUM2DBL(val);
	s_AtomAndMoleculeFromValue(self, &mp)->sigma.z = f;
	s_RegisterUndoForAtomAttrChange(self, s_SigmaZSym, val, oval);
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetV(VALUE self, VALUE val) {
	Vector v;
	Atom *ap;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetV(self);
	VectorFromValue(val, &v);
	s_AtomIndexFromValue(self, &ap, &mp);
	ap->v = v;
	s_RegisterUndoForAtomAttrChange(self, s_VSym, val, oval);
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetF(VALUE self, VALUE val) {
	Vector v;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetF(self);
	VectorFromValue(val, &v);
	VecScaleSelf(v, KCAL2INTERNAL);
	s_AtomAndMoleculeFromValue(self, &mp)->f = v;
	s_RegisterUndoForAtomAttrChange(self, s_FSym, val, oval);
	mp->needsMDCopyCoordinates = 1;
	return val;
}

static VALUE s_AtomRef_SetOccupancy(VALUE self, VALUE val) {
	VALUE oval = s_AtomRef_GetOccupancy(self);
	Molecule *mp;
	val = rb_Float(val);
	s_AtomAndMoleculeFromValue(self, &mp)->occupancy = NUM2DBL(val);
	s_RegisterUndoForAtomAttrChange(self, s_OccupancySym, val, oval);
	mp->needsMDCopyCoordinates = 1;  /*  Occupancy can be used to exclude particular atoms from MM/MD  */
	return val;
}

static VALUE s_AtomRef_SetTempFactor(VALUE self, VALUE val) {
	VALUE oval = s_AtomRef_GetTempFactor(self);
	val = rb_Float(val);
	s_AtomFromValue(self)->tempFactor = NUM2DBL(val);
	s_RegisterUndoForAtomAttrChange(self, s_TempFactorSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetAniso(VALUE self, VALUE val) {
	AtomRef *aref;
	int i, n, type;
	VALUE *valp;
	Double f[12];
	VALUE oval = s_AtomRef_GetAniso(self);
	Data_Get_Struct(self, AtomRef, aref);
	val = rb_funcall(val, rb_intern("to_a"), 0);
	n = RARRAY_LEN(val);
	valp = RARRAY_PTR(val);
	for (i = 0; i < 6; i++) {
		if (i < n)
			f[i] = NUM2DBL(rb_Float(valp[i]));
		else f[i] = 0.0;
	}
	if (n >= 7)
		type = NUM2INT(rb_Integer(valp[6]));
	else type = 0;
	if (n >= 13) {
		for (i = 0; i < 6; i++)
			f[i + 6] = NUM2DBL(rb_Float(valp[i + 7]));
	} else {
		for (i = 0; i < 6; i++)
			f[i + 6] = 0.0;
	}
	i = s_AtomIndexFromValue(self, NULL, NULL);
	MoleculeSetAniso(aref->mol, i, type, f[0], f[1], f[2], f[3], f[4], f[5], (n >= 13 ? f + 6 : NULL));
	s_RegisterUndoForAtomAttrChange(self, s_AnisoSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetSymop(VALUE self, VALUE val) {
	Molecule *mol;
	Atom *ap;
	int i, n;
	VALUE *valp;
	Int ival[5];
	VALUE oval = s_AtomRef_GetSymop(self);
	i = s_AtomIndexFromValue(self, &ap, &mol);
	if (val == Qnil) {
		ival[0] = ival[1] = ival[2] = ival[3] = ival[4] = 0;
	} else {
		val = rb_funcall(val, rb_intern("to_a"), 0);
		n = RARRAY_LEN(val);
		valp = RARRAY_PTR(val);
		for (i = 0; i < 5; i++) {
			if (i < n) {
				if (valp[i] == Qnil)
					ival[i] = -100000;
				else 
					ival[i] = NUM2INT(rb_Integer(valp[i]));
			} else ival[i] = -100000;
		}
	}
	if (ival[0] != -100000 && (ival[0] < 0 || (ival[0] != 0 && ival[0] >= mol->nsyms)))
		rb_raise(rb_eMolbyError, "index of symmetry (%d) is out of range (should be 0..%d)", ival[0], (mol->nsyms == 0 ? 0 : mol->nsyms - 1));
	if (ival[4] != -100000 && (ival[4] < 0 || ival[4] >= mol->natoms))
		rb_raise(rb_eMolbyError, "atom index number (%d) is out of range (should be 0..%d)", ival[4], mol->natoms - 1);
	if (ap->symop.sym == ival[0] && ap->symop.dx == ival[1] && ap->symop.dy == ival[2] && ap->symop.dz == ival[3])
		return val;  /*  No need to change  */
	if (ival[0] != -100000)
		ap->symop.sym = ival[0];
	if (ival[1] != -100000)
		ap->symop.dx = ival[1];
	if (ival[2] != -100000)
		ap->symop.dy = ival[2];
	if (ival[3] != -100000)
		ap->symop.dz = ival[3];
	ap->symop.alive = (SYMOP_ALIVE(ap->symop) != 0);
	if (ival[4] != -100000)
		ap->symbase = ival[4];
	if (ap->symop.alive && (ap->aniso != NULL || (ATOM_AT_INDEX(mol->atoms, ap->symbase))->aniso != NULL)) {
		/*  The anisotropic parameters should be recalculated  */
		VALUE oaval = s_AtomRef_GetAniso(self);
		MoleculeSetAnisoBySymop(mol, i);
		s_RegisterUndoForAtomAttrChange(self, s_AnisoSym, val, oaval);
	}
	s_RegisterUndoForAtomAttrChange(self, s_SymopSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetIntCharge(VALUE self, VALUE val) {
	VALUE oval = s_AtomRef_GetIntCharge(self);
	val = rb_Integer(val);
	s_AtomFromValue(self)->intCharge = NUM2INT(val);
	s_RegisterUndoForAtomAttrChange(self, s_IntChargeSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetFixForce(VALUE self, VALUE val) {
	Molecule *mp;
	VALUE oval = s_AtomRef_GetFixForce(self);
	val = rb_Float(val);
	s_AtomAndMoleculeFromValue(self, &mp)->fix_force = NUM2DBL(val) * KCAL2INTERNAL;
	s_RegisterUndoForAtomAttrChange(self, s_FixForceSym, val, oval);
	mp->needsMDRebuild = 1;
	return val;
}

static VALUE s_AtomRef_SetFixPos(VALUE self, VALUE val) {
	Vector v;
	Molecule *mp;
	VALUE oval = s_AtomRef_GetFixPos(self);
	VectorFromValue(val, &v);
	s_AtomAndMoleculeFromValue(self, &mp)->fix_pos = v;
	s_RegisterUndoForAtomAttrChange(self, s_FixPosSym, val, oval);
	mp->needsMDRebuild = 1;
	return val;
}

static VALUE s_AtomRef_SetExclusion(VALUE self, VALUE val) {
	rb_raise(rb_eMolbyError, "exclusion table is read-only.");
	return val; /* Not reached */
}

static VALUE s_AtomRef_SetMMExclude(VALUE self, VALUE val) {
	VALUE oval = s_AtomRef_GetIntCharge(self);
	val = rb_Integer(val);
	s_AtomFromValue(self)->mm_exclude = NUM2INT(val);
	s_RegisterUndoForAtomAttrChange(self, s_MMExcludeSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetPeriodicExclude(VALUE self, VALUE val) {
	VALUE oval = s_AtomRef_GetIntCharge(self);
	val = rb_Integer(val);
	s_AtomFromValue(self)->periodic_exclude = NUM2INT(val);
	s_RegisterUndoForAtomAttrChange(self, s_PeriodicExcludeSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetHidden(VALUE self, VALUE val) {
	Atom *ap = s_AtomFromValue(self);
	VALUE oval = ((ap->exflags & kAtomHiddenFlag) != 0 ? Qtrue : Qfalse);
	if (RTEST(val)) {
		ap->exflags |= kAtomHiddenFlag;
	} else {
		ap->exflags &= ~kAtomHiddenFlag;
	}
	s_RegisterUndoForAtomAttrChange(self, s_HiddenSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetAnchorList(VALUE self, VALUE val) {
	Int idx, i, j, k, n, *ip;
	Double *dp;
	Atom *ap;
	Molecule *mol;
	VALUE oval, v;
	AtomConnect ac;
	Int nUndoActions;
	MolAction **undoActions;
	memset(&ac, 0, sizeof(ac));
	idx = s_AtomIndexFromValue(self, &ap, &mol);
	oval = s_AtomRef_GetAnchorList(self);
	if (val != Qnil) {
		val = rb_ary_to_ary(val);
		n = RARRAY_LEN(val);
	} else n = 0;
	if (n == 0) {
		if (ap->anchor != NULL) {
			AtomConnectResize(&ap->anchor->connect, 0);
			free(ap->anchor->coeffs);
			free(ap->anchor);
			ap->anchor = NULL;
			s_RegisterUndoForAtomAttrChange(self, s_AnchorListSym, val, oval);
		}
		return val;
	}
	if (n < 2)
		rb_raise(rb_eMolbyError, "set_anchor_list: the argument should have at least two atom indices");
	if (ap->aname[0] == '_')
		rb_raise(rb_eMolbyError, "please avoid a name beginning with '_' for a pi anchor, because it causes some internal confusion.");
	ip = (Int *)malloc(sizeof(Int) * n);
	dp = NULL;
	for (i = 0; i < n; i++) {
		v = RARRAY_PTR(val)[i];
		if (rb_obj_is_kind_of(v, rb_cFloat))
			break;
		j = NUM2INT(rb_Integer(v));
		if (j < 0 || j >= mol->natoms)
			rb_raise(rb_eMolbyError, "set_anchor_list: atom index (%d) out of range", j);
		for (k = 0; k < i; k++) {
			if (ip[k] == j)
				rb_raise(rb_eMolbyError, "set_anchor_list: duplicate atom index (%d)", j);
		}
		ip[i] = j;
	}
	if (i < n) {
		if (i < 2)
			rb_raise(rb_eMolbyError, "set_anchor_list: the argument should have at least two atom indices");
		else if (i * 2 != n)
			rb_raise(rb_eMolbyError, "set_anchor_list: the weights should be given in the same number as the atom indices");
		dp = (Double *)malloc(sizeof(Double) * n / 2);
		for (i = 0; i < n / 2; i++) {
			dp[i] = NUM2DBL(rb_Float(RARRAY_PTR(val)[i + n / 2]));
			if (dp[i] <= 0.0)
				rb_raise(rb_eMolbyError, "set_anchor_list: the weights should be positive Floats");
		}
		n /= 2;
	}
	nUndoActions = 0;
	undoActions = NULL;
	i = MoleculeSetPiAnchorList(mol, idx, n, ip, dp, &nUndoActions, &undoActions);
	free(dp);
	free(ip);
	if (i != 0)
		rb_raise(rb_eMolbyError, "invalid argument");
	if (nUndoActions > 0) {
		for (i = 0; i < nUndoActions; i++) {
			MolActionCallback_registerUndo(mol, undoActions[i]);
			MolActionRelease(undoActions[i]);
		}
		free(undoActions);
	}
	s_RegisterUndoForAtomAttrChange(self, s_AnchorListSym, val, oval);
	return val;
}

static VALUE s_AtomRef_SetUFFType(VALUE self, VALUE val) {
	Atom *ap = s_AtomFromValue(self);
	char *p = StringValuePtr(val);
	VALUE oval = s_AtomRef_GetUFFType(self);
	strncpy(ap->uff_type, p, 5);
	ap->uff_type[5] = 0;
	s_RegisterUndoForAtomAttrChange(self, s_UFFTypeSym, val, oval);
	return val;
}

static struct s_AtomAttrDef {
	char *name;
	VALUE *symref;  /*  Address of s_IndexSymbol etc. */
	ID id;			/*  Will be set within InitMolby()  */
	VALUE (*getter)(VALUE);
	VALUE (*setter)(VALUE, VALUE);
} s_AtomAttrDefTable[] = {
	{"index",        &s_IndexSym,        0, s_AtomRef_GetIndex,        s_AtomRef_SetIndex},
	{"seg_seq",       &s_SegSeqSym,      0, s_AtomRef_GetSegSeq,       s_AtomRef_SetSegSeq},
	{"seg_name",      &s_SegNameSym,     0, s_AtomRef_GetSegName,      s_AtomRef_SetSegName},
	{"res_seq",       &s_ResSeqSym,      0, s_AtomRef_GetResSeq,       s_AtomRef_SetResSeqOrResName},
	{"res_name",      &s_ResNameSym,     0, s_AtomRef_GetResName,      s_AtomRef_SetResSeqOrResName},
	{"name",         &s_NameSym,         0, s_AtomRef_GetName,         s_AtomRef_SetName},
	{"atom_type",     &s_AtomTypeSym,    0, s_AtomRef_GetAtomType,     s_AtomRef_SetAtomType},
	{"charge",       &s_ChargeSym,       0, s_AtomRef_GetCharge,       s_AtomRef_SetCharge},
	{"weight",       &s_WeightSym,       0, s_AtomRef_GetWeight,       s_AtomRef_SetWeight},
	{"element",      &s_ElementSym,      0, s_AtomRef_GetElement,      s_AtomRef_SetElement},
	{"atomic_number", &s_AtomicNumberSym,0, s_AtomRef_GetAtomicNumber, s_AtomRef_SetAtomicNumber},
	{"connects",     &s_ConnectsSym,     0, s_AtomRef_GetConnects,     s_AtomRef_SetConnects},
	{"r",            &s_RSym,            0, s_AtomRef_GetR,            s_AtomRef_SetR},
	{"x",            &s_XSym,            0, s_AtomRef_GetX,            s_AtomRef_SetX},
	{"y",            &s_YSym,            0, s_AtomRef_GetY,            s_AtomRef_SetY},
    {"z",            &s_ZSym,            0, s_AtomRef_GetZ,            s_AtomRef_SetZ},
	{"fract_r",      &s_FractRSym,       0, s_AtomRef_GetFractionalR,  s_AtomRef_SetFractionalR},
	{"fract_x",      &s_FractXSym,       0, s_AtomRef_GetFractionalX,  s_AtomRef_SetFractionalX},
	{"fract_y",      &s_FractYSym,       0, s_AtomRef_GetFractionalY,  s_AtomRef_SetFractionalY},
	{"fract_z",      &s_FractZSym,       0, s_AtomRef_GetFractionalZ,  s_AtomRef_SetFractionalZ},
	{"sigma",        &s_SigmaSym,        0, s_AtomRef_GetSigma,        s_AtomRef_SetSigma},
	{"sigma_x",      &s_SigmaXSym,       0, s_AtomRef_GetSigmaX,       s_AtomRef_SetSigmaX},
	{"sigma_y",      &s_SigmaYSym,       0, s_AtomRef_GetSigmaY,       s_AtomRef_SetSigmaY},
	{"sigma_z",      &s_SigmaZSym,       0, s_AtomRef_GetSigmaZ,       s_AtomRef_SetSigmaZ},
	{"v",            &s_VSym,            0, s_AtomRef_GetV,            s_AtomRef_SetV},
	{"f",            &s_FSym,            0, s_AtomRef_GetF,            s_AtomRef_SetF},
	{"occupancy",    &s_OccupancySym,    0, s_AtomRef_GetOccupancy,    s_AtomRef_SetOccupancy},
	{"temp_factor",  &s_TempFactorSym,   0, s_AtomRef_GetTempFactor,   s_AtomRef_SetTempFactor},
	{"aniso",        &s_AnisoSym,        0, s_AtomRef_GetAniso,        s_AtomRef_SetAniso},
	{"symop",        &s_SymopSym,        0, s_AtomRef_GetSymop,        s_AtomRef_SetSymop},
	{"int_charge",   &s_IntChargeSym,    0, s_AtomRef_GetIntCharge,    s_AtomRef_SetIntCharge},
	{"fix_force",    &s_FixForceSym,     0, s_AtomRef_GetFixForce,     s_AtomRef_SetFixForce},
	{"fix_pos",      &s_FixPosSym,       0, s_AtomRef_GetFixPos,       s_AtomRef_SetFixPos},
	{"exclusion",    &s_ExclusionSym,    0, s_AtomRef_GetExclusion,    s_AtomRef_SetExclusion},
	{"mm_exclude",   &s_MMExcludeSym,    0, s_AtomRef_GetMMExclude,    s_AtomRef_SetMMExclude},
	{"periodic_exclude", &s_PeriodicExcludeSym, 0, s_AtomRef_GetPeriodicExclude, s_AtomRef_SetPeriodicExclude},
	{"hidden",       &s_HiddenSym,       0, s_AtomRef_GetHidden,       s_AtomRef_SetHidden},
	{"anchor_list",  &s_AnchorListSym,   0, s_AtomRef_GetAnchorList,   s_AtomRef_SetAnchorList},
	{"uff_type",     &s_UFFTypeSym,      0, s_AtomRef_GetUFFType,      s_AtomRef_SetUFFType},
	{NULL} /* Sentinel */
};

static VALUE
s_AtomRef_SetAttr(VALUE self, VALUE key, VALUE value)
{
	int i;
	ID kid;
	if (TYPE(key) != T_SYMBOL) {
		kid = rb_intern(StringValuePtr(key));
		key = ID2SYM(kid);
	} else kid = SYM2ID(key);
	for (i = 0; s_AtomAttrDefTable[i].name != NULL; i++) {
		if (s_AtomAttrDefTable[i].id == kid) {
			if (value == Qundef)
				return (*(s_AtomAttrDefTable[i].getter))(self);
			else
				return (*(s_AtomAttrDefTable[i].setter))(self, value);
		}
	}
	rb_raise(rb_eMolbyError, "unknown atom attribute \"%s\"", rb_id2name(kid));
	return Qnil; /* not reached */
}

static VALUE
s_AtomRef_GetAttr(VALUE self, VALUE key)
{
	return s_AtomRef_SetAttr(self, key, Qundef);
}

/*
 *  call-seq:
 *     self == atomRef -> boolean
 *
 *  True if the two references point to the same atom.
 */
static VALUE
s_AtomRef_Equal(VALUE self, VALUE val)
{
	if (rb_obj_is_kind_of(val, rb_cAtomRef)) {
		return (s_AtomFromValue(self) == s_AtomFromValue(val) ? Qtrue : Qfalse);
	} else return Qfalse;
}

#pragma mark ====== MolEnumerable Class ======

static int s_Molecule_AtomIndexFromValue(Molecule *, VALUE);

/*
 *  call-seq:
 *     self[idx] -> AtomRef or Array of Integers
 *  
 *  Get the idx-th atom, bond, etc. for the Molecule from which this MolEnuerable object is
 *  derived from. For the atom, the return value is AtomRef. For the residue, the return
 *  value is a String. Otherwise, the return value is an Array of Integers.
 */
static VALUE
s_MolEnumerable_Aref(VALUE self, VALUE arg1)
{
	MolEnumerable *mseq;
	Molecule *mol;
	int idx1, idx2;
    Data_Get_Struct(self, MolEnumerable, mseq);
	mol = mseq->mol;
	if (mseq->kind == kAtomKind) {
		return ValueFromMoleculeAndIndex(mol, s_Molecule_AtomIndexFromValue(mol, arg1));
	}
	idx1 = NUM2INT(arg1);
	switch (mseq->kind) {
		case kBondKind: {
			idx2 = (idx1 >= 0 ? idx1 : mol->nbonds + idx1);
			if (idx2 < 0 || idx2 >= mol->nbonds)
				rb_raise(rb_eIndexError, "bond index out of range (%d; should be %d..%d)", idx1, -mol->nbonds, mol->nbonds - 1);
			return rb_ary_new3(2, INT2NUM(mol->bonds[idx2 * 2]), INT2NUM(mol->bonds[idx2 * 2 + 1]));
		}
		case kAngleKind: {
			idx2 = (idx1 >= 0 ? idx1 : mol->nangles + idx1);
			if (idx2 < 0 || idx2 >= mol->nangles)
				rb_raise(rb_eIndexError, "angle index out of range (%d; should be %d..%d)", idx1, -mol->nangles, mol->nangles - 1);
			return rb_ary_new3(3, INT2NUM(mol->angles[idx2 * 3]), INT2NUM(mol->angles[idx2 * 3 + 1]), INT2NUM(mol->angles[idx2 * 3 + 2]));
		}
		case kDihedralKind: {
			idx2 = (idx1 >= 0 ? idx1 : mol->ndihedrals + idx1);
			if (idx2 < 0 || idx2 >= mol->ndihedrals)
				rb_raise(rb_eIndexError, "dihedral index out of range (%d; should be %d..%d)", idx1, -mol->ndihedrals, mol->ndihedrals - 1);
			return rb_ary_new3(4, INT2NUM(mol->dihedrals[idx2 * 4]), INT2NUM(mol->dihedrals[idx2 * 4 + 1]), INT2NUM(mol->dihedrals[idx2 * 4 + 2]), INT2NUM(mol->dihedrals[idx2 * 4 + 3]));
		}
		case kImproperKind: {
			idx2 = (idx1 >= 0 ? idx1 : mol->nimpropers + idx1);
			if (idx2 < 0 || idx2 >= mol->nimpropers)
				rb_raise(rb_eIndexError, "improper index out of range (%d; should be %d..%d)", idx1, -mol->nimpropers, mol->nimpropers - 1);
			return rb_ary_new3(4, INT2NUM(mol->impropers[idx2 * 4]), INT2NUM(mol->impropers[idx2 * 4 + 1]), INT2NUM(mol->impropers[idx2 * 4 + 2]), INT2NUM(mol->impropers[idx2 * 4 + 3]));
		}
		case kResidueKind: {
			char *p;
			idx2 = (idx1 >= 0 ? idx1 : mol->nresidues + idx1);
			if (idx2 < 0 || idx2 >= mol->nresidues)
				rb_raise(rb_eIndexError, "residue index out of range (%d; should be %d..%d)", idx1, -mol->nresidues, mol->nresidues - 1);
			p = mol->residues[idx2];
			return rb_str_new(p, strlen_limit(p, 4));
		}
	}
	return Qnil;
}

/*
 *  call-seq:
 *     length          -> Integer
 *  
 *  Returns the number of objects included in this enumerable.
 */
static VALUE
s_MolEnumerable_Length(VALUE self)
{
	MolEnumerable *mseq;
    Data_Get_Struct(self, MolEnumerable, mseq);
	switch (mseq->kind) {
		case kAtomKind:
			return INT2NUM(mseq->mol->natoms);
		case kBondKind:
			return INT2NUM(mseq->mol->nbonds);
		case kAngleKind:
			return INT2NUM(mseq->mol->nangles);
		case kDihedralKind:
			return INT2NUM(mseq->mol->ndihedrals);
		case kImproperKind:
			return INT2NUM(mseq->mol->nimpropers);
		case kResidueKind:
			return INT2NUM(mseq->mol->nresidues);
	}
	return INT2NUM(-1);
}

/*
 *  call-seq:
 *     each {|obj| ...}
 *  
 *  Call the block for each atom/bond/angle/dihedral/improper/residue. The block argument is
 *  an AtomRef for atoms, a String for residues, and an Array of Integers for others.
 *  For the atoms, a same AtomRef object is passed (with different internal information)
 *  for each invocation of block. Otherwise, a new Ruby object will be created and passed
 *  for each iteration.
 */
VALUE
s_MolEnumerable_Each(VALUE self)
{
	MolEnumerable *mseq;
	int i;
	int len = NUM2INT(s_MolEnumerable_Length(self));
    Data_Get_Struct(self, MolEnumerable, mseq);
	if (mseq->kind == kAtomKind) {
		/*  The same AtomRef object will be used during the loop  */
		AtomRef *aref = AtomRefNew(mseq->mol, 0);
		VALUE arval = Data_Wrap_Struct(rb_cAtomRef, 0, (void (*)(void *))AtomRefRelease, aref);
		for (i = 0; i < len; i++) {
			aref->idx = i;
			rb_yield(arval);
		}
    } else {
		/*  A new ruby object will be created at each iteration (not very efficient)  */
		for (i = 0; i < len; i++) {
			rb_yield(s_MolEnumerable_Aref(self, INT2NUM(i)));
		}
	}
    return self;
}

/*
 *  call-seq:
 *     self == molEnumerable -> boolean
 *
 *  True if the two arguments point to the same molecule and enumerable type.
 */
static VALUE
s_MolEnumerable_Equal(VALUE self, VALUE val)
{
	if (rb_obj_is_kind_of(val, rb_cMolEnumerable)) {
		MolEnumerable *mseq1, *mseq2;
		Data_Get_Struct(self, MolEnumerable, mseq1);
		Data_Get_Struct(val, MolEnumerable, mseq2);
		return ((mseq1->mol == mseq2->mol && mseq1->kind == mseq2->kind) ? Qtrue : Qfalse);
	} else return Qfalse;
}


#pragma mark ====== Molecule Class ======

/*  An malloc'ed string buffer. Retains the error/warning message from the last ***load/***save method.  */
/*  Accessible from Ruby as Molecule#error_message and Molecule#error_message=.  */
char *gLoadSaveErrorMessage = NULL;

#define MoleculeClearLoadSaveErrorMessage() (gLoadSaveErrorMessage != NULL ? (free(gLoadSaveErrorMessage), gLoadSaveErrorMessage = NULL) : NULL)

Molecule *
MoleculeFromValue(VALUE val)
{
	Molecule *mol;
	if (!rb_obj_is_kind_of(val, rb_cMolecule))
		rb_raise(rb_eMolbyError, "Molecule instance is expected");
    Data_Get_Struct(val, Molecule, mol);
	return mol;
}

static VALUE sMoleculeRetainArray = Qnil;

/*  The function is called from MoleculeRelease()  */
/*  The Ruby Molecule object is held as mol->exmolobj, and is also protected from */
/*  GC by registering in the sMoleculeRetainArray. This causes the same Ruby */
/*  object is always returned for the same Molecule.  */
/*  When the reference count of the Molecule becomes 1, then the Ruby object is */
/*  removed from sMoleculeRetainArray. In this situation, the Molecule is retained  */
/*  only by the currently alive Ruby containers.  When the Ruby Molecule object is */
/*  removed from all alive Ruby containers, the Ruby object will be collected by */
/*  the next GC invocation, and at that time the Molecule structure is properly released. */

/*  Register/unregister the exmolobj Ruby object  */
void
MoleculeReleaseExternalObj(Molecule *mol)
{
	if (mol != NULL && mol->base.refCount <= 2 && mol->exmolobjProtected == 1) {
		rb_ary_delete(sMoleculeRetainArray, (VALUE)mol->exmolobj);
		mol->exmolobjProtected = 0;
	}
}

void
MoleculeRetainExternalObj(Molecule *mol)
{
	if (mol != NULL && mol->base.refCount >= 2 && mol->exmolobj != NULL && mol->exmolobjProtected == 0) {
		if (sMoleculeRetainArray == Qnil) {
			rb_define_readonly_variable("molecules", &sMoleculeRetainArray);
			sMoleculeRetainArray = rb_ary_new();
		}
		
		rb_ary_push(sMoleculeRetainArray, (VALUE)mol->exmolobj);
		mol->exmolobjProtected = 1;
	}
}

/*  Release hook function for Ruby  */
void
MoleculeReleaseHook(Molecule *mol)
{
	if (mol->exmolobj != NULL) {
		/*  No need to remove from sMoleculeRetainArray  */
		mol->exmolobj = NULL;
		mol->exmolobjProtected = 0;
	}
	MoleculeRelease(mol);
}

VALUE
ValueFromMolecule(Molecule *mol)
{
	if (mol == NULL)
		return Qnil;
	if (mol->exmolobj != NULL)
		return (VALUE)mol->exmolobj;
	mol->exmolobj = (void *)Data_Wrap_Struct(rb_cMolecule, 0, (void (*)(void *))MoleculeReleaseHook, mol);
	MoleculeRetain(mol);  /*  MoleculeRetainExternalObj() is automatically called  */
	return (VALUE)mol->exmolobj;
}

/*  Allocator  */
static VALUE
s_Molecule_Alloc(VALUE klass)
{
	VALUE val;
	Molecule *mol = MoleculeNew();
	val = ValueFromMolecule(mol);
	MoleculeRelease(mol); /*  Only the returned Ruby object retains this molecule  */
	return val;
}

static int
s_Molecule_AtomIndexFromValue(Molecule *mol, VALUE val)
{
	int n;
	char *p = "";
	if (FIXNUM_P(val)) {
		n = FIX2INT(val);
		if (n >= 0 && n < mol->natoms)
			return n;
		n = -1; /*  No such atom  */
		val = rb_inspect(val);
	} else {
		n = MoleculeAtomIndexFromString(mol, StringValuePtr(val));
	}
	if (n >= 0 && n < mol->natoms)
		return n;
	p = StringValuePtr(val);
	if (n == -1)
		rb_raise(rb_eMolbyError, "no such atom: %s", p);
	else if (n == -2)
		rb_raise(rb_eMolbyError, "bad format of atom specification: %s", p);
	else
		rb_raise(rb_eMolbyError, "error in converting value to atom index: %s", p);
	return 0; /* Not reached */
}

static IntGroup *
s_Molecule_AtomGroupFromValue(VALUE self, VALUE val)
{
	IntGroup *ig;
	val = rb_funcall(self, rb_intern("atom_group"), 1, val);
	if (!rb_obj_is_kind_of(val, rb_cIntGroup))
		rb_raise(rb_eMolbyError, "IntGroup instance is expected");
	Data_Get_Struct(val, IntGroup, ig);
	IntGroupRetain(ig);
	return ig;
}

static void
s_Molecule_RaiseOnLoadSave(int status, const char *msg, const char *fname)
{
	if (gLoadSaveErrorMessage != NULL) {
		MyAppCallback_setConsoleColor(1);
		MyAppCallback_showScriptMessage("On loading %s:\n%s\n", fname, gLoadSaveErrorMessage);
		MyAppCallback_setConsoleColor(0);
	}
	if (status != 0)
		rb_raise(rb_eMolbyError, "%s %s", msg, fname);
}

/*
 *  call-seq:
 *     dup          -> Molecule
 *
 *  Duplicate a molecule. All entries are deep copied, so modifying the newly
 *  created object does not affect the old object in any sense.
 */
static VALUE
s_Molecule_InitCopy(VALUE self, VALUE arg)
{
	Molecule *mp1, *mp2;
	Data_Get_Struct(self, Molecule, mp1);
	mp2 = MoleculeFromValue(arg);
	if (MoleculeInitWithMolecule(mp1, mp2) == NULL)
		rb_raise(rb_eMolbyError, "Cannot duplicate molecule");
	return self;
}

/*
 *  call-seq:
 *     loadmbsf(file)       -> bool
 *
 *  Read a structure from a mbsf file.
 *  Return true if successful.
 */
static VALUE
s_Molecule_Loadmbsf(int argc, VALUE *argv, VALUE self)
{
	VALUE fname;
	char *fstr;
	Molecule *mol;
	int retval;
	MoleculeClearLoadSaveErrorMessage();
	Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "1", &fname);
	fstr = FileStringValuePtr(fname);
	retval = MoleculeLoadMbsfFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to load mbsf", fstr);
	return Qtrue;	
}

/*
 *  call-seq:
 *     loadpsf(file, pdbfile = nil)       -> bool
 *
 *  Read a structure from a psf file. molecule must be empty. The psf may be
 *  an "extended" version, which also contains coordinates. If pdbfile 
 *  is given, then atomic coordinates are read from that file.
 *  Return true if successful.
 */
static VALUE
s_Molecule_Loadpsf(int argc, VALUE *argv, VALUE self)
{
	VALUE fname, pdbname;
	char *fstr, *pdbstr;
	Molecule *mol;
	int retval;
	Data_Get_Struct(self, Molecule, mol);
	if (mol->natoms > 0)
		return Qnil;  /*  Must be a new molecule  */
	MoleculeClearLoadSaveErrorMessage();
	rb_scan_args(argc, argv, "11", &fname, &pdbname);
	fstr = FileStringValuePtr(fname);
	retval = MoleculeLoadPsfFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to load psf", fstr);
	pdbstr = NULL;
	if (!NIL_P(pdbname)) {
		pdbstr = strdup(FileStringValuePtr(pdbname));
		retval = MoleculeReadCoordinatesFromPdbFile(mol, pdbstr, &gLoadSaveErrorMessage);
		free(pdbstr);
		s_Molecule_RaiseOnLoadSave(retval, "Failed to load coordinates from pdb", pdbstr);
	}
	return Qtrue;
}

/*
 *  call-seq:
 *     loadpdb(file)       -> bool
 *
 *  Read coordinates from a pdb file. If molecule is empty, then structure is build
 *  by use of CONECT instructions. Otherwise, only the coordinates are read in.
 *  Return true if successful.
 */
static VALUE
s_Molecule_Loadpdb(int argc, VALUE *argv, VALUE self)
{
	VALUE fname;
	char *fstr;
	Molecule *mol;
	int retval;
	Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "1", &fname);
	MoleculeClearLoadSaveErrorMessage();
	fstr = FileStringValuePtr(fname);
	retval = MoleculeReadCoordinatesFromPdbFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to load pdb", fstr);
	return Qtrue;	
}

/*
 *  call-seq:
 *     loaddcd(file)       -> bool
 *
 *  Read coordinates from a dcd file. The molecule should not empty.
 *  Return true if successful.
 */
static VALUE
s_Molecule_Loaddcd(int argc, VALUE *argv, VALUE self)
{
	VALUE fname;
	char *fstr;
	Molecule *mol;
	int retval;
	Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "1", &fname);
	MoleculeClearLoadSaveErrorMessage();
	fstr = FileStringValuePtr(fname);
	retval = MoleculeReadCoordinatesFromDcdFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to load dcd", fstr);
	return Qtrue;	
}

/*
 *  call-seq:
 *     loadtep(file)       -> bool
 *
 *  Read coordinates from an ortep .tep file.
 *  Return true if successful.
 */
static VALUE
s_Molecule_Loadtep(int argc, VALUE *argv, VALUE self)
{
	VALUE fname;
	char *fstr;
	Molecule *mol;
	int retval;
	Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "1", &fname);
	MoleculeClearLoadSaveErrorMessage();
	fstr = FileStringValuePtr(fname);
	retval = MoleculeLoadTepFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to load ORTEP file", fstr);
	return Qtrue;	
}

/*
 *  call-seq:
 *     loadres(file)       -> bool
 *
 *  Read coordinates from a shelx .res file.
 *  Return true if successful.
 */
static VALUE
s_Molecule_Loadres(int argc, VALUE *argv, VALUE self)
{
	VALUE fname;
	char *fstr;
	Molecule *mol;
	int retval;
	Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "1", &fname);
	MoleculeClearLoadSaveErrorMessage();
	fstr = FileStringValuePtr(fname);
	retval = MoleculeLoadShelxFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to load SHELX res file", fstr);
	return Qtrue;	
}

/*
 *  call-seq:
 *     loadfchk(file)       -> bool
 *
 *  Read coordinates and MO information from a Gaussian fchk file. (TODO: implement this) 
 *  Return true if successful.
 */
static VALUE
s_Molecule_Loadfchk(int argc, VALUE *argv, VALUE self)
{
	VALUE fname;
	char *fstr;
	Molecule *mol;
	int retval;
	Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "1", &fname);
	MoleculeClearLoadSaveErrorMessage();
	fstr = FileStringValuePtr(fname);
	retval = MoleculeLoadGaussianFchkFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to load Gaussian fchk", fstr);
	return Qtrue;	
}

/*
 *  call-seq:
 *     loaddat(file)       -> bool
 *
 *  Read coordinates and ESP information from a GAMESS dat file. (TODO: read MO info as well) 
 *  Return true if successful.
 */
static VALUE
s_Molecule_Loaddat(int argc, VALUE *argv, VALUE self)
{
	VALUE fname;
	char *fstr;
	Molecule *mol;
	int retval;
	Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "1", &fname);
	MoleculeClearLoadSaveErrorMessage();
	fstr = FileStringValuePtr(fname);
	MyAppCallback_showProgressPanel("Loading GAMESS dat file...");
	retval = MoleculeLoadGamessDatFile(mol, fstr, &gLoadSaveErrorMessage);
	MyAppCallback_hideProgressPanel();
	s_Molecule_RaiseOnLoadSave(retval, "Failed to load GAMESS dat", fstr);
	return Qtrue;	
}

/*
 *  call-seq:
 *     savembsf(file)       -> bool
 *
 *  Write structure as a mbsf file. Returns true if successful.
 */
static VALUE
s_Molecule_Savembsf(VALUE self, VALUE fname)
{
	char *fstr;
    Molecule *mol;
	int retval;
    Data_Get_Struct(self, Molecule, mol);
	MoleculeClearLoadSaveErrorMessage();
	fstr = FileStringValuePtr(fname);
	retval = MoleculeWriteToMbsfFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to save mbsf", fstr);
	return Qtrue;
}

/*
 *  call-seq:
 *     savepsf(file)       -> bool
 *
 *  Write structure as a psf file. Returns true if successful.
 */
static VALUE
s_Molecule_Savepsf(VALUE self, VALUE fname)
{
	char *fstr;
    Molecule *mol;
	int retval;
    Data_Get_Struct(self, Molecule, mol);
	MoleculeClearLoadSaveErrorMessage();
	fstr = FileStringValuePtr(fname);
	retval = MoleculeWriteToPsfFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to save psf", fstr);
	return Qtrue;
}

/*
 *  call-seq:
 *     savepdb(file)       -> bool
 *
 *  Write coordinates as a pdb file. Returns true if successful.
 */
static VALUE
s_Molecule_Savepdb(VALUE self, VALUE fname)
{
	char *fstr;
    Molecule *mol;
	int retval;
    Data_Get_Struct(self, Molecule, mol);
	MoleculeClearLoadSaveErrorMessage();
	fstr = FileStringValuePtr(fname);
	retval = MoleculeWriteToPdbFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to save pdb", fstr);
	return Qtrue;
}

/*
 *  call-seq:
 *     savedcd(file)       -> bool
 *
 *  Write coordinates as a dcd file. Returns true if successful.
 */
static VALUE
s_Molecule_Savedcd(VALUE self, VALUE fname)
{
	char *fstr;
    Molecule *mol;
	int retval;
    Data_Get_Struct(self, Molecule, mol);
	MoleculeClearLoadSaveErrorMessage();
	fstr = FileStringValuePtr(fname);
	retval = MoleculeWriteToDcdFile(mol, fstr, &gLoadSaveErrorMessage);
	s_Molecule_RaiseOnLoadSave(retval, "Failed to save dcd", fstr);
	return Qtrue;
}

/*
 *  call-seq:
 *     savetep(file)       -> bool
 *
 *  Write coordinates as an ORTEP file. Returns true if successful.
 */
/*
static VALUE
s_Molecule_Savetep(VALUE self, VALUE fname)
{
	char *fstr;
    Molecule *mol;
	char errbuf[128];
    Data_Get_Struct(self, Molecule, mol);
	fstr = FileStringValuePtr(fname);
	if (MoleculeWriteToTepFile(mol, fstr, errbuf, sizeof errbuf) != 0)
		rb_raise(rb_eMolbyError, errbuf);
	return Qtrue;
}
*/

/*  load([ftype, ] fname, ...)  */
static VALUE
s_Molecule_LoadSave(int argc, VALUE *argv, VALUE self, int loadFlag)
{
	VALUE rval;
	char *argstr, *methname, *p, *type = "";
	ID mid = 0;
	int i;
	const char *ls = (loadFlag ? "load" : "save");
	int lslen = strlen(ls);

	if (argc == 0)
		return Qnil;

	if (argc == 0 || (argstr = StringValuePtr(argv[0])) == NULL)
		rb_raise(rb_eMolbyError, "the first argument must be either filename or \":filetype\"");
	if (argstr[0] == ':') {
		/*  Call "loadXXX" (or "saveXXX") for type ":XXX"  */
		methname = ALLOC_N(char, lslen + strlen(argstr));
		strcpy(methname, ls);
		strcat(methname, argstr + 1);
		type = argstr + 1;
		for (i = lslen; methname[i] != 0; i++)
			methname[i] = tolower(methname[i]);
		mid = rb_intern(methname);
		xfree(methname);
		argc--;
		argv++;
		rval = rb_funcall2(self, mid, argc, argv);
		if (rval == Qnil)
			goto failure;
		else
			goto success;
	}
	/*  Guess file type from extension  */
	p = strrchr(argstr, '.');
	if (p != NULL) {
		p++;
		type = p;
		for (methname = p; *methname != 0; methname++) {
			if (!isalpha(*methname))
				break;
		}
		if (*methname == 0) {
			methname = ALLOC_N(char, lslen + strlen(p) + 1);
			if (methname == NULL)
				rb_raise(rb_eMolbyError, "Low memory");
			strcpy(methname, ls);
			strcat(methname, p);
			for (i = lslen; methname[i] != 0; i++)
				methname[i] = tolower(methname[i]);
			mid = rb_intern(methname);
			xfree(methname);
			if (loadFlag) {
				if (rb_respond_to(self, mid)) {
					/*  Load: try to call the load procedure only if it is available  */
					rval = rb_funcall2(self, mid, argc, argv);
					if (rval != Qnil)
						goto success;
				}
			} else {
				/*  Save: call the save procedure, and if not found then call 'method_missing'  */
				rval = rb_funcall2(self, mid, argc, argv);
				if (rval != Qnil)
					goto success;
			}
		}
	}
failure:
	rval = rb_str_to_str(argv[0]);
	asprintf(&p, "Failed to %s file %s", (loadFlag ? "load" : "save"), type);
	s_Molecule_RaiseOnLoadSave(1, p, StringValuePtr(rval));
	return Qnil;  /*  Does not reach here  */

success:
	{
		/*  Register the path  */
		Molecule *mol;
	/*	Atom *ap; */
		Data_Get_Struct(self, Molecule, mol);
		MoleculeSetPath(mol, StringValuePtr(argv[0]));
		
		/*  Check if all occupancy factors are zero; if that is the case, then all set to 1.0  */
	/*	for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ap->occupancy != 0.0)
				break;
		}
		if (i == mol->natoms) {
			for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
				ap->occupancy = 1.0;
			}
		} */
	}
	return rval;
}

/*
 *  call-seq:
 *     molload(file, *args)       -> bool
 *
 *  Read a structure from the given file by calling the public method "loadXXX" (XXX is the
 *  file type given by the extension). If this method fails, then all defined (public)
 *  "loadXXX" methods are invoked, and raises an exception if none of them were successful.
 */
static VALUE
s_Molecule_Load(int argc, VALUE *argv, VALUE self)
{
	return s_Molecule_LoadSave(argc, argv, self, 1);
}

/*
 *  call-seq:
 *     molsave(file, *args)       -> bool
 *
 *  Write a structure/coordinate to the given file by calling the public method "saveXXX"
 *  (XXX is the file type given by the extension).
 */
static VALUE
s_Molecule_Save(int argc, VALUE *argv, VALUE self)
{
	return s_Molecule_LoadSave(argc, argv, self, 0);
}

/*
 *  call-seq:
 *     name       -> String
 *
 *  Returns the display name of the molecule. If the molecule has no associated
 *  document, then returns nil.
 */
static VALUE
s_Molecule_Name(VALUE self)
{
    Molecule *mol;
	char buf[1024];
    Data_Get_Struct(self, Molecule, mol);
	MoleculeCallback_displayName(mol, buf, sizeof buf);
	if (buf[0] == 0)
		return Qnil;
	else
		return rb_str_new2(buf);
}

/*
 *  call-seq:
 *     set_name(string) -> self
 *
 *  Set the name of an untitled molecule. If the molecule is not associated with window
 *  or it already has an associated file, then exception is thrown.
 */
static VALUE
s_Molecule_SetName(VALUE self, VALUE nval)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (MoleculeCallback_setDisplayName(mol, StringValuePtr(nval)))
		rb_raise(rb_eMolbyError, "Cannot change the window title");
	return self;
}


/*
 *  call-seq:
 *     path       -> String
 *
 *  Returns the full path name of the molecule, if it is associated with a file.
 *  If the molecule has no associated file, then returns nil.
 */
static VALUE
s_Molecule_Path(VALUE self)
{
    Molecule *mol;
	char buf[1024];
    Data_Get_Struct(self, Molecule, mol);
	MoleculeCallback_pathName(mol, buf, sizeof buf);
	if (buf[0] == 0)
		return Qnil;
	else
		return Ruby_NewFileStringValue(buf);
}

/*
 *  call-seq:
 *     dir       -> String
 *
 *  Returns the full path name of the directory in which the file associated with the
 *  molecule is located. If the molecule has no associated file, then returns nil.
 */
static VALUE
s_Molecule_Dir(VALUE self)
{
    Molecule *mol;
	char buf[1024], *p;
    Data_Get_Struct(self, Molecule, mol);
	MoleculeCallback_pathName(mol, buf, sizeof buf);
#if __WXMSW__
	translate_char(buf, '\\', '/');
#endif
	if (buf[0] == 0)
		return Qnil;
	else {
		p = strrchr(buf, '/');
		if (p != NULL)
			*p = 0;
		return rb_str_new2(buf);
	}
}

/*
 *  call-seq:
 *     inspect       -> String
 *
 *  Returns a string in the form "Molecule[name]" if the molecule has the associated
 *  document. Otherwise, a string "<Molecule:0x****>" (the address is the address of
 *  the Molecule structure) is returned.
 */
static VALUE
s_Molecule_Inspect(VALUE self)
{
    Molecule *mol;
	char buf[256];
    Data_Get_Struct(self, Molecule, mol);
	MoleculeCallback_displayName(mol, buf, sizeof buf);
	if (buf[0] == 0) {
		/*  No associated document  */
		snprintf(buf, sizeof buf, "#<Molecule:0x%lx>", self);
		return rb_str_new2(buf);
	} else {
		/*  Check whether the document name is duplicate  */
		char buf2[256];
		int idx, k, k2;
		Molecule *mol2;
		for (idx = 0, k = k2 = 0; (mol2 = MoleculeCallback_moleculeAtIndex(idx)) != NULL; idx++) {
			MoleculeCallback_displayName(mol2, buf2, sizeof buf2);
			if (strcmp(buf, buf2) == 0) {
				k++;
				if (mol == mol2)
					k2 = k;
			}
		}
		if (k > 1) {
			snprintf(buf2, sizeof buf2, "Molecule[\"%s\",%d]", buf, k2);
		} else {
			snprintf(buf2, sizeof buf2, "Molecule[\"%s\"]", buf);
		}
		return rb_str_new2(buf2);
	}
}

/*
 *  call-seq:
 *     open        -> Molecule
 *     open(file)  -> Molecule
 *
 *  Create a new molecule from file as a document. If file is not given, an untitled document is created.
 */
static VALUE
s_Molecule_Open(int argc, VALUE *argv, VALUE self)
{
	VALUE fname;
	const char *p;
	Molecule *mp;
	VALUE iflag;
	rb_scan_args(argc, argv, "01", &fname);
	if (NIL_P(fname))
		p = NULL;
	else
		p = FileStringValuePtr(fname);
	iflag = Ruby_SetInterruptFlag(Qfalse);
	mp = MoleculeCallback_openNewMolecule(p);
	Ruby_SetInterruptFlag(iflag);
	if (mp == NULL) {
		if (p == NULL)
			rb_raise(rb_eMolbyError, "Cannot create untitled document");
		else
			rb_raise(rb_eMolbyError, "Cannot open the file %s", p);
	}
	return ValueFromMolecule(mp);
}

/*
 *  call-seq:
 *     new  -> Molecule
 *     new(file, *args)  -> Molecule
 *
 *  Create a new molecule and call "load" method with the same arguments.
 */
static VALUE
s_Molecule_Initialize(int argc, VALUE *argv, VALUE self)
{
	if (argc > 0)
		return s_Molecule_Load(argc, argv, self);
	else return Qnil;  /*  An empty molecule (which is prepared in s_Molecule_Alloc()) is returned  */
}

static VALUE
s_Molecule_MolEnumerable(VALUE self, int kind)
{
    Molecule *mol;
	MolEnumerable *mseq;
    Data_Get_Struct(self, Molecule, mol);
	mseq = MolEnumerableNew(mol, kind);
	return Data_Wrap_Struct(rb_cMolEnumerable, 0, (void (*)(void *))MolEnumerableRelease, mseq);
}

/*
 *  call-seq:
 *     atoms       -> MolEnumerable
 *
 *  Returns a MolEnumerable object representing the array of atoms.
 */
static VALUE
s_Molecule_Atoms(VALUE self)
{
	return s_Molecule_MolEnumerable(self, kAtomKind);
}

/*
 *  call-seq:
 *     bonds       -> MolEnumerable
 *
 *  Returns a MolEnumerable object representing the array of bonds. A bond is represented
 *  by an array of two atom indices.
 */
static VALUE
s_Molecule_Bonds(VALUE self)
{
	return s_Molecule_MolEnumerable(self, kBondKind);
}

/*
 *  call-seq:
 *     angles       -> MolEnumerable
 *
 *  Returns a MolEnumerable object representing the array of angles. An angle is represented
 *  by an array of three atom indices.
 */
static VALUE
s_Molecule_Angles(VALUE self)
{
	return s_Molecule_MolEnumerable(self, kAngleKind);
}

/*
 *  call-seq:
 *     dihedrals       -> MolEnumerable
 *
 *  Returns a MolEnumerable object representing the array of dihedrals. A dihedral is represented
 *  by an array of four atom indices.
 */
static VALUE
s_Molecule_Dihedrals(VALUE self)
{
	return s_Molecule_MolEnumerable(self, kDihedralKind);
}

/*
 *  call-seq:
 *     impropers       -> MolEnumerable
 *
 *  Returns a MolEnumerable object representing the array of impropers. An improper is represented
 *  by an array of four atom indices.
 */
static VALUE
s_Molecule_Impropers(VALUE self)
{
	return s_Molecule_MolEnumerable(self, kImproperKind);
}

/*
 *  call-seq:
 *     residues       -> MolEnumerable
 *
 *  Returns a MolEnumerable object representing the array of residue names.
 */
static VALUE
s_Molecule_Residues(VALUE self)
{
	return s_Molecule_MolEnumerable(self, kResidueKind);
}

/*
 *  call-seq:
 *     natoms       -> Integer
 *
 *  Returns the number of atoms.
 */
static VALUE
s_Molecule_Natoms(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->natoms);
}

/*
 *  call-seq:
 *     nbonds       -> Integer
 *
 *  Returns the number of bonds.
 */
static VALUE
s_Molecule_Nbonds(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->nbonds);
}

/*
 *  call-seq:
 *     nangles       -> Integer
 *
 *  Returns the number of angles.
 */
static VALUE
s_Molecule_Nangles(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->nangles);
}

/*
 *  call-seq:
 *     ndihedrals       -> Integer
 *
 *  Returns the number of dihedrals.
 */
static VALUE
s_Molecule_Ndihedrals(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->ndihedrals);
}

/*
 *  call-seq:
 *     nimpropers       -> Integer
 *
 *  Returns the number of impropers.
 */
static VALUE
s_Molecule_Nimpropers(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->nimpropers);
}

/*
 *  call-seq:
 *     nresidues       -> Integer
 *
 *  Returns the number of residues.
 */
static VALUE
s_Molecule_Nresidues(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->nresidues);
}

static VALUE
s_Molecule_BondParIsObsolete(VALUE self, VALUE val)
{
	rb_raise(rb_eMolbyError, "Molecule#bond_par, angle_par, dihedral_par, improper_par, vdw_par are now obsolete. You can use MDArena#bond_par, angle_par, dihedral_par, improper_par, vdw_par instead, and probably these are what you really want.");
}

/*
 *  call-seq:
 *     bond_par(idx)    -> ParameterRef
 *
 *  Returns the MD parameter for the idx-th bond.
 */
/*
static VALUE
s_Molecule_BondPar(VALUE self, VALUE val)
{
    Molecule *mol;
	BondPar *bp;
	UInt t1, t2;
	Int i1, i2;
	Int ival;
    Data_Get_Struct(self, Molecule, mol);
	ival = NUM2INT(rb_Integer(val));
	if (ival < -mol->nbonds || ival >= mol->nbonds)
		rb_raise(rb_eMolbyError, "bond index (%d) out of range", ival);
	if (ival < 0)
		ival += mol->nbonds;
	s_RebuildMDParameterIfNecessary(self, Qtrue);
	i1 = mol->bonds[ival * 2];
	i2 = mol->bonds[ival * 2 + 1];
	t1 = ATOM_AT_INDEX(mol->atoms, i1)->type;
	t2 = ATOM_AT_INDEX(mol->atoms, i2)->type;
	bp = ParameterLookupBondPar(mol->par, t1, t2, i1, i2, 0);
	if (bp == NULL)
		return Qnil;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kBondParType, bp - mol->par->bondPars);
}
*/

/*
 *  call-seq:
 *     angle_par(idx)    -> ParameterRef
 *
 *  Returns the MD parameter for the idx-th angle.
 */
/*
static VALUE
s_Molecule_AnglePar(VALUE self, VALUE val)
{
    Molecule *mol;
	AnglePar *ap;
	UInt t1, t2, t3;
	Int i1, i2, i3;
	Int ival;
    Data_Get_Struct(self, Molecule, mol);
	ival = NUM2INT(rb_Integer(val));
	if (ival < -mol->nangles || ival >= mol->nangles)
		rb_raise(rb_eMolbyError, "angle index (%d) out of range", ival);
	if (ival < 0)
		ival += mol->nangles;
	s_RebuildMDParameterIfNecessary(self, Qtrue);
	i1 = mol->angles[ival * 3];
	i2 = mol->angles[ival * 3 + 1];
	i3 = mol->angles[ival * 3 + 2];
	t1 = ATOM_AT_INDEX(mol->atoms, i1)->type;
	t2 = ATOM_AT_INDEX(mol->atoms, i2)->type;
	t3 = ATOM_AT_INDEX(mol->atoms, i3)->type;
	ap = ParameterLookupAnglePar(mol->par, t1, t2, t3, i1, i2, i3, 0);
	if (ap == NULL)
		return Qnil;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kAngleParType, ap - mol->par->anglePars);
}
*/
/*
 *  call-seq:
 *     dihedral_par(idx)    -> ParameterRef
 *
 *  Returns the MD parameter for the idx-th dihedral.
 */
/*
static VALUE
s_Molecule_DihedralPar(VALUE self, VALUE val)
{
    Molecule *mol;
	Int ival;
	TorsionPar *tp;
	UInt t1, t2, t3, t4;
	Int i1, i2, i3, i4;
    Data_Get_Struct(self, Molecule, mol);
	ival = NUM2INT(rb_Integer(val));
	if (ival < -mol->ndihedrals || ival >= mol->ndihedrals)
		rb_raise(rb_eMolbyError, "dihedral index (%d) out of range", ival);
	if (ival < 0)
		ival += mol->ndihedrals;
	s_RebuildMDParameterIfNecessary(self, Qtrue);
	i1 = mol->dihedrals[ival * 4];
	i2 = mol->dihedrals[ival * 4 + 1];
	i3 = mol->dihedrals[ival * 4 + 2];
	i4 = mol->dihedrals[ival * 4 + 3];
	t1 = ATOM_AT_INDEX(mol->atoms, i1)->type;
	t2 = ATOM_AT_INDEX(mol->atoms, i2)->type;
	t3 = ATOM_AT_INDEX(mol->atoms, i3)->type;
	t4 = ATOM_AT_INDEX(mol->atoms, i4)->type;
	tp = ParameterLookupDihedralPar(mol->par, t1, t2, t3, t4, i1, i2, i3, i4, 0);
	if (tp == NULL)
		return Qnil;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kDihedralParType, tp - mol->par->dihedralPars);
}
*/
/*
 *  call-seq:
 *     improper_par(idx)    -> ParameterRef
 *
 *  Returns the MD parameter for the idx-th improper.
 */
/*
static VALUE
s_Molecule_ImproperPar(VALUE self, VALUE val)
{
    Molecule *mol;
	Int ival;
	TorsionPar *tp;
	UInt t1, t2, t3, t4;
	Int i1, i2, i3, i4;
    Data_Get_Struct(self, Molecule, mol);
	ival = NUM2INT(rb_Integer(val));
	if (ival < -mol->nimpropers || ival >= mol->nimpropers)
		rb_raise(rb_eMolbyError, "improper index (%d) out of range", ival);
	if (ival < 0)
		ival += mol->nimpropers;
	s_RebuildMDParameterIfNecessary(self, Qtrue);
	i1 = mol->impropers[ival * 4];
	i2 = mol->impropers[ival * 4 + 1];
	i3 = mol->impropers[ival * 4 + 2];
	i4 = mol->impropers[ival * 4 + 3];
	t1 = ATOM_AT_INDEX(mol->atoms, i1)->type;
	t2 = ATOM_AT_INDEX(mol->atoms, i2)->type;
	t3 = ATOM_AT_INDEX(mol->atoms, i3)->type;
	t4 = ATOM_AT_INDEX(mol->atoms, i4)->type;
	tp = ParameterLookupImproperPar(mol->par, t1, t2, t3, t4, i1, i2, i3, i4, 0);
	if (tp == NULL)
		return Qnil;
	return ValueFromMoleculeWithParameterTypeAndIndex(mol, kImproperParType, tp - mol->par->improperPars);
}
*/

/*
 *  call-seq:
 *     start_step       -> Integer
 *
 *  Returns the start step (defined by dcd format).
 */
static VALUE
s_Molecule_StartStep(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->startStep);
}

/*
 *  call-seq:
 *     start_step = Integer
 *
 *  Set the start step (defined by dcd format).
 */
static VALUE
s_Molecule_SetStartStep(VALUE self, VALUE val)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	mol->startStep = NUM2INT(rb_Integer(val));
	return val;
}

/*
 *  call-seq:
 *     steps_per_frame       -> Integer
 *
 *  Returns the number of steps between frames (defined by dcd format).
 */
static VALUE
s_Molecule_StepsPerFrame(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->stepsPerFrame);
}

/*
 *  call-seq:
 *     steps_per_frame = Integer
 *
 *  Set the number of steps between frames (defined by dcd format).
 */
static VALUE
s_Molecule_SetStepsPerFrame(VALUE self, VALUE val)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	mol->stepsPerFrame = NUM2INT(rb_Integer(val));
	return val;
}

/*
 *  call-seq:
 *     ps_per_step       -> Float
 *
 *  Returns the time increment (in picoseconds) for one step (defined by dcd format).
 */
static VALUE
s_Molecule_PsPerStep(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return rb_float_new(mol->psPerStep);
}

/*
 *  call-seq:
 *     ps_per_step = Float
 *
 *  Set the time increment (in picoseconds) for one step (defined by dcd format).
 */
static VALUE
s_Molecule_SetPsPerStep(VALUE self, VALUE val)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	mol->psPerStep = NUM2DBL(rb_Float(val));
	return val;
}

/*
 *  call-seq:
 *     find_angles     -> Integer
 *
 *  Find the angles from the bonds. Returns the number of angles newly created.
 */
/*
static VALUE
s_Molecule_FindAngles(VALUE self)
{
    Molecule *mol;
	Atom *ap;
	int n1, i, j, nc;
	Int *ip, nip, n[3];
    Data_Get_Struct(self, Molecule, mol);
	if (mol == NULL || mol->natoms == 0)
		return INT2NUM(0);
	ip = NULL;
	nip = 0;
	for (n1 = 0, ap = mol->atoms; n1 < mol->natoms; n1++, ap = ATOM_NEXT(ap)) {
		nc = ap->connect.count;
		n[1] = n1;
		for (i = 0; i < nc; i++) {
			n[0] = ap->connects[i];
			for (j = i + 1; j < nc; j++) {
				n[2] = ap->connects[j];
				if (MoleculeLookupAngle(mol, n[0], n[1], n[2]) < 0)
					AssignArray(&ip, &nip, sizeof(Int) * 3, nip, n);
			}
		}
	}
	if (nip > 0) {
		MolActionCreateAndPerform(mol, gMolActionAddAngles, nip * 3, ip, NULL);		
		free(ip);
	}
	return INT2NUM(nip);
}
*/
/*
 *  call-seq:
 *     find_dihedrals     -> Integer
 *
 *  Find the dihedrals from the bonds. Returns the number of dihedrals newly created.
 */
/*
static VALUE
s_Molecule_FindDihedrals(VALUE self)
{
    Molecule *mol;
	Atom *ap1, *ap2;
	int n1, i, j, k, nc1, nc2;
	Int *ip, nip, n[4];
    Data_Get_Struct(self, Molecule, mol);
	if (mol == NULL || mol->natoms == 0)
		return INT2NUM(0);
	ip = NULL;
	nip = 0;
	for (n1 = 0, ap1 = mol->atoms; n1 < mol->natoms; n1++, ap1 = ATOM_NEXT(ap1)) {
		nc1 = ap1->connect.count;
		n[1] = n1;
		for (i = 0; i < nc1; i++) {
			n[2] = ap1->connects[i];
			if (n[1] > n[2])
				continue;
			ap2 = ATOM_AT_INDEX(mol->atoms, n[2]);
			nc2 = ap2->connect.count;
			for (j = 0; j < nc1; j++) {
				n[0] = ap1->connects[j];
				if (n[0] == n[2])
					continue;
				for (k = 0; k < nc2; k++) {
					n[3] = ap2->connects[k];
					if (n[3] == n1 || n[3] == n[0])
						continue;
					if (MoleculeLookupDihedral(mol, n[0], n[1], n[2], n[3]) < 0)
						AssignArray(&ip, &nip, sizeof(Int) * 4, nip, n);
				}
			}
		}
	}
	if (nip > 0) {
		MolActionCreateAndPerform(mol, gMolActionAddDihedrals, nip * 4, ip, NULL);
		free(ip);
	}
	return INT2NUM(nip);
}
*/

/*
 *  call-seq:
 *     nresidues = Integer
 *
 *  Change the number of residues.
 */
static VALUE
s_Molecule_ChangeNresidues(VALUE self, VALUE val)
{
    Molecule *mol;
	int ival = NUM2INT(val);
    Data_Get_Struct(self, Molecule, mol);
	MolActionCreateAndPerform(mol, gMolActionChangeNumberOfResidues, ival);
	if (ival != mol->nresidues)
		rb_raise(rb_eMolbyError, "Cannot set the number of residues to %d (set to %d)", ival, mol->nresidues);
	return val;
}

/*
 *  call-seq:
 *     max_residue_number(atom_group = nil)     -> Integer
 *
 *  Returns the maximum residue number actually used. If an atom group is given, only
 *  these atoms are examined. If no atom is present, nil is returned.
 */
static VALUE
s_Molecule_MaxResSeq(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE gval;
	int maxSeq;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "01", &gval);
	ig = (gval == Qnil ? NULL : s_Molecule_AtomGroupFromValue(self, gval));
	maxSeq = MoleculeMaximumResidueNumber(mol, ig);
	return (maxSeq >= 0 ? INT2NUM(maxSeq) : Qnil);
}

/*
 *  call-seq:
 *     min_residue_number(atom_group = nil)     -> Integer
 *
 *  Returns the minimum residue number actually used. If an atom group is given, only
 *  these atoms are examined. If no atom is present, nil is returned.
 */
static VALUE
s_Molecule_MinResSeq(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE gval;
	int minSeq;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "01", &gval);
	ig = (gval == Qnil ? NULL : s_Molecule_AtomGroupFromValue(self, gval));
	minSeq = MoleculeMinimumResidueNumber(mol, ig);
	return (minSeq >= 0 ? INT2NUM(minSeq) : Qnil);
}

/*
 *  call-seq:
 *     each_atom(atom_group = nil) {|aref| ...}
 *
 *  Execute the block, with the AtomRef object for each atom as the argument. If an atom
 *  group is given, only these atoms are processed.
 *  If atom_group is nil, this is equivalent to self.atoms.each, except that the return value 
 *  is self (a Molecule object).
 */
static VALUE
s_Molecule_EachAtom(int argc, VALUE *argv, VALUE self)
{
	int i;
    Molecule *mol;
	AtomRef *aref;
	VALUE arval;
	VALUE gval;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "01", &gval);
	ig = (gval == Qnil ? NULL : s_Molecule_AtomGroupFromValue(self, gval));
	arval = ValueFromMoleculeAndIndex(mol, 0);
	Data_Get_Struct(arval, AtomRef, aref);
	for (i = 0; i < mol->natoms; i++) {
		aref->idx = i;
		if (ig == NULL || IntGroupLookup(ig, i, NULL))
			rb_yield(arval);
	}
	if (ig != NULL)
		IntGroupRelease(ig);
    return self;
}

/*
 *  call-seq:
 *     cell     -> [a, b, c, alpha, beta, gamma [, sig_a, sig_b, sig_c, sig_alpha, sig_beta, sig_gamma]]
 *
 *  Returns the unit cell parameters. If cell is not set, returns nil.
 */
static VALUE
s_Molecule_Cell(VALUE self)
{
    Molecule *mol;
	int i;
	VALUE val;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->cell == NULL)
		return Qnil;
	val = rb_ary_new2(6);
	for (i = 0; i < 6; i++)
		rb_ary_push(val, rb_float_new(mol->cell->cell[i]));
	if (mol->cell->has_sigma) {
		for (i = 0; i < 6; i++) {
			rb_ary_push(val, rb_float_new(mol->cell->cellsigma[i]));
		}
	}
	return val;
}

/*
 *  call-seq:
 *     cell = [a, b, c, alpha, beta, gamma [, sig_a, sig_b, sig_c, sig_alpha, sig_beta, sig_gamma]]
 *     set_cell([a, b, c, alpha, beta, gamma[, sig_a, sig_b, sig_c, sig_alpha, sig_beta, sig_gamma]], convert_coord = nil)
 *
 *  Set the unit cell parameters. If the cell value is nil, then clear the current cell.
    If the given argument has 12 or more members, then the second half of the parameters represents the sigma values.
    This operation is undoable.
    Convert_coord is a flag to specify that the coordinates should be transformed so that the fractional coordinates remain the same.
 */
static VALUE
s_Molecule_SetCell(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE val, cval;
	int i, convert_coord, n;
	double d[12];
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "11", &val, &cval);
	if (val == Qnil) {
		n = 0;
	} else {
		int len;
		val = rb_ary_to_ary(val);
		len = RARRAY_LEN(val);
		if (len >= 12) {
			n = 12;
		} else if (len >= 6) {
			n = 6;
		} else rb_raise(rb_eMolbyError, "too few members for cell parameters (6 or 12 required)");
		for (i = 0; i < n; i++)
			d[i] = NUM2DBL(rb_Float((RARRAY_PTR(val))[i]));
	}
	convert_coord = (RTEST(cval) ? 1 : 0);
	MolActionCreateAndPerform(mol, gMolActionSetCell, n, d, convert_coord);
	return val;
}

/*
 *  call-seq:
 *     box -> [avec, bvec, cvec, origin, flags]
 *
 *  Get the unit cell information in the form of a periodic bounding box.
 *  Avec, bvec, cvec, origin are Vector3D objects, and flags is a 3-member array of 
 *  Integers which define whether the system is periodic along the axis.
 *  If no unit cell is defined, nil is returned.
 */
static VALUE
s_Molecule_Box(VALUE self)
{
    Molecule *mol;
	VALUE v[5], val;
    Data_Get_Struct(self, Molecule, mol);
	if (mol == NULL || mol->cell == NULL)
		return Qnil;
	v[0] = ValueFromVector(&(mol->cell->axes[0]));
	v[1] = ValueFromVector(&(mol->cell->axes[1]));
	v[2] = ValueFromVector(&(mol->cell->axes[2]));
	v[3] = ValueFromVector(&(mol->cell->origin));
	v[4] = rb_ary_new3(3, INT2NUM(mol->cell->flags[0]), INT2NUM(mol->cell->flags[1]), INT2NUM(mol->cell->flags[2]));
	val = rb_ary_new4(5, v);
	return val;
}

/*
 *  call-seq:
 *     set_box(avec, bvec, cvec, origin = [0, 0, 0], flags = [1, 1, 1], convert_coordinates = nil)
 *     set_box(d, origin = [0, 0, 0])
 *     set_box
 *
 *  Set the unit cell parameters. Avec, bvec, and cvec can be either a Vector3D or a number.
 If it is a number, the x/y/z axis vector is multiplied with the given number and used
 as the box vector.
 Flags, if present, is a 3-member array of Integers defining whether the system is
 periodic along the axis.
 If convert_coordinates is true, then the coordinates are converted so that the fractional coordinates remain the same.
 In the second form, an isotropic box with cell-length d is set.
 In the third form, the existing box is cleared.
 Note: the sigma of the cell parameters is not cleared unless the periodic box itself is cleared.
 */
static VALUE
s_Molecule_SetBox(VALUE self, VALUE aval)
{
    Molecule *mol;
	VALUE v[6];
	static Vector ax[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	Vector vv[3];
	Vector origin = {0, 0, 0};
	char flags[3];
	Double d;
	int i, convertCoordinates = 0;
    Data_Get_Struct(self, Molecule, mol);
	if (aval == Qnil) {
		MolActionCreateAndPerform(mol, gMolActionClearBox);
		return self;
	}
	aval = rb_ary_to_ary(aval);
	for (i = 0; i < 6; i++) {
		if (i < RARRAY_LEN(aval))
			v[i] = (RARRAY_PTR(aval))[i];
		else v[i] = Qnil;
	}
	if (v[0] == Qnil) {
		MolActionCreateAndPerform(mol, gMolActionClearBox);
		return self;
	}
	if ((v[1] == Qnil || v[2] == Qnil) && rb_obj_is_kind_of(v[0], rb_cNumeric)) {
		d = NUM2DBL(rb_Float(v[0]));
		for (i = 0; i < 3; i++)
			VecScale(vv[i], ax[i], d);
		if (v[1] != Qnil)
			VectorFromValue(v[1], &origin);
		flags[0] = flags[1] = flags[2] = 1;
	} else {
		for (i = 0; i < 3; i++) {
			if (v[i] == Qnil) {
				VecZero(vv[i]);
			} else if (rb_obj_is_kind_of(v[i], rb_cNumeric)) {
				d = NUM2DBL(rb_Float(v[i]));
				VecScale(vv[i], ax[i], d);
			} else {
				VectorFromValue(v[i], &vv[i]);
			}
			flags[i] = (VecLength2(vv[i]) > 0.0);
		}
		if (v[3] != Qnil)
			VectorFromValue(v[3], &origin);
		if (v[4] != Qnil) {
			for (i = 0; i < 3; i++) {
				VALUE val = Ruby_ObjectAtIndex(v[4], i);
				flags[i] = (NUM2INT(rb_Integer(val)) != 0);
			}
		}
		if (RTEST(v[5]))
			convertCoordinates = 1;
	}
	MolActionCreateAndPerform(mol, gMolActionSetBox, &(vv[0]), &(vv[1]), &(vv[2]), &origin, (flags[0] * 4 + flags[1] * 2 + flags[2]), convertCoordinates);
	return self;
}

/*
 *  call-seq:
 *     cell_periodicity -> [n1, n2, n3]
 *
 *  Get flags denoting whether the cell is periodic along the a/b/c axes. If the cell is not defined
 *  nil is returned.
 */
static VALUE
s_Molecule_CellPeriodicity(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->cell == NULL)
		return Qnil;
	return rb_ary_new3(3, INT2FIX((int)mol->cell->flags[0]), INT2FIX((int)mol->cell->flags[1]), INT2FIX((int)mol->cell->flags[2]));
}

/*
 *  call-seq:
 *     self.cell_periodicity = [n1, n2, n3] or Integer or nil
 *     set_cell_periodicity = [n1, n2, n3] or Integer or nil
 *
 *  Set whether the cell is periodic along the a/b/c axes. If an integer is given as an argument,
 *  its bits 2/1/0 (from the lowest) correspond to the a/b/c axes. Nil is equivalent to [0, 0, 0].
 *  If cell is not defined, exception is raised.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_SetCellPeriodicity(VALUE self, VALUE arg)
{
    Molecule *mol;
	Int flag;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->cell == NULL)
		rb_raise(rb_eMolbyError, "periodic cell is not defined");
	if (arg == Qnil)
		flag = 0;
	else if (rb_obj_is_kind_of(arg, rb_cNumeric))
		flag = NUM2INT(rb_Integer(arg));
	else {
		Int i;
		VALUE arg0;
		arg = rb_ary_to_ary(arg);
		flag = 0;
		for (i = 0; i < 3 && i < RARRAY_LEN(arg); i++) {
			arg0 = RARRAY_PTR(arg)[i];
			if (arg0 != Qnil && arg0 != Qfalse && arg0 != INT2FIX(0))
				flag |= (1 << (2 - i));
		}
	}
	MolActionCreateAndPerform(mol, gMolActionSetCellPeriodicity, flag);
	return arg;
}

/*
 *  call-seq:
 *     cell_flexibility -> bool
 *
 *  Returns the unit cell is flexible or not
 */
static VALUE
s_Molecule_CellFlexibility(VALUE self)
{
	rb_warn("cell_flexibility is obsolete (unit cell is always frame dependent)");
	return Qtrue;
/*    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->cell == NULL)
		return Qfalse;
	if (mol->useFlexibleCell)
		return Qtrue;
	else return Qfalse; */
}

/*
 *  call-seq:
 *     self.cell_flexibility = bool
 *     set_cell_flexibility(bool)
 *
 *  Change the unit cell is flexible or not
 */
static VALUE
s_Molecule_SetCellFlexibility(VALUE self, VALUE arg)
{
	rb_warn("set_cell_flexibility is obsolete (unit cell is always frame dependent)");
	return self;
/*    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	MolActionCreateAndPerform(mol, gMolActionSetCellFlexibility, RTEST(arg) != 0);
	return self; */
}

/*
 *  call-seq:
 *     cell_transform -> Transform
 *
 *  Get the transform matrix that converts internal coordinates to cartesian coordinates.
 *  If cell is not defined, nil is returned.
 */
static VALUE
s_Molecule_CellTransform(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol == NULL || mol->cell == NULL)
		return Qnil;
	return ValueFromTransform(&(mol->cell->tr));
}

/*
 *  call-seq:
 *     symmetry -> Array of Transforms
 *     symmetries -> Array of Transforms
 *
 *  Get the currently defined symmetry operations. If no symmetry operation is defined,
 *  returns an empty array.
 */
static VALUE
s_Molecule_Symmetry(VALUE self)
{
    Molecule *mol;
	VALUE val;
	int i;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->nsyms <= 0)
		return rb_ary_new();
	val = rb_ary_new2(mol->nsyms);
	for (i = 0; i < mol->nsyms; i++) {
		rb_ary_push(val, ValueFromTransform(&mol->syms[i]));
	}
	return val;
}

/*
 *  call-seq:
 *     nsymmetries -> Integer
 *
 *  Get the number of currently defined symmetry operations.
 */
static VALUE
s_Molecule_Nsymmetries(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->nsyms);
}

/*
 *  call-seq:
 *     add_symmetry(Transform) -> Integer
 *
 *  Add a new symmetry operation. If no symmetry operation is defined and the
 *  given argument is not an identity transform, then also add an identity
 *  transform at the index 0.
 *  Returns the total number of symmetries after operation.
 */
static VALUE
s_Molecule_AddSymmetry(VALUE self, VALUE trans)
{
    Molecule *mol;
	Transform tr;
    Data_Get_Struct(self, Molecule, mol);
	TransformFromValue(trans, &tr);
	MolActionCreateAndPerform(mol, gMolActionAddSymmetryOperation, &tr);
	return INT2NUM(mol->nsyms);
}

/*
 *  call-seq:
 *     remove_symmetry(count = nil) -> Integer
 *     remove_symmetries(count = nil) -> Integer
 *
 *  Remove the specified number of symmetry operations. The last added ones are removed
 *  first. If count is nil, then all symmetry operations are removed. Returns the
 *  number of leftover symmetries.
 */
static VALUE
s_Molecule_RemoveSymmetry(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE cval;
	int i, n;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "01", &cval);
	if (cval == Qnil)
		n = mol->nsyms - 1;
	else {
		n = NUM2INT(rb_Integer(cval));
		if (n < 0 || n > mol->nsyms)
			rb_raise(rb_eMolbyError, "the given count of symops is out of range");
		if (n == mol->nsyms)
			n = mol->nsyms - 1;
	}
	for (i = 0; i < n; i++)
		MolActionCreateAndPerform(mol, gMolActionDeleteSymmetryOperation);
	return INT2NUM(mol->nsyms);
}

static VALUE
s_Molecule_AtomGroup_i(VALUE arg, VALUE values)
{
	Molecule *mol = (Molecule *)(((VALUE *)values)[0]);
	IntGroup *ig1 = (IntGroup *)(((VALUE *)values)[1]);
	int idx = s_Molecule_AtomIndexFromValue(mol, arg);
	IntGroup_RaiseIfError(IntGroupAdd(ig1, idx, 1));
	return Qnil;
}

/*
 *  call-seq:
 *     atom_group
 *     atom_group {|aref| ...}
 *     atom_group(arg1, arg2, ...)
 *     atom_group(arg1, arg2, ...) {|aref| ...}
 *
 *  Specify a group of atoms. If no arguments are given, IntGroup\[0...natoms] is the result.
 *  If arguments are given, then the atoms reprensented by the arguments are added to the
 *  group. For a conversion of a string to an atom index, see the description
 *  of Molecule#atom_index.
 *  If a block is given, it is evaluated with an AtomRef (not atom index integers)
 *  representing each atom, and the atoms are removed from the result if the block returns false.
 *
 */
static VALUE
s_Molecule_AtomGroup(int argc, VALUE *argv, VALUE self)
{
	IntGroup *ig1, *ig2;
    Molecule *mol;
	Int i, startPt, interval;
	VALUE retval = IntGroup_Alloc(rb_cIntGroup);
	Data_Get_Struct(retval, IntGroup, ig1);
    Data_Get_Struct(self, Molecule, mol);
	if (argc == 0) {
		IntGroup_RaiseIfError(IntGroupAdd(ig1, 0, mol->natoms));
	} else {
		while (argc > 0) {
			if (FIXNUM_P(*argv) || TYPE(*argv) == T_STRING) {
				i = s_Molecule_AtomIndexFromValue(mol, *argv);
				IntGroup_RaiseIfError(IntGroupAdd(ig1, i, 1));
			} else if (rb_obj_is_kind_of(*argv, rb_cIntGroup)) {
				ig2 = IntGroupFromValue(*argv);
				for (i = 0; (startPt = IntGroupGetStartPoint(ig2, i)) >= 0; i++) {
					interval = IntGroupGetInterval(ig2, i);
					IntGroup_RaiseIfError(IntGroupAdd(ig1, startPt, interval));
				}
				IntGroupRelease(ig2);
			} else if (rb_respond_to(*argv, rb_intern("each"))) {
				VALUE values[2];
				values[0] = (VALUE)mol;
				values[1] = (VALUE)ig1;
				rb_iterate(rb_each, *argv, s_Molecule_AtomGroup_i, (VALUE)values);
			} else
				IntGroup_RaiseIfError(IntGroupAdd(ig1, NUM2INT(*argv), 1));
			argc--;
			argv++;
		}
	}
	if (rb_block_given_p()) {
		/*  Evaluate the given block with an AtomRef as the argument, and delete
			the index if the block returns false  */
		AtomRef *aref = AtomRefNew(mol, 0);
		VALUE arval = Data_Wrap_Struct(rb_cAtomRef, 0, (void (*)(void *))AtomRefRelease, aref);
		ig2 = IntGroupNew();
		IntGroupCopy(ig2, ig1);
		for (i = 0; (startPt = IntGroupGetNthPoint(ig2, i)) >= 0; i++) {
			VALUE resval;
			if (startPt >= mol->natoms)
				break;
			aref->idx = startPt;
			resval = rb_yield(arval);
			if (!RTEST(resval))
				IntGroupRemove(ig1, startPt, 1);
		}
		IntGroupRelease(ig2);
	}
	
	/*  Remove points that are out of bounds */
	IntGroup_RaiseIfError(IntGroupRemove(ig1, mol->natoms, INT_MAX));

	return retval;			
}

/*
 *  call-seq:
 *     atom_index(val)       -> Integer
 *
 *  Returns the atom index represented by val. val can be either a non-negative integer
 *  (directly representing the atom index), a negative integer (representing <code>natoms - val</code>),
 *  a string of type "resname.resid:name" or "resname:name" or "resid:name" or "name", 
 *  where resname, resid, name are the residue name, residue id, and atom name respectively.
 *  If val is a string and multiple atoms match the description, the atom with the lowest index
 *  is returned.
 */
static VALUE
s_Molecule_AtomIndex(VALUE self, VALUE val)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(s_Molecule_AtomIndexFromValue(mol, val));
}

/*
 *  call-seq:
 *     extract(group, dummy_flag = nil)       -> Molecule
 *
 *  Extract the atoms given by group and return as a new molecule object.
 *  If dummy_flag is true, then the atoms that are not included in the group but are connected
 *  to any atoms in the group are converted to "dummy" atoms (i.e. with element "Du" and 
 *  names beginning with an underscore) and included in the new molecule object.
 */
static VALUE
s_Molecule_Extract(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol1, *mol2;
	IntGroup *ig;
	VALUE group, dummy_flag, retval;
    Data_Get_Struct(self, Molecule, mol1);
	rb_scan_args(argc, argv, "11", &group, &dummy_flag);
	ig = s_Molecule_AtomGroupFromValue(self, group);
	if (MoleculeExtract(mol1, &mol2, ig, (dummy_flag != Qnil && dummy_flag != Qfalse)) != 0) {
		retval = Qnil;
	} else {
		retval = ValueFromMolecule(mol2);
	}
	IntGroupRelease(ig);
	return retval;
}

/*
 *  call-seq:
 *     add(molecule2)       -> self
 *
 *  Combine two molecules. The residue numbers of the newly added atoms may be renumbered to avoid
    conflicts.
    This operation is undoable.
 */
static VALUE
s_Molecule_Add(VALUE self, VALUE val)
{
    Molecule *mol1, *mol2;
    Data_Get_Struct(self, Molecule, mol1);
	mol2 = MoleculeFromValue(val);
	MolActionCreateAndPerform(mol1, gMolActionMergeMolecule, mol2, NULL);
	return self; 
}

/*
 *  call-seq:
 *     remove(group)       -> Molecule
 *
 *  The atoms designated by the given group are removed from the molecule.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_Remove(VALUE self, VALUE group)
{
    Molecule *mol1;
	IntGroup *ig, *bg;
	Int i;
	IntGroupIterator iter;

    Data_Get_Struct(self, Molecule, mol1);
	group = rb_funcall(self, rb_intern("atom_group"), 1, group);
	if (!rb_obj_is_kind_of(group, rb_cIntGroup))
		rb_raise(rb_eMolbyError, "IntGroup instance is expected");
	Data_Get_Struct(group, IntGroup, ig);

	/*  Remove the bonds between the two fragments  */
	/*  (This is necessary for undo to work correctly)  */
	IntGroupIteratorInit(ig, &iter);
	bg = NULL;
	while ((i = IntGroupIteratorNext(&iter)) >= 0) {
		Atom *ap = ATOM_AT_INDEX(mol1->atoms, i);
		Int j, *cp;
		cp = AtomConnectData(&ap->connect);
		for (j = 0; j < ap->connect.count; j++) {
			int n = cp[j];
			if (!IntGroupLookup(ig, n, NULL)) {
				/*  bond i-n, i is in ig and n is not  */
				int k = MoleculeLookupBond(mol1, i, n);
				if (k >= 0) {
					if (bg == NULL)
						bg = IntGroupNew();
					IntGroupAdd(bg, k, 1);
				}
			}
		}
	}
	IntGroupIteratorRelease(&iter);
	if (bg != NULL) {
		/*  Remove bonds  */
		MolActionCreateAndPerform(mol1, gMolActionDeleteBonds, bg);
		IntGroupRelease(bg);
	}
	/*  Remove atoms  */
	if (MolActionCreateAndPerform(mol1, gMolActionUnmergeMolecule, ig) == 0)
		return Qnil;
	return self;
}

/*
 *  call-seq:
 *     create_atom(name, pos = -1)  -> AtomRef
 *
 *  Create a new atom with the specified name (may contain residue 
 *  information) and position (if position is out of range, the atom is appended at
 *  the end). Returns the reference to the new atom.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_CreateAnAtom(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
    Int i, pos;
	VALUE name, ival;
    Atom arec;
    AtomRef *aref;
	char *p, resName[6], atomName[6];
	int resSeq;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "02", &name, &ival);
	if (ival != Qnil)
		pos = NUM2INT(rb_Integer(ival));
	else pos = -1;
	if (name != Qnil) {
		p = StringValuePtr(name);
		if (p[0] != 0) {
			i = MoleculeAnalyzeAtomName(p, resName, &resSeq, atomName);
			if (atomName[0] == 0)
			  rb_raise(rb_eMolbyError, "bad atom name specification: %s", p);
		}
	} else p = NULL;
	if (p == NULL || p[0] == 0) {
		memset(atomName, 0, 4);
		resSeq = -1;
	}
    memset(&arec, 0, sizeof(arec));
    strncpy(arec.aname, atomName, 4);
    if (resSeq >= 0) {
      strncpy(arec.resName, resName, 4);
      arec.resSeq = resSeq;
    }
	arec.occupancy = 1.0;
//    i = MoleculeCreateAnAtom(mol, &arec);
	if (MolActionCreateAndPerform(mol, gMolActionAddAnAtom, &arec, pos, &pos) != 0)
		return Qnil;
    aref = AtomRefNew(mol, pos);
    return Data_Wrap_Struct(rb_cAtomRef, 0, (void (*)(void *))AtomRefRelease, aref);
}

/*
 *  call-seq:
 *     duplicate_atom(atomref, pos = -1)  -> AtomRef
 *
 *  Create a new atom with the same attributes (but no bonding information)
 *  with the specified atom. Returns the reference to the new atom.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_DuplicateAnAtom(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	const Atom *apsrc;
    Atom arec;
	AtomRef *aref;
	VALUE retval, aval, ival;
	Int pos;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "11", &aval, &ival);
	if (FIXNUM_P(aval)) {
		int idx = NUM2INT(aval);
		if (idx < 0 || idx >= mol->natoms)
			rb_raise(rb_eMolbyError, "atom index out of range: %d", idx);
		apsrc = ATOM_AT_INDEX(mol->atoms, idx);
	} else {
		apsrc = s_AtomFromValue(aval);
	}
	if (apsrc == NULL)
		rb_raise(rb_eMolbyError, "bad atom specification");
	if (ival != Qnil)
		pos = NUM2INT(rb_Integer(ival));
	else pos = -1;
	AtomDuplicate(&arec, apsrc);
	arec.connect.count = 0;
	if (MolActionCreateAndPerform(mol, gMolActionAddAnAtom, &arec, pos, &pos) != 0)
		retval = Qnil;
	else {
		aref = AtomRefNew(mol, pos);
		retval = Data_Wrap_Struct(rb_cAtomRef, 0, (void (*)(void *))AtomRefRelease, aref);
	}
	AtomClean(&arec);
	return retval;
}

/*
 *  call-seq:
 *     create_bond(n1, n2, ...)       -> Integer
 *
 *  Create bonds between atoms n1 and n2, n3 and n4, and so on. If the corresponding bond is already present for a particular pair,
 *  do nothing for that pair. Returns the number of bonds actually created.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_CreateBond(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	Int i, j, k, *ip, old_nbonds;
	if (argc == 0)
		rb_raise(rb_eMolbyError, "missing arguments");
	if (argc % 2 != 0)
		rb_raise(rb_eMolbyError, "bonds should be specified by pairs of atom indices");
    Data_Get_Struct(self, Molecule, mol);
	ip = ALLOC_N(Int, argc + 1);
	for (i = j = 0; i < argc; i++, j++) {
		ip[j] = s_Molecule_AtomIndexFromValue(mol, argv[i]);
		if (i % 2 == 1) {
			if (MoleculeLookupBond(mol, ip[j - 1], ip[j]) >= 0)
				j -= 2;  /*  This bond is already present: skip it  */
			else {
				for (k = 0; k < j - 1; k += 2) {
					if ((ip[k] == ip[j - 1] && ip[k + 1] == ip[j]) || (ip[k + 1] == ip[j - 1] && ip[k] == ip[j])) {
						j -= 2;   /*  The same entry is already in the argument  */
						break;
					}
				}
			}
		}
	}
	old_nbonds = mol->nbonds;
	if (j > 0) {
		ip[j] = kInvalidIndex;
		i = MolActionCreateAndPerform(mol, gMolActionAddBonds, j, ip, NULL);
	} else i = 0;
	xfree(ip);
	if (i == -1)
		rb_raise(rb_eMolbyError, "atom index out of range");
	else if (i == -2)
		rb_raise(rb_eMolbyError, "too many bonds");
	else if (i == -3)
		rb_raise(rb_eMolbyError, "duplicate bonds");
	else if (i == -5)
		rb_raise(rb_eMolbyError, "cannot create bond to itself");
	else if (i != 0)
		rb_raise(rb_eMolbyError, "error in creating bonds");
	return INT2NUM(mol->nbonds - old_nbonds);
}

/*
 *  call-seq:
 *     molecule.remove_bonds(n1, n2, ...)       -> Integer
 *
 *  Remove bonds between atoms n1 and n2, n3 and n4, and so on. If the corresponding bond is not present for
 *  a particular pair, do nothing for that pair. Returns the number of bonds actually removed.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_RemoveBond(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	Int i, j, n[2];
	IntGroup *bg;
	if (argc == 0)
		rb_raise(rb_eMolbyError, "missing arguments");
	if (argc % 2 != 0)
		rb_raise(rb_eMolbyError, "bonds should be specified by pairs of atom indices");
    Data_Get_Struct(self, Molecule, mol);
	bg = NULL;
	for (i = j = 0; i < argc; i++, j = 1 - j) {
		n[j] = s_Molecule_AtomIndexFromValue(mol, argv[i]);
		if (j == 1) {
			Int k = MoleculeLookupBond(mol, n[0], n[1]);
			if (k >= 0) {
				if (bg == NULL)
					bg = IntGroupNew();
				IntGroupAdd(bg, k, 1);
			}
		}
	}
	if (bg != NULL) {
		MolActionCreateAndPerform(mol, gMolActionDeleteBonds, bg);
		i = IntGroupGetCount(bg);
		IntGroupRelease(bg);
	} else i = 0;
	return INT2NUM(i);
}

/*
 *  call-seq:
 *     assign_bond_order(idx, d1)
 *     assign_bond_orders(group, [d1, d2, ...])
 *
 *  Assign bond order. In the first form, the bond order of the idx-th bond is set to d1 (a Float value).
 *  In the second form, the bond orders at the indices in the group are set to d1, d2, etc.
 *  At present, the bond orders are only used in UFF parameter guess, and not in the MM/MD calculations.
 *  (This may change in the future)
 *  This operation is undoable.
 */
static VALUE
s_Molecule_AssignBondOrder(VALUE self, VALUE idxval, VALUE dval)
{
    Molecule *mol;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	if (rb_obj_is_kind_of(idxval, rb_cNumeric)) {
		/*  The first form  */
		Int idx = NUM2INT(rb_Integer(idxval));
		Double d1 = NUM2DBL(rb_Float(dval));
		if (idx < 0 || idx >= mol->nbonds)
			rb_raise(rb_eMolbyError, "the bond index (%d) is out of bounds", idx);
		ig = IntGroupNewWithPoints(idx, 1, -1);
		MolActionCreateAndPerform(mol, gMolActionAssignBondOrders, 1, &d1, ig);
		IntGroupRelease(ig);
	} else {
		Int i, n;
		Double *dp;
		ig = IntGroupFromValue(idxval);
		n = IntGroupGetCount(ig);
		if (n == 0)
			rb_raise(rb_eMolbyError, "the bond index is empty");
		dval = rb_ary_to_ary(dval);
		dp = (Double *)calloc(sizeof(Double), n);
		for (i = 0; i < RARRAY_LEN(dval) && i < n; i++) {
			dp[i] = NUM2DBL(rb_Float(RARRAY_PTR(dval)[i]));
		}
		MolActionCreateAndPerform(mol, gMolActionAssignBondOrders, n, dp, ig);
		free(dp);
		IntGroupRelease(ig);
	}
	return self;
}

/*
 *  call-seq:
 *     get_bond_order(idx) -> Float
 *     get_bond_orders(group) -> Array
 *
 *  Get the bond order. In the first form, the bond order of the idx-th bond is returned.
 *  In the second form, the bond orders at the indices in the group are returned as an array.
 *  If no bond order information have been assigned, returns nil (the first form)
 *  or an empty array (the second form).
 */
static VALUE
s_Molecule_GetBondOrder(VALUE self, VALUE idxval)
{
    Molecule *mol;
	IntGroup *ig;
	Double *dp;
	VALUE retval;
	Int i, n, numericArg;
    Data_Get_Struct(self, Molecule, mol);
	if (rb_obj_is_kind_of(idxval, rb_cNumeric)) {
		/*  The first form  */
		Int idx = NUM2INT(rb_Integer(idxval));
		if (idx < 0 || idx >= mol->nbonds)
			rb_raise(rb_eMolbyError, "the bond index (%d) is out of bounds", idx);
		if (mol->bondOrders == NULL)
			return Qnil;
		ig = IntGroupNewWithPoints(idx, 1, -1);
		n = 1;
		numericArg = 1;
	} else {
		if (mol->bondOrders == NULL)
			return rb_ary_new();
		ig = IntGroupFromValue(idxval);
		n = IntGroupGetCount(ig);
		if (n == 0)
			rb_raise(rb_eMolbyError, "the bond index is empty");
		numericArg = 0;
	}
	dp = (Double *)calloc(sizeof(Double), n);
	MoleculeGetBondOrders(mol, dp, ig);
	if (numericArg)
		retval = rb_float_new(dp[0]);
	else {
		retval = rb_ary_new();
		for (i = 0; i < n; i++)
			rb_ary_push(retval, rb_float_new(dp[i]));
	}
	free(dp);
	IntGroupRelease(ig);
	return retval;
}

/*
 *  call-seq:
 *     bond_exist?(idx1, idx2) -> bool
 *
 *  Returns true if bond exists between atoms idx1 and idx2, otherwise returns false.
 *  Imaginary bonds between a pi-anchor and member atoms are not considered.
 */
static VALUE
s_Molecule_BondExist(VALUE self, VALUE ival1, VALUE ival2)
{
	Molecule *mol;
	Int idx1, idx2, i;
	Atom *ap;
	Int *cp;
    Data_Get_Struct(self, Molecule, mol);
	idx1 = NUM2INT(rb_Integer(ival1));
	idx2 = NUM2INT(rb_Integer(ival2));
	if (idx1 < 0 || idx1 >= mol->natoms || idx2 < 0 || idx2 >= mol->natoms)
		rb_raise(rb_eMolbyError, "Atom index (%d or %d) out of range", idx1, idx2);
	ap = ATOM_AT_INDEX(mol->atoms, idx1);
	cp = AtomConnectData(&ap->connect);
	for (i = 0; i < ap->connect.count; i++) {
		if (cp[i] == idx2)
			return Qtrue;
	}
	return Qfalse;
}

/*
 *  call-seq:
 *     add_angle(n1, n2, n3)       -> Molecule
 *
 *  Add angle n1-n2-n3. Returns self. Usually, angles are automatically added
 *  when a bond is created, so it is rarely necessary to use this method explicitly.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_AddAngle(VALUE self, VALUE v1, VALUE v2, VALUE v3)
{
	Int n[4];
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	n[0] = s_Molecule_AtomIndexFromValue(mol, v1);
	n[1] = s_Molecule_AtomIndexFromValue(mol, v2);
	n[2] = s_Molecule_AtomIndexFromValue(mol, v3);
	if (MoleculeLookupAngle(mol, n[0], n[1], n[2]) >= 0)
		rb_raise(rb_eMolbyError, "angle %d-%d-%d is already present", n[0], n[1], n[2]);
	n[3] = kInvalidIndex;
	MolActionCreateAndPerform(mol, gMolActionAddAngles, 3, n, NULL);
	return self;
}

/*
 *  call-seq:
 *     remove_angle(n1, n2, n3)       -> Molecule
 *
 *  Remove angle n1-n2-n3. Returns self. Usually, angles are automatically removed
 *  when a bond is removed, so it is rarely necessary to use this method explicitly.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_RemoveAngle(VALUE self, VALUE v1, VALUE v2, VALUE v3)
{
	Int n[4];
    Molecule *mol;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	n[0] = s_Molecule_AtomIndexFromValue(mol, v1);
	n[1] = s_Molecule_AtomIndexFromValue(mol, v2);
	n[2] = s_Molecule_AtomIndexFromValue(mol, v3);
	if ((n[3] = MoleculeLookupAngle(mol, n[0], n[1], n[2])) < 0)
		rb_raise(rb_eMolbyError, "angle %d-%d-%d is not present", n[0], n[1], n[2]);
	ig = IntGroupNewWithPoints(n[3], 1, -1);
	MolActionCreateAndPerform(mol, gMolActionDeleteAngles, ig);
	IntGroupRelease(ig);
	return self;
}

/*
 *  call-seq:
 *     add_dihedral(n1, n2, n3, n4)       -> Molecule
 *
 *  Add dihedral n1-n2-n3-n4. Returns self. Usually, dihedrals are automatically added
 *  when a bond is created, so it is rarely necessary to use this method explicitly.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_AddDihedral(VALUE self, VALUE v1, VALUE v2, VALUE v3, VALUE v4)
{
	Int n[5];
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	n[0] = s_Molecule_AtomIndexFromValue(mol, v1);
	n[1] = s_Molecule_AtomIndexFromValue(mol, v2);
	n[2] = s_Molecule_AtomIndexFromValue(mol, v3);
	n[3] = s_Molecule_AtomIndexFromValue(mol, v4);
	if (MoleculeLookupDihedral(mol, n[0], n[1], n[2], n[3]) >= 0)
		rb_raise(rb_eMolbyError, "dihedral %d-%d-%d-%d is already present", n[0], n[1], n[2], n[3]);
	n[4] = kInvalidIndex;
	MolActionCreateAndPerform(mol, gMolActionAddDihedrals, 4, n, NULL);
	return self;
}

/*
 *  call-seq:
 *     remove_dihedral(n1, n2, n3, n4)       -> Molecule
 *
 *  Remove dihedral n1-n2-n3-n4. Returns self. Usually, dihedrals are automatically removed
 *  when a bond is removed, so it is rarely necessary to use this method explicitly.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_RemoveDihedral(VALUE self, VALUE v1, VALUE v2, VALUE v3, VALUE v4)
{
	Int n[5];
    Molecule *mol;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	n[0] = s_Molecule_AtomIndexFromValue(mol, v1);
	n[1] = s_Molecule_AtomIndexFromValue(mol, v2);
	n[2] = s_Molecule_AtomIndexFromValue(mol, v3);
	n[3] = s_Molecule_AtomIndexFromValue(mol, v4);
	if ((n[4] = MoleculeLookupDihedral(mol, n[0], n[1], n[2], n[3])) < 0)
		rb_raise(rb_eMolbyError, "dihedral %d-%d-%d-%d is not present", n[0], n[1], n[2], n[3]);
	ig = IntGroupNewWithPoints(n[4], 1, -1);
	MolActionCreateAndPerform(mol, gMolActionDeleteDihedrals, ig);
	IntGroupRelease(ig);
	return self;
}

/*
 *  call-seq:
 *     add_improper(n1, n2, n3, n4)       -> Molecule
 *
 *  Add dihedral n1-n2-n3-n4. Returns self. Unlike angles and dihedrals, impropers are
 *  not automatically added when a new bond is created, so this method is more useful than
 *  the angle/dihedral counterpart.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_AddImproper(VALUE self, VALUE v1, VALUE v2, VALUE v3, VALUE v4)
{
	Int n[5];
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	n[0] = s_Molecule_AtomIndexFromValue(mol, v1);
	n[1] = s_Molecule_AtomIndexFromValue(mol, v2);
	n[2] = s_Molecule_AtomIndexFromValue(mol, v3);
	n[3] = s_Molecule_AtomIndexFromValue(mol, v4);
	if (MoleculeLookupImproper(mol, n[0], n[1], n[2], n[3]) >= 0)
		rb_raise(rb_eMolbyError, "improper %d-%d-%d-%d is already present", n[0], n[1], n[2], n[3]);
	n[4] = kInvalidIndex;
	MolActionCreateAndPerform(mol, gMolActionAddImpropers, 4, n, NULL);
	return self;
}

/*
 *  call-seq:
 *     remove_improper(n1, n2, n3, n4)       -> Molecule
 *     remove_improper(intgroup)             -> Molecule
 *
 *  Remove improper n1-n2-n3-n4, or the specified impropers (in indices) in IntGroup.
 *  Returns self. Unlike angles and dihedrals, impropers are
 *  not automatically added when a new bond is created, so this method is more useful than
 *  the angle/dihedral counterpart.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_RemoveImproper(int argc, VALUE *argv, VALUE self)
{
	Int n[5];
	VALUE v1, v2, v3, v4;
    Molecule *mol;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	if (argc == 1) {
		ig = IntGroupFromValue(argv[0]);
	} else {
		rb_scan_args(argc, argv, "40", &v1, &v2, &v3, &v4);
		n[0] = s_Molecule_AtomIndexFromValue(mol, v1);
		n[1] = s_Molecule_AtomIndexFromValue(mol, v2);
		n[2] = s_Molecule_AtomIndexFromValue(mol, v3);
		n[3] = s_Molecule_AtomIndexFromValue(mol, v4);
		if ((n[4] = MoleculeLookupImproper(mol, n[0], n[1], n[2], n[3])) < 0)
			rb_raise(rb_eMolbyError, "improper %d-%d-%d-%d is not present", n[0], n[1], n[2], n[3]);
		ig = IntGroupNewWithPoints(n[4], 1, -1);
	}
	MolActionCreateAndPerform(mol, gMolActionDeleteImpropers, ig);
	IntGroupRelease(ig);
	return self;
}

/*
 *  call-seq:
 *     assign_residue(group, res)       -> Molecule
 *
 *  Assign the specified atoms as the given residue. res can either be an integer, "resname"
 *  or "resname.resno". When the residue number is not specified, the residue number of
 *  the first atom in the group is used.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_AssignResidue(VALUE self, VALUE range, VALUE res)
{
    Molecule *mol;
	IntGroup *ig;
	char *p, *pp, buf[16];
	Int resid, n;
	Atom *ap;
    Data_Get_Struct(self, Molecule, mol);
	
	/*  Parse the argument res  */
	if (FIXNUM_P(res)) {
		/*  We can assume Fixnum here because Bignum is non-realistic as residue numbers  */
		resid = NUM2INT(res);
		buf[0] = 0;
	} else {
		p = StringValuePtr(res);
		pp = strchr(p, '.');
		if (pp != NULL) {
			resid = atoi(pp + 1);
			n = pp - p;
		} else {
			resid = -1;
			n = strlen(p);
		}
		if (n > sizeof buf - 1)
			n = sizeof buf - 1;
		strncpy(buf, p, n);
		buf[n] = 0;
	}
	ig = s_Molecule_AtomGroupFromValue(self, range);
	if (ig == NULL || IntGroupGetCount(ig) == 0)
		return Qnil;

	if (resid < 0) {
		/*  Use the residue number of the first specified atom  */
		n = IntGroupGetNthPoint(ig, 0);
		if (n >= mol->natoms)
			rb_raise(rb_eMolbyError, "Atom index (%d) out of range", n);
		ap = ATOM_AT_INDEX(mol->atoms, n);
		resid = ap->resSeq;
	}
	/*  Change the residue number  */
	MolActionCreateAndPerform(mol, gMolActionChangeResidueNumber, ig, resid);
	/*  Change the residue name if necessary  */
	if (buf[0] != 0) {
	/*	Int seqs[2];
		seqs[0] = resid;
		seqs[1] = kInvalidIndex; */
		MolActionCreateAndPerform(mol, gMolActionChangeResidueNames, 1, &resid, 4, buf);
	}
	IntGroupRelease(ig);
	return self;
}

/*
 *  call-seq:
 *     offset_residue(group, offset)       -> Molecule
 *
 *  Offset the residue number of the specified atoms. If any of the residue number gets
 *  negative, then exception is thrown.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_OffsetResidue(VALUE self, VALUE range, VALUE offset)
{
    Molecule *mol;
	IntGroup *ig;
	int ofs, result;
    Data_Get_Struct(self, Molecule, mol);
	ig = s_Molecule_AtomGroupFromValue(self, range);
	ofs = NUM2INT(offset);
	result = MolActionCreateAndPerform(mol, gMolActionOffsetResidueNumbers, ig, ofs, -1);
	if (result > 0)
		rb_raise(rb_eMolbyError, "residue number of atom %d becomes negative", result - 1);
	IntGroupRelease(ig);
	return self;
}

/*
 *  call-seq:
 *     renumber_atoms(array)       -> IntGroup
 *
 *  Change the order of atoms so that the atoms specified in the array argument appear
 *  in this order from the top of the molecule. The atoms that are not included in array
 *  are placed after these atoms, and these atoms are returned as an intGroup.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_RenumberAtoms(VALUE self, VALUE array)
{
    Molecule *mol;
	Int *new2old;
	IntGroup *ig;
	int i, n;
	VALUE *valp, retval;
    Data_Get_Struct(self, Molecule, mol);
	if (TYPE(array) != T_ARRAY)
		array = rb_funcall(array, rb_intern("to_a"), 0);
	n = RARRAY_LEN(array);
	valp = RARRAY_PTR(array);
	new2old = ALLOC_N(Int, n + 1);
	for (i = 0; i < n; i++)
		new2old[i] = s_Molecule_AtomIndexFromValue(mol, valp[i]);
	new2old[i] = kInvalidIndex;
	i = MolActionCreateAndPerform(mol, gMolActionRenumberAtoms, i, new2old);
	if (i == 1)
		rb_raise(rb_eMolbyError, "Atom index out of range");
	else if (i == 2)
		rb_raise(rb_eMolbyError, "Duplicate entry");
	else if (i == 3)
		rb_raise(rb_eMolbyError, "Internal inconsistency during atom renumbering");
	retval = IntGroup_Alloc(rb_cIntGroup);
	Data_Get_Struct(retval, IntGroup, ig);
	if (mol->natoms > n)
		IntGroup_RaiseIfError(IntGroupAdd(ig, n, mol->natoms - n));
	xfree(new2old);
	return retval;
}

/*
 *  call-seq:
 *     find_close_atoms(atom, limit = 1.2, radius = 0.77)   -> array of Integers (atom indices)
 *
 *  Find atoms that are within the threshold distance from the given atom.
 *  (The atom argument can also be a vector, representing a cartesian coordinate. In that case, the van der Waals of the atom can also be specified.)
 *  If limit is a positive number, the threshold distance is the sum of the vdw radii times limit.
 *  If limit is a negative number, its absolute value is used for the threshold distance in angstrom.
 *  If limit is not given, a default value of 1.2 is used.
 *  An array of atom indices is returned. If no atoms are found, an empty array is returned.
 */
static VALUE
s_Molecule_FindCloseAtoms(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE aval, limval, radval;
	double limit, radius;
	Int n1, nbonds, *bonds, an;
	Vector v;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "12", &aval, &limval, &radval);
	if (rb_obj_is_kind_of(aval, rb_cVector3D) || rb_obj_is_kind_of(aval, rb_cLAMatrix) || rb_obj_is_kind_of(aval, rb_mEnumerable)) {
		VectorFromValue(aval, &v);
		if (radval == Qnil)
			radius = gElementParameters[6].radius;
		else
			radius = NUM2DBL(rb_Float(radval));
		n1 = 0;
	} else {
		n1 = s_Molecule_AtomIndexFromValue(mol, aval);
		v = ATOM_AT_INDEX(mol->atoms, n1)->r;
		an = ATOM_AT_INDEX(mol->atoms, n1)->atomicNumber;
		if (an >= 0 && an < gCountElementParameters)
			radius = gElementParameters[an].radius;
		else radius = gElementParameters[6].radius;
	}
	if (limval == Qnil)
		limit = 1.2;
	else
		limit = NUM2DBL(rb_Float(limval));
	nbonds = 0;  /*  This initialization is necessary: see comments in MoleculeFindCloseAtoms()  */
	bonds = NULL;
	MoleculeFindCloseAtoms(mol, &v, radius, limit, &nbonds, &bonds, n1);
	aval = rb_ary_new();
	if (nbonds > 0) {
		for (n1 = 0; n1 < nbonds; n1++)
			rb_ary_push(aval, INT2NUM(bonds[n1 * 2 + 1]));
		free(bonds);
	}
	return aval;
}

/*
 *  call-seq:
 *     guess_bonds(limit = 1.2)       -> Integer
 *
 *  Create bonds between atoms that are within the threshold distance.
 *  If limit is a positive number, the threshold distance is the sum of the vdw radii times limit.
 *  If limit is a negative number, its absolute value is used for the threshold distance in angstrom.
 *  If limit is not given, a default value of 1.2 is used.
 *  The number of the newly created bonds is returned.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_GuessBonds(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE limval;
	double limit;
	Int nbonds, *bonds;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "01", &limval);
	if (limval == Qnil)
		limit = 1.2;
	else
		limit = NUM2DBL(rb_Float(limval));
	MoleculeGuessBonds(mol, limit, &nbonds, &bonds);
	if (nbonds > 0) {
		MolActionCreateAndPerform(mol, gMolActionAddBonds, nbonds * 2, bonds, NULL);
		free(bonds);
	}
	return INT2NUM(nbonds);
}
	
/*
 *  call-seq:
 *     register_undo(script, *args)
 *
 *  Register an undo operation with the current molecule.
 */
static VALUE
s_Molecule_RegisterUndo(int argc, VALUE *argv, VALUE self)
{
	Molecule *mol;
	VALUE script, args;
	MolAction *act;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "1*", &script, &args);
	act = MolActionNew(SCRIPT_ACTION("R"), StringValuePtr(script), args);
	MolActionCallback_registerUndo(mol, act);
	return script;
}

/*
 *  call-seq:
 *     undo_enabled? -> bool
 *
 *  Returns true if undo is enabled for this molecule; otherwise no.
 */
static VALUE
s_Molecule_UndoEnabled(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (MolActionCallback_isUndoRegistrationEnabled(mol))
		return Qtrue;
	else return Qfalse;
}

/*
 *  call-seq:
 *     undo_enabled = bool
 *
 *  Enable or disable undo.
 */
static VALUE
s_Molecule_SetUndoEnabled(VALUE self, VALUE val)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	MolActionCallback_setUndoRegistrationEnabled(mol, (val != Qfalse && val != Qnil));
	return val;
}

/*
 *  call-seq:
 *     selection       -> IntGroup
 *
 *  Returns the current selection.
 */
static VALUE
s_Molecule_Selection(VALUE self)
{
    Molecule *mol;
	IntGroup *ig;
	VALUE val;
    Data_Get_Struct(self, Molecule, mol);
	if (mol != NULL && (ig = MoleculeGetSelection(mol)) != NULL) {
		ig = IntGroupNewFromIntGroup(ig);  /*  Duplicate, so that the change from GUI does not affect the value  */
		val = ValueFromIntGroup(ig);
		IntGroupRelease(ig);
	} else {
		val = IntGroup_Alloc(rb_cIntGroup);
	}
	return val;
}

static VALUE
s_Molecule_SetSelectionSub(VALUE self, VALUE val, int undoable)
{
    Molecule *mol;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	if (val == Qnil)
		ig = NULL;
	else
		ig = s_Molecule_AtomGroupFromValue(self, val);
	if (undoable)
		MolActionCreateAndPerform(mol, gMolActionSetSelection, ig);
	else
		MoleculeSetSelection(mol, ig);
	if (ig != NULL)
		IntGroupRelease(ig);
	return val;
}

/*
 *  call-seq:
 *     selection = IntGroup
 *
 *  Set the current selection. The right-hand operand may be nil.
 *  This operation is _not_ undoable. If you need undo, use set_undoable_selection instead.
 */
static VALUE
s_Molecule_SetSelection(VALUE self, VALUE val)
{
	return s_Molecule_SetSelectionSub(self, val, 0);
}

/*
 *  call-seq:
 *     set_undoable_selection(IntGroup)
 *
 *  Set the current selection with undo registration. The right-hand operand may be nil.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_SetUndoableSelection(VALUE self, VALUE val)
{
	return s_Molecule_SetSelectionSub(self, val, 1);
}

/*
 *  call-seq:
 *     hidden_atoms       -> IntGroup
 *
 *  Returns the currently hidden atoms.
 */
static VALUE
s_Molecule_HiddenAtoms(VALUE self)
{
	rb_raise(rb_eMolbyError, "set_hidden_atoms is now obsolete. Try using Molecule#is_atom_visible or AtomRef#hidden.");
	return Qnil;  /*  Not reached  */
/*    Molecule *mol;
	IntGroup *ig;
	VALUE val;
    Data_Get_Struct(self, Molecule, mol);
	if (mol != NULL) {
		Atom *ap;
		int i;
		ig = IntGroupNew();
		for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ap->exflags & kAtomHiddenFlag)
				IntGroupAdd(ig, i, 1);
		}
		val = ValueFromIntGroup(ig);
		IntGroupRelease(ig);
		rb_obj_freeze(val);
		return val;
	} else return Qnil; */
}

/*
 *  call-seq:
 *     set_hidden_atoms(IntGroup)
 *     self.hidden_atoms = IntGroup
 *
 *  Hide the specified atoms. This operation is _not_ undoable.
 */
static VALUE
s_Molecule_SetHiddenAtoms(VALUE self, VALUE val)
{
	rb_raise(rb_eMolbyError, "set_hidden_atoms is now obsolete. Try using Molecule#is_atom_visible or AtomRef#hidden.");
	return Qnil;  /*  Not reached  */
/*
	Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol != NULL) {
		Atom *ap;
		int i;
		IntGroup *ig;
		if (val == Qnil)
			ig = NULL;
		else
			ig = s_Molecule_AtomGroupFromValue(self, val);
		for (i = 0, ap = mol->atoms; i < mol->natoms; i++, ap = ATOM_NEXT(ap)) {
			if (ig != NULL && IntGroupLookup(ig, i, NULL)) {
				ap->exflags |= kAtomHiddenFlag;
			} else {
				ap->exflags &= kAtomHiddenFlag;
			}
		}
		if (ig != NULL)
			IntGroupRelease(ig);
		MoleculeCallback_notifyModification(mol, 0);
	}
	return val; */
}

/*
 *  call-seq:
 *     select_frame(index)
 *     frame = index
 *
 *  Select the specified frame. If successful, returns true, otherwise returns false.
 */
static VALUE
s_Molecule_SelectFrame(VALUE self, VALUE val)
{
    Molecule *mol;
	int ival = NUM2INT(val);
    Data_Get_Struct(self, Molecule, mol);
	ival = MoleculeSelectFrame(mol, ival, 1);
	if (ival >= 0)
		return Qtrue;
	else return Qfalse;
}

/*
 *  call-seq:
 *     frame -> Integer
 *
 *  Get the current frame.
 */
static VALUE
s_Molecule_Frame(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->cframe);
}

/*
 *  call-seq:
 *     nframes -> Integer
 *
 *  Get the number of frames.
 */
static VALUE
s_Molecule_Nframes(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(MoleculeGetNumberOfFrames(mol));
}

/*
 *  call-seq:
 *     insert_frame(integer, coordinates = nil, cell_axes = nil) -> bool
 *     insert_frames(intGroup = nil, coordinates = nil, cell_axes = nil) -> bool
 *
 *  Insert new frames at the indices specified by the intGroup. If the first argument is
 *  an integer, a single new frame is inserted at that index. If the first argument is 
 *  nil, a new frame is inserted at the last. If non-nil coordinates is given, it
 *  should be an array of arrays of Vector3Ds, then those coordinates are set 
 *  to the new frame. Otherwise, the coordinates of current molecule are copied 
 *  to the new frame.
 *  Returns an intGroup representing the inserted frames if successful, nil if not.
 */
static VALUE
s_Molecule_InsertFrames(int argc, VALUE *argv, VALUE self)
{
	VALUE val, coords, cells;
    Molecule *mol;
	IntGroup *ig;
	int count, ival, i, j, len, len_c, len2, nframes;
	VALUE *ptr, *ptr2;
	Vector *vp, *vp2;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "12", &val, &coords, &cells);
	if (coords != Qnil) {
		if (TYPE(coords) != T_ARRAY)
			rb_raise(rb_eTypeError, "the coordinates should be given as an array of Vector3D");
		len = RARRAY_LEN(coords);
	} else len = 0;
	if (cells != Qnil) {
		if (mol->cell == NULL)
			rb_raise(rb_eTypeError, "the unit cell is not defined but the cell axes are given");
		if (TYPE(cells) != T_ARRAY)
			rb_raise(rb_eTypeError, "the cell axes should be given as an array of Vector3D");
		len_c = RARRAY_LEN(cells);
	} else len_c = 0;
	count = (len > len_c ? len : len_c);  /*  May be zero; will be updated later  */
	nframes = MoleculeGetNumberOfFrames(mol);
	if (val == Qnil) {
		ig = IntGroupNewWithPoints(nframes, (count > 0 ? count : 1), -1);
		val = ValueFromIntGroup(ig);
	} else {
		ig = IntGroupFromValue(val);
	}
	count = IntGroupGetCount(ig);  /*  Count is updated here  */
	vp = ALLOC_N(Vector, mol->natoms * count);
	if (cells != Qnil)
		vp2 = ALLOC_N(Vector, 4 * count);
	else vp2 = NULL;
	if (len > 0) {
		if (len < count)
			rb_raise(rb_eMolbyError, "the coordinates should contain no less than %d arrays (for frames)", count);
		ptr = RARRAY_PTR(coords);
		for (i = 0; i < count; i++) {
			if (TYPE(ptr[i]) != T_ARRAY)
				rb_raise(rb_eTypeError, "the coordinate array contains non-array object at index %d", i);
			len2 = RARRAY_LEN(ptr[i]);
			if (len2 < mol->natoms)
				rb_raise(rb_eMolbyError, "the array at index %d contains less than %d elements", i, mol->natoms);
			ptr2 = RARRAY_PTR(ptr[i]);
			for (j = 0; j < mol->natoms; j++)
				VectorFromValue(ptr2[j], &vp[i * mol->natoms + j]);
		}
	} else {
		Atom *ap;
		for (i = 0; i < count; i++) {
			for (j = 0, ap = mol->atoms; j < mol->natoms; j++, ap = ATOM_NEXT(ap)) {
				vp[i * mol->natoms + j] = ap->r;
			}
		}
	}
	if (len_c > 0) {
		if (len_c < count)
			rb_raise(rb_eMolbyError, "the cell vectors should contain no less than %d arrays (for frames)", count);
		ptr = RARRAY_PTR(cells);
		for (i = 0; i < count; i++) {
			if (TYPE(ptr[i]) != T_ARRAY)
				rb_raise(rb_eTypeError, "the cell parameter array contains non-array object at index %d", i);
			len2 = RARRAY_LEN(ptr[i]);
			if (len2 < 4)
				rb_raise(rb_eMolbyError, "the cell parameter should contain 4 vectors");
			ptr2 = RARRAY_PTR(ptr[i]);
			for (j = 0; j < 4; j++)
				VectorFromValue(ptr2[j], &vp2[i * 4 + j]);
		}
	}
	ival = MolActionCreateAndPerform(mol, gMolActionInsertFrames, ig, mol->natoms * count, vp, (vp2 != NULL ? 4 * count : 0), vp2);
	IntGroupRelease(ig);
	xfree(vp);
	if (vp2 != NULL)
		xfree(vp2);
	return (ival >= 0 ? val : Qnil);
}

/*
 *  call-seq:
 *     create_frame(coordinates = nil) -> Integer
 *     create_frames(coordinates = nil) -> Integer
 *
 *  Same as molecule.insert_frames(nil, coordinates).
 */
static VALUE
s_Molecule_CreateFrames(int argc, VALUE *argv, VALUE self)
{
	VALUE vals[3];
	rb_scan_args(argc, argv, "02", &vals[1], &vals[2]);
	vals[0] = Qnil;
	return s_Molecule_InsertFrames(3, vals, self);
}

/*
 *  call-seq:
 *     remove_frames(IntGroup, wantCoordinates = false)
 *
 *  Remove the frames at group. If wantsCoordinates is false (default), returns true if successful
 *  and nil otherwise. If wantsCoordinates is true, an array of arrays of the coordinates in the
 *  removed frames is returned if operation is successful.
 */
static VALUE
s_Molecule_RemoveFrames(int argc, VALUE *argv, VALUE self)
{
	VALUE val, flag;
	VALUE retval;
    Molecule *mol;
	IntGroup *ig;
	int count;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "11", &val, &flag);
	ig = IntGroupFromValue(val);
	count = IntGroupGetCount(ig);
	if (RTEST(flag)) {
		/*  Create return value before removing frames  */
		VALUE coords;
		int i, j, n;
		Atom *ap;
		Vector v;
		retval = rb_ary_new2(count);
		for (i = 0; i < count; i++) {
			n = IntGroupGetNthPoint(ig, i);
			coords = rb_ary_new2(mol->natoms);
			for (j = 0, ap = mol->atoms; j < mol->natoms; j++, ap = ATOM_NEXT(ap)) {
				if (n < ap->nframes && n != mol->cframe)
					v = ap->frames[n];
				else v = ap->r;
				rb_ary_push(coords, ValueFromVector(&v));
			}
			rb_ary_push(retval, coords);
		}
	} else retval = Qtrue;
	if (MolActionCreateAndPerform(mol, gMolActionRemoveFrames, ig) >= 0)
		return retval;
	else return Qnil;
}

/*
 *  call-seq:
 *     each_frame {|n| ...}
 *
 *  Set the frame number from 0 to nframes-1 and execute the block. The block argument is
 *  the frame number. After completion, the original frame number is restored.
 */
static VALUE
s_Molecule_EachFrame(VALUE self)
{
	int i, cframe, nframes;
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	cframe = mol->cframe;
	nframes = MoleculeGetNumberOfFrames(mol);
	if (nframes > 0) {
		for (i = 0; i < nframes; i++) {
			MoleculeSelectFrame(mol, i, 1);
			rb_yield(INT2NUM(i));
		}
		MoleculeSelectFrame(mol, cframe, 1);
	}
    return self;
}

/*
 *  call-seq:
 *     set_atom_attr(index, key, value)
 *
 *  Set the atom attribute for the specified atom.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_SetAtomAttr(VALUE self, VALUE idx, VALUE key, VALUE val)
{
	Molecule *mol;
	VALUE aref, oldval;
    Data_Get_Struct(self, Molecule, mol);
	aref = ValueFromMoleculeAndIndex(mol, s_Molecule_AtomIndexFromValue(mol, idx));
	oldval = s_AtomRef_GetAttr(aref, key);
	if (val == Qundef)
		return oldval;
	s_AtomRef_SetAttr(aref, key, val);
	return val;
}

/*
 *  call-seq:
 *     get_atom_attr(index, key)
 *
 *  Get the atom attribute for the specified atom.
 */
static VALUE
s_Molecule_GetAtomAttr(VALUE self, VALUE idx, VALUE key)
{
	return s_Molecule_SetAtomAttr(self, idx, key, Qundef);
}

/*
 *  call-seq:
 *     get_coord_from_frame(index, group = nil)
 *
 *  Copy the coordinates from the indicated frame. If group is specified, only the specified atoms
 *  are modified. Third argument (cflag) is now obsolete (it used to specify whether the cell parameters are to be
 *  copied; now they are always copied)
 */
static VALUE
s_Molecule_GetCoordFromFrame(int argc, VALUE *argv, VALUE self)
{
	Molecule *mol;
	VALUE ival, gval, cval;
	Int index, i, j, n, nn;
	IntGroup *ig;
	IntGroupIterator iter;
	Atom *ap;
	Vector *vp;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "12", &ival, &gval, &cval);
	if (argc == 3)
		rb_warn("The 3rd argument to get_coord_from_frame() is now obsolete");
	index = NUM2INT(rb_Integer(ival));
	if (index < 0 || index >= (n = MoleculeGetNumberOfFrames(mol))) {
		if (n == 0)
			rb_raise(rb_eMolbyError, "No frame is present");
		else
			rb_raise(rb_eMolbyError, "Frame index (%d) out of range (should be 0..%d)", index, n - 1);
	}
	if (gval == Qnil) {
		ig = IntGroupNewWithPoints(0, mol->natoms, -1);
	} else {
		ig = s_Molecule_AtomGroupFromValue(self, gval);
	}
	n = IntGroupGetCount(ig);
	if (n > 0) {
		vp = (Vector *)calloc(sizeof(Vector), n);
		IntGroupIteratorInit(ig, &iter);
		j = 0;
		nn = 0;
		while ((i = IntGroupIteratorNext(&iter)) >= 0) {
			ap = ATOM_AT_INDEX(mol->atoms, i);
			if (index < ap->nframes) {
				vp[j] = ap->frames[index];
				nn++;
			} else {
				vp[j] = ap->r;
			}
			j++;
		}
		if (nn > 0)
			MolActionCreateAndPerform(mol, gMolActionSetAtomPositions, ig, n, vp);
		free(vp);
		if (mol->cell != NULL && mol->frame_cells != NULL && index < mol->nframe_cells) {
			vp = mol->frame_cells + index * 4;
			MolActionCreateAndPerform(mol, gMolActionSetBox, vp, vp + 1, vp + 2, vp + 3, -1, 0);
		}
		IntGroupIteratorRelease(&iter);
	}
	IntGroupRelease(ig);
	return self;
}

/*
 *  call-seq:
 *     fragment(n1, *exatoms)  -> IntGroup
 *     fragment(group, *exatoms)  -> IntGroup
 *
 *  Get the fragment including the atom n1 or the atom group. If additional arguments are given,
 *  those atoms will not be counted during the search.
 */
static VALUE
s_Molecule_Fragment(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	IntGroup *baseg, *ig, *exatoms;
	int n;
	volatile VALUE nval, exval;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "1*", &nval, &exval);
	if (rb_obj_is_kind_of(nval, rb_cNumeric) || rb_obj_is_kind_of(nval, rb_cString)) {
		baseg = NULL;
		n = NUM2INT(s_Molecule_AtomIndex(self, nval));
	} else {
		baseg = s_Molecule_AtomGroupFromValue(self, nval);
	}
	if (RARRAY_LEN(exval) == 0) {
		exatoms = NULL;
	} else {
		exval = s_Molecule_AtomGroup(RARRAY_LEN(exval), RARRAY_PTR(exval), self);
		Data_Get_Struct(exval, IntGroup, exatoms);
	}
	if (baseg == NULL) {
		ig = MoleculeFragmentExcludingAtomGroup(mol, n, exatoms);
	} else {
		IntGroupIterator iter;
		IntGroupIteratorInit(baseg, &iter);
		if ((n = IntGroupIteratorNext(&iter)) < 0) {
			ig = IntGroupNew();
		} else {
			ig = MoleculeFragmentExcludingAtomGroup(mol, n, exatoms);
			if (ig != NULL) {
				while ((n = IntGroupIteratorNext(&iter)) >= 0) {
					IntGroup *subg;
					subg = MoleculeFragmentExcludingAtomGroup(mol, n, exatoms);
					if (subg != NULL) {
						IntGroupAddIntGroup(ig, subg);
						IntGroupRelease(subg);
					}
				}
			}
		}
		IntGroupIteratorRelease(&iter);
		IntGroupRelease(baseg);
	}
	if (ig == NULL)
		rb_raise(rb_eMolbyError, "invalid specification of molecular fragment");
	nval = ValueFromIntGroup(ig);
	IntGroupRelease(ig);
	return nval;
}

/*
 *  call-seq:
 *     fragments(exclude = nil)
 *
 *  Returns the fragments as an array of IntGroups. 
 *  If exclude is given (as an array or an IntGroup), then those atoms are excluded
 *  in defining the fragment.
 */
static VALUE
s_Molecule_Fragments(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	IntGroup *ag, *fg, *eg;
	VALUE gval, exval, retval;
    Data_Get_Struct(self, Molecule, mol);
	if (mol == NULL)
		return Qnil;
	if (mol->natoms == 0)
		return rb_ary_new();
	rb_scan_args(argc, argv, "01", &exval);
	if (exval == Qnil)
		eg = NULL;
	else
		eg = IntGroupFromValue(exval);
	ag = IntGroupNewWithPoints(0, mol->natoms, -1);
	if (eg != NULL)
		IntGroupRemoveIntGroup(ag, eg);
	retval = rb_ary_new();
	while (IntGroupGetCount(ag) > 0) {
		int n = IntGroupGetNthPoint(ag, 0);
		fg = MoleculeFragmentExcludingAtomGroup(mol, n, eg);
		if (fg == NULL)
			rb_raise(rb_eMolbyError, "internal error during each_fragment");
		gval = ValueFromIntGroup(fg);
		rb_ary_push(retval, gval);
		IntGroupRemoveIntGroup(ag, fg);
		IntGroupRelease(fg);
	}
	IntGroupRelease(ag);
	if (eg != NULL)
		IntGroupRelease(eg);
	return retval;
}

/*
 *  call-seq:
 *     each_fragment(exclude = nil) {|group| ...}
 *
 *  Execute the block, with the IntGroup object for each fragment as the argument.
 *  Atoms or bonds should not be added or removed during the execution of the block.
 *  If exclude is given (as an array or an IntGroup), then those atoms are excluded
 *  in defining the fragment.
 */
static VALUE
s_Molecule_EachFragment(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	IntGroup *ag, *fg, *eg;
	VALUE gval, exval;
    Data_Get_Struct(self, Molecule, mol);
	if (mol == NULL || mol->natoms == 0)
		return self;
	rb_scan_args(argc, argv, "01", &exval);
	if (exval == Qnil)
		eg = NULL;
	else
		eg = IntGroupFromValue(exval);
	ag = IntGroupNewWithPoints(0, mol->natoms, -1);
	if (eg != NULL)
		IntGroupRemoveIntGroup(ag, eg);
	while (IntGroupGetCount(ag) > 0) {
		int n = IntGroupGetNthPoint(ag, 0);
		fg = MoleculeFragmentExcludingAtomGroup(mol, n, eg);
		if (fg == NULL)
			rb_raise(rb_eMolbyError, "internal error during each_fragment");
		gval = ValueFromIntGroup(fg);
		rb_yield(gval);
		IntGroupRemoveIntGroup(ag, fg);
		IntGroupRelease(fg);
	}
	IntGroupRelease(ag);
	if (eg != NULL)
		IntGroupRelease(eg);
	return self;
}

/*
 *  call-seq:
 *     detachable?(group)  -> [n1, n2]
 *
 *  Check whether the group is 'detachable', i.e. the group is bound to the rest 
 *  of the molecule via only one bond. If it is, then the indices of the atoms
 *  belonging to the bond is returned, the first element being the atom included
 *  in the fragment. Otherwise, Qnil is returned.
 */
static VALUE
s_Molecule_Detachable_P(VALUE self, VALUE gval)
{
	Molecule *mol;
	IntGroup *ig;
	int n1, n2;
	VALUE retval;
    Data_Get_Struct(self, Molecule, mol);
	ig = s_Molecule_AtomGroupFromValue(self, gval);
	if (MoleculeIsFragmentDetachable(mol, ig, &n1, &n2)) {
		retval = rb_ary_new3(2, INT2NUM(n1), INT2NUM(n2));
	} else retval = Qnil;
	IntGroupRelease(ig);
	return retval;
}

/*
 *  call-seq:
 *     bonds_on_border(group = selection)  -> Array of Array of two Integers
 *
 *  Returns an array of bonds that connect an atom in the group and an atom out
 *  of the group. The first atom in the bond always belongs to the group. If no
 *  such bonds are present, an empty array is returned.
 */
static VALUE
s_Molecule_BondsOnBorder(int argc, VALUE *argv, VALUE self)
{
	Molecule *mol;
	IntGroup *ig, *bg;
	VALUE gval, retval;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "01", &gval);
	if (gval == Qnil) {
		ig = MoleculeGetSelection(mol);
		if (ig != NULL)
			IntGroupRetain(ig);
	} else {
		ig = s_Molecule_AtomGroupFromValue(self, gval);
	}
	retval = rb_ary_new();
	if (ig == NULL)
		return retval;
	bg = MoleculeSearchBondsAcrossAtomGroup(mol, ig);
	if (bg != NULL) {
		IntGroupIterator iter;
		Int i;
		IntGroupIteratorInit(bg, &iter);
		while ((i = IntGroupIteratorNext(&iter)) >= 0) {
			/*  The atoms at the border  */
			Int n1, n2;
			n1 = mol->bonds[i * 2];
			n2 = mol->bonds[i * 2 + 1];
			if (IntGroupLookupPoint(ig, n1) < 0) {
				int w = n1;
				n1 = n2;
				n2 = w;
				if (IntGroupLookupPoint(ig, n1) < 0)
					continue;  /*  Actually this is an internal error  */
			}
			rb_ary_push(retval, rb_ary_new3(2, INT2NUM(n1), INT2NUM(n2)));
		}
		IntGroupIteratorRelease(&iter);
	}
	IntGroupRelease(bg);
	IntGroupRelease(ig);
	return retval;
}

/*
 *  call-seq:
 *     translate(vec, group = nil)       -> Molecule
 *
 *  Translate the molecule by vec. If group is given, only atoms in the group are moved.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_Translate(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE vec, group;
	Vector v;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "11", &vec, &group);
	ig = (NIL_P(group) ? NULL : s_Molecule_AtomGroupFromValue(self, group));
	VectorFromValue(vec, &v);
//	MoleculeTranslate(mol, &v, ig);
	MolActionCreateAndPerform(mol, gMolActionTranslateAtoms, &v, ig);
	if (ig != NULL)
		IntGroupRelease(ig);
	return self;
}

/*
 *  call-seq:
 *     rotate(axis, angle, center = [0,0,0], group = nil)       -> Molecule
 *
 *  Rotate the molecule. The axis must not a zero vector. angle is given in degree.
 *  If group is given, only atoms in the group are moved.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_Rotate(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	volatile VALUE aval, anval, cval, gval;
	Double angle;
	Vector av, cv;
	Transform tr;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "22", &aval, &anval, &cval, &gval);
	ig = (NIL_P(gval) ? NULL : s_Molecule_AtomGroupFromValue(self, gval));
	angle = NUM2DBL(rb_Float(anval)) * kDeg2Rad;
	VectorFromValue(aval, &av);
	if (NIL_P(cval))
		cv.x = cv.y = cv.z = 0.0;
	else
		VectorFromValue(cval, &cv);
	if (TransformForRotation(tr, &av, angle, &cv))
		rb_raise(rb_eMolbyError, "rotation axis cannot be a zero vector");
	MolActionCreateAndPerform(mol, gMolActionTransformAtoms, &tr, ig);
	if (ig != NULL)
		IntGroupRelease(ig);
	return self;
}

/*
 *  call-seq:
 *     reflect(axis, center = [0,0,0], group = nil)       -> Molecule
 *
 *  Reflect the molecule by the plane which is perpendicular to axis and including center. 
 *  axis must not be a zero vector.
 *  If group is given, only atoms in the group are moved.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_Reflect(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	volatile VALUE aval, cval, gval;
	Vector av, cv;
	Transform tr;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "12", &aval, &cval, &gval);
	ig = (NIL_P(gval) ? NULL : s_Molecule_AtomGroupFromValue(self, gval));
	VectorFromValue(aval, &av);
	if (NIL_P(cval))
		cv.x = cv.y = cv.z = 0.0;
	else
		VectorFromValue(cval, &cv);
	if (TransformForReflection(tr, &av, &cv))
		rb_raise(rb_eMolbyError, "reflection axis cannot be a zero vector");
	MolActionCreateAndPerform(mol, gMolActionTransformAtoms, &tr, ig);
	if (ig != NULL)
		IntGroupRelease(ig);
	return self;
}

/*
 *  call-seq:
 *     invert(center = [0,0,0], group = nil)       -> Molecule
 *
 *  Invert the molecule with the given center.
 *  If group is given, only atoms in the group are moved.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_Invert(int argc, VALUE *argv, VALUE self)
{
	Molecule *mol;
	volatile VALUE cval, gval;
	Vector cv;
	Transform tr;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "02", &cval, &gval);
	ig = (NIL_P(gval) ? NULL : s_Molecule_AtomGroupFromValue(self, gval));
	if (NIL_P(cval))
		cv.x = cv.y = cv.z = 0.0;
	else
		VectorFromValue(cval, &cv);
	TransformForInversion(tr, &cv);
	MolActionCreateAndPerform(mol, gMolActionTransformAtoms, &tr, ig);
	if (ig != NULL)
		IntGroupRelease(ig);
	return self;
}

/*
 *  call-seq:
 *     transform(transform, group = nil)       -> Molecule
 *
 *  Transform the molecule by the given Transform object.
 *  If group is given, only atoms in the group are moved.
 *  This operation is undoable.
 */
static VALUE
s_Molecule_Transform(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE trans, group;
	Transform tr;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "11", &trans, &group);
	ig = (NIL_P(group) ? NULL : s_Molecule_AtomGroupFromValue(self, group));
	TransformFromValue(trans, &tr);
/*	MoleculeTransform(mol, tr, ig); */
	MolActionCreateAndPerform(mol, gMolActionTransformAtoms, &tr, ig);
	if (ig != NULL)
		IntGroupRelease(ig);
	return self;
}

static void
s_Molecule_DoCenterOfMass(Molecule *mol, Vector *outv, IntGroup *ig)
{
	switch (MoleculeCenterOfMass(mol, outv, ig)) {
		case 2: rb_raise(rb_eMolbyError, "atom group is empty"); break;
		case 3: rb_raise(rb_eMolbyError, "weight is zero --- atomic weights are not defined?"); break;
		case 0: break;
		default: rb_raise(rb_eMolbyError, "cannot calculate center of mass"); break;
	}
}

/*
 *  call-seq:
 *     center_of_mass(group = nil)       -> Vector3D
 *
 *  Calculate the center of mass for the given set of atoms. The argument
 *  group is null, then all atoms are considered.
 */
static VALUE
s_Molecule_CenterOfMass(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE group;
	IntGroup *ig;
	Vector v;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "01", &group);
	ig = (NIL_P(group) ? NULL : s_Molecule_AtomGroupFromValue(self, group));
	s_Molecule_DoCenterOfMass(mol, &v, ig);
	if (ig != NULL)
		IntGroupRelease(ig);
	return ValueFromVector(&v);
}

/*
 *  call-seq:
 *     centralize(group = nil)       -> self
 *
 *  Translate the molecule so that the center of mass of the given group is located
 *  at (0, 0, 0). Equivalent to molecule.translate(molecule.center_of_mass(group) * -1).
 */
static VALUE
s_Molecule_Centralize(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE group;
	IntGroup *ig;
	Vector v;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "01", &group);
	ig = (NIL_P(group) ? NULL : s_Molecule_AtomGroupFromValue(self, group));
	s_Molecule_DoCenterOfMass(mol, &v, ig);
	if (ig != NULL)
		IntGroupRelease(ig);
	v.x = -v.x;
	v.y = -v.y;
	v.z = -v.z;
	MolActionCreateAndPerform(mol, gMolActionTranslateAtoms, &v, NULL);
	return self;
}

/*
 *  call-seq:
 *     bounds(group = nil)       -> [min, max]
 *
 *  Calculate the boundary. The return value is an array of two Vector3D objects.
 */
static VALUE
s_Molecule_Bounds(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE group;
	IntGroup *ig;
	Vector vmin, vmax;
	int n;
	Atom *ap;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "01", &group);
	ig = (NIL_P(group) ? NULL : s_Molecule_AtomGroupFromValue(self, group));
	if (ig != NULL && IntGroupGetCount(ig) == 0)
		rb_raise(rb_eMolbyError, "atom group is empty");
	vmin.x = vmin.y = vmin.z = 1e30;
	vmax.x = vmax.y = vmax.z = -1e30;
	for (n = 0, ap = mol->atoms; n < mol->natoms; n++, ap = ATOM_NEXT(ap)) {
		Vector r;
		if (ig != NULL && IntGroupLookup(ig, n, NULL) == 0)
			continue;
		r = ap->r;
		if (r.x < vmin.x)
			vmin.x = r.x;
		if (r.y < vmin.y)
			vmin.y = r.y;
		if (r.z < vmin.z)
			vmin.z = r.z;
		if (r.x > vmax.x)
			vmax.x = r.x;
		if (r.y > vmax.y)
			vmax.y = r.y;
		if (r.z > vmax.z)
			vmax.z = r.z;
	}
	return rb_ary_new3(2, ValueFromVector(&vmin), ValueFromVector(&vmax));
}

/*  Get atom position or a vector  */
static void
s_Molecule_GetVectorFromArg(Molecule *mol, VALUE val, Vector *vp)
{
	if (rb_obj_is_kind_of(val, rb_cInteger) || rb_obj_is_kind_of(val, rb_cString)) {
		int n1 = s_Molecule_AtomIndexFromValue(mol, val);
		*vp = ATOM_AT_INDEX(mol->atoms, n1)->r;
	} else {
		VectorFromValue(val, vp);
	}
}

/*
 *  call-seq:
 *     measure_bond(n1, n2)       -> Float
 *
 *  Calculate the bond length. The arguments can either be atom indices, the "residue:name" representation, 
 *  or Vector3D values.
 *  If the crystallographic cell is defined, then the internal coordinates are convereted to the cartesian.
 */
static VALUE
s_Molecule_MeasureBond(VALUE self, VALUE nval1, VALUE nval2)
{
    Molecule *mol;
	Vector v1, v2;
    Data_Get_Struct(self, Molecule, mol);
	s_Molecule_GetVectorFromArg(mol, nval1, &v1);
	s_Molecule_GetVectorFromArg(mol, nval2, &v2);
	return rb_float_new(MoleculeMeasureBond(mol, &v1, &v2));
}

/*
 *  call-seq:
 *     measure_angle(n1, n2, n3)       -> Float
 *
 *  Calculate the bond angle. The arguments can either be atom indices, the "residue:name" representation, 
 *  or Vector3D values. The return value is in degree.
 *  If the crystallographic cell is defined, then the internal coordinates are convereted to the cartesian.
 */
static VALUE
s_Molecule_MeasureAngle(VALUE self, VALUE nval1, VALUE nval2, VALUE nval3)
{
    Molecule *mol;
	Vector v1, v2, v3;
	Double d;
    Data_Get_Struct(self, Molecule, mol);
	s_Molecule_GetVectorFromArg(mol, nval1, &v1);
	s_Molecule_GetVectorFromArg(mol, nval2, &v2);
	s_Molecule_GetVectorFromArg(mol, nval3, &v3);	
	d = MoleculeMeasureAngle(mol, &v1, &v2, &v3);
	if (isnan(d))
		return Qnil;  /*  Cannot define  */
	else return rb_float_new(d);
}

/*
 *  call-seq:
 *     measure_dihedral(n1, n2, n3, n4)       -> Float
 *
 *  Calculate the dihedral angle. The arguments can either be atom indices, the "residue:name" representation, 
 *  or Vector3D values. The return value is in degree.
 *  If the crystallographic cell is defined, then the internal coordinates are convereted to the cartesian.
 */
static VALUE
s_Molecule_MeasureDihedral(VALUE self, VALUE nval1, VALUE nval2, VALUE nval3, VALUE nval4)
{
    Molecule *mol;
	Vector v1, v2, v3, v4;
	Double d;
    Data_Get_Struct(self, Molecule, mol);
	s_Molecule_GetVectorFromArg(mol, nval1, &v1);
	s_Molecule_GetVectorFromArg(mol, nval2, &v2);
	s_Molecule_GetVectorFromArg(mol, nval3, &v3);	
	s_Molecule_GetVectorFromArg(mol, nval4, &v4);	
	d = MoleculeMeasureDihedral(mol, &v1, &v2, &v3, &v4);
	if (isnan(d))
		return Qnil;  /*  Cannot define  */
	else return rb_float_new(d);
}

/*
 *  call-seq:
 *     expand_by_symmetry(group, sym, dx=0, dy=0, dz=0) -> Array
 *
 *  Expand the specified part of the molecule by the given symmetry operation.
 *  Returns the array of atom indices corresponding to the expanded atoms.
 */
static VALUE
s_Molecule_ExpandBySymmetry(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE gval, sval, xval, yval, zval, rval;
	IntGroup *ig;
	Int n[4];
	Int natoms;
	Int nidx, *idx;

    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "23", &gval, &sval, &xval, &yval, &zval);
	n[0] = NUM2INT(rb_Integer(sval));
	n[1] = (xval == Qnil ? 0 : NUM2INT(rb_Integer(xval)));
	n[2] = (yval == Qnil ? 0 : NUM2INT(rb_Integer(yval)));
	n[3] = (zval == Qnil ? 0 : NUM2INT(rb_Integer(zval)));
	ig = s_Molecule_AtomGroupFromValue(self, gval);
	if (n[0] < 0 || (n[0] > 0 && n[0] >= mol->nsyms))
		rb_raise(rb_eMolbyError, "symmetry index is out of bounds");
	natoms = mol->natoms;
	
	MolActionCreateAndPerform(mol, gMolActionExpandBySymmetry, ig, n[1], n[2], n[3], n[0], &nidx, &idx);

	rval = rb_ary_new2(nidx);
	while (--nidx >= 0) {
		rb_ary_store(rval, nidx, INT2NUM(idx[nidx]));
	}
/*	if (natoms == mol->natoms)
		rval = Qnil;
	else {
		rval = IntGroup_Alloc(rb_cIntGroup);
		Data_Get_Struct(rval, IntGroup, ig);
		IntGroup_RaiseIfError(IntGroupAdd(ig, natoms, mol->natoms - natoms));
	} */
	return rval;
}

/*
 *  call-seq:
 *     amend_by_symmetry(group = nil) -> IntGroup
 *
 *  Expand the specified part of the molecule by the given symmetry operation.
 *  Returns an IntGroup containing the added atoms.
 */
static VALUE
s_Molecule_AmendBySymmetry(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	IntGroup *ig, *ig2;
	VALUE rval, gval;
    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "01", &gval);
	if (gval != Qnil)
		ig = s_Molecule_AtomGroupFromValue(self, gval);
	else ig = NULL;
	MolActionCreateAndPerform(mol, gMolActionAmendBySymmetry, ig, &ig2);
	rval = ValueFromIntGroup(ig2);
	IntGroupRelease(ig2);
	return rval;
}

/*
 *  call-seq:
 *     transform_for_symop(symop, is_cartesian = nil) -> Transform
 *
 *  Get the transform corresponding to the symmetry operation. The symop can either be
 *  an integer (index of symmetry operation) or [sym, dx, dy, dz].
 *  If is_cartesian is true, the returned transform is for cartesian coordinates.
 *  Otherwise, the returned transform is for fractional coordinates.
 *  Raises exception when no cell or no transform are defined.
 */
static VALUE
s_Molecule_TransformForSymop(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE sval, fval;
	Symop symop;
	Transform tr;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->cell == NULL)
		rb_raise(rb_eMolbyError, "no unit cell is defined");
	if (mol->nsyms == 0)
		rb_raise(rb_eMolbyError, "no symmetry operation is defined");
	rb_scan_args(argc, argv, "11", &sval, &fval);
	if (rb_obj_is_kind_of(sval, rb_cNumeric)) {
		symop.sym = NUM2INT(rb_Integer(sval));
		symop.dx = symop.dy = symop.dz = 0;
	} else {
		sval = rb_ary_to_ary(sval);
		if (RARRAY_LEN(sval) < 4)
			rb_raise(rb_eMolbyError, "missing arguments as symop; at least four integers should be given");
		symop.sym = NUM2INT(rb_Integer(RARRAY_PTR(sval)[0]));
		symop.dx = NUM2INT(rb_Integer(RARRAY_PTR(sval)[1]));
		symop.dy = NUM2INT(rb_Integer(RARRAY_PTR(sval)[2]));
		symop.dz = NUM2INT(rb_Integer(RARRAY_PTR(sval)[3]));
	}
	if (symop.sym >= mol->nsyms)
		rb_raise(rb_eMolbyError, "index of symmetry operation (%d) is out of range", symop.sym);
	MoleculeGetTransformForSymop(mol, symop, &tr, (RTEST(fval) != 0));
	return ValueFromTransform(&tr);
}
	
/*
 *  call-seq:
 *     symop_for_transform(transform, is_cartesian = nil) -> [sym, dx, dy, dz]
 *
 *  Get the symmetry operation corresponding to the given transform.
 *  If is_cartesian is true, the given transform is for cartesian coordinates.
 *  Otherwise, the given transform is for fractional coordinates.
 *  Raises exception when no cell or no transform are defined.
 */
static VALUE
s_Molecule_SymopForTransform(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE tval, fval;
	Symop symop;
	Transform tr;
	int n;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->cell == NULL)
		rb_raise(rb_eMolbyError, "no unit cell is defined");
	if (mol->nsyms == 0)
		rb_raise(rb_eMolbyError, "no symmetry operation is defined");
	rb_scan_args(argc, argv, "11", &tval, &fval);
	TransformFromValue(tval, &tr);
	n = MoleculeGetSymopForTransform(mol, tr, &symop, (RTEST(fval) != 0));
	if (n == 0) {
		return rb_ary_new3(4, INT2NUM(symop.sym), INT2NUM(symop.dx), INT2NUM(symop.dy), INT2NUM(symop.dz));
	} else {
		return Qnil;  /*  Not found  */
	}
}

/*
 *  call-seq:
 *     wrap_unit_cell(group) -> Vector3D
 *
 *  Move the specified group so that the center of mass of the group is within the
 *  unit cell. The offset vector is returned. If no periodic box is defined, 
 *  exception is raised.
 */
static VALUE
s_Molecule_WrapUnitCell(VALUE self, VALUE gval)
{
    Molecule *mol;
	IntGroup *ig;
	Vector v, cv, dv;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->cell == NULL)
		rb_raise(rb_eMolbyError, "no unit cell is defined");
	ig = s_Molecule_AtomGroupFromValue(self, gval);
	s_Molecule_DoCenterOfMass(mol, &cv, ig);
	TransformVec(&v, mol->cell->rtr, &cv);
	if (mol->cell->flags[0])
		v.x -= floor(v.x);
	if (mol->cell->flags[1])
		v.y -= floor(v.y);
	if (mol->cell->flags[2])
		v.z -= floor(v.z);
	TransformVec(&dv, mol->cell->tr, &v);
	VecDec(dv, cv);
	MolActionCreateAndPerform(mol, gMolActionTranslateAtoms, &dv, ig);
	IntGroupRelease(ig);
	return ValueFromVector(&dv);
}

/*
 *  call-seq:
 *     find_conflicts(limit[, group1[, group2 [, ignore_exclusion]]]) -> [[n1, n2], [n3, n4], ...]
 *
 *  Find pairs of atoms that are within the limit distance. If group1 and group2 are given, the
 *  first and second atom in the pair should belong to group1 and group2, respectively.
 *  If ignore_exclusion is true, then 1-2 (bonded), 1-3, 1-4 pairs are also considered.
 */
static VALUE
s_Molecule_FindConflicts(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE limval, gval1, gval2, rval, igval;
	IntGroup *ig1, *ig2;
	IntGroupIterator iter1, iter2;
	Int npairs, *pairs;
	Int n[2], i;
	Double lim;
	Vector r1;
	Atom *ap1, *ap2;
	MDExclusion *exinfo;
	Int *exlist;

    Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "13", &limval, &gval1, &gval2, &igval);
	lim = NUM2DBL(rb_Float(limval));
	if (lim <= 0.0)
		rb_raise(rb_eMolbyError, "the limit (%g) should be positive", lim);
	if (gval1 != Qnil)
		ig1 = s_Molecule_AtomGroupFromValue(self, gval1);
	else
		ig1 = IntGroupNewWithPoints(0, mol->natoms, -1);
	if (gval2 != Qnil)
		ig2 = s_Molecule_AtomGroupFromValue(self, gval2);
	else
		ig2 = IntGroupNewWithPoints(0, mol->natoms, -1);
	
	if (!RTEST(igval)) {
		/*  Use the exclusion table in MDArena  */
		if (mol->par == NULL || mol->arena == NULL || mol->arena->is_initialized == 0 || mol->needsMDRebuild) {
			VALUE mval = ValueFromMolecule(mol);
			s_RebuildMDParameterIfNecessary(mval, Qnil);
		}
		exinfo = mol->arena->exinfo;  /*  May be NULL  */
		exlist = mol->arena->exlist;	
	} else {
		exinfo = NULL;
		exlist = NULL;
	}
	IntGroupIteratorInit(ig1, &iter1);
	IntGroupIteratorInit(ig2, &iter2);
	npairs = 0;
	pairs = NULL;
	while ((n[0] = IntGroupIteratorNext(&iter1)) >= 0) {
		Int exn1, exn2;
		ap1 = ATOM_AT_INDEX(mol->atoms, n[0]);
		r1 = ap1->r;
		if (exinfo != NULL) {
			exn1 = exinfo[n[0]].index1;
			exn2 = exinfo[n[0] + 1].index1;
		} else exn1 = exn2 = -1;
		IntGroupIteratorReset(&iter2);
		while ((n[1] = IntGroupIteratorNext(&iter2)) >= 0) {
			ap2 = ATOM_AT_INDEX(mol->atoms, n[1]);
			if (n[0] == n[1])
				continue;  /*  Same atom  */
			if (exinfo != NULL) {
				/*  Look up exclusion table to exclude 1-2, 1-3, and 1-4 pairs  */
				for (i = exn1; i < exn2; i++) {
					if (exlist[i] == n[1])
						break;
				}
				if (i < exn2)
					continue;  /*  Should be excluded  */
			}
			if (MoleculeMeasureBond(mol, &r1, &(ap2->r)) < lim) {
				/*  Is this pair already registered?  */
				Int *ip;
				for (i = 0, ip = pairs; i < npairs; i++, ip += 2) {
					if ((ip[0] == n[0] && ip[1] == n[1]) || (ip[0] == n[1] && ip[1] == n[0]))
						break;
				}
				if (i >= npairs) {
					/*  Not registered yet  */
					AssignArray(&pairs, &npairs, sizeof(Int) * 2, npairs, n);
				}
			}
		}
	}
	IntGroupIteratorRelease(&iter2);
	IntGroupIteratorRelease(&iter1);
	IntGroupRelease(ig2);
	IntGroupRelease(ig1);
	rval = rb_ary_new2(npairs);
	if (pairs != NULL) {
		for (i = 0; i < npairs; i++) {
			rb_ary_push(rval, rb_ary_new3(2, INT2NUM(pairs[i * 2]), INT2NUM(pairs[i * 2 + 1])));
		}
		free(pairs);
	}
	return rval;
}

/*  Calculate the transform that moves the current coordinates to the reference
 coordinates with least displacements.   */
static Double
s_Molecule_FitCoordinates_Sub(Molecule *mol, IntGroup *ig, Vector *ref, Double *weights, Transform trans)
{
	Atom *ap, *ap1;
	Int natoms, nn;
	Vector org1, org2;
	Int i, in, j, k;
	Double w, w1;
	Mat33 r, q, u;
	Double eigen_val[3];
	Vector eigen_vec[3];
	Vector s[3];
	IntGroupIterator iter;

	natoms = mol->natoms;
	ap = mol->atoms;
	IntGroupIteratorInit(ig, &iter);
	
	/*  Calculate the weighted center  */
	VecZero(org1);
	VecZero(org2);
	w = 0.0;
	for (i = 0; (in = IntGroupIteratorNext(&iter)) >= 0; i++) {
		ap1 = ATOM_AT_INDEX(ap, in);
		w1 = (weights != NULL ? weights[i] : ap1->weight);
		VecScaleInc(org1, ap1->r, w1);
		VecScaleInc(org2, ref[i], w1);
		w += w1;
	}
	w = 1.0 / w;
	VecScaleSelf(org1, w);
	VecScaleSelf(org2, w);

    /*  R = sum(weight[n]^2 * x[n] * t(y[n]));  */
    /*  Matrix to diagonalize = R * tR    */
	memset(r, 0, sizeof(Mat33));
	memset(q, 0, sizeof(Mat33));
	memset(u, 0, sizeof(Mat33));
	nn = 0;
	IntGroupIteratorReset(&iter);
	for (i = 0; (in = IntGroupIteratorNext(&iter)) >= 0; i++) {
		Vector v1, v2;
		ap1 = ATOM_AT_INDEX(ap, in);
		w1 = (weights != NULL ? weights[i] : ap1->weight);
		w1 *= w1;
		VecSub(v1, ap1->r, org1);
		VecSub(v2, ref[i], org2);
		r[0] += w1 * v1.x * v2.x;
		r[1] += w1 * v1.y * v2.x;
		r[2] += w1 * v1.z * v2.x;
		r[3] += w1 * v1.x * v2.y;
		r[4] += w1 * v1.y * v2.y;
		r[5] += w1 * v1.z * v2.y;
		r[6] += w1 * v1.x * v2.z;
		r[7] += w1 * v1.y * v2.z;
		r[8] += w1 * v1.z * v2.z;
		nn++;
	}
	for (i = 0; i < 9; i++)
		r[i] /= (nn * nn);
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				q[i+j*3] += r[i+k*3] * r[j+k*3];
			}
		}
	}
	
	if (MatrixSymDiagonalize(q, eigen_val, eigen_vec) != 0) {
		IntGroupIteratorRelease(&iter);
		return -1.0;  /*  Cannot determine the eigenvector  */
	}

    /*  s[i] = tR * v[i] / sqrt(eigenval[i])  */
    /*  U = s0*t(v0) + s1*t(v1) + s2*t(v2)  */
	MatrixTranspose(r, r);
	for (i = 0; i < 3; i++) {
		MatrixVec(&s[i], r, &eigen_vec[i]);
		w1 = 1.0 / sqrt(eigen_val[i]);
		VecScaleSelf(s[i], w1);
	}
	for (k = 0; k < 3; k++) {
		u[0] += s[k].x * eigen_vec[k].x;
		u[1] += s[k].y * eigen_vec[k].x;
		u[2] += s[k].z * eigen_vec[k].x;
		u[3] += s[k].x * eigen_vec[k].y;
		u[4] += s[k].y * eigen_vec[k].y;
		u[5] += s[k].z * eigen_vec[k].y;
		u[6] += s[k].x * eigen_vec[k].z;
		u[7] += s[k].y * eigen_vec[k].z;
		u[8] += s[k].z * eigen_vec[k].z;
	}
	
	/*  y = U*(x - org1) + org2 = U*x + (org2 - U*org1)  */
	MatrixVec(&org1, u, &org1);
	VecDec(org2, org1);
	for (i = 0; i < 9; i++)
		trans[i] = u[i];
	trans[9] = org2.x;
	trans[10] = org2.y;
	trans[11] = org2.z;
	
	/*  Calculate rmsd  */
	IntGroupIteratorReset(&iter);
	w = 0.0;
	for (i = 0; (in = IntGroupIteratorNext(&iter)) >= 0; i++) {
		Vector tv;
		ap1 = ATOM_AT_INDEX(ap, in);
		TransformVec(&tv, trans, &ap1->r);
		VecDec(tv, ref[i]);
		w += VecLength2(tv);
	}
	w = sqrt(w / nn);
	IntGroupIteratorRelease(&iter);
	return w;
}

/*
 *  call-seq:
 *     fit_coordinates(group, ref, weight = nil) -> [transform, rmsd]
 *
 *  Calculate the transform to fit the given group to the set of reference coordinates.
 *  The reference coordinates ref is given as either a frame number, an array of
 *  Vector3Ds or arrays, or an LAMatrix. Weight can be optionally given as an array
 *  of numbers or an LAMatrix. If weight is not given, the atomic weights are used.
 *  Return values are the transform (that converts the present coordinates to the
 *  target coordinates) and root mean square deviation (without weight).
 */
static VALUE
s_Molecule_FitCoordinates(int argc, VALUE *argv, VALUE self)
{
	Molecule *mol;
	Atom *ap;
	VALUE gval, rval, wval;
	IntGroup *ig;
	IntGroupIterator iter;
	int nn, errno, i, j, in, status;
	Vector *ref;
	Double *weights, dval[3];
	Transform tr;

	Data_Get_Struct(self, Molecule, mol);
	rb_scan_args(argc, argv, "21", &gval, &rval, &wval);
	if (gval == Qnil)
		ig = IntGroupNewWithPoints(0, mol->natoms, -1);
	else
		ig = IntGroupFromValue(gval);
	if (ig == NULL || (nn = IntGroupGetCount(ig)) == 0) {
		IntGroupRelease(ig);
		rb_raise(rb_eMolbyError, "atom group is not given correctly");
	}
	ref = (Vector *)calloc(sizeof(Vector), nn);
	weights = (Double *)calloc(sizeof(Double), nn);
	IntGroupIteratorInit(ig, &iter);
	if (rb_obj_is_kind_of(rval, rb_cNumeric)) {
		int fn = NUM2INT(rb_Integer(rval));
		if (fn < 0 || fn >= MoleculeGetNumberOfFrames(mol)) {
			errno = 1;
			status = fn;
			goto err;
		}
		for (i = 0; (in = IntGroupIteratorNext(&iter)) >= 0; i++) {
			ap = ATOM_AT_INDEX(mol->atoms, in);
			if (fn < ap->nframes)
				ref[i] = ap->frames[fn];
			else ref[i] = ap->r;
		}
	} else if (rb_obj_is_kind_of(rval, rb_cLAMatrix)) {
		LAMatrix *m = LAMatrixFromValue(rval, NULL, 0, 0);
		if (m->row * m->column < nn * 3) {
			errno = 2;
			goto err;
		}
		for (i = 0; i < nn; i++) {
			ref[i].x = m->data[i * 3];
			ref[i].y = m->data[i * 3 + 1];
			ref[i].z = m->data[i * 3 + 2];
		}
	} else {
		VALUE aval;
		rval = rb_protect(rb_ary_to_ary, rval, &status);
		if (status != 0) {
			errno = 3;
			goto err;
		}
		if (RARRAY_LEN(rval) < nn) {
			errno = 2;
			goto err;
		}
		if (rb_obj_is_kind_of((RARRAY_PTR(rval))[0], rb_cNumeric)) {
			/*  Array of 3*nn numbers  */
			if (RARRAY_LEN(rval) < nn * 3) {
				errno = 2;
				goto err;
			}
			for (i = 0; i < nn; i++) {
				for (j = 0; j < 3; j++) {
					aval = rb_protect(rb_Float, (RARRAY_PTR(rval))[i * 3 + j], &status);
					if (status != 0) {
						errno = 3;
						goto err;
					}
					dval[j] = NUM2DBL(aval);
				}
				ref[i].x = dval[0];
				ref[i].y = dval[1];
				ref[i].z = dval[2];
			}
		} else {
			/*  Array of nn Vector3Ds or Arrays  */
			for (i = 0; i < nn; i++) {
				aval = (RARRAY_PTR(rval))[i];
				if (rb_obj_is_kind_of(aval, rb_cVector3D)) {
					VectorFromValue(aval, &ref[i]);
				} else {
					aval = rb_protect(rb_ary_to_ary, aval, &status);
					if (status != 0) {
						errno = 3;
						goto err;
					}
					if (RARRAY_LEN(aval) < 3) {
						errno = 4;
						status = i;
						goto err;
					}
					for (j = 0; j < 3; j++) {
						VALUE aaval = rb_protect(rb_Float, (RARRAY_PTR(aval))[j], &status);
						if (status != 0) {
							errno = 3;
							goto err;
						}
						dval[j] = NUM2DBL(aaval);
					}
					ref[i].x = dval[0];
					ref[i].y = dval[1];
					ref[i].z = dval[2];
				}
			}
		}
	}
	if (wval == Qnil) {
		/*  Use atomic weights  */
		IntGroupIteratorReset(&iter);
		for (i = 0; (in = IntGroupIteratorNext(&iter)) >= 0; i++) {
			ap = ATOM_AT_INDEX(mol->atoms, in);
			weights[i] = ap->weight;
		}
	} else {
		wval = rb_protect(rb_ary_to_ary, wval, &status);
		if (status != 0) {
			errno = 3;
			goto err;
		}
		if (RARRAY_LEN(wval) < nn) {
			errno = 5;
			goto err;
		}
		for (i = 0; i < nn; i++) {
			VALUE wwval = rb_protect(rb_Float, (RARRAY_PTR(wval))[i], &status);
			if (status != 0) {
				errno = 3;
				goto err;
			}
			weights[i] = NUM2DBL(wwval);
		}
	}
	dval[0] = s_Molecule_FitCoordinates_Sub(mol, ig, ref, weights, tr);
	if (dval[0] < 0) {
		errno = 6;
		goto err;
	}
	errno = 0;
err:
	IntGroupIteratorRelease(&iter);
	free(ref);
	free(weights);
	if (errno == 0) {
		return rb_ary_new3(2, ValueFromTransform(&tr), rb_float_new(dval[0]));
	} else if (errno == 1) {
		rb_raise(rb_eMolbyError, "frame index (%d) is out of range", status);
	} else if (errno == 2) {
		rb_raise(rb_eMolbyError, "insufficient number of reference coordinates");
	} else if (errno == 3) {
		rb_jump_tag(status);
	} else if (errno == 4) {
		rb_raise(rb_eMolbyError, "less than 3 elements for index %d of reference coordinates", status);
	} else if (errno == 5) {
		rb_raise(rb_eMolbyError, "insufficient number of weight values");
	} else if (errno == 6) {
		rb_raise(rb_eMolbyError, "matrix calculation failed during coordinate fitting");
	}
	return Qnil;  /*  Not reached  */
}

/*
 *  call-seq:
 *     display
 *
 *  Refresh the display if this molecule is bound to a view. Otherwise do nothing.
 */
static VALUE
s_Molecule_Display(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview != NULL)
		MainViewCallback_display(mol->mview);
	return Qnil;
}

/*
 *  call-seq:
 *     make_front
 *
 *  Make the window frontmost if this molecule is bound to a view. Otherwise do nothing.
 */
static VALUE
s_Molecule_MakeFront(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview != NULL)
		MainViewCallback_makeFront(mol->mview);
	return Qnil;
}

/*
 *  call-seq:
 *     update_enabled? -> bool
 *
 *  Returns true if screen update is enabled; otherwise no.
 */
static VALUE
s_Molecule_UpdateEnabled(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview != NULL && !mol->mview->freezeScreen)
		return Qtrue;
	else return Qfalse;
}

/*
 *  call-seq:
 *     update_enabled = bool
 *
 *  Enable or disable screen update. This is effective for automatic update on modification.
 *  Explicit call to molecule.display() always updates the screen.
 */
static VALUE
s_Molecule_SetUpdateEnabled(VALUE self, VALUE val)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	val = ((val != Qfalse && val != Qnil) ? Qtrue : Qfalse);
	if (mol->mview != NULL)
		mol->mview->freezeScreen = (val == Qfalse);
	else val = Qfalse;
	return val;
}

/*
 *  call-seq:
 *     show_unitcell
 *     show_unitcell(bool)
 *     show_unitcell = bool
 *
 *  Set the flag whether to show the unit cell. If no argument is given, the
 *  current flag is returned.
 */
static VALUE
s_Molecule_ShowUnitCell(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	if (argc > 0) {
		mol->mview->showUnitCell = (RTEST(argv[0]) != 0);
		MainViewCallback_setNeedsDisplay(mol->mview, 1);
	}
	return (mol->mview->showUnitCell ? Qtrue : Qfalse);
}

/*
 *  call-seq:
 *     show_hydrogens
 *     show_hydrogens(bool)
 *     show_hydrogens = bool
 *
 *  Set the flag whether to show the hydrogen atoms. If no argument is given, the
 *  current flag is returned.
 */
static VALUE
s_Molecule_ShowHydrogens(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	if (argc > 0) {
		mol->mview->showHydrogens = (RTEST(argv[0]) != 0);
		MainViewCallback_setNeedsDisplay(mol->mview, 1);
	}
	return (mol->mview->showHydrogens ? Qtrue : Qfalse);
}

/*
 *  call-seq:
 *     show_dummy_atoms
 *     show_dummy_atoms(bool)
 *     show_dummy_atoms = bool
 *
 *  Set the flag whether to show the dummy atoms. If no argument is given, the
 *  current flag is returned.
 */
static VALUE
s_Molecule_ShowDummyAtoms(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	if (argc > 0) {
		mol->mview->showDummyAtoms = (RTEST(argv[0]) != 0);
		MainViewCallback_setNeedsDisplay(mol->mview, 1);
	}
	return (mol->mview->showDummyAtoms ? Qtrue : Qfalse);
}

/*
 *  call-seq:
 *     show_expanded
 *     show_expanded(bool)
 *     show_expanded = bool
 *
 *  Set the flag whether to show the expanded atoms. If no argument is given, the
 *  current flag is returned.
 */
static VALUE
s_Molecule_ShowExpanded(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	if (argc > 0) {
		mol->mview->showExpandedAtoms = (RTEST(argv[0]) != 0);
		MainViewCallback_setNeedsDisplay(mol->mview, 1);
	}
	return (mol->mview->showExpandedAtoms ? Qtrue : Qfalse);
}

/*
 *  call-seq:
 *     show_ellipsoids
 *     show_ellipsoids(bool)
 *     show_ellipsoids = bool
 *
 *  Set the flag whether to show the thermal ellipsoids. If no argument is given, the
 *  current flag is returned.
 */
static VALUE
s_Molecule_ShowEllipsoids(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	if (argc > 0) {
		mol->mview->showEllipsoids = (RTEST(argv[0]) != 0);
		MainViewCallback_setNeedsDisplay(mol->mview, 1);
	}
	return (mol->mview->showEllipsoids ? Qtrue : Qfalse);
}

/*
 *  call-seq:
 *     is_atom_visible(index)  -> Boolean
 *
 *  Check is an atom is visible. It examines the atom attribute (ap->exflags & kAtomHiddenFlag)
 *  as well as the molecule attributes (showHydrogens, etc.)
 */
static VALUE
s_Molecule_IsAtomVisible(VALUE self, VALUE ival)
{
	Molecule *mol;
	Int idx;
	Atom *ap;
    Data_Get_Struct(self, Molecule, mol);
	idx = s_Molecule_AtomIndexFromValue(mol, ival);
	if (idx < 0 || idx >= mol->natoms)
		return Qnil;
	ap = ATOM_AT_INDEX(mol->atoms, idx);
	if (mol->mview != NULL) {
		if (mol->mview->showHydrogens == 0 && ap->atomicNumber == 1)
			return Qfalse;
		if (mol->mview->showExpandedAtoms == 0 && SYMOP_ALIVE(ap->symop))
			return Qfalse;
		if (mol->mview->showDummyAtoms == 0 && ap->atomicNumber == 0)
			return Qfalse;
	}
	return ((ap->exflags & kAtomHiddenFlag) != 0 ? Qfalse : Qtrue);
}

/*
 *  call-seq:
 *     show_graphite -> Integer
 *     show_graphite = Integer
 *     show_graphite = boolean
 *
 *  Set whether to show the graphite plane. If the argument is positive, it also indicates the
 *  number of rings to display for each direction.
 *  If the argument is boolean, only the show/hide flag is set.
 */
static VALUE
s_Molecule_ShowGraphite(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	if (argc > 0) {
		if (argv[0] == Qnil || argv[0] == Qfalse)
			mol->mview->showGraphiteFlag = 0;
		else if (argv[0] == Qtrue)
			mol->mview->showGraphiteFlag = 1;
		else {
			int n = NUM2INT(rb_Integer(argv[0]));
			if (n < 0)
				rb_raise(rb_eMolbyError, "The argument must be non-negative integer");
			mol->mview->showGraphite = n;
		}
		MainViewCallback_setNeedsDisplay(mol->mview, 1);
	}
	return INT2NUM(mol->mview->showGraphite);
}

/*
 *  call-seq:
 *     show_graphite? -> boolean
 *
 *  Return whether the graphite is set visible or not.
*/
static VALUE
s_Molecule_ShowGraphiteFlag(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	return (mol->mview->showGraphiteFlag ? Qtrue : Qfalse);
}
	
/*
 *  call-seq:
 *     show_periodic_image -> [amin, amax, bmin, bmax, cmin, cmax]
 *     show_periodic_image = [amin, amax, bmin, bmax, cmin, cmax]
 *     show_periodic_image = boolean
 *
 *  Set to show the periodic image of the atoms. If the unit cell is not defined, the values are
 *  set but no visual effects are observed.
 *  If the argument is boolean, only the show/hide flag is modified.
 */
static VALUE
s_Molecule_ShowPeriodicImage(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	VALUE val;
	int ival[6];
	int i;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	rb_scan_args(argc, argv, "01", &val);
	if (argc > 0) {
		/*  Change current settings  */
		if (val == Qnil || val == Qfalse)
			mol->mview->showPeriodicImageFlag = 0;
		else if (val == Qtrue)
			mol->mview->showPeriodicImageFlag = 1;
		else {
			val = rb_ary_to_ary(val);
			for (i = 0; i < 6; i++) {
				if (i < RARRAY_LEN(val))
					ival[i] = NUM2INT(rb_Integer(RARRAY_PTR(val)[i]));
			}
			if (ival[0] > 0 || ival[1] < 0 || ival[2] > 0 || ival[3] < 0 || ival[4] > 0 || ival[5] < 0)
				rb_raise(rb_eMolbyError, "bad arguments");
			for (i = 0; i < 6; i++)
				mol->mview->showPeriodicImage[i] = ival[i];
		}
		MainViewCallback_setNeedsDisplay(mol->mview, 1);
	}
	val = rb_ary_new();
	for (i = 0; i < 6; i++)
		rb_ary_push(val, INT2NUM(mol->mview->showPeriodicImage[i]));
	return val;
}

/*
 *  call-seq:
 *     show_periodic_image? -> boolean
 *
 *  Return whether the periodic images are set to visible or not. This flag is
 *  independent from the show_periodic_image settings.
 */
static VALUE
s_Molecule_ShowPeriodicImageFlag(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	return (mol->mview->showPeriodicImageFlag ? Qtrue : Qfalse);
}

/*
 *  call-seq:
 *     line_mode
 *     line_mode(bool)
 *     line_mode = bool
 *
 *  Set the flag whether to draw the model in line mode. If no argument is given, the
 *  current flag is returned.
 */
static VALUE
s_Molecule_LineMode(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	if (argc > 0) {
		mol->mview->lineMode = (RTEST(argv[0]) != 0);
		MainViewCallback_setNeedsDisplay(mol->mview, 1);
	}
	return (mol->mview->lineMode ? Qtrue : Qfalse);
}

/*
 *  call-seq:
 *     resize_to_fit
 *
 *  Resize the model drawing to fit in the window.
 */
static VALUE
s_Molecule_ResizeToFit(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview != NULL)
		MainView_resizeToFit(mol->mview);
	return self;	
}

/*
 *  call-seq:
 *     get_view_rotation -> [[ax, ay, az], angle]
 *
 *  Get the current rotation for the view. Angle is in degree, not radian.
 */
static VALUE
s_Molecule_GetViewRotation(VALUE self)
{
    Molecule *mol;
	float f[4];
	Vector v;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	TrackballGetRotate(mol->mview->track, f);
	f[0] = -f[0];  /*  Convert to left-handed screw (to be consistent with Transform)  */
	v.x = f[1];
	v.y = f[2];
	v.z = f[3];
	return rb_ary_new3(2, ValueFromVector(&v), rb_float_new(f[0]));
}

/*
 *  call-seq:
 *     get_view_scale -> float
 *
 *  Get the current scale for the view.
 */
static VALUE
s_Molecule_GetViewScale(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	return rb_float_new(TrackballGetScale(mol->mview->track));
}

/*
 *  call-seq:
 *     get_view_center -> Vector
 *
 *  Get the current center point of the view.
 */
static VALUE
s_Molecule_GetViewCenter(VALUE self)
{
    Molecule *mol;
	float f[4];
	Vector v;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	TrackballGetTranslate(mol->mview->track, f);
	v.x = -f[0] * mol->mview->dimension;
	v.y = -f[1] * mol->mview->dimension;
	v.z = -f[2] * mol->mview->dimension;
	return ValueFromVector(&v);
}

/*
 *  call-seq:
 *     set_view_rotation([ax, ay, az], angle) -> self
 *
 *  Set the current rotation for the view. Angle is in degree, not radian.
 */
static VALUE
s_Molecule_SetViewRotation(VALUE self, VALUE aval, VALUE angval)
{
    Molecule *mol;
	float f[4];
	Vector v;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	VectorFromValue(aval, &v);
	if (NormalizeVec(&v, &v))
		rb_raise(rb_eMolbyError, "Cannot normalize nearly zero vector");
	f[1] = v.x;
	f[2] = v.y;
	f[3] = v.z;
	f[0] = -NUM2DBL(rb_Float(angval));
	TrackballSetRotate(mol->mview->track, f);
	MainViewCallback_setNeedsDisplay(mol->mview, 1);
	return self;
}

/*
 *  call-seq:
 *     set_view_scale(scale) -> self
 *
 *  Set the current scale for the view.
 */
static VALUE
s_Molecule_SetViewScale(VALUE self, VALUE aval)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	TrackballSetScale(mol->mview->track, NUM2DBL(rb_Float(aval)));
	MainViewCallback_setNeedsDisplay(mol->mview, 1);
	return self;
}

/*
 *  call-seq:
 *     set_view_center(vec) -> self
 *
 *  Set the current center point of the view.
 */
static VALUE
s_Molecule_SetViewCenter(VALUE self, VALUE aval)
{
    Molecule *mol;
	Vector v;
	float f[4];
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	VectorFromValue(aval, &v);
	f[0] = -v.x / mol->mview->dimension;
	f[1] = -v.y / mol->mview->dimension;
	f[2] = -v.z / mol->mview->dimension;
	TrackballSetTranslate(mol->mview->track, f);
	MainViewCallback_setNeedsDisplay(mol->mview, 1);
	return self;
}

/*
 *  call-seq:
 *     set_background_color(red, green, blue)
 *
 *  Set the background color of the model window.
 */
static VALUE
s_Molecule_SetBackgroundColor(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview != NULL) {
		VALUE rval, gval, bval;
		rb_scan_args(argc, argv, "30", &rval, &gval, &bval);
		MainView_setBackgroundColor(mol->mview, NUM2DBL(rb_Float(rval)), NUM2DBL(rb_Float(gval)), NUM2DBL(rb_Float(bval)));
	}
	return self;	
}

/*
 *  call-seq:
 *     create_graphic(kind, color, points, fill = nil) -> integer
 *
 *  Create a new graphic object.
 *   kind: a symbol representing the kind of the graphic. :line, :poly, :cylinder, :cone, :ellipsoid
 *   color: an array of 3 (rgb) or 4 (rgba) floating numbers
 *   points: an array of Vectors
 *   
 */
static VALUE
s_Molecule_CreateGraphic(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	MainViewGraphic g;
	int i, n, ni;
	const char *p;
	VALUE kval, cval, pval, fval;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		rb_raise(rb_eMolbyError, "this molecule has no associated graphic view");
	rb_scan_args(argc, argv, "31", &kval, &cval, &pval, &fval);
	kval = rb_obj_as_string(kval);
	memset(&g, 0, sizeof(g));
	g.visible = 1;
	p = RSTRING_PTR(kval);
	if (strcmp(p, "line") == 0)
		g.kind = kMainViewGraphicLine;
	else if (strcmp(p, "poly") == 0)
		g.kind = kMainViewGraphicPoly;
	else if (strcmp(p, "cylinder") == 0)
		g.kind = kMainViewGraphicCylinder;
	else if (strcmp(p, "cone") == 0)
		g.kind = kMainViewGraphicCone;
	else if (strcmp(p, "ellipsoid") == 0)
		g.kind = kMainViewGraphicEllipsoid;
	else rb_raise(rb_eMolbyError, "unknown graphic object type: %s", p);
	g.closed = (RTEST(fval) ? 1 : 0);
	cval = rb_ary_to_ary(cval);
	n = RARRAY_LEN(cval);
	if (n < 3 || n >= 5)
		rb_raise(rb_eArgError, "the color should have 3 or 4 elements");
	if (n == 3)
		g.rgba[3] = 1.0;
	for (i = 0; i < n; i++)
		g.rgba[i] = NUM2DBL(rb_Float(RARRAY_PTR(cval)[i]));
	pval = rb_ary_to_ary(pval);
	n = RARRAY_LEN(pval);
	ni = -1;  /*  If this is non-negative, then ni-th control point is [number, 0, 0] */
	if (n <= 0)
		rb_raise(rb_eArgError, "no control points are given");
	switch (g.kind) {
		case kMainViewGraphicLine:
			if (n < 2)
				rb_raise(rb_eArgError, "the line object must have at least two control points");
			break;
		case kMainViewGraphicPoly:
			if (n < 3)
				rb_raise(rb_eArgError, "the polygon object must have at least three control points");
			break;
		case kMainViewGraphicCylinder:
		case kMainViewGraphicCone:
			if (n != 3)
				rb_raise(rb_eArgError, "the %s object must have two control points and one number (radius)", (g.kind == kMainViewGraphicCylinder ? "cylinder" : "cone"));
			ni = 2;
			break;
		case kMainViewGraphicEllipsoid:
			if (n == 2) {
				ni = 1;
			} else if (n != 4)
				rb_raise(rb_eArgError, "the ellipsoid object must have either one point and one number (radius) or four points (center and three main axes)");
			break;
	}
	NewArray(&g.points, &g.npoints, sizeof(GLfloat) * 3, n);
	for (i = 0; i < n; i++) {
		Vector v;
		if (i == ni) {
			v.x = NUM2DBL(rb_Float(RARRAY_PTR(pval)[i]));
			v.y = v.z = 0;
		} else {
			VectorFromValue(RARRAY_PTR(pval)[i], &v);
		}
		g.points[i * 3] = v.x;
		g.points[i * 3 + 1] = v.y;
		g.points[i * 3 + 2] = v.z;
	}
	if (g.kind == kMainViewGraphicEllipsoid && ni == 1) {
		/*  Sphere  */
		AssignArray(&g.points, &g.npoints, sizeof(GLfloat) * 3, 4, NULL);
		g.points[6] = g.points[8] = g.points[9] = g.points[10] = 0;
		g.points[7] = g.points[11] = g.points[3];
	}
	MainView_insertGraphic(mol->mview, -1, &g);
	return INT2NUM(mol->mview->ngraphics - 1);	
}

/*
 *  call-seq:
 *     remove_graphic(index) -> integer
 *
 *  Remove a graphic object.
 */
static VALUE
s_Molecule_RemoveGraphic(VALUE self, VALUE ival)
{
    Molecule *mol;
	int i;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		rb_raise(rb_eMolbyError, "this molecule has no associated graphic view");
	i = NUM2INT(rb_Integer(ival));
	if (i < 0 || i >= mol->mview->ngraphics)
		rb_raise(rb_eArgError, "graphic index is out of range");
	MainView_removeGraphic(mol->mview, i);
	return ival;
}

/*
 *  call-seq:
 *     ngraphics -> integer
 *
 *  Get the number of graphic objects.
 */
static VALUE
s_Molecule_NGraphics(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		rb_raise(rb_eMolbyError, "this molecule has no associated graphic view");
	return INT2NUM(mol->mview->ngraphics);
}
	
/*
 *  call-seq:
 *     set_graphic_point(graphic_index, point_index, new_value) -> new_value
 *
 *  Change the point_index-th control point of graphic_index-th graphic object
 *   
 */
static VALUE
s_Molecule_SetGraphicPoint(VALUE self, VALUE gval, VALUE pval, VALUE nval)
{
	MainViewGraphic *gp;
    Molecule *mol;
	int index;
	Vector v;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		rb_raise(rb_eMolbyError, "this molecule has no associated graphic view");
	index = NUM2INT(rb_Integer(gval));
	if (index < 0 || index >= mol->mview->ngraphics)
		rb_raise(rb_eArgError, "the graphic index is out of range");
	gp = mol->mview->graphics + index;
	index = NUM2INT(rb_Integer(pval));
	if (index < 0 || index >= gp->npoints)
		rb_raise(rb_eArgError, "the point index is out of range");
	if (rb_obj_is_kind_of(nval, rb_cNumeric)) {
		if ((gp->kind == kMainViewGraphicCylinder || gp->kind == kMainViewGraphicCone) && index == 2) {
			v.x = NUM2DBL(rb_Float(nval));
			v.y = v.z = 0;
		} else if (gp->kind == kMainViewGraphicEllipsoid && index == 1) {
			gp->points[3] = gp->points[7] = gp->points[11] = NUM2DBL(rb_Float(nval));
			gp->points[4] = gp->points[5] = gp->points[6] = gp->points[8] = gp->points[9] = gp->points[10] = 0;
			return nval;
		} else rb_raise(rb_eArgError, "the argument must be an array-like object");
	} else {
		if (nval == Qnil) {
			v.x = kInvalidFloat;
			v.y = v.z = 0.0;
		} else VectorFromValue(nval, &v);
	}
	gp->points[index * 3] = v.x;
	gp->points[index * 3 + 1] = v.y;
	gp->points[index * 3 + 2] = v.z;
	MoleculeCallback_notifyModification(mol, 0);
	return nval;
}

/*
 *  call-seq:
 *     set_graphic_color(graphic_index, new_value) -> new_value
 *
 *  Change the color of graphic_index-th graphic object
 *   
 */
static VALUE
s_Molecule_SetGraphicColor(VALUE self, VALUE gval, VALUE cval)
{
	MainViewGraphic *gp;
    Molecule *mol;
	int index, n;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		rb_raise(rb_eMolbyError, "this molecule has no associated graphic view");
	index = NUM2INT(rb_Integer(gval));
	if (index < 0 || index >= mol->mview->ngraphics)
		rb_raise(rb_eArgError, "the graphic index is out of range");
	gp = mol->mview->graphics + index;
	cval = rb_ary_to_ary(cval);
	n = RARRAY_LEN(cval);
	if (n != 3 && n != 4)
		rb_raise(rb_eArgError, "the color argument must have 3 or 4 numbers");
	for (index = 0; index < n; index++) {
		gp->rgba[index] = NUM2DBL(rb_Float(RARRAY_PTR(cval)[index]));
	}
	if (n == 3)
		gp->rgba[3] = 1.0;
	MoleculeCallback_notifyModification(mol, 0);
	return cval;
}

/*
 *  call-seq:
 *     show_graphic(graphic_index) -> self
 *
 *  Enable the visible flag of the graphic_index-th graphic object
 *   
 */
static VALUE
s_Molecule_ShowGraphic(VALUE self, VALUE gval)
{
	MainViewGraphic *gp;
    Molecule *mol;
	int index;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		rb_raise(rb_eMolbyError, "this molecule has no associated graphic view");
	index = NUM2INT(rb_Integer(gval));
	if (index < 0 || index >= mol->mview->ngraphics)
		rb_raise(rb_eArgError, "the graphic index is out of range");
	gp = mol->mview->graphics + index;
	gp->visible = 1;
	MoleculeCallback_notifyModification(mol, 0);
	return self;
}

/*
 *  call-seq:
 *     hide_graphic(graphic_index) -> self
 *
 *  Disable the visible flag of the graphic_index-th graphic object
 *   
 */
static VALUE
s_Molecule_HideGraphic(VALUE self, VALUE gval)
{
	MainViewGraphic *gp;
    Molecule *mol;
	int index;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		rb_raise(rb_eMolbyError, "this molecule has no associated graphic view");
	index = NUM2INT(rb_Integer(gval));
	if (index < 0 || index >= mol->mview->ngraphics)
		rb_raise(rb_eArgError, "the graphic index is out of range");
	gp = mol->mview->graphics + index;
	gp->visible = 0;
	MoleculeCallback_notifyModification(mol, 0);
	return self;
}

/*
 *  call-seq:
 *     show_text(string)
 *
 *  Show the string in the info text box.
 */
static VALUE
s_Molecule_ShowText(VALUE self, VALUE arg)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview != NULL)
		MainViewCallback_drawInfoText(mol->mview, StringValuePtr(arg));
	return Qnil;
}

/*
 *  call-seq:
 *     md_arena -> MDArena
 *
 *  Returns the MDArena object associated to this molecule. If no MDArena is associated to
 *  this molecule, a new arena is created.
 */
static VALUE
s_Molecule_MDArena(VALUE self)
{
    Molecule *mol;
	VALUE retval;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->arena == NULL)
		md_arena_new(mol);
	retval = ValueFromMDArena(mol->arena);
	return retval;
}

/*
 *  call-seq:
 *     set_parameter_attr(type, index, key, value, src) -> value
 *
 *  This method is used only internally.
 */
static VALUE
s_Molecule_SetParameterAttr(VALUE self, VALUE tval, VALUE ival, VALUE kval, VALUE vval, VALUE sval)
{
	/*  This method is called from MolAction to change a MM parameter attribute.  */
    Molecule *mol;
	VALUE pval;
	ParameterRef *pref;
	UnionPar *up;
    Data_Get_Struct(self, Molecule, mol);
	pval = ValueFromMoleculeWithParameterTypeAndIndex(mol, FIX2INT(tval), NUM2INT(ival));
	vval = s_ParameterRef_SetAttr(pval, kval, vval);
	
	/*  This is the special part of this method; it allows modification of the src field. */
	/*  (ParameterRef#set_attr sets 0 to the src field)  */
	Data_Get_Struct(pval, ParameterRef, pref);
	up = ParameterRefGetPar(pref);
	up->bond.src = FIX2INT(sval);
	
	return vval;
}

/*
 *  call-seq:
 *     parameter -> Parameter
 *
 *  Get the local parameter of this molecule. If not defined, returns nil.
 */
static VALUE
s_Molecule_Parameter(VALUE self)
{
    Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
/*	if (mol->par == NULL)
		return Qnil; */
	return s_NewParameterValueFromValue(self);
}

/*
 *  call-seq:
 *     selectedMO -> IntGroup
 *
 *  Returns a group of selected mo in the "MO Info" table. If the MO info table
 *  is not selected, returns nil. If the MO info table is selected but no MOs 
 *  are selected, returns an empty IntGroup. The numbers in the table are 1-based.
 */
static VALUE
s_Molecule_SelectedMO(VALUE self)
{
    Molecule *mol;
	IntGroup *ig;
	VALUE val;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->mview == NULL)
		return Qnil;
	ig = MainView_selectedMO(mol->mview);
	if (ig == NULL)
		return Qnil;
	IntGroupOffset(ig, 1);
	val = ValueFromIntGroup(ig);
	IntGroupRelease(ig);
	return val;
}

/*
 *  call-seq:
 *     default_MO_grid(npoints = 80*80*80) -> [origin, dx, dy, dz, nx, ny, nz]
 *
 *  Returns a default MO grid for cube file generation. Origin: Vector, dx, dy, dz: float, nx, ny, nz: integer.
 *  If the molecule does not contain a basis set information, then returns nil.
 */
static VALUE
s_Molecule_GetDefaultMOGrid(int argc, VALUE *argv, VALUE self)
{
    Molecule *mol;
	Vector o, dx, dy, dz;
	Int nx, ny, nz;
	VALUE nval;
	Int npoints = 80 * 80 * 80;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->bset == NULL)
		return Qnil;
	rb_scan_args(argc, argv, "01", &nval);
	if (nval != Qnil)
		npoints = NUM2INT(rb_Integer(nval));
	if (MoleculeGetDefaultMOGrid(mol, npoints, &o, &dx, &dy, &dz, &nx, &ny, &nz) != 0)
		return Qnil;
	return rb_ary_new3(7, ValueFromVector(&o), rb_float_new(VecLength(dx)), rb_float_new(VecLength(dy)), rb_float_new(VecLength(dz)), INT2NUM(nx), INT2NUM(ny), INT2NUM(nz));
}

static int
s_Cubegen_callback(double progress, void *ref)
{
	MyAppCallback_setProgressValue(progress);
	if (MyAppCallback_checkInterrupt())
		return 1;
	else return 0;
}

/*
 *  call-seq:
 *     cubegen(fname, mo, npoints=1000000 [, iflag [, beta]])
 *     cubegen(fname, mo, origin, dx, dy, dz, nx, ny, nz [, iflag [, beta]])
 *
 *  Calculate the molecular orbital with number mo and create a 'cube' file.
 *  In the first form, the cube size is estimated from the atomic coordinates. In the
 *  second form, the cube dimension is explicitly given.
 *  Returns fname when successful, nil otherwise.
 *  If iflag is non-false, then MyAppCallback_checkInterrupt() is periodically called during calculation.
 *  If beta is non-false, then the beta MO is calculated (valid only for UHF MO)
 *  (The interrupt monitoring thread is disabled, so as not to enter Ruby interpreter during calculation)
 */
static VALUE
s_Molecule_Cubegen(int argc, VALUE *argv, VALUE self)
{
	VALUE fval, mval, oval, dxval, dyval, dzval, nxval, nyval, nzval, ival, bval;
    Molecule *mol;
	Int mono, nx, ny, nz, npoints;
	Vector o, dx, dy, dz;
	int index, n;
	char buf[1024];
    Data_Get_Struct(self, Molecule, mol);
	if (mol->bset == NULL)
		rb_raise(rb_eMolbyError, "The molecule does not contain MO information");
	rb_scan_args(argc, argv, "29", &fval, &mval, &oval, &dxval, &dyval, &dzval, &nxval, &nyval, &nzval, &ival, &bval);
	
	/*  Set up parameters  */
	mono = NUM2INT(rb_Integer(mval));
	if (mono <= 0 || mono > mol->bset->ncomps)
		rb_raise(rb_eMolbyError, "The MO number (%d) is out of range (should be 1..%d)", mono, mol->bset->ncomps);
	if (RTEST(bval)) {
		if (mol->bset->rflag != 0)
			rb_raise(rb_eMolbyError, "Beta MO is requested but not present");
		mono += mol->bset->ncomps;
	}
		
	if (oval == Qnil || dxval == Qnil || dyval == Qnil) {
		/*  Automatic grid formation  */
		if (oval != Qnil)
			npoints = NUM2INT(rb_Integer(oval));
		else npoints = 0;
		if (npoints == 0)
			npoints = 1000000;
		else if (npoints < 8)
			rb_raise(rb_eMolbyError, "The number of points (%d) should be at least 8", npoints);
		if (MoleculeGetDefaultMOGrid(mol, npoints, &o, &dx, &dy, &dz, &nx, &ny, &nz) != 0)
			rb_raise(rb_eMolbyError, "Cannot determine cube grids");
		ival = dxval;
		bval = dyval;
	} else {
		VectorFromValue(oval, &o);
		if (TYPE(dxval) == T_ARRAY || rb_obj_is_kind_of(dxval, rb_cVector3D))
			VectorFromValue(dxval, &dx);
		else {
			dx.x = NUM2DBL(rb_Float(dxval));
			dx.y = dx.z = 0.0;
		}
		if (TYPE(dyval) == T_ARRAY || rb_obj_is_kind_of(dyval, rb_cVector3D))
			VectorFromValue(dyval, &dy);
		else {
			dy.y = NUM2DBL(rb_Float(dyval));
			dy.x = dy.z = 0.0;
		}
		if (TYPE(dzval) == T_ARRAY || rb_obj_is_kind_of(dzval, rb_cVector3D))
			VectorFromValue(dzval, &dz);
		else {
			dz.z = NUM2DBL(rb_Float(dzval));
			dz.x = dz.y = 0.0;
		}
		nx = NUM2INT(rb_Integer(nxval));
		ny = NUM2INT(rb_Integer(nyval));
		nz = NUM2INT(rb_Integer(nzval));
		if (nx <= 0 || ny <= 0 || nz <= 0)
			rb_raise(rb_eMolbyError, "The number of points (%d,%d,%d) should be positive integers", nx, ny, nz);
		else if (nx > 10000 || ny > 10000 || nz > 10000 || nx * ny * nz > 100000000)
			rb_raise(rb_eMolbyError, "The number of points (%d,%d,%d) seem to be too large; please keep them less than 10000 and nx*ny*nz < 100000000", nx, ny, nz);
	}
	
	/*  Calc MO  */
	index = MoleculeCalcMO(mol, mono, &o, &dx, &dy, &dz, nx, ny, nz, (RTEST(ival) ? s_Cubegen_callback : NULL), NULL);
	if (index == -2)
		rb_interrupt();
	else if (index < 0)
		rb_raise(rb_eMolbyError, "Cannot calculate MO %d", mono);
	
	/*  Output to file  */
	MoleculeCallback_displayName(mol, buf, sizeof buf);
	n = MoleculeOutputCube(mol, index, FileStringValuePtr(fval), buf);
	if (n != 0)
		rb_raise(rb_eMolbyError, "Cannot output MO cube to file %p", FileStringValuePtr(fval));
	
	/*  Discard the cube  */
	MoleculeClearCubeAtIndex(mol, index);
	return fval;
}

/*
 *  call-seq:
 *     nelpots
 *
 *  Get the number of electrostatic potential info.
 */
static VALUE
s_Molecule_NElpots(VALUE self)
{
	Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	return INT2NUM(mol->nelpots);
}

/*
 *  call-seq:
 *     elpot(idx)
 *
 *  Get the electrostatic potential info at the given index. If present, then the
 *  return value is [Vector, Float] (position and potential). If not present, then
 *  returns nil.
 */
static VALUE
s_Molecule_Elpot(VALUE self, VALUE ival)
{
	Molecule *mol;
	int idx;
    Data_Get_Struct(self, Molecule, mol);
	idx = NUM2INT(rb_Integer(ival));
	if (idx < 0 || idx >= mol->nelpots)
		return Qnil;
	return rb_ary_new3(2, ValueFromVector(&mol->elpots[idx].pos), rb_float_new(mol->elpots[idx].esp));
}

/*
 *  call-seq:
 *     add_gaussian_orbital_shell(sym, nprims, atom_index)
 *
 *  To be used internally. Add a gaussian orbital shell with symmetry code, number of primitives,
 *  and the corresponding atom index. Symmetry code: 0, S-type; 1, P-type; -1, SP-type; 2, D-type;
 *  -2, D5-type.
 */
static VALUE
s_Molecule_AddGaussianOrbitalShell(VALUE self, VALUE symval, VALUE npval, VALUE aval)
{
	Molecule *mol;
	int sym, nprims, a_idx, n;
    Data_Get_Struct(self, Molecule, mol);
	sym = NUM2INT(rb_Integer(symval));
	nprims = NUM2INT(rb_Integer(npval));
	a_idx = NUM2INT(rb_Integer(aval));
	n = MoleculeAddGaussianOrbitalShell(mol, sym, nprims, a_idx);
	if (n == -1)
		rb_raise(rb_eMolbyError, "Molecule is emptry");
	else if (n == -2)
		rb_raise(rb_eMolbyError, "Low memory");
	else if (n == -3)
		rb_raise(rb_eMolbyError, "Unknown orbital type");
	else if (n != 0)
		rb_raise(rb_eMolbyError, "Unknown error");
	return self;
}

/*
 *  call-seq:
 *     add_gaussian_primitive_coefficients(exponent, contraction, contraction_sp)
 *
 *  To be used internally. Add a gaussian primitive coefficients.
 */
static VALUE
s_Molecule_AddGaussianPrimitiveCoefficients(VALUE self, VALUE expval, VALUE cval, VALUE cspval)
{
	Molecule *mol;
	Int n;
	Double exponent, contraction, contraction_sp;
    Data_Get_Struct(self, Molecule, mol);
	exponent = NUM2DBL(rb_Float(expval));
	contraction = NUM2DBL(rb_Float(cval));
	contraction_sp = NUM2DBL(rb_Float(cspval));
	n = MoleculeAddGaussianPrimitiveCoefficients(mol, exponent, contraction, contraction_sp);
	if (n == -1)
		rb_raise(rb_eMolbyError, "Molecule is emptry");
	else if (n == -2)
		rb_raise(rb_eMolbyError, "Low memory");
	else if (n != 0)
		rb_raise(rb_eMolbyError, "Unknown error");
	return self;
}

/*
 *  call-seq:
 *     mo_type
 *
 *  Returns either "RHF", "UHF", or "ROHF". If no MO info is present, returns nil.
 */
static VALUE
s_Molecule_MOType(VALUE self)
{
	Molecule *mol;
    Data_Get_Struct(self, Molecule, mol);
	if (mol != NULL && mol->bset != NULL) {
		const char *s;
		int rflag = mol->bset->rflag;
		if (rflag == 0)
			s = "UHF";
		else if (rflag == 2)
			s = "ROHF";
		else s = "RHF";
		return rb_str_new2(s);
	} else return Qnil;
}

/*
 *  call-seq:
 *     set_mo_coefficients(idx, energy, coefficients)
 *
 *  To be used internally. Add a MO coefficients. Idx is the MO index (for open shell system, 
 *  beta MOs comes after all alpha MOs), energy is the MO energy, coefficients is an array
 *  of MO coefficients.
 */
static VALUE
s_Molecule_SetMOCoefficients(VALUE self, VALUE ival, VALUE eval, VALUE aval)
{
	Molecule *mol;
	Int idx, ncomps, i;
	Double energy;
	Double *coeffs;
    Data_Get_Struct(self, Molecule, mol);
	idx = NUM2INT(rb_Integer(ival));
	energy = NUM2DBL(rb_Float(eval));
	aval = rb_ary_to_ary(aval);
	ncomps = RARRAY_LEN(aval);
	coeffs = (Double *)calloc(sizeof(Double), ncomps);
	if (coeffs == NULL) {
		i = -2;
		goto end;
	}
	for (i = 0; i < ncomps; i++)
		coeffs[i] = NUM2DBL(rb_Float(RARRAY_PTR(aval)[i]));
	i = MoleculeSetMOCoefficients(mol, idx, energy, ncomps, coeffs);
end:
	if (i == -1)
		rb_raise(rb_eMolbyError, "Molecule is emptry");
	else if (i == -2)
		rb_raise(rb_eMolbyError, "Low memory");
	else if (i == -3)
		rb_raise(rb_eMolbyError, "Bad or inconsistent number of MOs");
	else if (i == -4)
		rb_raise(rb_eMolbyError, "Bad MO index");
	else if (i == -5)
		rb_raise(rb_eMolbyError, "Insufficient number of coefficients are given");
	else if (i != 0)
		rb_raise(rb_eMolbyError, "Unknown error");
	return self;
}

/*
 *  call-seq:
 *     allocate_basis_set_record(rflag, ne_alpha, ne_beta)
 *
 *  To be used internally. Allocate a basis set record. rflag: 0, unrestricted; 1, restricted.
 *  ne_alpha, ne_beta: number of alpha/beta electrons.
 */
static VALUE
s_Molecule_AllocateBasisSetRecord(VALUE self, VALUE rval, VALUE naval, VALUE nbval)
{
	Molecule *mol;
	Int rflag, na, nb, n;
    Data_Get_Struct(self, Molecule, mol);
	rflag = NUM2INT(rb_Integer(rval));
	na = NUM2INT(rb_Integer(naval));
	nb = NUM2INT(rb_Integer(nbval));
	n = MoleculeAllocateBasisSetRecord(mol, rflag, na, nb);
	if (n == -1)
		rb_raise(rb_eMolbyError, "Molecule is emptry");
	else if (n == -2)
		rb_raise(rb_eMolbyError, "Low memory");
	else if (n != 0)
		rb_raise(rb_eMolbyError, "Unknown error");
	return self;
}

/*
 *  call-seq:
 *     search_equivalent_atoms(ig = nil)
 *
 *  Search equivalent atoms (within the atom group if given). Returns an array of integers.
 */
static VALUE
s_Molecule_SearchEquivalentAtoms(int argc, VALUE *argv, VALUE self)
{
	Molecule *mol;
	Int *result, i;
	VALUE val;
	IntGroup *ig;
    Data_Get_Struct(self, Molecule, mol);
	if (mol->natoms == 0)
		return Qnil;
	rb_scan_args(argc, argv, "01", &val);
	if (val != Qnil)
		ig = IntGroupFromValue(val);
	else ig = NULL;
	result = MoleculeSearchEquivalentAtoms(mol, ig);
	if (result == NULL)
		rb_raise(rb_eMolbyError, "Failed to search equivalent atoms (out of memory?)");
	if (ig != NULL)
		IntGroupRelease(ig);
	val = rb_ary_new2(mol->natoms);
	for (i = 0; i < mol->natoms; i++)
		rb_ary_push(val, INT2NUM(result[i]));
	free(result);
	return val;
}

/*
 *  call-seq:
 *     create_pi_anchor(name, group [, type] [, weights] [, index]) -> AtomRef
 *
 *  Create a pi anchor, which is a "dummy" atom to represent pi-metal bonds.
 *  Name is the name of the new pi anchor, and group is the atoms that define
 *  the pi system. Type (a String) is an atom type for MM implementation.
 *  Weights represent the relative significance of the component atoms; if omitted, then
 *  1.0/n (n is the number of component atoms) is assumed for all atoms.
 *  The weight values will be normalized so that the sum of the weights is 1.0.
 *  The weight values must be positive.
 *  Index is the atom index where the created pi-anchor is inserted in the 
 *  atoms array; if omitted, the pi-anchor is inserted after the component atom
 *  having the largest index.
 *  Pi anchors are appear in the atom list along with other ordinary atoms. The list
 *  of component atoms (and their weights) can be retrived by AtomRef#anchor_list.
 */
static VALUE
s_Molecule_CreatePiAnchor(int argc, VALUE *argv, VALUE self)
{
	Molecule *mol;
	VALUE nval, gval;
	IntGroup *ig;
	Int i, n, idx, last_component;
	Atom a, *ap;
	PiAnchor an;
	AtomRef *aref;
	if (argc < 2 || argc >= 6)
		rb_raise(rb_eMolbyError, "too %s arguments (should be 2..5)", (argc < 2 ? "few" : "many"));
	nval = *argv++;
	gval = *argv++;
	argc -= 2;
    Data_Get_Struct(self, Molecule, mol);
	ig = IntGroupFromValue(gval);
	memset(&a, 0, sizeof(a));
	memset(&an, 0, sizeof(an));
	strncpy(a.aname, StringValuePtr(nval), 4);
	if (a.aname[0] == '_')
		rb_raise(rb_eMolbyError, "please avoid a name beginning with '_' for a pi anchor, because it causes some internal confusion.");
	a.type = AtomTypeEncodeToUInt("##");  /*  Default atom type for pi_anchor  */
	for (i = 0; (n = IntGroupGetNthPoint(ig, i)) >= 0; i++) {
		if (n >= mol->natoms) {
			AtomConnectResize(&an.connect, 0);
			rb_raise(rb_eMolbyError, "atom index (%d) out of range", n);
		}
		AtomConnectInsertEntry(&an.connect, an.connect.count, n);
		last_component = n;
	}
	if (an.connect.count == 0)
		rb_raise(rb_eMolbyError, "no atoms are specified");
	NewArray(&an.coeffs, &an.ncoeffs, sizeof(Double), an.connect.count);
	for (i = 0; i < an.connect.count; i++) {
		an.coeffs[i] = 1.0 / an.connect.count;
	}
	if (argc > 0 && (argv[0] == Qnil || rb_obj_is_kind_of(argv[0], rb_cString))) {
		/*  Atom type  */
		if (argv[0] != Qnil)
			a.type = AtomTypeEncodeToUInt(StringValuePtr(argv[0]));
		argc--;
		argv++;
	}
	if (argc > 0 && (argv[0] == Qnil || rb_obj_is_kind_of(argv[0], rb_mEnumerable))) {
		if (argv[0] != Qnil) {
			VALUE aval = rb_ary_to_ary(argv[0]);
			Double d, sum;
			if (RARRAY_LEN(aval) != an.connect.count)
				rb_raise(rb_eMolbyError, "the number of weight values does not match the number of atoms");
			for (i = 0, sum = 0.0; i < an.connect.count; i++) {
				d = NUM2DBL(rb_Float(RARRAY_PTR(aval)[i]));
				if (d <= 0.0)
					rb_raise(rb_eMolbyError, "the weight value must be positive");
				sum += d;
				an.coeffs[i] = d;
			}
			for (i = 0; i < an.connect.count; i++)
				an.coeffs[i] /= sum;
		}
		argc--;
		argv++;
	}
	if (argc > 0 && argv[0] != Qnil) {
		/*  Index  */
		idx = NUM2INT(rb_Integer(argv[0]));
	} else idx = -1;
	if (idx < 0 || idx > mol->natoms) {
		/*  Immediately after the last specified atom  */
		idx = last_component + 1;
	}
	a.anchor = (PiAnchor *)malloc(sizeof(PiAnchor));
	memmove(a.anchor, &an, sizeof(PiAnchor));
	/*  Use residue information of the last specified atom  */
	ap = ATOM_AT_INDEX(mol->atoms, last_component);
	a.resSeq = ap->resSeq;
	strncpy(a.resName, ap->resName, 4);
	if (MolActionCreateAndPerform(mol, gMolActionAddAnAtom, &a, idx, &idx) != 0)
		return Qnil;
	MoleculeCalculatePiAnchorPosition(mol, idx);
    aref = AtomRefNew(mol, idx);
    return Data_Wrap_Struct(rb_cAtomRef, 0, (void (*)(void *))AtomRefRelease, aref);
}

/*
 *  call-seq:
 *     current       -> Molecule
 *
 *  Get the currently "active" molecule.
 */
static VALUE
s_Molecule_Current(VALUE klass)
{
	return ValueFromMolecule(MoleculeCallback_currentMolecule());
}

/*
 *  call-seq:
 *     Molecule[]          -> Molecule
 *     Molecule[n]         -> Molecule
 *     Molecule[name]      -> Molecule
 *     Molecule[name, k]   -> Molecule
 *     Molecule[regex]     -> Molecule
 *     Molecule[regex, k]  -> Molecule
 *
 *  Molecule[] is equivalent to Molecule.current.
 *  Molecule[n] (n is an integer) is equivalent to Molecule.list[n].
 *  Molecule[name] gives the first document (in the order of creation time) that has
 *  the given name. If a second argument (k) is given, the k-th document that has the
 *  given name is returned.
 *  Molecule[regex] gives the first document (in the order of creation time) that
 *  has a name matching the regular expression. If a second argument (k) is given, 
 *  the k-th document that has a name matching the re is returned.
 */
static VALUE
s_Molecule_MoleculeAtIndex(int argc, VALUE *argv, VALUE klass)
{
	VALUE val, kval;
	int idx, k;
	Molecule *mol;
	char buf[1024];
	rb_scan_args(argc, argv, "02", &val, &kval);
	if (val == Qnil)
		return s_Molecule_Current(klass);
	if (rb_obj_is_kind_of(val, rb_cInteger)) {
		idx = NUM2INT(val);
		mol = MoleculeCallback_moleculeAtIndex(idx);
	} else if (rb_obj_is_kind_of(val, rb_cString)) {
		char *p = StringValuePtr(val);
		k = (kval == Qnil ? 1 : NUM2INT(rb_Integer(kval)));
		for (idx = 0; (mol = MoleculeCallback_moleculeAtIndex(idx)) != NULL; idx++) {
			MoleculeCallback_displayName(mol, buf, sizeof buf);
			if (strcmp(buf, p) == 0 && --k == 0)
				break;
		}
	} else if (rb_obj_is_kind_of(val, rb_cRegexp)) {
		k = (kval == Qnil ? 1 : NUM2INT(rb_Integer(kval)));
		for (idx = 0; (mol = MoleculeCallback_moleculeAtIndex(idx)) != NULL; idx++) {
			VALUE name;
			MoleculeCallback_displayName(mol, buf, sizeof buf);
			name = rb_str_new2(buf);
			if (rb_reg_match(val, name) != Qnil && --k == 0)
				break;
		}	
	} else rb_raise(rb_eTypeError, "the argument should be either empty, an integer, a string, or a regular expression");
	
	if (mol == NULL)
		return Qnil;
	else return ValueFromMolecule(mol);
}

/*
 *  call-seq:
 *     list         -> array of Molecules
 *
 *  Get the list of molecules associated to the documents, in the order of creation
 *  time of the document. If no document is open, returns an empry array.
 */
static VALUE
s_Molecule_List(VALUE klass)
{
	Molecule *mol;
	int i;
	VALUE ary;
	i = 0;
	ary = rb_ary_new();
	while ((mol = MoleculeCallback_moleculeAtIndex(i)) != NULL) {
		rb_ary_push(ary, ValueFromMolecule(mol));
		i++;
	}
	return ary;
}

/*
 *  call-seq:
 *     ordered_list         -> array of Molecules
 *
 *  Get the list of molecules associated to the documents, in the order of front-to-back
 *  ordering of the associated window. If no document is open, returns an empry array.
 */
static VALUE
s_Molecule_OrderedList(VALUE klass)
{
	Molecule *mol;
	int i;
	VALUE ary;
	i = 0;
	ary = rb_ary_new();
	while ((mol = MoleculeCallback_moleculeAtOrderedIndex(i)) != NULL) {
		rb_ary_push(ary, ValueFromMolecule(mol));
		i++;
	}
	return ary;
}

/*
 *  call-seq:
 *     error_message       -> String
 *
 *  Get the error_message from the last load/save method. If no error, returns nil.
 */
static VALUE
s_Molecule_ErrorMessage(VALUE klass)
{
	if (gLoadSaveErrorMessage == NULL)
		return Qnil;
	else return rb_str_new2(gLoadSaveErrorMessage);
}

/*
 *  call-seq:
 *     set_error_message(String)
 *     Molecule.error_message = String
 *
 *  Get the error_message from the last load/save method. If no error, returns nil.
 */
static VALUE
s_Molecule_SetErrorMessage(VALUE klass, VALUE sval)
{
	if (gLoadSaveErrorMessage != NULL) {
		free(gLoadSaveErrorMessage);
		gLoadSaveErrorMessage = NULL;
	}
	if (sval != Qnil) {
		sval = rb_str_to_str(sval);
		gLoadSaveErrorMessage = strdup(StringValuePtr(sval));
	}
	return sval;
}

/*
 *  call-seq:
 *     self == Molecule -> boolean
 *
 *  True if the two arguments point to the same molecule.
 */
static VALUE
s_Molecule_Equal(VALUE self, VALUE val)
{
	if (rb_obj_is_kind_of(val, rb_cMolecule)) {
		Molecule *mol1, *mol2;
		Data_Get_Struct(self, Molecule, mol1);
		Data_Get_Struct(val, Molecule, mol2);
		return (mol1 == mol2 ? Qtrue : Qfalse);
	} else return Qfalse;
}

/*  The callback functions for call_subprocess_async  */
static int
s_Molecule_CallSubProcessAsync_EndCallback(Molecule *mol, int status)
{
	int ruby_status;
	VALUE procval, retval, args[2];
	args[0] = ValueFromMolecule(mol);
	args[1] = INT2NUM(status);
	procval = rb_ivar_get(args[0], rb_intern("end_proc"));
	if (procval != Qnil) {
		retval = Ruby_funcall2_protect(procval, rb_intern("call"), 2, args, &ruby_status);
		if (ruby_status != 0 || retval == Qnil || retval == Qfalse)
			return 1;
	}
	return 0;
}

static int
s_Molecule_CallSubProcessAsync_TimerCallback(Molecule *mol, int tcount)
{
	int ruby_status;
	VALUE procval, retval, args[2];
	args[0] = ValueFromMolecule(mol);
	args[1] = INT2NUM(tcount);
	procval = rb_ivar_get(args[0], rb_intern("timer_proc"));
	if (procval != Qnil) {
		retval = Ruby_funcall2_protect(procval, rb_intern("call"), 2, args, &ruby_status);
		if (ruby_status != 0 || retval == Qnil || retval == Qfalse)
			return 1;
	}
	return 0;
}

/*
 *  call-seq:
 *     call_subprocess_async(cmd [, end_callback [, timer_callback [, standard_output_file [, error_output_file]]]])
 *
 *  Call subprocess asynchronically.
 *  If end_callback is given, it will be called (with two arguments self and termination status)
 *  when the subprocess terminated.
 *  If timer_callback is given, it will be called (also with two arguments, self and timer count).
 *  If timer_callback returns nil or false, then the subprocess will be interrupted.
 *  If stdout_file or stderr_file is a filename, then the message will be sent to the file; if the
 *  filename begins with ">>", then the message will be appended to the file.
 *  If the filename is "/dev/null" or "NUL", then the message will be lost.
 *  If the argument is nil, then the message will be sent to the Ruby console.
 *  Returns the process ID as an integer.
 */
static VALUE
s_Molecule_CallSubProcessAsync(int argc, VALUE *argv, VALUE self)
{
	VALUE cmd, end_proc, timer_proc, stdout_val, stderr_val;
	Molecule *mol;
	char *sout, *serr;
	int n;
	FILE *fpout, *fperr;
	rb_scan_args(argc, argv, "14", &cmd, &end_proc, &timer_proc, &stdout_val, &stderr_val);
	Data_Get_Struct(self, Molecule, mol);

	if (stdout_val == Qnil) {
		fpout = (FILE *)1;
	} else {
		sout = StringValuePtr(stdout_val);
		if (strcmp(sout, "/dev/null") == 0 || strcmp(sout, "NUL") == 0)
			fpout = NULL;
		else {
			if (strncmp(sout, ">>", 2) == 0) {
				sout += 2;
				fpout = fopen(sout, "a");
			} else {
				if (*sout == '>')
					sout++;
				fpout = fopen(sout, "w");
			}
			if (fpout == NULL)
				rb_raise(rb_eMolbyError, "Cannot open file for standard output (%s)", sout);
		}
	}
	if (stderr_val == Qnil) {
		fperr = (FILE *)1;
	} else {
		serr = StringValuePtr(stderr_val);
		if (strcmp(serr, "/dev/null") == 0 || strcmp(serr, "NUL") == 0)
			fperr = NULL;
		else {
			if (strncmp(serr, ">>", 2) == 0) {
				serr += 2;
				fpout = fopen(serr, "a");
			} else {
				if (*serr == '>')
					serr++;
				fperr = fopen(serr, "w");
			}
			if (fperr == NULL)
				rb_raise(rb_eMolbyError, "Cannot open file for standard output (%s)", serr);
		}
	}
	
	/*  Register procs as instance variables  */
	rb_ivar_set(self, rb_intern("end_proc"), end_proc);
	rb_ivar_set(self, rb_intern("timer_proc"), timer_proc);
	n = MoleculeCallback_callSubProcessAsync(mol, StringValuePtr(cmd), s_Molecule_CallSubProcessAsync_EndCallback, (timer_proc == Qnil ? NULL : s_Molecule_CallSubProcessAsync_TimerCallback), fpout, fperr);
	if (fpout != NULL && fpout != (FILE *)1)
		fclose(fpout);
	if (fperr != NULL && fperr != (FILE *)1)
		fclose(fperr);
	return INT2NUM(n);
}

void
Init_Molby(void)
{
	int i;
	
	/*  Define module Molby  */
	rb_mMolby = rb_define_module("Molby");
	
	/*  Define Vector3D, Transform, IntGroup  */
	Init_MolbyTypes();
	
	/*  Define MDArena  */
	Init_MolbyMDTypes();

	/*  class Molecule  */
	rb_cMolecule = rb_define_class_under(rb_mMolby, "Molecule", rb_cObject);
	rb_define_alloc_func(rb_cMolecule, s_Molecule_Alloc);
    rb_define_private_method(rb_cMolecule, "initialize", s_Molecule_Initialize, -1);
    rb_define_private_method(rb_cMolecule, "initialize_copy", s_Molecule_InitCopy, 1);
    rb_define_method(rb_cMolecule, "loadmbsf", s_Molecule_Loadmbsf, -1);
    rb_define_method(rb_cMolecule, "loadpsf", s_Molecule_Loadpsf, -1);
    rb_define_alias(rb_cMolecule, "loadpsfx", "loadpsf");
    rb_define_method(rb_cMolecule, "loadpdb", s_Molecule_Loadpdb, -1);
    rb_define_method(rb_cMolecule, "loaddcd", s_Molecule_Loaddcd, -1);
    rb_define_method(rb_cMolecule, "loadtep", s_Molecule_Loadtep, -1);
    rb_define_method(rb_cMolecule, "loadres", s_Molecule_Loadres, -1);
    rb_define_method(rb_cMolecule, "loadfchk", s_Molecule_Loadfchk, -1);
    rb_define_alias(rb_cMolecule, "loadfch", "loadfchk");
    rb_define_method(rb_cMolecule, "loaddat", s_Molecule_Loaddat, -1);
    rb_define_method(rb_cMolecule, "molload", s_Molecule_Load, -1);
    rb_define_method(rb_cMolecule, "molsave", s_Molecule_Save, -1);
	rb_define_method(rb_cMolecule, "savembsf", s_Molecule_Savembsf, 1);
    rb_define_method(rb_cMolecule, "savepsf", s_Molecule_Savepsf, 1);
    rb_define_alias(rb_cMolecule, "savepsfx", "savepsf");
    rb_define_method(rb_cMolecule, "savepdb", s_Molecule_Savepdb, 1);
    rb_define_method(rb_cMolecule, "savedcd", s_Molecule_Savedcd, 1);
/*    rb_define_method(rb_cMolecule, "savetep", s_Molecule_Savetep, 1); */
    rb_define_method(rb_cMolecule, "name", s_Molecule_Name, 0);
	rb_define_method(rb_cMolecule, "set_name", s_Molecule_SetName, 1);
    rb_define_method(rb_cMolecule, "path", s_Molecule_Path, 0);
    rb_define_method(rb_cMolecule, "dir", s_Molecule_Dir, 0);
    rb_define_method(rb_cMolecule, "inspect", s_Molecule_Inspect, 0);
    rb_define_method(rb_cMolecule, "atoms", s_Molecule_Atoms, 0);
    rb_define_method(rb_cMolecule, "bonds", s_Molecule_Bonds, 0);
    rb_define_method(rb_cMolecule, "angles", s_Molecule_Angles, 0);
    rb_define_method(rb_cMolecule, "dihedrals", s_Molecule_Dihedrals, 0);
    rb_define_method(rb_cMolecule, "impropers", s_Molecule_Impropers, 0);
    rb_define_method(rb_cMolecule, "residues", s_Molecule_Residues, 0);
	rb_define_method(rb_cMolecule, "natoms", s_Molecule_Natoms, 0);
	rb_define_method(rb_cMolecule, "nbonds", s_Molecule_Nbonds, 0);
	rb_define_method(rb_cMolecule, "nangles", s_Molecule_Nangles, 0);
	rb_define_method(rb_cMolecule, "ndihedrals", s_Molecule_Ndihedrals, 0);
	rb_define_method(rb_cMolecule, "nimpropers", s_Molecule_Nimpropers, 0);

	rb_define_method(rb_cMolecule, "bond_par", s_Molecule_BondParIsObsolete, 1);
	rb_define_method(rb_cMolecule, "angle_par", s_Molecule_BondParIsObsolete, 1);
	rb_define_method(rb_cMolecule, "dihedral_par", s_Molecule_BondParIsObsolete, 1);
	rb_define_method(rb_cMolecule, "improper_par", s_Molecule_BondParIsObsolete, 1);
	rb_define_method(rb_cMolecule, "vdw_par", s_Molecule_BondParIsObsolete, 1);

	rb_define_method(rb_cMolecule, "start_step", s_Molecule_StartStep, 0);
	rb_define_method(rb_cMolecule, "start_step=", s_Molecule_SetStartStep, 1);
	rb_define_method(rb_cMolecule, "steps_per_frame", s_Molecule_StepsPerFrame, 0);
	rb_define_method(rb_cMolecule, "steps_per_frame=", s_Molecule_SetStepsPerFrame, 1);
	rb_define_method(rb_cMolecule, "ps_per_step", s_Molecule_PsPerStep, 0);
	rb_define_method(rb_cMolecule, "ps_per_step=", s_Molecule_SetPsPerStep, 1);
	
	rb_define_method(rb_cMolecule, "nresidues", s_Molecule_Nresidues, 0);
	rb_define_method(rb_cMolecule, "nresidues=", s_Molecule_ChangeNresidues, 1);
	rb_define_method(rb_cMolecule, "max_residue_number", s_Molecule_MaxResSeq, -1);
	rb_define_method(rb_cMolecule, "min_residue_number", s_Molecule_MinResSeq, -1);
	rb_define_method(rb_cMolecule, "each_atom", s_Molecule_EachAtom, -1);
	rb_define_method(rb_cMolecule, "cell", s_Molecule_Cell, 0);
	rb_define_method(rb_cMolecule, "cell=", s_Molecule_SetCell, -1);
	rb_define_alias(rb_cMolecule, "set_cell", "cell=");
	rb_define_method(rb_cMolecule, "box", s_Molecule_Box, 0);
	rb_define_method(rb_cMolecule, "box=", s_Molecule_SetBox, 1);
	rb_define_method(rb_cMolecule, "set_box", s_Molecule_SetBox, -2);
	rb_define_method(rb_cMolecule, "cell_transform", s_Molecule_CellTransform, 0);
	rb_define_method(rb_cMolecule, "cell_periodicity", s_Molecule_CellPeriodicity, 0);
	rb_define_method(rb_cMolecule, "cell_periodicity=", s_Molecule_SetCellPeriodicity, 1);
	rb_define_alias(rb_cMolecule, "set_cell_periodicity", "cell_periodicity=");
	rb_define_method(rb_cMolecule, "cell_flexibility", s_Molecule_CellFlexibility, 0);
	rb_define_method(rb_cMolecule, "cell_flexibility=", s_Molecule_SetCellFlexibility, 1);
	rb_define_alias(rb_cMolecule, "set_cell_flexibility", "cell_flexibility=");
	rb_define_method(rb_cMolecule, "symmetry", s_Molecule_Symmetry, 0);
	rb_define_alias(rb_cMolecule, "symmetries", "symmetry");
	rb_define_method(rb_cMolecule, "nsymmetries", s_Molecule_Nsymmetries, 0);
	rb_define_method(rb_cMolecule, "add_symmetry", s_Molecule_AddSymmetry, 1);
	rb_define_method(rb_cMolecule, "remove_symmetry", s_Molecule_RemoveSymmetry, -1);
	rb_define_alias(rb_cMolecule, "remove_symmetries", "remove_symmetry");
	rb_define_method(rb_cMolecule, "extract", s_Molecule_Extract, -1);
    rb_define_method(rb_cMolecule, "add", s_Molecule_Add, 1);
	rb_define_alias(rb_cMolecule, "+", "add");
    rb_define_method(rb_cMolecule, "remove", s_Molecule_Remove, 1);
	rb_define_method(rb_cMolecule, "atom_group", s_Molecule_AtomGroup, -1);
	rb_define_method(rb_cMolecule, "atom_index", s_Molecule_AtomIndex, 1);
	rb_define_method(rb_cMolecule, "create_atom", s_Molecule_CreateAnAtom, -1);
	rb_define_method(rb_cMolecule, "duplicate_atom", s_Molecule_DuplicateAnAtom, -1);
	rb_define_method(rb_cMolecule, "create_bond", s_Molecule_CreateBond, -1);
	rb_define_alias(rb_cMolecule, "create_bonds", "create_bond");
	rb_define_method(rb_cMolecule, "remove_bond", s_Molecule_RemoveBond, -1);
	rb_define_alias(rb_cMolecule, "remove_bonds", "remove_bond");
	rb_define_method(rb_cMolecule, "assign_bond_order", s_Molecule_AssignBondOrder, 2);
	rb_define_alias(rb_cMolecule, "assign_bond_orders", "assign_bond_order");
	rb_define_method(rb_cMolecule, "get_bond_order", s_Molecule_GetBondOrder, 1);
	rb_define_alias(rb_cMolecule, "get_bond_orders", "get_bond_order");
	rb_define_method(rb_cMolecule, "bond_exist?", s_Molecule_BondExist, 2);
	rb_define_method(rb_cMolecule, "add_angle", s_Molecule_AddAngle, 3);
	rb_define_method(rb_cMolecule, "remove_angle", s_Molecule_RemoveAngle, 3);
	rb_define_method(rb_cMolecule, "add_dihedral", s_Molecule_AddDihedral, 4);
	rb_define_method(rb_cMolecule, "remove_dihedral", s_Molecule_RemoveDihedral, 4);
	rb_define_method(rb_cMolecule, "add_improper", s_Molecule_AddImproper, 4);
	rb_define_method(rb_cMolecule, "remove_improper", s_Molecule_RemoveImproper, -1);
	rb_define_method(rb_cMolecule, "assign_residue", s_Molecule_AssignResidue, 2);
	rb_define_method(rb_cMolecule, "offset_residue", s_Molecule_OffsetResidue, 2);
	rb_define_method(rb_cMolecule, "renumber_atoms", s_Molecule_RenumberAtoms, 1);
	rb_define_method(rb_cMolecule, "find_close_atoms", s_Molecule_FindCloseAtoms, -1);
	rb_define_method(rb_cMolecule, "guess_bonds", s_Molecule_GuessBonds, -1);
	rb_define_method(rb_cMolecule, "selection", s_Molecule_Selection, 0);
	rb_define_method(rb_cMolecule, "selection=", s_Molecule_SetSelection, 1);
	rb_define_method(rb_cMolecule, "set_undoable_selection", s_Molecule_SetUndoableSelection, 1);
	rb_define_method(rb_cMolecule, "hidden_atoms", s_Molecule_HiddenAtoms, 0);  /*  obsolete  */
	rb_define_method(rb_cMolecule, "hidden_atoms=", s_Molecule_SetHiddenAtoms, 1);	/*  obsolete  */
	rb_define_alias(rb_cMolecule, "set_hidden_atoms", "hidden_atoms=");  /*  obsolete  */
	rb_define_method(rb_cMolecule, "frame", s_Molecule_Frame, 0);
	rb_define_method(rb_cMolecule, "frame=", s_Molecule_SelectFrame, 1);
	rb_define_alias(rb_cMolecule, "select_frame", "frame=");
	rb_define_method(rb_cMolecule, "nframes", s_Molecule_Nframes, 0);
	rb_define_method(rb_cMolecule, "create_frame", s_Molecule_CreateFrames, -1);
	rb_define_method(rb_cMolecule, "insert_frame", s_Molecule_InsertFrames, -1);
	rb_define_method(rb_cMolecule, "remove_frame", s_Molecule_RemoveFrames, -1);
	rb_define_alias(rb_cMolecule, "create_frames", "create_frame");
	rb_define_alias(rb_cMolecule, "insert_frames", "insert_frame");
	rb_define_alias(rb_cMolecule, "remove_frames", "remove_frame");
	rb_define_method(rb_cMolecule, "each_frame", s_Molecule_EachFrame, 0);
	rb_define_method(rb_cMolecule, "get_coord_from_frame", s_Molecule_GetCoordFromFrame, -1);
	rb_define_method(rb_cMolecule, "register_undo", s_Molecule_RegisterUndo, -1);
	rb_define_method(rb_cMolecule, "undo_enabled?", s_Molecule_UndoEnabled, 0);
	rb_define_method(rb_cMolecule, "undo_enabled=", s_Molecule_SetUndoEnabled, 1);
	rb_define_method(rb_cMolecule, "set_atom_attr", s_Molecule_SetAtomAttr, 3);
	rb_define_method(rb_cMolecule, "get_atom_attr", s_Molecule_GetAtomAttr, 2);
	rb_define_method(rb_cMolecule, "fragment", s_Molecule_Fragment, -1);
	rb_define_method(rb_cMolecule, "fragments", s_Molecule_Fragments, -1);
	rb_define_method(rb_cMolecule, "each_fragment", s_Molecule_EachFragment, -1);
	rb_define_method(rb_cMolecule, "detachable?", s_Molecule_Detachable_P, 1);
	rb_define_method(rb_cMolecule, "bonds_on_border", s_Molecule_BondsOnBorder, -1);
	rb_define_method(rb_cMolecule, "translate", s_Molecule_Translate, -1);
	rb_define_method(rb_cMolecule, "rotate", s_Molecule_Rotate, -1);
	rb_define_method(rb_cMolecule, "reflect", s_Molecule_Reflect, -1);
	rb_define_method(rb_cMolecule, "invert", s_Molecule_Invert, -1);
	rb_define_method(rb_cMolecule, "transform", s_Molecule_Transform, -1);
	rb_define_method(rb_cMolecule, "center_of_mass", s_Molecule_CenterOfMass, -1);
	rb_define_method(rb_cMolecule, "centralize", s_Molecule_Centralize, -1);
	rb_define_method(rb_cMolecule, "bounds", s_Molecule_Bounds, -1);
	rb_define_method(rb_cMolecule, "measure_bond", s_Molecule_MeasureBond, 2);
	rb_define_method(rb_cMolecule, "measure_angle", s_Molecule_MeasureAngle, 3);
	rb_define_method(rb_cMolecule, "measure_dihedral", s_Molecule_MeasureDihedral, 4);
	rb_define_method(rb_cMolecule, "expand_by_symmetry", s_Molecule_ExpandBySymmetry, -1);
	rb_define_method(rb_cMolecule, "amend_by_symmetry", s_Molecule_AmendBySymmetry, -1);
	rb_define_method(rb_cMolecule, "transform_for_symop", s_Molecule_TransformForSymop, -1);
	rb_define_method(rb_cMolecule, "symop_for_transform", s_Molecule_SymopForTransform, -1);
	rb_define_method(rb_cMolecule, "wrap_unit_cell", s_Molecule_WrapUnitCell, 1);
	rb_define_method(rb_cMolecule, "find_conflicts", s_Molecule_FindConflicts, -1);
	rb_define_method(rb_cMolecule, "fit_coordinates", s_Molecule_FitCoordinates, -1);
	rb_define_method(rb_cMolecule, "display", s_Molecule_Display, 0);
	rb_define_method(rb_cMolecule, "make_front", s_Molecule_MakeFront, 0);
	rb_define_method(rb_cMolecule, "update_enabled?", s_Molecule_UpdateEnabled, 0);
	rb_define_method(rb_cMolecule, "update_enabled=", s_Molecule_SetUpdateEnabled, 1);	
	rb_define_method(rb_cMolecule, "show_unitcell", s_Molecule_ShowUnitCell, -1);
	rb_define_method(rb_cMolecule, "show_hydrogens", s_Molecule_ShowHydrogens, -1);
	rb_define_method(rb_cMolecule, "show_hydrogens=", s_Molecule_ShowHydrogens, -1);
	rb_define_method(rb_cMolecule, "show_dummy_atoms", s_Molecule_ShowDummyAtoms, -1);
	rb_define_method(rb_cMolecule, "show_dummy_atoms=", s_Molecule_ShowDummyAtoms, -1);
	rb_define_method(rb_cMolecule, "show_expanded", s_Molecule_ShowExpanded, -1);
	rb_define_method(rb_cMolecule, "show_expanded=", s_Molecule_ShowExpanded, -1);
	rb_define_method(rb_cMolecule, "is_atom_visible", s_Molecule_IsAtomVisible, 1);
	rb_define_method(rb_cMolecule, "show_ellipsoids", s_Molecule_ShowEllipsoids, -1);
	rb_define_method(rb_cMolecule, "show_ellipsoids=", s_Molecule_ShowEllipsoids, -1);
	rb_define_method(rb_cMolecule, "show_graphite", s_Molecule_ShowGraphite, -1);
	rb_define_method(rb_cMolecule, "show_graphite=", s_Molecule_ShowGraphite, -1);
	rb_define_method(rb_cMolecule, "show_graphite?", s_Molecule_ShowGraphiteFlag, 0);
	rb_define_method(rb_cMolecule, "show_periodic_image", s_Molecule_ShowPeriodicImage, -1);
	rb_define_method(rb_cMolecule, "show_periodic_image=", s_Molecule_ShowPeriodicImage, -1);
	rb_define_method(rb_cMolecule, "show_periodic_image?", s_Molecule_ShowPeriodicImageFlag, 0);
	rb_define_alias(rb_cMolecule, "show_unitcell=", "show_unitcell");
	rb_define_alias(rb_cMolecule, "show_unit_cell", "show_unitcell");
	rb_define_alias(rb_cMolecule, "show_unit_cell=", "show_unitcell");
	rb_define_alias(rb_cMolecule, "show_hydrogens=", "show_hydrogens");
	rb_define_alias(rb_cMolecule, "show_dummy_atoms=", "show_dummy_atoms");
	rb_define_alias(rb_cMolecule, "show_expanded=", "show_expanded");
	rb_define_alias(rb_cMolecule, "show_ellipsoids=", "show_ellipsoids");
	rb_define_alias(rb_cMolecule, "show_graphite=", "show_graphite");
	rb_define_alias(rb_cMolecule, "show_periodic_image=", "show_periodic_image");
	rb_define_method(rb_cMolecule, "line_mode", s_Molecule_LineMode, -1);
	rb_define_alias(rb_cMolecule, "line_mode=", "line_mode");
	rb_define_method(rb_cMolecule, "resize_to_fit", s_Molecule_ResizeToFit, 0);
#if 1 || !defined(__CMDMAC__)
	rb_define_method(rb_cMolecule, "get_view_rotation", s_Molecule_GetViewRotation, 0);
	rb_define_method(rb_cMolecule, "get_view_scale", s_Molecule_GetViewScale, 0);
	rb_define_method(rb_cMolecule, "get_view_center", s_Molecule_GetViewCenter, 0);
	rb_define_method(rb_cMolecule, "set_view_rotation", s_Molecule_SetViewRotation, 2);
	rb_define_method(rb_cMolecule, "set_view_scale", s_Molecule_SetViewScale, 1);
	rb_define_method(rb_cMolecule, "set_view_center", s_Molecule_SetViewCenter, 1);
	rb_define_method(rb_cMolecule, "set_background_color", s_Molecule_SetBackgroundColor, -1);
	rb_define_method(rb_cMolecule, "create_graphic", s_Molecule_CreateGraphic, -1);
	rb_define_method(rb_cMolecule, "remove_graphic", s_Molecule_RemoveGraphic, 1);
	rb_define_method(rb_cMolecule, "ngraphics", s_Molecule_NGraphics, 0);
	rb_define_method(rb_cMolecule, "set_graphic_point", s_Molecule_SetGraphicPoint, 3);
	rb_define_method(rb_cMolecule, "set_graphic_color", s_Molecule_SetGraphicColor, 2);
	rb_define_method(rb_cMolecule, "show_graphic", s_Molecule_ShowGraphic, 1);
	rb_define_method(rb_cMolecule, "hide_graphic", s_Molecule_HideGraphic, 1);
#endif
	rb_define_method(rb_cMolecule, "show_text", s_Molecule_ShowText, 1);
	rb_define_method(rb_cMolecule, "md_arena", s_Molecule_MDArena, 0);
	rb_define_method(rb_cMolecule, "set_parameter_attr", s_Molecule_SetParameterAttr, 5);
	rb_define_method(rb_cMolecule, "parameter", s_Molecule_Parameter, 0);
	rb_define_method(rb_cMolecule, "selected_MO", s_Molecule_SelectedMO, 0);
	rb_define_method(rb_cMolecule, "default_MO_grid", s_Molecule_GetDefaultMOGrid, -1);
	rb_define_method(rb_cMolecule, "cubegen", s_Molecule_Cubegen, -1);
	rb_define_method(rb_cMolecule, "nelpots", s_Molecule_NElpots, 0);
	rb_define_method(rb_cMolecule, "elpot", s_Molecule_Elpot, 1);
	rb_define_method(rb_cMolecule, "add_gaussian_orbital_shell", s_Molecule_AddGaussianOrbitalShell, 3);
	rb_define_method(rb_cMolecule, "add_gaussian_primitive_coefficients", s_Molecule_AddGaussianPrimitiveCoefficients, 3);
	rb_define_method(rb_cMolecule, "mo_type", s_Molecule_MOType, 0);
	rb_define_method(rb_cMolecule, "set_mo_coefficients", s_Molecule_SetMOCoefficients, 3);
	rb_define_method(rb_cMolecule, "allocate_basis_set_record", s_Molecule_AllocateBasisSetRecord, 3);
	rb_define_method(rb_cMolecule, "search_equivalent_atoms", s_Molecule_SearchEquivalentAtoms, -1);
	
	rb_define_method(rb_cMolecule, "create_pi_anchor", s_Molecule_CreatePiAnchor, -1);
	rb_define_method(rb_cMolecule, "==", s_Molecule_Equal, 1);
	rb_define_method(rb_cMolecule, "call_subprocess_async", s_Molecule_CallSubProcessAsync, -1);

	rb_define_singleton_method(rb_cMolecule, "current", s_Molecule_Current, 0);
	rb_define_singleton_method(rb_cMolecule, "[]", s_Molecule_MoleculeAtIndex, -1);
	rb_define_singleton_method(rb_cMolecule, "open", s_Molecule_Open, -1);
	rb_define_singleton_method(rb_cMolecule, "list", s_Molecule_List, 0);
	rb_define_singleton_method(rb_cMolecule, "ordered_list", s_Molecule_OrderedList, 0);
	rb_define_singleton_method(rb_cMolecule, "error_message", s_Molecule_ErrorMessage, 0);
	rb_define_singleton_method(rb_cMolecule, "set_error_message", s_Molecule_SetErrorMessage, 1);
	rb_define_singleton_method(rb_cMolecule, "error_message=", s_Molecule_SetErrorMessage, 1);
	
	/*  class MolEnumerable  */
	rb_cMolEnumerable = rb_define_class_under(rb_mMolby, "MolEnumerable", rb_cObject);
    rb_include_module(rb_cMolEnumerable, rb_mEnumerable);
	rb_define_method(rb_cMolEnumerable, "[]", s_MolEnumerable_Aref, 1);
	rb_define_method(rb_cMolEnumerable, "length", s_MolEnumerable_Length, 0);
    rb_define_alias(rb_cMolEnumerable, "size", "length");
	rb_define_method(rb_cMolEnumerable, "each", s_MolEnumerable_Each, 0);
	rb_define_method(rb_cMolEnumerable, "==", s_MolEnumerable_Equal, 1);

	/*  class AtomRef  */
	rb_cAtomRef = rb_define_class_under(rb_mMolby, "AtomRef", rb_cObject);
	for (i = 0; s_AtomAttrDefTable[i].name != NULL; i++) {
		char buf[64];
		snprintf(buf, sizeof buf - 1, "%s", s_AtomAttrDefTable[i].name);
		rb_define_method(rb_cAtomRef, buf, s_AtomAttrDefTable[i].getter, 0);
		s_AtomAttrDefTable[i].id = rb_intern(buf);
		*(s_AtomAttrDefTable[i].symref) = ID2SYM(s_AtomAttrDefTable[i].id);
		strcat(buf, "=");
		rb_define_method(rb_cAtomRef, buf, s_AtomAttrDefTable[i].setter, 1);
	}
	rb_define_method(rb_cAtomRef, "[]=", s_AtomRef_SetAttr, 2);
	rb_define_alias(rb_cAtomRef, "set_attr", "[]=");
	rb_define_method(rb_cAtomRef, "[]", s_AtomRef_GetAttr, 1);
	rb_define_alias(rb_cAtomRef, "get_attr", "[]");
	s_SetAtomAttrString = rb_str_new2("set_atom_attr");
	rb_global_variable(&s_SetAtomAttrString);
	rb_define_method(rb_cAtomRef, "molecule", s_AtomRef_GetMolecule, 0);
	rb_define_method(rb_cAtomRef, "==", s_AtomRef_Equal, 1);

	/*  class Parameter  */
	rb_cParameter = rb_define_class_under(rb_mMolby, "Parameter", rb_cObject);
	rb_define_method(rb_cParameter, "bond", s_Parameter_Bond, 1);
	rb_define_method(rb_cParameter, "angle", s_Parameter_Angle, 1);
	rb_define_method(rb_cParameter, "dihedral", s_Parameter_Dihedral, 1);
	rb_define_method(rb_cParameter, "improper", s_Parameter_Improper, 1);
	rb_define_method(rb_cParameter, "vdw", s_Parameter_Vdw, 1);
	rb_define_method(rb_cParameter, "vdw_pair", s_Parameter_VdwPair, 1);
	rb_define_method(rb_cParameter, "vdw_cutoff", s_Parameter_VdwCutoff, 1);
	rb_define_method(rb_cParameter, "element", s_Parameter_Element, 1);
	rb_define_method(rb_cParameter, "nbonds", s_Parameter_Nbonds, 0);
	rb_define_method(rb_cParameter, "nangles", s_Parameter_Nangles, 0);
	rb_define_method(rb_cParameter, "ndihedrals", s_Parameter_Ndihedrals, 0);
	rb_define_method(rb_cParameter, "nimpropers", s_Parameter_Nimpropers, 0);
	rb_define_method(rb_cParameter, "nvdws", s_Parameter_Nvdws, 0);
	rb_define_method(rb_cParameter, "nvdw_pairs", s_Parameter_NvdwPairs, 0);
	rb_define_method(rb_cParameter, "nvdw_cutoffs", s_Parameter_NvdwCutoffs, 0);
	rb_define_method(rb_cParameter, "nelements", s_Parameter_Nelements, 0);
	rb_define_method(rb_cParameter, "bonds", s_Parameter_Bonds, 0);
	rb_define_method(rb_cParameter, "angles", s_Parameter_Angles, 0);
	rb_define_method(rb_cParameter, "dihedrals", s_Parameter_Dihedrals, 0);
	rb_define_method(rb_cParameter, "impropers", s_Parameter_Impropers, 0);
	rb_define_method(rb_cParameter, "vdws", s_Parameter_Vdws, 0);
	rb_define_method(rb_cParameter, "vdw_pairs", s_Parameter_VdwPairs, 0);
	rb_define_method(rb_cParameter, "vdw_cutoffs", s_Parameter_VdwCutoffs, 0);
	rb_define_method(rb_cParameter, "elements", s_Parameter_Elements, 0);
	rb_define_method(rb_cParameter, "lookup", s_Parameter_LookUp, -1);
	rb_define_method(rb_cParameter, "==", s_Parameter_Equal, 1);
	rb_define_singleton_method(rb_cParameter, "builtin", s_Parameter_Builtin, 0);
	rb_define_singleton_method(rb_cParameter, "bond", s_Parameter_Bond, 1);
	rb_define_singleton_method(rb_cParameter, "angle", s_Parameter_Angle, 1);
	rb_define_singleton_method(rb_cParameter, "dihedral", s_Parameter_Dihedral, 1);
	rb_define_singleton_method(rb_cParameter, "improper", s_Parameter_Improper, 1);
	rb_define_singleton_method(rb_cParameter, "vdw", s_Parameter_Vdw, 1);
	rb_define_singleton_method(rb_cParameter, "vdw_pair", s_Parameter_VdwPair, 1);
	rb_define_singleton_method(rb_cParameter, "vdw_cutoff", s_Parameter_VdwCutoff, 1);
	rb_define_singleton_method(rb_cParameter, "element", s_Parameter_Element, 1);
	rb_define_singleton_method(rb_cParameter, "nbonds", s_Parameter_Nbonds, 0);
	rb_define_singleton_method(rb_cParameter, "nangles", s_Parameter_Nangles, 0);
	rb_define_singleton_method(rb_cParameter, "ndihedrals", s_Parameter_Ndihedrals, 0);
	rb_define_singleton_method(rb_cParameter, "nimpropers", s_Parameter_Nimpropers, 0);
	rb_define_singleton_method(rb_cParameter, "nvdws", s_Parameter_Nvdws, 0);
	rb_define_singleton_method(rb_cParameter, "nvdw_pairs", s_Parameter_NvdwPairs, 0);
	rb_define_singleton_method(rb_cParameter, "nvdw_cutoffs", s_Parameter_NvdwCutoffs, 0);
	rb_define_singleton_method(rb_cParameter, "nelements", s_Parameter_Nelements, 0);
	rb_define_singleton_method(rb_cParameter, "bonds", s_Parameter_Bonds, 0);
	rb_define_singleton_method(rb_cParameter, "angles", s_Parameter_Angles, 0);
	rb_define_singleton_method(rb_cParameter, "dihedrals", s_Parameter_Dihedrals, 0);
	rb_define_singleton_method(rb_cParameter, "impropers", s_Parameter_Impropers, 0);
	rb_define_singleton_method(rb_cParameter, "vdws", s_Parameter_Vdws, 0);
	rb_define_singleton_method(rb_cParameter, "vdw_pairs", s_Parameter_VdwPairs, 0);
	rb_define_singleton_method(rb_cParameter, "vdw_cutoffs", s_Parameter_VdwCutoffs, 0);
	rb_define_singleton_method(rb_cParameter, "elements", s_Parameter_Elements, 0);
	rb_define_singleton_method(rb_cParameter, "lookup", s_Parameter_LookUp, -1);
	rb_define_const(rb_cParameter, "Builtin", Data_Wrap_Struct(rb_cParameter, 0, NULL, NULL));

	/*  class ParEnumerable  */
	rb_cParEnumerable = rb_define_class_under(rb_mMolby, "ParEnumerable", rb_cObject);
    rb_include_module(rb_cParEnumerable, rb_mEnumerable);
	rb_define_method(rb_cParEnumerable, "par_type", s_ParEnumerable_ParType, 0);
	rb_define_method(rb_cParEnumerable, "[]", s_ParEnumerable_Aref, 1);
	rb_define_method(rb_cParEnumerable, "length", s_ParEnumerable_Length, 0);
	rb_define_method(rb_cParEnumerable, "each", s_ParEnumerable_Each, 0);
	rb_define_method(rb_cParEnumerable, "reverse_each", s_ParEnumerable_ReverseEach, 0);
	rb_define_method(rb_cParEnumerable, "insert", s_ParEnumerable_Insert, -1);
	rb_define_method(rb_cParEnumerable, "delete", s_ParEnumerable_Delete, 1);
	rb_define_method(rb_cParEnumerable, "lookup", s_ParEnumerable_LookUp, -1);
	rb_define_method(rb_cParEnumerable, "==", s_ParEnumerable_Equal, 1);
	
	/*  class ParameterRef  */
	rb_cParameterRef = rb_define_class_under(rb_mMolby, "ParameterRef", rb_cObject);
	for (i = 0; s_ParameterAttrDefTable[i].name != NULL; i++) {
		char buf[64];
		snprintf(buf, sizeof buf - 1, "%s", s_ParameterAttrDefTable[i].name);
		rb_define_method(rb_cParameterRef, buf, s_ParameterAttrDefTable[i].getter, 0);
		s_ParameterAttrDefTable[i].id = rb_intern(buf);
		if (s_ParameterAttrDefTable[i].symref != NULL)
			*(s_ParameterAttrDefTable[i].symref) = ID2SYM(s_ParameterAttrDefTable[i].id);
		if (s_ParameterAttrDefTable[i].setter != NULL) {
			strcat(buf, "=");
			rb_define_method(rb_cParameterRef, buf, s_ParameterAttrDefTable[i].setter, 1);
		}
	}
	rb_define_method(rb_cParameterRef, "[]=", s_ParameterRef_SetAttr, 2);
	rb_define_alias(rb_cParameterRef, "set_attr", "[]=");
	rb_define_method(rb_cParameterRef, "[]", s_ParameterRef_GetAttr, 1);
	rb_define_alias(rb_cParameterRef, "get_attr", "[]");
	rb_define_method(rb_cParameterRef, "to_hash", s_ParameterRef_ToHash, 0);
	rb_define_method(rb_cParameterRef, "to_s", s_ParameterRef_ToString, 0);
	rb_define_method(rb_cParameterRef, "keys", s_ParameterRef_Keys, 0);
	rb_define_method(rb_cParameterRef, "==", s_ParameterRef_Equal, 1);

	/*  class MolbyError  */
	rb_eMolbyError = rb_define_class_under(rb_mMolby, "MolbyError", rb_eStandardError);

	/*  module Kernel  */
	rb_define_method(rb_mKernel, "check_interrupt", s_Kernel_CheckInterrupt, 0);
	rb_define_method(rb_mKernel, "get_interrupt_flag", s_GetInterruptFlag, 0);
	rb_define_method(rb_mKernel, "set_interrupt_flag", s_SetInterruptFlag, 1);
	rb_define_method(rb_mKernel, "show_progress_panel", s_ShowProgressPanel, -1);
	rb_define_method(rb_mKernel, "hide_progress_panel", s_HideProgressPanel, 0);
	rb_define_method(rb_mKernel, "set_progress_value", s_SetProgressValue, 1);
	rb_define_method(rb_mKernel, "set_progress_message", s_SetProgressMessage, 1);
	rb_define_method(rb_mKernel, "ask", s_Kernel_Ask, -1);
	rb_define_method(rb_mKernel, "register_menu", s_Kernel_RegisterMenu, 2);
	rb_define_method(rb_mKernel, "lookup_menu", s_Kernel_LookupMenu, 1);
	rb_define_method(rb_mKernel, "get_global_settings", s_Kernel_GetGlobalSettings, 1);
	rb_define_method(rb_mKernel, "set_global_settings", s_Kernel_SetGlobalSettings, 2);
	rb_define_method(rb_mKernel, "execute_script", s_Kernel_ExecuteScript, 1);
	rb_define_method(rb_mKernel, "document_home", s_Kernel_DocumentHome, 0);
	rb_define_method(rb_mKernel, "call_subprocess", s_Kernel_CallSubProcess, -1);
	rb_define_method(rb_mKernel, "backquote", s_Kernel_Backquote, 1);
	rb_define_method(rb_mKernel, "message_box", s_Kernel_MessageBox, -1);
	rb_define_method(rb_mKernel, "error_message_box", s_Kernel_ErrorMessageBox, 1);
	rb_define_method(rb_mKernel, "show_console_window", s_Kernel_ShowConsoleWindow, 0);
	rb_define_method(rb_mKernel, "hide_console_window", s_Kernel_HideConsoleWindow, 0);
	rb_define_method(rb_mKernel, "bell", s_Kernel_Bell, 0);
	rb_define_method(rb_mKernel, "play_sound", s_Kernel_PlaySound, -1);
	rb_define_method(rb_mKernel, "stop_sound", s_Kernel_StopSound, 0);
	rb_define_method(rb_mKernel, "export_to_clipboard", s_Kernel_ExportToClipboard, 1);

	s_ID_equal = rb_intern("==");
	g_RubyID_call = rb_intern("call");
}

#pragma mark ====== External functions ======

static VALUE s_ruby_top_self = Qfalse;

static VALUE
s_evalRubyScriptOnMoleculeSub(VALUE val)
{
	void **ptr = (void **)val;
	Molecule *mol = (Molecule *)ptr[1];
	VALUE sval = rb_str_new2((char *)ptr[0]);
	VALUE fnval;
	if (s_ruby_top_self == Qfalse) {
		s_ruby_top_self = rb_eval_string("eval(\"self\",TOPLEVEL_BINDING)");
	}
	if (ptr[2] == NULL) {
		fnval = Qnil;
	} else {
		fnval = Ruby_NewFileStringValue((char *)ptr[2]);
	}
	if (mol == NULL) {
		if (fnval == Qnil)
			return rb_funcall(s_ruby_top_self, rb_intern("eval"), 1, sval);
		else
			return rb_funcall(s_ruby_top_self, rb_intern("eval"), 4, sval, Qnil, fnval, INT2FIX(1));
	} else {
		VALUE mval = ValueFromMolecule(mol);
		if (fnval == Qnil)
			return rb_funcall(mval, rb_intern("instance_eval"), 1, sval);
		else return rb_funcall(mval, rb_intern("instance_eval"), 3, sval, fnval, INT2FIX(1));
	}
}

RubyValue
Molby_evalRubyScriptOnMolecule(const char *script, Molecule *mol, const char *fname, int *status)
{
	RubyValue retval;
	void *args[3];
	VALUE save_interrupt_flag;
	char *save_ruby_sourcefile;
	int save_ruby_sourceline;
	if (gMolbyIsCheckingInterrupt) {
		MolActionAlertRubyIsRunning();
		*status = -1;
		return (RubyValue)Qnil;
	}
	gMolbyRunLevel++;
	args[0] = (void *)script;
	args[1] = (void *)mol;
	args[2] = (void *)fname;
	save_interrupt_flag = s_SetInterruptFlag(Qnil, Qtrue);
	save_ruby_sourcefile = ruby_sourcefile;
	save_ruby_sourceline = ruby_sourceline;
	retval = (RubyValue)rb_protect(s_evalRubyScriptOnMoleculeSub, (VALUE)args, status);
	if (*status != 0) {
		/*  Is this 'exit' exception?  */
		VALUE last_exception = rb_gv_get("$!");
		if (rb_obj_is_kind_of(last_exception, rb_eSystemExit)) {
			/*  Capture exit and return the status value  */
			retval = (RubyValue)rb_funcall(last_exception, rb_intern("status"), 0);
			*status = 0;
		}
	}
	s_SetInterruptFlag(Qnil, save_interrupt_flag);
	ruby_sourcefile = save_ruby_sourcefile;
	ruby_sourceline = save_ruby_sourceline;
	gMolbyRunLevel--;
	return retval;
}

void
Molby_showRubyValue(RubyValue value, char **outValueString)
{
	VALUE val = (VALUE)value;
	if (gMolbyIsCheckingInterrupt) {
		MolActionAlertRubyIsRunning();
		return;
	}
	if (val != Qnil) {
		int status;
		char *str;
		gMolbyRunLevel++;
		val = rb_protect(rb_inspect, val, &status);
		gMolbyRunLevel--;
		str = StringValuePtr(val);
		if (outValueString != NULL)
			*outValueString = strdup(str);
		MyAppCallback_showScriptMessage("%s", str);
	}
}

void
Molby_showError(int status)
{
	static const int tag_raise = 6;
	char *msg = NULL, *msg2;
	VALUE val, backtrace;
	int interrupted = 0;
	if (status == tag_raise) {
		VALUE eclass = CLASS_OF(ruby_errinfo);
		if (eclass == rb_eInterrupt) {
			msg = "Interrupt";
			interrupted = 1;
		}
	}
	gMolbyRunLevel++;
	backtrace = rb_eval_string_protect("$backtrace = $!.backtrace.join(\"\\n\")", &status);
	if (msg == NULL) {
		val = rb_eval_string_protect("$!.to_s", &status);
		if (status == 0)
			msg = RSTRING_PTR(val);
		else msg = "(message not available)";
	}
	asprintf(&msg2, "%s\n%s", msg, RSTRING_PTR(backtrace));
	MyAppCallback_messageBox(msg2, (interrupted == 0 ? "Molby script error" : "Molby script interrupted"), 0, 3);
	free(msg2);
	gMolbyRunLevel--;
}

char *
Molby_getDescription(void)
{
	extern const char *gVersionString, *gCopyrightString;
	extern int gRevisionNumber;
	extern char *gLastBuildString;
	char *s;
	char *revisionString;
	if (gRevisionNumber > 0) {
		asprintf(&revisionString, ", revision %d", gRevisionNumber);
	} else revisionString = "";
	asprintf(&s, 
			 "Molby %s%s\n%s\nLast compile: %s\n"
#if !defined(__CMDMAC__)
			 "\nIncluding:\n"
			 "%s"
#else
			 "Including "
#endif
			 "ruby %s, http://www.ruby-lang.org/\n"
			 "%s\n"
			 "FFTW 3.3.2, http://www.fftw.org/\n"
			 "  Copyright (C) 2003, 2007-11 Matteo Frigo\n"
			 "  and Massachusetts Institute of Technology",
			 gVersionString, revisionString, gCopyrightString, gLastBuildString,
#if !defined(__CMDMAC__)
			 MyAppCallback_getGUIDescriptionString(),
#endif
			 gRubyVersion, gRubyCopyright);
	if (revisionString[0] != 0)
		free(revisionString);
	return s;
}

void
Molby_startup(const char *script, const char *dir)
{
	VALUE val;
	int status;
	char *libpath;
	char *respath, *p, *wbuf;

	/*  Get version/copyright string from Ruby interpreter  */
	{
		gRubyVersion = strdup(ruby_version);
		asprintf(&gRubyCopyright, "%sCopyright (C) %d-%d %s",
#if defined(__CMDMAC__)
				 "",
#else
				 "  ",  /*  Indent for displaying in About dialog  */
#endif
				 RUBY_BIRTH_YEAR, RUBY_RELEASE_YEAR, RUBY_AUTHOR);
	}
	
	/*  Read build and revision information for Molby  */
/*	{
		char buf[200];
		extern int gRevisionNumber;
		extern char *gLastBuildString;
		FILE *fp = fopen("../buildInfo.txt", "r");
		gLastBuildString = "";
		if (fp != NULL) {
			if (fgets(buf, sizeof(buf), fp) != NULL) {
				char *p1 = strchr(buf, '\"');
				char *p2 = strrchr(buf, '\"');
				if (p1 != NULL && p2 != NULL && p2 - p1 > 1) {
					memmove(buf, p1 + 1, p2 - p1 - 1);
					buf[p2 - p1 - 1] = 0;
					asprintf(&gLastBuildString, "Last compile: %s\n", buf);
				}
			}
			fclose(fp);
		}
		fp = fopen("../revisionInfo.txt", "r");
		gRevisionNumber = 0;
		if (fp != NULL) {
			if (fgets(buf, sizeof(buf), fp) != NULL) {
				gRevisionNumber = strtol(buf, NULL, 0);
			}
			fclose(fp);
		}
    } */

#if defined(__CMDMAC__)
	wbuf = Molby_getDescription();
	printf("%s\n", wbuf);
	free(wbuf);
#endif
	
	/*  Read atom display parameters  */
	if (ElementParameterInitialize("element.par", &wbuf) != 0) {
#if defined(__CMDMAC__)
		fprintf(stderr, "%s\n", wbuf);
#else
		MyAppCallback_setConsoleColor(1);
		MyAppCallback_showScriptMessage("%s", wbuf);
		MyAppCallback_setConsoleColor(0);
#endif
		free(wbuf);
	}
	
	/*  Read default parameters  */
	ParameterReadFromFile(gBuiltinParameters, "default.par", &wbuf, NULL);
	if (wbuf != NULL) {
#if defined(__CMDMAC__)
		fprintf(stderr, "%s\n", wbuf);
#else
		MyAppCallback_setConsoleColor(1);
		MyAppCallback_showScriptMessage("%s", wbuf);
		MyAppCallback_setConsoleColor(0);
#endif
		free(wbuf);
	}
		
	/*  Initialize Ruby interpreter  */
	ruby_init();
	
	/*  Initialize loadpath; the specified directory, "lib" subdirectory, and "."  */
	ruby_incpush(".");
	asprintf(&libpath, "%s%clib", dir, PATH_SEPARATOR);
	ruby_incpush(libpath);
	free(libpath);
	ruby_incpush(dir);

	ruby_script("Molby");
	
	/*  Find the resource path (the parent directory of the given directory)  */
	respath = strdup(dir);
	p = strrchr(respath, '/');
	if (p == NULL && PATH_SEPARATOR != '/')
		p = strrchr(respath, PATH_SEPARATOR);
	if (p != NULL)
		*p = 0;
	val = Ruby_NewFileStringValue(respath);
	rb_define_global_const("MolbyResourcePath", val);
	free(respath);

	/*  Define Molby classes  */
	Init_Molby();
	RubyDialogInitClass();

	rb_define_const(rb_mMolby, "ResourcePath", val);
	val = Ruby_NewFileStringValue(dir);
	rb_define_const(rb_mMolby, "ScriptPath", val);
	asprintf(&p, "%s%c%s", dir, PATH_SEPARATOR, "mbsf");
	val = Ruby_NewFileStringValue(p);
	rb_define_const(rb_mMolby, "MbsfPath", val);	
	free(p);
	
#if defined(__CMDMAC__)
	rb_define_const(rb_mMolby, "HasGUI", Qfalse);
#else
	rb_define_const(rb_mMolby, "HasGUI", Qtrue);
#endif

#if !__CMDMAC__
	
	/*  Create objects for stdout and stderr  */
	val = rb_funcall(rb_cObject, rb_intern("new"), 0);
	rb_define_singleton_method(val, "write", s_StandardOutput, 1);
	rb_gv_set("$stdout", val);
	val = rb_funcall(rb_cObject, rb_intern("new"), 0);
	rb_define_singleton_method(val, "write", s_StandardErrorOutput, 1);
	rb_gv_set("$stderr", val);

	/*  Create objects for stdin  */
	val = rb_funcall(rb_cObject, rb_intern("new"), 0);
	rb_define_singleton_method(val, "gets", s_StandardInputGets, -1);
	rb_define_singleton_method(val, "readline", s_StandardInputGets, -1);
	rb_define_singleton_method(val, "method_missing", s_StandardInputMethodMissing, -1);
	rb_gv_set("$stdin", val);
	
#endif
	
	/*  Global variable to hold backtrace  */
	rb_define_variable("$backtrace", &gMolbyBacktrace);

#if !__CMDMAC__
	/*  Register interrupt check code  */
	rb_add_event_hook(s_Event_Callback, RUBY_EVENT_ALL);
#endif
	
#if !__CMDMAC__
	/*  Start interval timer (for periodic polling of interrupt); firing every 50 msec  */
	s_SetIntervalTimer(0, 50);
#endif
	
	/*  Read the startup script  */
	if (script != NULL && script[0] != 0) {
		MyAppCallback_showScriptMessage("Evaluating %s...\n", script);
		gMolbyRunLevel++;
		rb_load_protect(rb_str_new2(script), 0, &status);
		gMolbyRunLevel--;
		if (status != 0)
			Molby_showError(status);
		else
			MyAppCallback_showScriptMessage("Done.\n");
	}
}

void
Molby_buildARGV(int argc, const char **argv)
{
	int i;
    rb_ary_clear(rb_argv);
    for (i = 0; i < argc; i++) {
		VALUE arg = rb_tainted_str_new2(argv[i]);
		OBJ_FREEZE(arg);
		rb_ary_push(rb_argv, arg);
    }
}
