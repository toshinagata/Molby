/*
 *  MyThread.cpp
 *  Molby
 *
 *  Created by Toshi Nagata on 09/10/21.
 *  Copyright 2009 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#include "MyThread.h"

wxThreadError
MyThread::DetachNewThread(int (*entry_func)(void *, int), int (*exit_func)(void *, int), void *argptr, int argnum)
{
	wxThreadError err;
	MyThread *thread = new MyThread();
	if (thread == NULL)
		return wxTHREAD_NO_RESOURCE;
	thread->m_entry_func = entry_func;
	thread->m_exit_func = exit_func;
	thread->m_argptr = argptr;
	thread->m_argnum = argnum;
	if ((err = thread->Create()) != wxTHREAD_NO_ERROR)
		return err;
	return thread->Run();
}

MyThread::ExitCode
MyThread::Entry()
{
	ExitCode code = (ExitCode)((*m_entry_func)(m_argptr, m_argnum));
	if (m_exit_func)
		(*m_exit_func)(m_argptr, m_argnum);
	Exit(code);
	return 0;
}
