/*
 *  MyThread.h
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

#ifndef __MyThread_h__
#define __MyThread_h__

#include "wx/thread.h"

class MyThread: public wxThread
{
public:
	static wxThreadError DetachNewThread(int (*entry_func)(void *, int), int (*exit_func)(void *, int), void *argptr, int argnum);
	virtual ExitCode Entry();
	int (*m_entry_func)(void *, int);
	int (*m_exit_func)(void *, int);
	void *m_argptr;
	int m_argnum;
};

#endif /* __MyThread_h__ */
