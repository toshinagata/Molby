/*
 *  modalwindow_osx.mm
 *  Molby
 *
 *  Created by Toshi Nagata on 14/09/19.
 *  Copyright 2014 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#import <Cocoa/Cocoa.h>

void
MacRunModalForWindow(void *w)
{
	[NSApp runModalForWindow:(NSWindow *)w];
}

void
MacStopModal(void)
{
	[NSApp stopModal];
}

void *
MacGetActiveWindow(void)
{
	return [NSApp mainWindow];
}

/*
int
MacCheckEscapeKeyPressed(void)
{
	NSEvent *ep = [NSApp nextEventMatchingMask:NSKeyDownMask untilDate:nil inMode:NSEventTrackingRunLoopMode dequeue:YES];
	if (ep != nil && [[ep charactersIgnoringModifiers] characterAtIndex:0] == 0x1e) {
		return 1;
	} else return 0;
}
*/
