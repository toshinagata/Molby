/*
 *  MyMBConv.h
 *  Molby
 *
 *  Created by Toshi Nagata on 11/09/30.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __MyMBConv_h__
#define __MyMBConv_h__

#if defined(__WXMAC__)
#define WX_DEFAULT_CONV wxConvUTF8
#elif defined(__WXMSW__)
#define WX_DEFAULT_CONV *wxConvCurrent
#endif

#endif  /* __MyMBConv_h__ */
