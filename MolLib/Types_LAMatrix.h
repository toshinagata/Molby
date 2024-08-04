/*
 *  Types_LAMatrix.h
 *
 *  Created by Toshi Nagata on 2024/07/31.
 *  Copyright 2024-2024 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#ifndef __Types_LAMatrix_h__
#define __Types_LAMatrix_h__

/*  CBLAS/CLAPACK headers  */
#if defined(__WXMAC__) || defined(__CMDMAC__)
/*  On Mac OS X, CLAPACK is in Accelerate.framework  */
//#include <cblas.h>
//#include <clapack.h>
#include <Accelerate/Accelerate.h>
#else
#include <f2c.h>
#include <blaswrap.h>
#include <clapack.h>
#endif

/*  Get the eigenvalue/eigenvector for a real symmetric matrix (3x3)  */
#if !defined(__WXMAC__) && !defined(__CMDMAC__)
typedef integer        __CLPK_integer;
typedef logical        __CLPK_logical;
typedef real           __CLPK_real;
typedef doublereal     __CLPK_doublereal;
#endif

/*  Wrapper struct for CLAPACK routines  */
struct LAMatrix {
  __CLPK_integer row, column;
  __CLPK_doublereal data[1];
};

#endif /* __Types_LAMatrix_h__ */
