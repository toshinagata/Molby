/*
 *  mollib.h
 *
 *  Created by Toshi Nagata on 2005/06/03.
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

#ifndef __MolLib_h__
#define __MolLib_h__

#include "Types.h"
#include "Parameter.h"
#include "Molecule.h"
#include "MainView.h"
#include "Trackball.h"
#include "MolAction.h"

#ifdef __cplusplus
extern "C" {
#endif
	
typedef struct MolArena {
	Molecule *mp;       /*  Current molecule  */
	struct MainView *mview;    /*  Active view  */
	Int npars;          /*  Number of parameters in pars[] */
	Parameter **pars;   /*  Parameters. pars[npars - 1] is always the 'custom' parameters, that are defined individually from TCL command. Other parameters are those from files, and they are registered by the order they are read from files.  */
	
	/*  These two fields are used by T_MolNew, T_MolRelease  */
	Int nmollist;
	Molecule **mollist;

} MolArena;

#ifdef __cplusplus
}
#endif
		
#endif /* __MolLib_h__ */
