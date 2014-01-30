/*
 *  MDForce.h
 *
 *  Created by Toshi Nagata on 2005/06/07.
 *  Copyright 2005 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#ifndef __MDFORCE_H__
#define __MDFORCE_H__

#include "MDCore.h"

#ifdef __cplusplus
extern "C" {
#endif
	
void calc_force(MDArena *arena);
void calc_pair_interaction(MDArena *arena, const MDGroupFlags *group_flags_1, const MDGroupFlags *group_flags_2);
void calc_surface_force(MDArena *arena);
	
#ifdef __cplusplus
}
#endif
		
#endif /* __MD_FORCE_H__ */
