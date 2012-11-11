/*
 *  MDEwald.h
 *  Molby
 *
 *  Created by Toshi Nagata on 2012/11/07.
 *  Copyright 2012 Toshi Nagata. All rights reserved.
 *
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation version 2 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 */

#ifndef __MDEWALD_H__
#define __MDEWALD_H__

#include "MDCore.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_DIM_SPLINE  10

typedef struct MDPME MDPME;

void calc_ewald_force(MDArena *arena);
void pme_release(MDArena *arena);
void pme_init(MDArena *arena);

#ifdef __cplusplus
}
#endif
		
#endif /* __MDEWALD_H__ */
