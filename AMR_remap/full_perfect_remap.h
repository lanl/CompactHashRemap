/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */


#ifndef FULL_PERFECT_REMAP_H
#define FULL_PERFECT_REMAP_H

#include "meshgen/meshgen.h"

void full_perfect_remap (cell_list icells, cell_list ocells);
#ifdef _OPENMP
void full_perfect_remap_openMP (cell_list icells, cell_list ocells);
#endif

#endif
