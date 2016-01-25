/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */

#ifndef HIERARCHICAL_REMAP_H
#define HIERARCHICAL_REMAP_H

#include "meshgen/meshgen.h"
#include "HashFactory/HashFactory.h"

void h_remap (cell_list icells, cell_list ocells);
void h_remap_compact (cell_list icells, cell_list ocells, intintHash_Factory *factory);
void h_remap_openMP (cell_list icells, cell_list ocells);
void h_remap_compact_openMP (cell_list icells, cell_list ocells, intintHash_Factory *factory);

#endif
