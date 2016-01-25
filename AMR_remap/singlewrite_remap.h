/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */


#ifndef SINGLEWRITE_REMAP_H
#define SINGLEWRITE_REMAP_H

#include "meshgen/meshgen.h"

void singlewrite_remap (cell_list icells, cell_list ocells);
void singlewrite_remap_openMP (cell_list icells, cell_list ocells);
void singlewrite_remap_compact (cell_list icells, cell_list ocells);
void singlewrite_remap_compact_openMP (cell_list icells, cell_list ocells);

#endif
