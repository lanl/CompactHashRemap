/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */

#ifndef H_REMAP_GPU_H
#define H_REMAP_GPU_H

#include "HashFactory/HashFactory.h"

void setup_cl ();
void cleanup_cl();

double cl_full_perfect_remap (cell_list icells, cell_list ocells, int run_tests);
double cl_singlewrite_remap (cell_list icells, cell_list ocells, int run_tests);
double cl_compact_singlewrite_remap (cell_list icells, cell_list ocells, int run_tests);
double cl_hierarchical_remap (cell_list icells, cell_list ocells, int run_tests);
double cl_compact_hierarchical_remap (cell_list icells, cell_list ocells,
     intintHash_Factory *CLFactory, int run_tests);

#endif
