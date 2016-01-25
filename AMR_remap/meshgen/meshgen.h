/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */

#ifndef MESHGEN_H
#define MESHGEN_H

#include <stdlib.h>

typedef unsigned int uint;

#define USE_MACROS
#define USE_ASSERT

#ifdef USE_MACROS
#define two_to_the(ishift)       (1<<(ishift) )
#define four_to_the(ishift)      (1 << ( (ishift)*2 ) )
#define key_to_i(key, lev)       ( (key) % two_to_the(lev) )
#define key_to_j(key, lev)       ( (key) / two_to_the(lev) )
#define truncate_base(val, val2) ( ((val)/(val2)) +val2 );
#endif

typedef struct {
    uint ncells;    // number of cells in the mesh
    uint ibasesize; // number of coarse cells across the x dimension for the minimum level of the mesh
    uint jbasesize; // number of coarse cells across the y dimension for the minimum level of the mesh
    uint levmax;    // number of refinement levels in addition to the base mesh
    uint *i;
    uint *j;
    uint *level;
    double *values;
} cell_list;

cell_list new_cell_list(uint *x, uint *y, uint *lev, double *val);
cell_list create_cell_list(cell_list a, int length);
void destroy(cell_list a);

cell_list mesh_maker_level (cell_list clist, uint levels_diff, uint *length,
    uint *max_level, uint *min_level);
cell_list mesh_maker_sparsity (cell_list clist, uint levels_diff, uint *length,
    uint *max_level, uint *min_level, double sparsity);
cell_list adaptiveMeshConstructorWij(cell_list icells, const int n, const int levmax, float threshold,
    int target_ncells);
void divide_cell (uint super_i, uint super_j, uint super_level, cell_list cells, 
    uint cell_count, uint cell_id);
void print_cell_list (cell_list cells, uint length);

uint translate_cell (uint i, uint j, uint lev, uint new_lev);
uint translate_cell (uint i, uint j, uint lev, uint new_lev, int ibasesize);

#ifndef USE_MACROS
int two_to_the (int val);
uint four_to_the (int val);
uint truncate_base (uint val, uint val2);
uint key_to_i (uint key, uint lev);
uint key_to_j (uint key, uint lev);
#endif

#endif
