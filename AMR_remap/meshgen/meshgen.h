/* Copyright 2015-16.  Los Alamos National Security, LLC. This material was produced
 * under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
 * Laboratory (LANL), which is operated by Los Alamos National Security, LLC
 * for the U.S. Department of Energy. The U.S. Government has rights to use,
 * reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS
 * ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
 * to produce derivative works, such modified software should be clearly marked,
 * so as not to confuse it with the version available from LANL.   
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy
 * of the License at 
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed
 * under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
 * CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License.
 *
 * Under this license, it is required to include a reference to this work.
 *
 * This is LANL Copyright Disclosure C16017/LA-CC-15-102
 *
 * Authors: Bob Robey         XCP-2   brobey@lanl.gov
 *          Gerald Collom     XCP-2   gcollom@lanl.gov
 *          Colin Redman      XCP-2   credman@lanl.gov 
 */

#ifndef MESHGEN_H
#define MESHGEN_H

#include <stdlib.h>

typedef unsigned int uint;

#define USE_MACROS
#define USE_ASSERT

#ifdef USE_MACROS
#define two_to_the(ishift)       (1u <<(ishift) )
#define four_to_the(ishift)      (1u << ( (ishift)*2 ) )
#define key_to_i(key, lev)       ( (key) % two_to_the(lev) )
#define key_to_j(key, lev)       ( (key) / two_to_the(lev) )
#define truncate_base(val, val2) ( ((val)/(val2)) +val2 );
#endif

typedef struct {
    uint ncells;    // number of cells in the mesh
    uint ibasesize; // number of coarse cells across the x dimension for the minimum level of the mesh
    uint levmax;    // number of refinement levels in addition to the base mesh
    uint *dist;     // distribution of cells across levels of refinemnt
    uint *i;
    uint *j;
    uint *level;
    double *values;
} cell_list;

cell_list new_cell_list(uint *x, uint *y, uint *lev, double *val);
cell_list create_cell_list(cell_list a, int length);
void destroy(cell_list a);

cell_list mesh_maker (cell_list clist, uint levels_diff, uint *length,
    uint *max_level, double sparsity, uint min_base_size);
cell_list adaptiveMeshConstructorWij(cell_list icells, const uint n, const uint levmax, float threshold,
    uint target_ncells);
void divide_cell (uint super_i, uint super_j, uint super_level, cell_list cells, 
    uint cell_count, uint cell_id);
void print_cell_list (cell_list cells, uint length);
cell_list shuffle_cell_list(cell_list clist, uint num);

#ifndef USE_MACROS
int two_to_the (int val);
uint four_to_the (int val);
uint truncate_base (uint val, uint val2);
uint key_to_i (uint key, uint lev);
uint key_to_j (uint key, uint lev);
#endif

#endif
