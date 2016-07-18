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

#include "meshgen/meshgen.h"
#include "genmalloc/genmalloc.h"
#include "full_perfect_remap.h"
#include "stdio.h"
#include <assert.h>


double avg_sub_cells (cell_list icells, uint ji, uint ii, uint level, uint *hash) {

    uint key, i_max, jump;
    double sum = 0.0;
    i_max = icells.ibasesize*two_to_the(icells.levmax);
    jump = two_to_the(icells.levmax - level - 1);
    
    for (uint j = 0; j < 2; j++) {
        for (uint i = 0; i < 2; i++) {
            key = ((ji + (j*jump)) * i_max) + (ii + (i*jump));
            uint probe = hash[key];
            if (icells.level[probe] == (level + 1)) {
                sum += icells.values[probe];
            } else {
                sum += avg_sub_cells(icells, ji + (j*jump), ii + (i*jump), level + 1, hash);
            }
        }
    }
    
    return sum/4.0;
}

void full_perfect_remap (cell_list icells, cell_list ocells) {

    // Allocate a hash table the size of the finest level of the grid
    size_t hash_size = icells.ibasesize*icells.ibasesize*four_to_the(icells.levmax);

    uint *hash = (uint *) malloc(hash_size * sizeof(uint));
    uint lev_mod;
    // levmax+1?
    uint i_max = icells.ibasesize*two_to_the(icells.levmax);
    
    //for (uint i = 0; i < icells.ncells; i++) {
    //    printf("%u\t%u\t%u\n", icells.i[i], icells.j[i], icells.level[i]);
    //} printf("\n%u\t%u\n", icells.ibasesize, icells.levmax);

    // Fill Hash Table from Input mesh
    for (uint ic = 0; ic < icells.ncells; ic++){
        uint lev = icells.level[ic];
        uint i = icells.i[ic];
        uint j = icells.j[ic];
        // If at the maximum level just set the one cell
        if (lev == icells.levmax) {
            hash[(j*i_max)+i] = ic;
        } else {
            // Set the square block of cells at the finest level
            // to the index number
            uint lev_mod = two_to_the(icells.levmax - lev);
            for (uint jj = j*lev_mod; jj < (j+1)*lev_mod; jj++) {
                for (uint ii = i*lev_mod; ii < (i+1)*lev_mod; ii++) {
                    hash[(jj*i_max)+ii] = ic;
                }
            }
        }
    }

    // Use Hash Table to Perform Remap
    for (uint ic = 0; ic < ocells.ncells; ic++){
        uint ii, jj;
        uint i = ocells.i[ic];
        uint j = ocells.j[ic];
        uint lev = ocells.level[ic];

        if (lev < ocells.levmax) {
            lev_mod = two_to_the(ocells.levmax - lev);
            ii = i*lev_mod;
            jj = j*lev_mod;
        } else {
            lev_mod = two_to_the(lev - ocells.levmax);
            ii = i/lev_mod;
            jj = j/lev_mod;
        }

        uint key = hash[(jj*i_max)+ii];
        
        if (lev >= icells.level[key]) {
            ocells.values[ic] = icells.values[key];
        } else {
            ocells.values[ic] = avg_sub_cells(icells, jj, ii, lev, hash);
        }
    }

    // Deallocate hash table
    free(hash);
}

#ifdef _OPENMP
void full_perfect_remap_openMP (cell_list icells, cell_list ocells) {

    // Allocate a hash table the size of the finest level of the grid
    uint i_max = icells.ibasesize*two_to_the(icells.levmax);
    uint j_max = icells.ibasesize*two_to_the(icells.levmax);
    uint *hash = (uint *)malloc(i_max*j_max*sizeof(uint));

#pragma omp parallel default(none) shared(icells, ocells, hash, i_max)
    {
        uint ilength = icells.ncells;
        uint olength = ocells.ncells;
        uint max_lev = icells.levmax;
        uint lev_mod;

        // Fill Hash Table from Input mesh
#pragma omp for
        for (uint ic = 0; ic < ilength; ic++){
            uint lev = icells.level[ic];
            uint i = icells.i[ic];
            uint j = icells.j[ic];
            // If at the maximum level just set the one cell
            if (lev == max_lev) {
                //printf("%u\t%u\n", i, j);
                hash[j*i_max + i] = ic;
            } else {
                // Set the square block of cells at the finest level
                // to the index number
                lev_mod = two_to_the(max_lev - lev);
                //printf("%u\t%u\n", i*lev_mod, j*lev_mod);
                for (uint jj = j*lev_mod; jj < (j+1)*lev_mod; jj++) {
                    for (uint ii = i*lev_mod; ii < (i+1)*lev_mod; ii++) {
                        hash[jj*i_max + ii] = ic;
                    }
                }
            }
        }

    // Use Hash Table to Perform Remap
#pragma omp for
        for (uint ic = 0; ic < olength; ic++){
            uint lev = ocells.level[ic];
            uint i = ocells.i[ic];
            uint j = ocells.j[ic];
            uint ii, jj;
            
            if (lev < ocells.levmax) {
                lev_mod = two_to_the(ocells.levmax - lev);
                ii = i*lev_mod;
                jj = j*lev_mod;
            } else {
                lev_mod = two_to_the(lev - ocells.levmax);
                ii = i/lev_mod;
                jj = j/lev_mod;
            }
            
            uint key = hash[(jj*i_max)+ii];
        
            if (lev >= icells.level[key]) {
                ocells.values[ic] = icells.values[key];
            } else {
                ocells.values[ic] = avg_sub_cells(icells, jj, ii, lev, hash);
            }
            
        }
    }

    // Deallocate hash table
    free(hash);
}
#endif
