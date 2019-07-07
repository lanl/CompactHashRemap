/* Copyright 2015-19.  Triad National Security, LLC. This material was produced
 * under U.S. Government contract 89233218CNA000001 for Los Alamos National 
 * Laboratory (LANL), which is operated by Triad National Security, LLC
 * for the U.S. Department of Energy. The U.S. Government has rights to use,
 * reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
 * TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
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

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "meshgen/meshgen.h"
#include "simplehash/simplehash.h"
#include "singlewrite_remap.h"

// These subroutines are private to this file
double avg_sub_cells (cell_list icells, uint jo, uint io, uint lev, int *hash);
double avg_sub_cells_compact (cell_list icells, uint ji, uint ii, uint level, int *hash);
double avg_sub_cells_compact_openMP (cell_list icells, uint ji, uint ii, uint level, int *hash, uint max_lev);

double avg_sub_cells (cell_list icells, uint ji, uint ii, uint level, int *hash) {

    uint key, i_max, jump;
    double sum = 0.0;
    i_max = icells.ibasesize*two_to_the(icells.levmax);
    jump = two_to_the(icells.levmax - level - 1);
    
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            key = ((ji + (j*jump)) * i_max) + (ii + (i*jump));
            int ic = hash[key];
            // Getting sub averages failed
            assert(ic >= 0);
            if (icells.level[ic] == (level + 1)) {
                sum += icells.values[ic];
            } else {
                sum += avg_sub_cells(icells, ji + (j*jump), ii + (i*jump), level + 1, hash);
            }
        }
    }
    
    return sum/4.0;
}

double avg_sub_cells_compact (cell_list icells, uint ji, uint ii, uint level, int *hash) {

    uint key, i_max, jump;
    double sum = 0.0;
    i_max = icells.ibasesize*two_to_the(icells.levmax);
    jump = two_to_the(icells.levmax - level - 1);
    
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            key = ((ji + (j*jump)) * i_max) + (ii + (i*jump));
            int ic = read_hash(key, hash);
            // Getting sub averages failed
            assert(ic >= 0);
            if (icells.level[ic] == (level + 1)) {
                sum += icells.values[ic];
            } else {
                sum += avg_sub_cells_compact(icells, ji + (j*jump), ii + (i*jump), level + 1, hash);
            }
        }
    }
    
    return sum/4.0;
}

double avg_sub_cells_compact_openMP (cell_list icells, uint ji, uint ii, uint level, int *hash, uint max_lev) {

    uint key, i_max, jump;
    double sum = 0;
    i_max = two_to_the(max_lev);
    jump = two_to_the(max_lev - level - 1);
    
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
            key = ((ji + (j*jump)) * i_max) + (ii + (i*jump));
            int ic = read_hash(key, hash);
            // Getting sub averages failed
            assert(ic >= 0);
            if (icells.level[ic] == (level + 1)) {
                sum += icells.values[ic];
            } else {
                sum += avg_sub_cells_compact_openMP(icells, ji + (j*jump), ii + (i*jump), level + 1, 
                    hash, max_lev);
            }
        }
    }
    
    return sum/4;
}


void singlewrite_remap (cell_list icells, cell_list ocells) {
    
    size_t hash_size = icells.ibasesize*two_to_the(icells.levmax)*
                       icells.ibasesize*two_to_the(icells.levmax);
    int *hash = (int *) malloc(hash_size * sizeof(int));
    uint i_max = icells.ibasesize*two_to_the(icells.levmax);
    
    memset(hash, 0xFFFFFFFF, hash_size*sizeof(uint));
    
    for (uint i = 0; i < icells.ncells; i++) {
        uint lev_mod = two_to_the(icells.levmax - icells.level[i]);
        hash[((icells.j[i] * lev_mod) * i_max) + (icells.i[i] * lev_mod)] = i;
    }
    
    for (uint i = 0; i < ocells.ncells; i++) {
        uint io = ocells.i[i];
        uint jo = ocells.j[i];
        uint lev = ocells.level[i];
        
        uint lev_mod = two_to_the(ocells.levmax - lev);
        uint ii = io*lev_mod;
        uint ji = jo*lev_mod;
        
        uint key = ji*i_max + ii;
        int probe = hash[key];

        if (lev > ocells.levmax){lev = ocells.levmax;}
        
        while(probe < 0 && lev > 0) {
            lev--;
            uint lev_diff = ocells.levmax - lev;
            ii >>= lev_diff;
            ii <<= lev_diff;
            ji >>= lev_diff;
            ji <<= lev_diff;
            key = ji*i_max + ii;
            probe = hash[key];
        }
        if (lev >= icells.level[probe]) {
            ocells.values[i] = icells.values[probe];
        } else {
            ocells.values[i] = avg_sub_cells(icells, ji, ii, lev, hash);
        }
    }
    free(hash);
}

void singlewrite_remap_compact (cell_list icells, cell_list ocells) {
    
    uint i_max = icells.ibasesize*two_to_the(icells.levmax);
    uint j_max = icells.ibasesize*two_to_the(icells.levmax);
    int *hash = compact_hash_init(icells.ncells, i_max, j_max, 1, 0);

//  compact_hash_initializes the keys to -1 (not the values)
    
    for (uint i = 0; i < icells.ncells; i++) {
        uint lev_mod = two_to_the(icells.levmax - icells.level[i]);
        write_hash(i, ((icells.j[i] * lev_mod) * i_max) + (icells.i[i] * lev_mod), hash);
    }
    
    i_max = ocells.ibasesize*two_to_the(ocells.levmax);
    for (uint i = 0; i < ocells.ncells; i++) {
        uint ii, ji;
        uint io = ocells.i[i];
        uint jo = ocells.j[i];
        int lev = ocells.level[i];
        
        uint lev_mod = two_to_the(ocells.levmax - lev);
        ii = io*lev_mod;
        ji = jo*lev_mod;
        
        uint key = ji*i_max + ii;
        int ic = read_hash(key, hash);

        if (lev > (int)ocells.levmax) lev = ocells.levmax;
        while (ic < 0 && lev > 0) {
            lev--;
            uint lev_diff = ocells.levmax - lev;
            ii >>= lev_diff;
            ii <<= lev_diff;
            ji >>= lev_diff;
            ji <<= lev_diff;
            key = ji*i_max + ii;
            ic = read_hash(key, hash);
        }
        if (lev >= (int)icells.level[ic]) {
            ocells.values[i] = icells.values[ic];
        } else {
            ocells.values[i] = avg_sub_cells_compact(icells, ji, ii, lev, hash);
        }
    }
    compact_hash_delete(hash);
}

#ifdef _OPENMP
void singlewrite_remap_openMP (cell_list icells, cell_list ocells) {

    size_t hash_size = icells.ibasesize*two_to_the(icells.levmax)*
                       icells.ibasesize*two_to_the(icells.levmax);
    int *hash = (int *) malloc(hash_size * sizeof(int));

#pragma omp parallel default(none) firstprivate(hash_size) shared(ocells, icells, hash)
    {
        uint ilength = icells.ncells;
        uint olength = ocells.ncells;
        uint max_lev = icells.levmax;

        uint i_max = icells.ibasesize*two_to_the(max_lev);
    
#pragma omp for
        for (uint i = 0; i < hash_size; i++) {
            hash[i] = -1;
        }
    
#pragma omp for
        for (uint i = 0; i < ilength; i++) {
            uint lev_mod = two_to_the(max_lev - icells.level[i]);
            hash[((icells.j[i] * lev_mod) * i_max) + (icells.i[i] * lev_mod)] = i;
        }
    
#pragma omp for
        for (uint i = 0; i < olength; i++) {
            uint ii, ji, lev_mod;
            uint io = ocells.i[i];
            uint jo = ocells.j[i];
            int lev = ocells.level[i];
        
            if (lev < (int)max_lev) {
                lev_mod = two_to_the(max_lev - lev);
                ii = io*lev_mod;
                ji = jo*lev_mod;
            } else {
                lev_mod = two_to_the(lev - max_lev);
                ii = io/lev_mod;
                ji = jo/lev_mod;
            }

            uint key = ji*i_max + ii;
            int ic = hash[key];

            if (lev > (int)max_lev) lev = max_lev;
            while(ic < 0 && lev > 0) {
                lev--;
                uint lev_diff = max_lev - lev;
                ii >>= lev_diff;
                ii <<= lev_diff;
                ji >>= lev_diff;
                ji <<= lev_diff;
                key = ji*i_max + ii;
                ic = hash[key];
            }
            if (lev >= (int)icells.level[ic]) {
                ocells.values[i] = icells.values[ic];
            } else {
                ocells.values[i] = avg_sub_cells(icells, ji, ii, lev, hash);
            }
        }
    }
    free(hash);
}

void singlewrite_remap_compact_openMP (cell_list icells, cell_list ocells) {
    
    uint i_max = icells.ibasesize*two_to_the(icells.levmax);
    uint j_max = icells.ibasesize*two_to_the(icells.levmax);
#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
    int *hash = compact_hash_init_openmp(icells.ncells, i_max, j_max, 1, 0);
#else
    omp_lock_t *lock = NULL;
    int *hash = compact_hash_init_openmp(icells.ncells, i_max, j_max, 1, 0, &lock);
#endif

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
#pragma omp parallel default(none) firstprivate(i_max, j_max) shared(write_hash_openmp, read_hash) shared(ocells, icells, hash)
#else
#pragma omp parallel default(none) firstprivate(i_max, j_max) shared(write_hash_openmp, read_hash) shared(ocells, icells, hash, lock)
#endif
    {
        uint ilength = icells.ncells;
        uint olength = ocells.ncells;
        uint max_lev = icells.levmax;
        //size_t hash_size = i_max*j_max;

        uint i_max = icells.ibasesize*two_to_the(max_lev);
    
// Done in compact hash init
//#pragma omp for
//        for (uint i = 0; i < hash_size; i++) {
//            hash[i] = -1;
//        }
    
#pragma omp for
        for (uint i = 0; i < ilength; i++) {
            uint lev_mod = two_to_the(max_lev - icells.level[i]);
#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
            write_hash_openmp(i, ((icells.j[i] * lev_mod) * i_max) + (icells.i[i] * lev_mod), hash);
#else
            write_hash_openmp(i, ((icells.j[i] * lev_mod) * i_max) + (icells.i[i] * lev_mod), hash, lock);
#endif
        }
    
        /*for (int j = i_max-1; j >= 0; j--) {
            for (int i = 0; i < i_max; i ++) {
                printf("%i\t", hash[j*(i_max) + i]);
            }
            printf("\n");
        }
        printf("\n");*/
    
#pragma omp for
        for (uint i = 0; i < olength; i++) {
            uint ii, ji;
            uint io = ocells.i[i];
            uint jo = ocells.j[i];
            int lev = ocells.level[i];
        
            if (lev < (int)max_lev) {
                uint lev_mod = two_to_the(max_lev - lev);
                ii = io*lev_mod;
                ji = jo*lev_mod;
            } else {
                uint lev_mod = two_to_the(lev - max_lev);
                ii = io/lev_mod;
                ji = jo/lev_mod;
            }
        
            uint key = ji*i_max + ii;
            int ic = read_hash(key, hash);

            if (lev > (int)max_lev) lev = max_lev;
            while (ic < 0 && lev > 0) {
                lev--;
                uint lev_diff = max_lev - lev;
                ii >>= lev_diff;
                ii <<= lev_diff;
                ji >>= lev_diff;
                ji <<= lev_diff;
                key = ji*i_max + ii;
                ic = read_hash(key, hash);
            }
            if (lev >= (int)icells.level[ic]) {
                ocells.values[i] = icells.values[ic];
            } else {
                ocells.values[i] = avg_sub_cells_compact(icells, ji, ii, lev, hash);
            }
            //printf("%i\t%i\t%i\t%f\n", ocells[i].i, ocells[i].j, ocells[i].lev, ocells[i].values);
            //print_cell(ocells[i]);
        }
        //printf("\n");
    }

#ifdef __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4
    compact_hash_delete_openmp(hash);
#else
    compact_hash_delete_openmp(hash, lock);
#endif
}
#endif
