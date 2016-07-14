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
#include "brute_force_remap.h"

uint translate_cell (uint i, uint j, uint lev, uint new_lev, uint ibasesize) {
    uint j_comp, i_comp;
    if (new_lev < lev) {
        //j_comp = j / two_to_the(lev - new_lev);
        //i_comp = i / two_to_the(lev - new_lev);
        j_comp = j >> (lev - new_lev);
        i_comp = i >> (lev - new_lev);
    } else {
        //j_comp = j * two_to_the(new_lev - lev);
        //i_comp = i * two_to_the(new_lev - lev);
        j_comp = j << (new_lev - lev);
        i_comp = i << (new_lev - lev);
    }

    //printf("j_comp: %u\tx_comp: %u\n", j_comp, i_comp);
    uint key = (j_comp * ibasesize*two_to_the(new_lev)) + i_comp;
    //printf("newKey: %d\n", newKey);
    return key;
}

void brute_force_remap (cell_list icells, cell_list ocells) {
    
    for(uint o = 0; o<ocells.ncells; o++){
        ocells.values[o] = 0;
    }
    
    for (uint o = 0; o < ocells.ncells; o++){
        for (uint i = 0; i < icells.ncells; i++){
            uint lev_out = ocells.level[o];
            uint lev_in = icells.level[i];
            uint i_out = ocells.i[o];
            uint j_out = ocells.j[o];
            uint i_in = icells.i[i];
            uint j_in = icells.j[i];
            if (lev_out > lev_in) {
                uint key_in = j_in * icells.ibasesize*two_to_the(lev_in) + i_in;
                uint key_out = translate_cell(i_out, j_out, lev_out, lev_in, ocells.ibasesize);
                if (key_in == key_out){
                    ocells.values[o] = icells.values[i];
                    continue;
                }
            }else{
                uint key_out = j_out * ocells.ibasesize*two_to_the(lev_out) + i_out;
                uint key_in = translate_cell(i_in, j_in, lev_in, lev_out, icells.ibasesize);
                if (key_out == key_in){
                    ocells.values[o] += (icells.values[i] / four_to_the(lev_in - lev_out));
                }
            }
        }
    }
}

