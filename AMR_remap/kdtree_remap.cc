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

#include "meshgen/meshgen.h"
#include "KDTree/KDTree2d.h"
#include "kdtree_remap.h"

void remap_kDtree2d(cell_list icells, cell_list ocells) {

    int num;
    int index_list[ocells.ncells];
    TKDTree2d tree;

    KDTree_Initialize2d(&tree);

    TBounds2d box;

    for(uint ic = 0; ic < icells.ncells; ic++) {
        uint ifactor = two_to_the(icells.levmax - icells.level[ic]);
        box.min.x = icells.i[ic] * ifactor;
        box.max.x = box.min.x    + ifactor - 0.1;
        box.min.y = icells.j[ic] * ifactor;
        box.max.y = box.min.y    + ifactor - 0.1;
        KDTree_AddElement2d(&tree, &box);
    }

    for(uint ic = 0; ic < ocells.ncells; ic++) {
        
        uint ifactor = two_to_the(ocells.levmax - ocells.level[ic]);
        box.min.x = ocells.i[ic] * ifactor;
        box.max.x = box.min.x    + ifactor - 0.1;
        box.min.y = ocells.j[ic] * ifactor;
        box.max.y = box.min.y    + ifactor - 0.1;

        KDTree_QueryBoxIntersect2d(&tree, &num, &(index_list[0]), &box);

        for(int jc = 0; jc < num; jc++) {
            
            uint ifactor = two_to_the(icells.levmax - icells.level[index_list[jc]]);
            double xmin_a = icells.i[index_list[jc]] * ifactor;
            double xmax_a = xmin_a                   + ifactor - 0.1;
            double ymin_a = icells.j[index_list[jc]] * ifactor;
            double ymax_a = ymin_a                   + ifactor - 0.1;

            double overlap_x = MIN(xmax_a, box.max.x) - MAX(xmin_a, box.min.x);
            double overlap_y = MIN(ymax_a, box.max.y) - MAX(ymin_a, box.min.y);

            if(overlap_x > 0 && overlap_y > 0) {
                if (icells.level[index_list[jc]] <= ocells.level[ic]) {
                    //printf("setting b to %f\n", icells.values[index_list[jc]]);
                    ocells.values[ic] = icells.values[index_list[jc]];
                } else {
                    //printf("adding %f to b\n", icells.values[index_list[jc]]
                    //    / four_to_the(levmx - icells.level[index_list[jc]]));
                    ocells.values[ic] += icells.values[index_list[jc]] 
                        / four_to_the(icells.level[index_list[jc]]
                        - ocells.level[ic]);
                }
            }
        }
    }
    KDTree_Destroy2d(&tree);

   return;

}
