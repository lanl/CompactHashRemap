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

#include "unstructured_types.h"

void get_bounding_box(uint cell_id, double bin_size, cell_list cells, node_list nodes,
    face_list faces, uint *min_x, uint *max_x, uint *min_y, uint *max_y,
    uint *dx, uint *dy);
void get_local_hash (int *local_hash, uint cell_id, double bin_size, cell_list cells,
    face_list faces, node_list nodes, uint num_faces, uint min_x, uint dx, uint min_y, 
    uint dy);
void write_edge (int *hash, uint hash_width, double bin_size, node_list nodes, 
    face_list faces,uint edge, uint cell, uint min_x, uint min_y);
void fill_polygon (uint cell_id, int *hash, uint hash_width, uint hash_height);
uint get_key(uint i, uint j, uint hash_width);
uint key_to_i(uint key, uint hash_width);
uint key_to_j(uint key, uint hash_width);
