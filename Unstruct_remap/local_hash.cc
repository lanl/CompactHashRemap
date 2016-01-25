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

#include <limits.h>
#include <stdio.h>
#include "local_hash.h"

void get_bounding_box(uint cell_id, double bin_size, cell_list cells, node_list nodes,
    face_list faces, uint *min_x, uint *max_x, uint *min_y, uint *max_y,
    uint *dx, uint *dy) {
    
    *min_x = INT_MAX;
    *max_x = 0;
    *min_y = INT_MAX;
    *max_y = 0;
    uint num_faces = cells.num_faces[cell_id];
    //I don't want to write a new filling algorithm that is generic for dimensionality
    //so I am starting it here but hard-coding 2d for now.
    uint num_nodes_per_face = 2;
    
    for (uint n = 0; n < num_faces; n++) {
        for (uint m = 0; m < num_nodes_per_face; m++) {
            //printf("c: %u\t", cells.faces[cell_id][n]);
            //printf("f: %u\t", faces.nodes[cells.faces[cell_id][n]][m]);
            //fflush(stdout);
            //printf("n: %f\t", nodes.i[faces.nodes[cells.faces[cell_id][n]][m]]);
            if ((uint) (nodes.i[faces.nodes[cells.faces[cell_id][n]][m]] / bin_size) 
            > *max_x) {
            *max_x = (uint) (nodes.i[faces.nodes[cells.faces[cell_id][n]][m]] /
            bin_size);
            }
        if ((uint) (nodes.i[faces.nodes[cells.faces[cell_id][n]][m]] / bin_size) 
            < *min_x) {
            *min_x = (uint) (nodes.i[faces.nodes[cells.faces[cell_id][n]][m]] /
            bin_size);
            }
        if ((uint) (nodes.j[faces.nodes[cells.faces[cell_id][n]][m]] / bin_size) 
            > *max_y) {
            *max_y = (uint) (nodes.j[faces.nodes[cells.faces[cell_id][n]][m]] /
            bin_size);
            }
        if ((uint) (nodes.j[faces.nodes[cells.faces[cell_id][n]][m]] / bin_size) 
            < *min_y) {
            *min_y = (uint) (nodes.j[faces.nodes[cells.faces[cell_id][n]][m]] /
            bin_size);
            }
        }
    }
    
    //printf("minx: %i, miny %i\n", *min_x, *min_y);
    
    *dx = (*max_x - *min_x) + 1;
    *dy = (*max_y - *min_y) + 1;
    
    return;
}

void get_local_hash (int *local_hash, uint cell_id, double bin_size, cell_list cells,
    face_list faces, node_list nodes, uint num_faces, uint min_x, uint dx, uint min_y, 
    uint dy) {
    
    for (uint n = 0; n < (dx * dy); n++) {
        local_hash[n] = -1;
    }
    
    for (uint n = 0; n < num_faces; n++) {
        //printf("n: %u\n", n);
        write_edge (local_hash, dx, bin_size, nodes, faces,
            cells.faces[cell_id][n], cell_id, min_x, min_y);
    }
    
    fill_polygon (cell_id, local_hash, dx, dy);
    
    /*for (int i = 0; i < dx*dy; i++) {
        printf("%i\t", local_hash[i]);
    }
    printf("\n");*/
}

void write_edge (int *hash, uint hash_width, double bin_size, node_list nodes, 
    face_list faces, uint edge, uint cell, uint min_x, uint min_y) {
    
    double node_a_x, node_a_y, node_b_x, node_b_y, y_0, y_f, x_0, x_f, m, y_left, y_right;
    uint x, x_end, y, y_end;
    
    node_a_x = nodes.i[faces.nodes[edge][0]];
    node_a_y = nodes.j[faces.nodes[edge][0]];
    node_b_x = nodes.i[faces.nodes[edge][1]];
    node_b_y = nodes.j[faces.nodes[edge][1]];
    
    if (node_a_x < node_b_x) {
        x_0 = node_a_x / bin_size - min_x;
        y_0 = node_a_y / bin_size - min_y;
        x_f = node_b_x / bin_size - min_x;
        y_f = node_b_y / bin_size - min_y;
    } else if (node_a_x > node_b_x) {
        x_0 = node_b_x / bin_size - min_x;
        y_0 = node_b_y / bin_size - min_y;
        x_f = node_a_x / bin_size - min_x;
        y_f = node_a_y / bin_size - min_y;
    } else {
        x_0 = node_a_x / bin_size - min_x;
        y_0 = node_a_y / bin_size - min_y;
        y_f = node_b_y / bin_size - min_y;
        x = (uint) x_0;
        y = (uint) y_0;
        y_end = (uint) y_f;
        if (y_end > y) {
            for (uint i = y; i < y_end + 1; i++) {
                hash[get_key(x, i, hash_width)] = cell;
            }
        } else {
            for (int i = y; i >= (int) y_end && i > -1; i--) {
                hash[get_key(x, i, hash_width)] = cell;
            }
        }
        return;
    }
    
    x = (uint) x_0;
    x_end = (uint) x_f;

    m = (y_f - y_0) / (x_f - x_0);
    
    y_left = m * (x - x_0) + y_0;
    y_right = y_left + m;
    y = (uint) y_0;
    y_end = (uint) y_right;
    y_left = y_0;
    
    //skip the first bin so the end of another face doesn't redundantly write
    //y++;
    
    if (m > 0) {
        for (uint i = y; i < y_end + 1; i++) {
            hash[get_key(x, i, hash_width)] = cell;
        }
    } else {
        for (int i = y; i >= (int) y_end && i > -1; i--) {
            hash[get_key(x, i, hash_width)] = cell;
        }
    }
    
    uint key1;
    for (uint i = x + 1; i < x_end + 1; i++) {
        y_left = y_right;
        y_right = y_left + m;
        y = (uint) y_left;
        y_end = (uint) y_right;
        if (m > 0) {
            for (uint j = y; j < y_end + 1 && j <= (uint) y_f; j++) {
                key1 = get_key(i, j, hash_width);
                hash[key1] = cell;
            }
        } else {
            for (int j = y; j >= (int) y_end && j >= (int) y_f && j > -1; j--) {
                key1 = get_key(i, j, hash_width);
                hash[key1] = cell;
            }
        }
    }
    return;
    
}

void fill_polygon (uint cell_id, int *hash, uint hash_width, uint hash_height) {
    uint key;
    
    for (uint j = 1; j < hash_height - 1; j++) {
        int flag = 0;
        for (uint i = 0; i < hash_width - 1; i++) {
            key = get_key(i, j, hash_width);
            if (hash[key] == (int) cell_id) {
                if (hash[key + 1] != (int) cell_id){
                    if (flag == 0) {
                        flag = 1;
                    } else {
                        break;
                    }
                }
                continue;
            }
            if (flag == 1) {
                hash[key] = cell_id;
            }
        }
    }
    //go left to right until hit first 
}

uint get_key(uint i, uint j, uint hash_width) {
    return ((j * hash_width) + i);
}

uint key_to_i(uint key, uint hash_width) {
    return (key % hash_width);
}

uint key_to_j(uint key, uint hash_width) {
    return (key / hash_width);
}
