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

#include <list>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "local_hash.h"
#include "unstructured_read.h"

void unstructured_read_list (uint num_cells_target, uint hash_width, 
    double bin_size, node_list nodes, face_list faces, cell_list cells, 
    int **hash, uint *histogram, uint **candidates) {
    
    int probe;
    uint key1, key2;
    
    for (uint n = 0; n < num_cells_target; n++) {
        uint dx, dy, min_x, min_y, max_x, max_y;
        get_bounding_box(n, bin_size, cells, nodes, faces, &min_x, &max_x, &min_y, 
            &max_y, &dx, &dy);
            
        int *local_read = (int *) malloc(dx * dy * sizeof(int));
        
        get_local_hash(local_read, n, bin_size, cells, faces, nodes, cells.num_faces[n],
            min_x, dx, min_y, dy);
        
        std::list<uint> candidates_list;
            
        for (uint m = 0; m < dx * dy; m++) {
            if (local_read[m] > -1) {
                key1 = get_key(key_to_i(m, dx) + min_x, key_to_j(m, dx) + min_y, 
                    hash_width);
                
                for (key2 = 0; key2 < histogram[key1]; key2++) {
                    probe = hash[key1][key2];
                    candidates_list.push_back(probe);
                }
            }
        }
        candidates_list.sort();
        candidates_list.unique();
        candidates[n] = (uint *) malloc(candidates_list.size() * sizeof(uint));
        
        uint m = 0;
        for (std::list<uint>::iterator it=candidates_list.begin(); 
            it != candidates_list.end(); ++it) {
            
            candidates[n][m] = *it;
            m++;
        }
    }
    return;    
}

/*void unstructured_read_hash (uint num_cells_target, uint hash_width, 
    double bin_size, node_list nodes, face_list faces, cell_list cells, 
    int **hash, uint *histogram, uint **candidates) {
    
    int probe;
    uint key1, key2;
    
    for (uint n = 0; n < num_cells_target; n++) {
        uint dx, dy, min_x, min_y, max_x, max_y, num_collisions = 1;
        get_bounding_box(n, bin_size, cells, nodes, faces, &min_x, &max_x, &min_y, 
            &max_y, &dx, &dy);
            
        int *local_read = (int *) malloc(dx * dy * sizeof(int));
        
        get_local_hash(local_read, n, bin_size, cells, faces, nodes, cells.num_faces[n],
            min_x, dx, min_y, dy);
        
        std::list<uint> candidates_list;
            
        for (uint m = 0; m < dx * dy; m++) {
            if (local_read[m] > -1) {
                key1 = get_key(key_to_i(m, dx) + min_x, key_to_j(m, dx) + min_y, 
                    hash_width);
                
                for (key2 = 0; key2 < histogram[key1]; key2++) {
                    probe = hash[key1][key2];
                    candidates_list.push_back(probe);
                }
            }
        }
        candidates_list.sort();
        candidates_list.unique();
        candidates[n] = (uint *) malloc(candidates_list.size() * sizeof(uint));
        
        uint m = 0;
        for (std::list<uint>::iterator it=candidates_list.begin(); 
            it != candidates_list.end(); ++it) {
            
            candidates[n][m] = *it;
            m++;
        }
    }
    return;    
}*/

/*uint **try_to_hash (uint id, uint *candidates, uint *length, uint *length, 
    uint *max, uint *min, uint *range, uint *bin_size) {
    if (*length > 1) {
        uint key = (id - *min) / *bin_size;
        if (candidates[key] != id) {
                
        }
    }
}*/

int test_read() {
    int success = 1;
    
    uint hash_width = 7;
    uint hash_height = 8;
    double bin_size = 1;
    
    int **hash = (int **) malloc (hash_width * hash_height * sizeof(uint *));
    uint *histogram = (uint *) malloc(hash_width * hash_height * sizeof(uint));
    
    
    uint num_cells = 7;
    uint num_nodes = 12;
    uint num_faces = 27;
    
    uint **candidates = (uint **) malloc(num_cells * sizeof(uint *));    
    
    node_list nodes;
    nodes.i = (double *) malloc(num_nodes * sizeof(double));
    nodes.j = (double *) malloc(num_nodes * sizeof(double));
    nodes.k = (double *) malloc(num_nodes * sizeof(double));
    
    face_list faces;
    faces.num_nodes = (uint *) malloc(num_faces * sizeof(uint));
    faces.nodes = (uint **) malloc(num_faces * sizeof(uint *));
    
    cell_list cells;
    cells.num_faces = (uint *) malloc(num_cells * sizeof(uint));
    cells.faces = (uint **) malloc(num_cells * sizeof(uint *)); 
    
    nodes.i[0] = 0.5;
    nodes.i[1] = 1.2;
    nodes.i[2] = 4.2;
    nodes.i[3] = 3.5;
    nodes.i[4] = 5.3;
    nodes.i[5] = 2.5;
    nodes.i[6] = 1.7;
    nodes.i[7] = 3.2;
    nodes.i[8] = 5.7;
    nodes.i[9] = 3.3;
    nodes.i[10] = 4.8;
    nodes.i[11] = 6.5;
    
    nodes.j[0] = 3.5;
    nodes.j[1] = 6.5;
    nodes.j[2] = 7.1;
    nodes.j[3] = 5.65;
    nodes.j[4] = 5.5;
    nodes.j[5] = 4.5;
    nodes.j[6] = 1.8;
    nodes.j[7] = 3.5;
    nodes.j[8] = 3.5;
    nodes.j[9] = 0.3;
    nodes.j[10] = 1.5;
    nodes.j[11] = 1.5;
    
    for (uint n = 0; n < num_faces; n++) {
        faces.num_nodes[n] = 2;
    }
    
    for (uint n = 0; n < num_faces; n++) {
        faces.nodes[n] = (uint *) malloc(faces.num_nodes[n] * sizeof(uint));
    }
    
    faces.nodes[0][0] = 1;
    faces.nodes[0][1] = 0;
    faces.nodes[1][0] = 1;
    faces.nodes[1][1] = 3;
    faces.nodes[2][0] = 5;
    faces.nodes[2][1] = 3;
    faces.nodes[3][0] = 0;
    faces.nodes[3][1] = 5;
    faces.nodes[4][0] = 2;
    faces.nodes[4][1] = 1;
    faces.nodes[5][0] = 3;
    faces.nodes[5][1] = 2;
    faces.nodes[6][0] = 1;
    faces.nodes[6][1] = 3;
    faces.nodes[7][0] = 2;
    faces.nodes[7][1] = 3;
    faces.nodes[8][0] = 2;
    faces.nodes[8][1] = 4;
    faces.nodes[9][0] = 4;
    faces.nodes[9][1] = 8;
    faces.nodes[10][0] = 3;
    faces.nodes[10][1] = 8;
    faces.nodes[11][0] = 3;
    faces.nodes[11][1] = 5;
    faces.nodes[12][0] = 3;
    faces.nodes[12][1] = 8;
    faces.nodes[13][0] = 5;
    faces.nodes[13][1] = 7;
    faces.nodes[14][0] = 8;
    faces.nodes[14][1] = 7;
    faces.nodes[15][0] = 0;
    faces.nodes[15][1] = 6;
    faces.nodes[16][0] = 5;
    faces.nodes[16][1] = 0;
    faces.nodes[17][0] = 7;
    faces.nodes[17][1] = 5;
    faces.nodes[18][0] = 7;
    faces.nodes[18][1] = 6;
    faces.nodes[19][0] = 7;
    faces.nodes[19][1] = 6;
    faces.nodes[20][0] = 7;
    faces.nodes[20][1] = 10;
    faces.nodes[21][0] = 9;
    faces.nodes[21][1] = 6;
    faces.nodes[22][0] = 9;
    faces.nodes[22][1] = 10;
    faces.nodes[23][0] = 7;
    faces.nodes[23][1] = 10;
    faces.nodes[24][0] = 7;
    faces.nodes[24][1] = 8;
    faces.nodes[25][0] = 8;
    faces.nodes[25][1] = 11;
    faces.nodes[26][0] = 10;
    faces.nodes[26][1] = 11;
    
    cells.num_faces[0] = 4;
    cells.num_faces[1] = 3;
    cells.num_faces[2] = 4;
    cells.num_faces[3] = 4;
    cells.num_faces[4] = 4;
    cells.num_faces[5] = 4;
    cells.num_faces[6] = 4;
    
    for (uint n = 0; n < num_cells; n++) {
        cells.faces[n] = (uint *) malloc(cells.num_faces[n] * sizeof(uint));
    }
    
    cells.faces[0][0] = 0;
    cells.faces[0][1] = 1;
    cells.faces[0][2] = 2;
    cells.faces[0][4] = 3;
    cells.faces[1][0] = 4;
    cells.faces[1][1] = 5;
    cells.faces[1][2] = 6;
    cells.faces[2][0] = 7;
    cells.faces[2][1] = 8;
    cells.faces[2][2] = 9;
    cells.faces[2][3] = 10;
    cells.faces[3][0] = 11;
    cells.faces[3][1] = 12;
    cells.faces[3][2] = 13;
    cells.faces[3][3] = 14;
    cells.faces[4][0] = 15;
    cells.faces[4][1] = 16;
    cells.faces[4][2] = 17;
    cells.faces[4][3] = 18;
    cells.faces[5][0] = 19;
    cells.faces[5][1] = 20;
    cells.faces[5][2] = 21;
    cells.faces[5][3] = 22;
    cells.faces[6][0] = 23;
    cells.faces[6][1] = 24;
    cells.faces[6][2] = 25;
    cells.faces[6][3] = 26;
    
    for (uint n = 0; n < hash_width * hash_height; n++) {
        histogram[n] = 0;
    }
    
    histogram[3] = 1;
    histogram[4] = 1;
    histogram[8] = 1;
    histogram[9] = 1;
    histogram[10] = 1;
    histogram[11] = 1;
    histogram[12] = 1;
    histogram[15] = 2;
    histogram[16] = 2;
    histogram[17] = 2;
    histogram[18] = 3;
    histogram[19] = 2;
    histogram[20] = 2;
    histogram[22] = 1;
    histogram[23] = 1;
    histogram[24] = 1;
    histogram[25] = 2;
    histogram[26] = 1;
    histogram[27] = 1;
    histogram[29] = 1;
    histogram[30] = 1;
    histogram[31] = 1;
    histogram[32] = 2;
    histogram[33] = 1;
    histogram[34] = 1;
    histogram[36] = 2;
    histogram[37] = 2;
    histogram[38] = 3;
    histogram[39] = 4;
    histogram[40] = 1;
    histogram[41] = 1;
    histogram[43] = 1;
    histogram[44] = 2;
    histogram[45] = 2;
    histogram[46] = 2;
    histogram[47] = 2;
    histogram[50] = 1;
    histogram[51] = 2;
    histogram[52] = 2;
    histogram[53] = 1;
    histogram[54] = 2;
    histogram[55] = 2;
    
    for (uint n = 0; n < hash_width * hash_height; n++) {
        hash[n] = (int *) malloc(histogram[n] * sizeof(int));
    }
    
    for (uint n = 0; n < hash_width * hash_height; n++) {
        for (uint m = 0; m < histogram[n]; m++) {
            hash[n][m] = -1;
        }
    }
    

    hash[3][0] = 4;
    hash[4][0] = 4;
    hash[8][0] = 4;
    hash[9][0] = 4;
    hash[10][0] = 4;
    hash[11][0] = 4;
    hash[12][0] = 4;
    hash[15][0] = 2;
    hash[15][1] = 4;
    hash[16][0] = 2;
    hash[16][1] = 4;
    hash[17][0] = 2;
    hash[17][1] = 4;
    hash[18][0] = 2;
    hash[18][1] = 3;
    hash[18][2] = 4;
    hash[19][0] = 3;
    hash[19][1] = 4;
    hash[20][0] = 3;
    hash[20][1] = 4;
    hash[22][0] = 2;
    hash[23][0] = 2;
    hash[24][0] = 2;
    hash[25][0] = 2;
    hash[25][1] = 3;
    hash[26][0] = 3;
    hash[27][0] = 3;
    hash[29][0] = 2;
    hash[30][0] = 2;
    hash[31][0] = 2;
    hash[32][0] = 2;
    hash[32][1] = 3;
    hash[33][0] = 3;
    hash[34][0] = 3;
    hash[36][0] = 0;
    hash[36][1] = 2;
    hash[37][0] = 0;
    hash[37][1] = 2;
    hash[38][0] = 0;
    hash[38][1] = 1;
    hash[38][2] = 2;
    hash[39][0] = 0;
    hash[39][1] = 1;
    hash[39][2] = 2;
    hash[39][3] = 3;
    hash[40][0] = 3;
    hash[41][0] = 3;
    hash[43][0] = 0;
    hash[44][0] = 0;
    hash[44][1] = 1;
    hash[45][0] = 0;
    hash[45][1] = 1;
    hash[46][0] = 1;
    hash[46][1] = 3;
    hash[47][0] = 1;
    hash[47][1] = 3;
    hash[48][0] = 3;
    hash[51][0] = 0;
    hash[51][1] = 1;
    hash[52][0] = 0;
    hash[52][1] = 1;
    hash[53][0] = 1;
    hash[54][0] = 1;
    hash[54][1] = 3;
    hash[55][0] = 1;
    hash[55][1] = 3;
    
    unstructured_read_list (num_cells, hash_width, bin_size, nodes, faces, 
        cells, hash, histogram, candidates);
    
    uint **candidates_true = (uint **) malloc(num_cells * sizeof(uint *));
    
    candidates_true[0] = (uint *) malloc(3 * sizeof(uint));
    candidates_true[1] = (uint *) malloc(4 * sizeof(uint));
    candidates_true[2] = (uint *) malloc(4 * sizeof(uint));
    candidates_true[3] = (uint *) malloc(4 * sizeof(uint));
    candidates_true[4] = (uint *) malloc(2 * sizeof(uint));
    candidates_true[5] = (uint *) malloc(3 * sizeof(uint));
    candidates_true[6] = (uint *) malloc(3 * sizeof(uint));
    
    candidates_true[0][0] = 0;
    candidates_true[0][1] = 1;
    candidates_true[0][2] = 2;
    candidates_true[1][0] = 0;
    candidates_true[1][1] = 1;
    candidates_true[1][2] = 2;
    candidates_true[1][3] = 3;
    candidates_true[2][0] = 0;
    candidates_true[2][1] = 1;
    candidates_true[2][2] = 2;
    candidates_true[2][3] = 3;
    candidates_true[3][0] = 0;
    candidates_true[3][1] = 1;
    candidates_true[3][2] = 2;
    candidates_true[3][3] = 3;
    candidates_true[4][0] = 2;
    candidates_true[4][1] = 4;
    candidates_true[5][0] = 2;
    candidates_true[5][1] = 3;
    candidates_true[5][2] = 4;
    candidates_true[6][0] = 2;
    candidates_true[6][1] = 3;
    candidates_true[6][2] = 4;
    
    for (uint n = 0; n < 3; n++) {
        if (candidates[0][n] != candidates_true[0][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 4; n++) {
        if (candidates[1][n] != candidates_true[1][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 4; n++) {
        if (candidates[2][n] != candidates_true[2][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 4; n++) {
        if (candidates[3][n] != candidates_true[3][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 2; n++) {
        if (candidates[4][n] != candidates_true[4][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 3; n++) {
        if (candidates[5][n] != candidates_true[5][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 3; n++) {
        if (candidates[6][n] != candidates_true[6][n]) {
            success = 0;
        }
    }
    return success;
}

void swap_uint(uint* a, uint* b) {
   uint c = *a;
   *a = *b;
   *b = c;
}

void bubblesort(uint length, uint *arr) {
	uint n = length;
	bool swapped;
	do {
		swapped = false;
		for (uint i = 1; i < n; i++) {
			if (arr[i-1] > arr[i]) {
				swap_uint(&arr[i-1], &arr[i]);
				swapped = true;
			}
		}
		n--;
	} while (swapped == true);
}

uint* procrastisort (uint length, uint *arr) {
	/*stops the recursion*/
    if (length > 1) {
		uint min = arr[0];
		uint max = arr[0];
		
		/*finds the min and max. demanding procrastisort has min and max as an argument*/
		for (uint i = 1; i < length; i++) {
			if (arr[i] < min) {
				min = arr[i];
			}
			if (arr[i] > max) {
				max = arr[i];
			}
		}
		    
    	uint range = max - min;
    	
		/*if range is 0, the numbers are all the same and don't need to be sorted*/
    	if (range > 0) {
    		
    		uint *sorted = (uint *) malloc(length*sizeof(uint));
    	
    		uint binIndex;
    		double binSize = (range)/(length-1);

			/*subtract min out*/
    		if (min > 0) {
    			//#pragma omp parallel for
    			for (uint i = 0; i < length; i++) {
    				arr[i] -= min;
    			}
    		}
    
			/*get histogram and prefix scan*/
	    	uint *num = (uint*) malloc(length*sizeof(uint));
	    	uint *sum = (uint*) malloc((length+1)*sizeof(uint));
    		
			/*initialise to zero*/
			memset(num, 0, (length)*sizeof(uint));
    
    		//#pragma omp parallel for
	    	for(uint i = 0; i < length; i++) {
				binIndex = (uint) ((arr[i])/binSize);
				num[binIndex]++;
	    	}
	    	
			/*prefix scanning*/
			sum[0] = 0;
	    	for (uint i = 1; i < length+1; i++) {
				sum[i] = sum[i-1] + num[i-1];
	    	}
	    	free(num);

	    	/*set all elements of hash array to -1*/
	    	int *hash = (int*)malloc((length)*sizeof(int));
	        //memset(hash, -1, (length)*sizeof(int));

			/*hashing usig the prefix scan*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = sum[(int) ((arr[i])/binSize)];
				hash[binIndex] = i; 
				sum[binIndex]++;
	    	}
    		
    		/*finsih by using the keys in the hash*/
			if (min > 0) {
				for (uint i = 0; i < length; i++) {
	    			arr[hash[i]] += min;
	    		    sorted[i] = arr[hash[i]];
	    		}
	    	} else {
	    		for (uint i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]];
	    		}
	    	}
    		
    		/*handle collisions*/
	    	bubblesort(length, sorted);
			
			
	    	free(sum);
	    	free(hash);
	    	return sorted;
	    } else {
	    	return arr;
    	}
	} else {
		return arr;
	}
}
