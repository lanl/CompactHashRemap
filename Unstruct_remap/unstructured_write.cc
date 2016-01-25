/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */

#include <stdlib.h>
#include <stdio.h>
#include "local_hash.h"
#include "unstructured_write.h"

void write_to_hash (uint num_cells, uint hash_width, 
    uint hash_height, double bin_size, node_list nodes, face_list faces, cell_list cells, 
    int **hash, uint *histogram) {
    
    uint *dx_list = (uint *) malloc(num_cells * sizeof(uint));
    uint *dy_list = (uint *) malloc(num_cells * sizeof(uint));
    uint *min_x_list = (uint *) malloc(num_cells * sizeof(uint));
    uint *min_y_list = (uint *) malloc(num_cells * sizeof(uint));
    int **local_hashes = (int **) malloc(num_cells * sizeof(int *));
    
    uint max_x, max_y;
    int probe;
    uint key1, key2;
    
    for (uint i = 0; i < hash_width * hash_height; i++) {
            histogram[i] = 0;
    }
    
    for (uint n = 0; n < num_cells; n++) {
        
        //printf("Cell\t%i:\n", n);
        
        get_bounding_box(n, bin_size, cells, nodes, faces, &min_x_list[n], &max_x, 
            &min_y_list[n], &max_y, &dx_list[n], &dy_list[n]);
        
        //printf("dx * dy:\t%i\n", dx_list[n]*dy_list[n]);
        
        local_hashes[n] = (int *) malloc(dx_list[n] * dy_list[n] * sizeof(int));
        
        get_local_hash(local_hashes[n], n, bin_size, cells, faces, nodes,
            cells.num_faces[n], min_x_list[n], dx_list[n], min_y_list[n], dy_list[n]);
        
        for (uint m = 0; m < dx_list[n] * dy_list[n]; m++) {
            //printf("%i\t", local_hashes[n][m]);
            if (local_hashes[n][m] > -1) {
                histogram[get_key(key_to_i(m, dx_list[n]) + min_x_list[n], 
                    key_to_j(m, dx_list[n]) + min_y_list[n], hash_width)]++;
            }
        }
        //printf("\n");
    }
    
//    for (uint n = 0; n < num_cells; n++) {
//        for (int j = dy_list[n] - 1; j >= 0; j--) {
//            for (uint i = 0; i < dx_list[n]; i++) {
//                int key = get_key(i, j, dx_list[n]);
//                printf("%i\t", local_hashes[n][key]);
//            }
//            printf("\n");
//        }
//        printf("\n");
//    }
    
    for (uint n = 0; n < hash_width * hash_height; n++) {
        hash[n] = (int *) malloc(histogram[n] * sizeof(int));
        for (uint m = 0; m < histogram[n]; m++) {
            hash[n][m] = -1;
        }
    }
    
    for (uint n = 0; n < num_cells; n++) {
        for (uint m = 0; m < dx_list[n] * dy_list[n]; m++) {
            if (local_hashes[n][m] > -1) {
                key1 = get_key(key_to_i(m, dx_list[n]) + min_x_list[n], 
                    key_to_j(m, dx_list[n]) + min_y_list[n], hash_width);
                key2 = 0;
                probe = hash[key1][key2];
                while (probe > -1) {
                    key2++;
                    probe = hash[key1][key2];
                }
                hash[key1][key2] = n;
            }
        }
        free(local_hashes[n]);
    }
}

int test_write() {

    int success = 1;
    
    uint hash_width = 7;
    uint hash_height = 8;
    double bin_size = 1;
    int **hash_true = (int **) malloc (hash_width * hash_height * sizeof(int *));
    uint *histogram_true = (uint *) malloc(hash_width * hash_height * sizeof(uint));
    
    uint num_cells = 5;
    uint num_nodes = 8;
    uint num_faces = 18;
    
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
    
    nodes.i[0] = 2.8;
    nodes.i[1] = 6.5;
    nodes.i[2] = 1.5;
    nodes.i[3] = 4;
    nodes.i[4] = 1.5;
    nodes.i[5] = 4;
    nodes.i[6] = 6.5;
    nodes.i[7] = 4;
    
    nodes.j[0] = 7.99;
    nodes.j[1] = 7.99;
    nodes.j[2] = 5.5;
    nodes.j[3] = 5.5;
    nodes.j[4] = 2.5;
    nodes.j[5] = 2.5;
    nodes.j[6] = 2.5;
    nodes.j[7] = 0;
    
    for (uint n = 0; n < num_faces; n++) {
        faces.num_nodes[n] = 2;
    }
    
    for (uint n = 0; n < num_faces; n++) {
        faces.nodes[n] = (uint *)malloc(faces.num_nodes[n] * sizeof(uint));
    }
    
    faces.nodes[0][0] = 2;
    faces.nodes[0][1] = 0;
    faces.nodes[1][0] = 0;
    faces.nodes[1][1] = 3;
    faces.nodes[2][0] = 2;
    faces.nodes[2][1] = 3;
    faces.nodes[3][0] = 0;
    faces.nodes[3][1] = 3;
    faces.nodes[4][0] = 0;
    faces.nodes[4][1] = 1;
    faces.nodes[5][0] = 3;
    faces.nodes[5][1] = 1;
    faces.nodes[6][0] = 2;
    faces.nodes[6][1] = 4;
    faces.nodes[7][0] = 2;
    faces.nodes[7][1] = 3;
    faces.nodes[8][0] = 3;
    faces.nodes[8][1] = 5;
    faces.nodes[9][0] = 4;
    faces.nodes[9][1] = 5;
    faces.nodes[10][0] = 3;
    faces.nodes[10][1] = 5;
    faces.nodes[11][0] = 3;
    faces.nodes[11][1] = 1;
    faces.nodes[12][0] = 1;
    faces.nodes[12][1] = 6;
    faces.nodes[13][0] = 5;
    faces.nodes[13][1] = 6;
    faces.nodes[14][0] = 4;
    faces.nodes[14][1] = 7;
    faces.nodes[15][0] = 4;
    faces.nodes[15][1] = 5;
    faces.nodes[16][0] = 5;
    faces.nodes[16][1] = 6;
    faces.nodes[17][0] = 7;
    faces.nodes[17][1] = 6;
    
    cells.num_faces[0] = 3;
    cells.num_faces[1] = 3;
    cells.num_faces[2] = 4;
    cells.num_faces[3] = 4;
    cells.num_faces[4] = 4;
    
    for (uint n = 0; n < num_cells; n++) {
        cells.faces[n] = (uint*) malloc(cells.num_faces[n] * sizeof(uint));
    }
    
    cells.faces[0][0] = 0;
    cells.faces[0][1] = 1;
    cells.faces[0][2] = 2;
    cells.faces[1][0] = 3;
    cells.faces[1][1] = 4;
    cells.faces[1][2] = 5;
    cells.faces[2][0] = 6;
    cells.faces[2][1] = 7;
    cells.faces[2][2] = 8;
    cells.faces[2][3] = 9;
    cells.faces[3][0] = 10;
    cells.faces[3][1] = 11;
    cells.faces[3][2] = 12;
    cells.faces[3][3] = 13;
    cells.faces[4][0] = 14;
    cells.faces[4][1] = 15;
    cells.faces[4][2] = 16;
    cells.faces[4][3] = 17;
    
    int **hash = (int **) malloc(hash_width * hash_height * sizeof(int *));
    uint *histogram = (uint *) malloc(hash_width* hash_height * sizeof(uint)); 
    
    write_to_hash (num_cells, hash_width, hash_height, 
        bin_size, nodes, faces, cells, hash, histogram);
    
    for (uint n = 0; n < hash_width * hash_height; n++) {
        histogram_true[n] = 0;
    }

    histogram_true[3] = 1;
    histogram_true[4] = 1;
    histogram_true[8] = 1;
    histogram_true[9] = 1;
    histogram_true[10] = 1;
    histogram_true[11] = 1;
    histogram_true[12] = 1;
    histogram_true[15] = 2;
    histogram_true[16] = 2;
    histogram_true[17] = 2;
    histogram_true[18] = 3;
    histogram_true[19] = 2;
    histogram_true[20] = 2;
    histogram_true[22] = 1;
    histogram_true[23] = 1;
    histogram_true[24] = 1;
    histogram_true[25] = 2;
    histogram_true[26] = 1;
    histogram_true[27] = 1;
    histogram_true[29] = 1;
    histogram_true[30] = 1;
    histogram_true[31] = 1;
    histogram_true[32] = 2;
    histogram_true[33] = 1;
    histogram_true[34] = 1;
    histogram_true[36] = 2;
    histogram_true[37] = 2;
    histogram_true[38] = 3;
    histogram_true[39] = 4;
    histogram_true[40] = 1;
    histogram_true[41] = 1;
    histogram_true[43] = 1;
    histogram_true[44] = 2;
    histogram_true[45] = 2;
    histogram_true[46] = 2;
    histogram_true[47] = 2;
    histogram_true[50] = 1;
    histogram_true[51] = 2;
    histogram_true[52] = 2;
    histogram_true[53] = 1;
    histogram_true[54] = 2;
    histogram_true[55] = 2;
    
    for (uint n = 0; n < hash_width * hash_height; n++) {
        hash_true[n] = (int *) malloc(histogram_true[n] * sizeof(int));
    }

    hash_true[3][0] = 4;
    hash_true[4][0] = 4;
    hash_true[9][0] = 4;
    hash_true[10][0] = 4;
    hash_true[11][0] = 4;
    hash_true[12][0] = 4;
    hash_true[15][0] = 2;
    hash_true[15][1] = 4;
    hash_true[16][0] = 2;
    hash_true[16][1] = 4;
    hash_true[17][0] = 2;
    hash_true[17][1] = 4;
    hash_true[18][0] = 2;
    hash_true[18][1] = 3;
    hash_true[18][2] = 4;
    hash_true[19][0] = 3;
    hash_true[19][1] = 4;
    hash_true[20][0] = 3;
    hash_true[20][1] = 4;
    hash_true[22][0] = 2;
    hash_true[23][0] = 2;
    hash_true[24][0] = 2;
    hash_true[25][0] = 2;
    hash_true[25][1] = 3;
    hash_true[26][0] = 3;
    hash_true[27][0] = 3;
    hash_true[29][0] = 2;
    hash_true[30][0] = 2;
    hash_true[31][0] = 2;
    hash_true[32][0] = 2;
    hash_true[32][1] = 3;
    hash_true[33][0] = 3;
    hash_true[34][0] = 3;
    hash_true[36][0] = 0;
    hash_true[36][1] = 2;
    hash_true[37][0] = 0;
    hash_true[37][1] = 2;
    hash_true[38][0] = 0;
    hash_true[38][1] = 1;
    hash_true[38][2] = 2;
    hash_true[39][0] = 0;
    hash_true[39][1] = 1;
    hash_true[39][2] = 2;
    hash_true[39][3] = 3;
    hash_true[40][0] = 3;
    hash_true[41][0] = 3;
    hash_true[43][0] = 0;
    hash_true[44][0] = 0;
    hash_true[45][0] = 0;
    hash_true[45][1] = 1;
    hash_true[46][0] = 1;
    hash_true[46][1] = 3;
    hash_true[47][0] = 1;
    hash_true[47][1] = 3;
    hash_true[48][0] = 3;
    hash_true[51][0] = 0;
    hash_true[51][1] = 1;
    hash_true[52][0] = 0;
    hash_true[52][1] = 1;
    hash_true[53][0] = 1;
    hash_true[54][0] = 1;
    hash_true[54][1] = 3;
    hash_true[55][0] = 1;
    hash_true[55][1] = 3;
    
//    for (int j = hash_height - 1; j >= 0; j--) {
//        for (uint i = 0; i < hash_width; i++) {
//            printf("%i\t\t", histogram[get_key(i, j, hash_width)]);
//        }
//        printf("\n");
//    }
//    printf("\n");
//    
//    for (int j = hash_height - 1; j >= 0; j--) {
//        for (uint i = 0; i < hash_width; i++) {
//            if (histogram[get_key(i, j, hash_width)] == 0) {
//                printf("x\t\t");
//                continue;
//            }
//            printf("%i", hash[get_key(i, j, hash_width)][0]);
//            for(uint k = 1; k < histogram[get_key(i, j, hash_width)]; k++) {
//                printf(",%i", hash[get_key(i, j, hash_width)][k]);
//            }
//            printf("\t\t");
//        }
//        printf("\n");
//    }
//    printf("\n");
    
    for (uint n = 0; n < hash_width * hash_height; n++) {
        for (uint m = 0; m < histogram_true[n]; m++) {
            if (hash[n][m] != hash_true[n][m]) {
                success = 0;
                printf("Error in cell (%u, %u) key (%u, %u): expected %i, but found %i.\n",
                    key_to_i(n, hash_width), key_to_j(n, hash_width), n, m, hash_true[n][m], 
                    hash[n][m]);
            }
        }
    }
    
    free(hash_true);
    free(histogram_true);
    free(hash);
    free(histogram);
    
    return success;
}
