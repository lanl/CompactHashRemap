/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */

#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include "unstructured_search.h"

uint **unstructured_hash_remap (const char *cells_in_file, const char *cells_out_file,
    double bin_size) {
    uint num_nodes_in, num_faces_in, num_cells_in;
    node_list nodes_in;
    face_list faces_in;
    cell_list cells_in;
    
    load_mesh(&num_nodes_in, &nodes_in, &num_faces_in, &faces_in, &num_cells_in,
        &cells_in, cells_in_file);
        
//    for (uint n = 0; n < num_nodes_in; n++) {
//        printf("node #%u\tx: %f\ty: %f\n", n, nodes_in.i[n], nodes_in.j[n]);
//    }
//    
//    for (uint n = 0; n < num_faces_in; n++) {
//        printf("face #%u\tn1: %u\tn2: %u\n", n, faces_in.nodes[n][0], 
//            faces_in.nodes[n][1]);
//    }
//    
//    for (uint n = 0; n < num_cells_in; n++) {
//        printf("cell #%u\tf1: %u\tf2: %u\tf3: %u", n, cells_in.faces[n][0], 
//            cells_in.faces[n][1], cells_in.faces[n][2]);
//        if (cells_in.num_faces[n] == 4) {
//            printf("\tf4: %u", cells_in.faces[n][3]);
//        }
//        printf("\n");
//    }
    
    uint num_nodes_out, num_faces_out, num_cells_out;
    node_list nodes_out;
    face_list faces_out;
    cell_list cells_out;
    
    load_mesh(&num_nodes_out, &nodes_out, &num_faces_out, &faces_out, 
        &num_cells_out, &cells_out, cells_out_file);
        
    double dx_in, dy_in;
    get_hash_range(nodes_in, num_nodes_in, &dx_in, &dy_in);
    uint hash_width_in = ((uint) (dx_in / bin_size)) + 1;
    uint hash_height_in = ((uint) (dy_in / bin_size)) + 1;
        
    uint *histogram = (uint *) malloc(hash_width_in * hash_height_in * sizeof(uint));
    int **hash = (int **) malloc(hash_width_in * hash_height_in * sizeof(int *));
    
    write_to_hash(num_cells_in, hash_width_in, 
        hash_height_in, bin_size, nodes_in, faces_in, cells_in, hash, histogram);
        
    uint **results = (uint **) malloc(num_cells_out * sizeof(uint *));
    unstructured_read_list(num_cells_out, hash_width_in, bin_size, nodes_out, 
        faces_out, cells_out, hash, histogram, results);
    return results;
}

void get_hash_range (node_list nodes, uint num_nodes, double *dx, double *dy) {
    
    double max_x = 0;
    double max_y = 0;
    
    for (uint n = 0; n < num_nodes; n++) {
        if (nodes.i[n] > max_x) {
            max_x = nodes.i[n];
        }
        if (nodes.j[n] > max_y) {
            max_y = nodes.j[n];
        }
    }
    
    //printf("max_y: %f\tmax_x: %f\n", max_y, max_x);
    *dx = max_x;
    *dy = max_y;
    
    return;
}

int test_unstructured_hash_remap() {
    int success = 1;
    uint **test_result;
    uint num_cells = 7;
    
    test_result = unstructured_hash_remap("test_in.x3d.00001", 
        "test_out.x3d.00001", 1.0);
    
    uint **true_result = (uint **) malloc(num_cells * sizeof(uint *));
    
    true_result[0] = (uint *) malloc(3 * sizeof(uint));
    true_result[1] = (uint *) malloc(4 * sizeof(uint));
    true_result[2] = (uint *) malloc(4 * sizeof(uint));
    true_result[3] = (uint *) malloc(4 * sizeof(uint));
    true_result[4] = (uint *) malloc(2 * sizeof(uint));
    true_result[5] = (uint *) malloc(3 * sizeof(uint));
    true_result[6] = (uint *) malloc(3 * sizeof(uint));
    
    true_result[0][0] = 0;
    true_result[0][1] = 1;
    true_result[0][2] = 2;
    true_result[1][0] = 0;
    true_result[1][1] = 1;
    true_result[1][2] = 2;
    true_result[1][3] = 3;
    true_result[2][0] = 0;
    true_result[2][1] = 1;
    true_result[2][2] = 2;
    true_result[2][3] = 3;
    true_result[3][0] = 0;
    true_result[3][1] = 1;
    true_result[3][2] = 2;
    true_result[3][3] = 3;
    true_result[4][0] = 2;
    true_result[4][1] = 4;
    true_result[5][0] = 2;
    true_result[5][1] = 3;
    true_result[5][2] = 4;
    true_result[6][0] = 2;
    true_result[6][1] = 3;
    true_result[6][2] = 4;
    
    for (uint n = 0; n < 3; n++) {
        if (test_result[0][n] != true_result[0][n]) {
            success = 0;
        }
    }
    for (uint n = 0; n < 4; n++) {
        if (test_result[1][n] != true_result[1][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 4; n++) {
        if (test_result[2][n] != true_result[2][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 4; n++) {
        if (test_result[3][n] != true_result[3][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 2; n++) {
        if (test_result[4][n] != true_result[4][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 3; n++) {
        if (test_result[5][n] != true_result[5][n]) {
            success = 0;
        }
    }
    
    for (uint n = 0; n < 3; n++) {
        if (test_result[6][n] != true_result[6][n]) {
            success = 0;
        }
    }
    return success;
    
}
