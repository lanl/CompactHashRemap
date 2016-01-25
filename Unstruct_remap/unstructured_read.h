
#include "unstructured_types.h"

void unstructured_read_list (uint num_cells_target, uint hash_width, 
    double bin_size, node_list nodes, face_list faces, cell_list cells, 
    int **hash, uint *histogram, uint **candidates);
void get_bounding_box(uint cell_id, double bin_size, cell_list cells, node_list nodes,
    face_list faces, uint *min_x, uint *max_x, uint *min_y, uint *max_y,
    uint *dx, uint *dy);
void get_local_hash (int *local_hash, uint cell_id, double bin_size, cell_list cells,
    face_list faces, node_list nodes, uint num_faces, uint min_x, uint dx, uint min_y, 
    uint dy);
void write_edge (int *hash, uint hash_width, double bin_size, node_list nodes, 
    face_list faces, uint edge, uint cell, uint min_x, uint min_y);
void fill_polygon (uint cell_id, int *hash, uint hash_width, uint hash_height);
uint get_key(uint i, uint j, uint hash_width);
uint key_to_i(uint key, uint hash_width);
uint key_to_j(uint key, uint hash_width);
int test_read();
