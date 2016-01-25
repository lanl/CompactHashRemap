
#include "unstructured_types.h"

int load_mesh(uint *num_nodes, node_list *nodes, uint *num_faces, face_list *faces, uint *num_cells, cell_list *cells, const char *file_name);
int test_mesh_loader ();
void destroy_nodes (node_list nodes);
void destroy_faces (face_list faces, uint num_faces);
void destroy_cells (cell_list cells, uint num_cells);
void destroy_mesh (node_list nodes, face_list faces, uint num_faces, cell_list cells, uint num_cells);
