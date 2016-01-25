#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "unstructured_parse.h"

int main(int argc, char** argv){
    int repetitions = 1;
    char test = 0;
    const char *mesh_file_name;
    char file_found = 0;
    for (int i = 1; i < argc; i++){
        if (!strcmp(argv[i], "-c")){
            i++;
            repetitions = atoi(argv[i]);
        } else
        if (!strcmp(argv[i], "-test")){
            test = 1;
        } else {
            file_found = 1;
            mesh_file_name = argv[i];
        }
    }
    if (!file_found){
        printf ("No file specified - using example\n");
        mesh_file_name =  "example.x3d.00001";
    }    
    
    
    for (int i = 0; i < repetitions; i++){
        node_list nodes;
        face_list faces;
        cell_list cells;
        
        uint num_nodes, num_faces, num_cells;
        
        load_mesh(&num_nodes, &nodes, &num_faces, &faces, &num_cells, &cells, mesh_file_name);
        
        if (test){
            int errors = test_mesh_loader();
            printf ("Errors: %i\n", errors);
        }
        destroy_mesh (nodes, faces, num_faces, cells, num_cells);
    }
    
    
    return 0;
}
