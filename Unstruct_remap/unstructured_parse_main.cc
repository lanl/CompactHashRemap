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
