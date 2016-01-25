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

#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include "unstructured_parse.h"

int load_mesh(uint *num_nodes, node_list *nodes, uint *num_faces, face_list *faces, uint *num_cells, cell_list *cells, const char *file_name){
    FILE *file = fopen(file_name, "rt");
    // XXX if lines are longer than 256 characters something strange might happen
    char line[256];
    
    // 0 means nothing to read. 1 means header. 2 means nodes. 3 means faces. 4 means cells.
    char element_flag = 0;
    
    // Keeps track of how far after the start of the most recent section we are
    int index = 0;
    
    while(fgets(line, 256, file) != NULL){
    	if (!strncmp(line,"header",6)){
    	    element_flag = 1;
    	    index = 0;
    	}else if (!strncmp(line,"nodes",5)){
    	    element_flag = 2;
    	    index = 0;
    	}else if (!strncmp(line,"faces",5)){
    	    element_flag = 3;
    	    index = 0;
    	}else if (!strncmp(line,"cells",5)){
    	    element_flag = 4;
    	    index = 0;
    	}else if (!strncmp(line,"end_",4)){
    	    element_flag = 0;
    	}else if (element_flag){
    	    index++;
    	    // 1 means dealing with the header
    	    if (element_flag==1){
    	        // header line 4 is the nodes
    	        if (index == 4){
    	            // The first call is a setup and returns the first value
            	    strtok (line, " ");
            	    
            	    // this will be the second value on the line
            	    *num_nodes = atoi(strtok(NULL, " "));
            	    
            	    nodes->i = (double*) malloc (*num_nodes*sizeof(double));
            	    nodes->j = (double*) malloc (*num_nodes*sizeof(double));
            	    nodes->k = (double*) malloc (*num_nodes*sizeof(double));
    	        } else if (index == 5) {
    	            // The first call is a setup and returns the first value
            	    strtok (line, " ");
            	    
            	    // this will be the second value on the line
            	    *num_faces = atoi(strtok(NULL, " "));
            	    
            	    faces->num_nodes = (uint *) malloc (*num_faces*sizeof(uint));
            	    faces->nodes = (uint **) malloc (*num_faces*sizeof(uint*));
            	    faces->processor = (uint *) malloc (*num_faces*sizeof(uint));
            	    faces->neighbor_processor = (uint *) malloc (*num_faces*sizeof(uint));
            	    faces->neighbor_face = (uint *) malloc (*num_faces*sizeof(uint));
            	    
    	        } else
    	        // header line 6 is the cells
    	        if (index == 6){
    	            // The first call is a setup and returns the first value
            	    strtok (line, " ");
            	    
            	    // this will be the second value on the line
            	    *num_cells = atoi(strtok(NULL, " "));
            	    
            	    cells->num_faces = (uint *) malloc (*num_cells*sizeof(uint));
            	    cells->faces = (uint **) malloc (*num_cells*sizeof(uint*));
    	        }
    	    } else
    	    // 2 = nodes
    	    if (element_flag == 2){
    	        // The first call is a setup and returns the first value
    	        // This is the index and we should already know
            	strtok (line, " ");
            	
            	char *i_string = strtok(NULL, " ");
            	// The strings occasionally have a return character at the end.
            	// We need to remove it because it will crash the parser, if it exists.
            	if (i_string[strlen(i_string)-1] == '\n'){
            	    i_string[strlen(i_string)-1] = '\0';
            	}
            	char *j_string = strtok(NULL, " ");
            	// The strings occasionally have a return character at the end.
            	// We need to remove it because it will crash the parser, if it exists.
            	if (j_string[strlen(j_string)-1] == '\n'){
            	    j_string[strlen(j_string)-1] = '\0';
            	}
            	char *k_string = strtok(NULL, " ");
            	// The strings occasionally have a return character at the end.
            	// We need to remove it because it will crash the parser, if it exists.
            	if (k_string[strlen(k_string)-1] == '\n'){
            	    k_string[strlen(k_string)-1] = '\0';
            	}
            	
            	// subtract 1 from index because the format starts counting at 1 not 0
            	sscanf(i_string,"%lg",&(nodes->i[index-1]));
                sscanf(j_string,"%lg",&(nodes->j[index-1]));
            	sscanf(k_string,"%lg",&(nodes->k[index-1]));
            	
       	    } else
    	    // 3 = faces
    	    if (element_flag == 3){
    	        // The first call is a setup and returns the first value
    	        // This is the index and we should already know
            	strtok (line, " ");
            	
            	int n_nodes = atoi(strtok(NULL, " "));
            	// subtract 1 from index because the format starts counting at 1 not 0
            	faces->num_nodes[index-1] = n_nodes;
            	faces->nodes[index-1] = (uint *) malloc (n_nodes*sizeof(uint));
            	for (int i = 0; i < n_nodes; i++){
            	    faces->nodes[index-1][i] = atoi(strtok(NULL, " ")) - 1;
            	}
            	faces->processor[index-1] = atoi(strtok(NULL, " "));
            	faces->neighbor_processor[index-1] = atoi(strtok(NULL, " "));
            	faces->neighbor_face[index-1] = atoi(strtok(NULL, " "));
    	    } else if (element_flag == 4) {
    	        // The first call is a setup and returns the first value
    	        // This is the index and we should already know
            	strtok (line, " ");
            	int n_faces = atoi(strtok(NULL, " "));
            	// subtract 1 from index because the format starts counting at 1 not 0
            	cells->num_faces[index-1] = n_faces;
            	cells->faces[index-1] = (uint *) malloc (n_faces*sizeof(uint));
            	for (int i = 0; i < n_faces; i++) {
            	    cells->faces[index-1][i] = atoi(strtok(NULL, " ")) - 1;
            	}
    	    }
    	}
    }
    fclose (file);
    return 0;
}

// The test compares what is parsed from the example file to what is actually hand copied from that file
int test_mesh_loader(){
    // stores the number of errors
    int error = 0;


    /* make the test nodes */
    node_list test_nodes;
    uint real_num_nodes = 20;
    
    double test_nodes_i[20] = { 0, 7.499999999999999e-01, 0, 0, 4.999999999999999e-01, 2.500000000000000e-01, 0, 9.999999999999999e-01, 4.999999999999999e-01, 0, 3.749999999999999e-01, 0, 1.666666666666666e-01, 3.333333333333333e-01, 0, 3.333333333333333e-01, 6.666666666666665e-01, 2.207790270515294e-01, 2.748913874363922e-01, 5.221529395202574e-01 };
    test_nodes.i = &test_nodes_i[0];
    
    double test_nodes_j[20] = { 1.300000000000000e+00, 1.299038105676658e+00, 1.156666666666667e+00, 1.443333333333333e+00, 8.660254037844387e-01, 8.680127018922194e-01, 8.700000000000000e-01, 1.732050807568877e+00, 1.731025403784439e+00, 1.730000000000000e+00, 1.299519052838329e+00, 8.700000000000000e-01, 8.686751345948129e-01, 8.673502691896258e-01, 1.730000000000000e+00, 1.730683602522959e+00, 1.731367205045918e+00, 1.152977306437331e+00, 1.437279478279848e+00, 1.367410615103820e+00};
    test_nodes.j = &test_nodes_j[0];
    
    double test_nodes_k[20] = {0};
    test_nodes.k = &test_nodes_k[0];
    
    
    /* make the test faces list */
    face_list test_faces;
    uint real_num_faces = 45;
    
    test_faces.num_nodes = (uint *) malloc (real_num_faces*sizeof(uint));
    for (uint i = 0; i < real_num_faces; i++){
        test_faces.num_nodes[i] = 2;
    }
    
    uint test_faces_nodes[45][2] = { {12,13}, {13,18}, {18,3}, {3,12}, {3,18}, {18,19}, {19,4}, {4,3}, {4,19}, {19,16}, {16,15}, {15,4}, {13,14}, {14,20}, {20,19}, {19,18}, {18,13}, {19,20}, {20,17}, {17,16}, {16,19}, {14,5}, {5,2}, {2,20}, {20,14}, {20,2}, {2,8}, {8,17}, {17,20}, {5,2}, {2,11}, {11,6}, {6,5}, {2,8}, {8,9}, {9,11}, {11,2}, {6,11}, {11,1}, {1,7}, {7,6}, {11,9}, {9,10}, {10,1}, {1,11} };
    
    test_faces.nodes = (uint **) malloc (45*sizeof(uint*));
    for (int i = 0; i < 45; i++){
        test_faces.nodes[i] = test_faces_nodes[i];
    }                                           
    test_faces.nodes[0] = test_faces_nodes[0];
    test_faces.nodes[1] = test_faces_nodes[1];
    
    test_faces.processor = (uint *) malloc (real_num_faces*sizeof(uint));
    // all of our faces use processor 1
    for (int i = 0; i < 45; i++){
        test_faces.processor[i] = 1;
    }
    
    uint test_neighbor_processor[45] = {0,1,1,0,1,1,1,0,1,1,0,0,0,1,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,1,1,1,0,1,0,1,1,1,1,0,0,1,0,0,1};
    
    test_faces.neighbor_processor = &test_neighbor_processor[0];
    
    uint test_neighbor_face[45] = {0,17,5,0,3,16,9,0,7,21,0,0,0,25,18,6,2,15,29,0,10,0,30,26,14,24,34,0,19,23,37,38,0,27,0,42,31,32,45,0,0,36,0,0,39};
    
    test_faces.neighbor_face = &test_neighbor_face[0];
    
    /* make the test cells list */
    cell_list test_cells;
    uint real_num_cells = 11;
    
    test_cells.num_faces = (uint *) malloc (real_num_faces*sizeof(uint));
    // all of the cells except one have four faces, so we set them all then go back and fix it
    for (int i = 0; i<11; i++){
        test_cells.num_faces[i] = 4;
    }
    test_cells.num_faces[3] = 5;
    
    // all of the cells' faces except for the 4 length one; same deal as above
    uint test_cell_faces[11][4] = { {1,2,3,4}, {5,6,7,8}, {9,10,11,12}, {0,0,0,0}, {18,19,20,21}, {22,23,24,25}, {26,27,28,29}, {30,31,32,33}, {34,35,36,37}, {38,39,40,41}, {42,43,44,45} };
    test_cells.faces = (uint **) malloc (11*sizeof(uint*));
    for (int i = 0; i < 11; i++){
        if (i!=3){
            test_cells.faces[i] = test_cell_faces[i];
        } else {
            test_cells.faces[i] = (uint *) malloc (5 * sizeof(uint));
            test_cells.faces[i][0] = 13;
            test_cells.faces[i][1] = 14;
            test_cells.faces[i][2] = 15;
            test_cells.faces[i][3] = 16;
            test_cells.faces[i][4] = 17;
        }
    }
    
    node_list nodes;
    face_list faces;
    cell_list cells;
    uint num_nodes, num_faces, num_cells;
    
    load_mesh(&num_nodes, &nodes, &num_faces, &faces, &num_cells, &cells, "example.x3d.00001");
    
    if (num_nodes!=real_num_nodes){
        printf("Loader test failed - wrong num nodes: expected %i, but got %i\n", real_num_nodes, num_nodes);
        error += 1;
    }
    
    if (num_faces!=real_num_faces){
        printf("Loader test failed - wrong num faces: expected %i, but got %i\n", real_num_faces, num_faces);
        error += 1;
    }
    
    if (num_cells!=real_num_cells){
        printf("Loader test failed - wrong num cells: expected %i, but got %i\n", real_num_cells, num_cells);
        error += 1;
    }
    
    for (uint ni = 0; ni < num_nodes; ni ++){
        if (nodes.i[ni]!=test_nodes.i[ni]){
            printf ("Loader test failed - node list i incorrect at %i: expected %.15f, but got %.15f\n",ni + 1,test_nodes.i[ni],nodes.i[ni]);
            error += 1;
        }
        
        if (nodes.j[ni]!=test_nodes.j[ni]){
            printf ("Loader test failed - node list j incorrect at %i: expected %.15f, but got %.15f\n",ni + 1,test_nodes.j[ni],nodes.j[ni]);
            error += 1;
        }
        
        if (nodes.k[ni]!=test_nodes.k[ni]){
            printf ("Loader test failed - node list k incorrect at %i: expected %.15f, but got %.15f\n",ni + 1,test_nodes.k[ni],nodes.k[ni]);
            error += 1;
        }
    }
    for (uint fi = 0; fi < num_faces; fi ++){
        uint n_nodes = faces.num_nodes[fi];
        if (faces.num_nodes[fi]!=test_faces.num_nodes[fi]){
            printf ("Loader test failed - face list num nodes incorrect at %i: expected %i, but got %i\n",fi + 1,test_faces.num_nodes[fi],faces.num_nodes[fi]);
            error += 1;
        }
        for (uint node = 0; node < n_nodes; node++){
            if (faces.nodes[fi][node] != test_faces.nodes[fi][node]){
                printf ("Loader test failed - face list node incorrect at %i, #%i: expected %i, but got %i\n",fi + 1,node,test_faces.nodes[fi][node],faces.nodes[fi][node]);
                error += 1;
            }
        }
        if (faces.neighbor_processor[fi]!=test_faces.neighbor_processor[fi]){
            printf ("Loader test failed - face list neighbor processor incorrect at %i: expected %i, but got %i\n",fi + 1,test_faces.neighbor_processor[fi],faces.neighbor_processor[fi]);
            error += 1;
        }
        if (faces.neighbor_face[fi]!=test_faces.neighbor_face[fi]){
            printf ("Loader test failed - face list neighbor face incorrect at %i: expected %i, but got %i\n",fi + 1,test_faces.neighbor_face[fi],faces.neighbor_face[fi]);
            error += 1;
        }
    }
    for (uint ci = 0; ci < num_cells; ci ++){
        uint n_faces = cells.num_faces[ci];
        if (cells.num_faces[ci]!=test_cells.num_faces[ci]){
            printf ("Loader test failed - cell list num faces incorrect at %i: expected %i, but got %i\n",ci + 1,test_cells.num_faces[ci],cells.num_faces[ci]);
            error += 1;
        }
        for (uint face = 0; face < n_faces; face++){
            if (cells.faces[ci][face] != test_cells.faces[ci][face]){
                printf ("Loader test failed - cell list face incorrect at %i, #%i: expected %i, but got %i\n",ci + 1,face,test_cells.faces[ci][face],cells.faces[ci][face]);
                error += 1;
            }
        }
    }
    
    destroy_mesh (nodes, faces, num_faces, cells, num_cells);
    
    // many of the arrays for the tests were auto-alloc'd, so we have to manually free the others (cant use destroy);
    free(test_faces.num_nodes);
    free(test_faces.nodes);
    free(test_faces.processor);
    free(test_cells.num_faces);
    free(test_cells.faces[3]);
    free(test_cells.faces);
    
    return error;
}

void destroy_nodes(node_list nodes){
    free(nodes.i);
    free(nodes.j);
    free(nodes.k);
}

void destroy_faces(face_list faces, uint num_faces){
    free(faces.num_nodes);
    for (uint i = 0; i < num_faces; i++){
        free(faces.nodes[i]);
    }
    free(faces.nodes);
    free(faces.processor);
    free(faces.neighbor_processor);
    free(faces.neighbor_face);
}

void destroy_cells (cell_list cells, uint num_cells){
    free(cells.num_faces);
    for (uint i = 0; i < num_cells; i++){
        free (cells.faces[i]);
    }
    free(cells.faces);
}

void destroy_mesh(node_list nodes, face_list faces, uint num_faces, cell_list cells, uint num_cells){
    destroy_nodes (nodes);
    destroy_faces (faces, num_faces);
    destroy_cells (cells, num_cells);
}

