#ifndef UNSTRUCTURED_TYPES_H
#define UNSTRUCTURED_TYPES_H

typedef unsigned int uint;

typedef struct {
    double* i;
    double* j;
    double* k;
} node_list;

typedef struct {
    uint *num_nodes;
    uint **nodes;
    uint *processor;
    uint *neighbor_processor;
    uint *neighbor_face;
} face_list;

typedef struct {
    uint *num_faces;
    uint **faces;
} cell_list;

#endif
