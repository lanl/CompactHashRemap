
#include "unstructured_types.h"
#include "unstructured_parse.h"
#include "unstructured_write.h"
#include "unstructured_read.h"

void get_hash_range (node_list nodes, uint num_nodes, double *dx, double *dy);
uint **unstructured_hash_remap (const char *cells_in_file, const char *cells_out_file,
    double bin_size);
int test_unstructured_hash_remap();
