#include <stdio.h>
#include "unstructured_write.h"

int main () {
    
    int local_hash_success = test_write();
    printf("Local hash success: %i\n", local_hash_success);
    
    return 0;
}
