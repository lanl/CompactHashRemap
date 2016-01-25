#include <stdio.h>
#include "unstructured_read.h"

int main () {
    int read_success = test_read();
    printf("Success: %i\n", read_success);
    return 0;
}
