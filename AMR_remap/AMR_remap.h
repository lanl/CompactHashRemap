/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */


#ifndef AMR_REMAP_H
#define AMR_REMAP_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdbool.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <unistd.h>
#include <string.h>

#include "meshgen/meshgen.h"

#ifndef _HASH_H

typedef unsigned int uint;
#endif

void check_output(const char *string, uint olength, double *output_val, double *val_test_answer);

#endif

