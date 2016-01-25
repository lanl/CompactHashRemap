/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */

#include <stdio.h>

#include "simplehash/simplehash.h"
#include "hierarchical_remap.h"
#include "meshgen/meshgen.h"
//
//#define DETAILED_TIMING
#ifdef DETAILED_TIMING
#include "timer.h"
#endif

#ifndef DEBUG
#define DEBUG 0
#endif

// Private functions to this routine
double avg_sub_cells_h (cell_list icells, uint i, uint j, uint lev, int **h_hash, uint ibasesize);
double avg_sub_cells_h_compact (cell_list icells, uint i, uint j, uint lev, intintHash_Table** h_hashTable, uint ibasesize);

double avg_sub_cells_h (cell_list icells, uint i, uint j, uint lev, int **h_hash, uint ibasesize) {

    int probe;
    double sum = 0.0;

    lev++;
    i *= 2;
    j *= 2;
    uint istride = ibasesize*two_to_the(lev);
    uint key = j*istride + i;

    uint key_new[4];
    key_new[0] = key;
    key_new[1] = key + 1;
    key_new[2] = key + istride;
    key_new[3] = key + istride + 1;
    
    for (int ic = 0; ic < 4; ic++){
        key = key_new[ic];
        probe = h_hash[lev][key];
        if (probe >= 0) {
            sum += icells.values[probe];
        } else {
            i = key % istride;
            j = key / istride;
            sum += avg_sub_cells_h(icells, i, j, lev, h_hash, ibasesize);
        }
    }

    return sum/4.0;
}

double avg_sub_cells_h_compact (cell_list icells, uint i, uint j, uint lev, intintHash_Table** h_hashTable, uint ibasesize) {

    int probe;
    double sum = 0.0;

    lev++;
    i *= 2;
    j *= 2;
    uint istride = ibasesize*two_to_the(lev);
    uint key = j*istride + i;
    
    uint key_new[4];
    key_new[0] = key;
    key_new[1] = key + 1;
    key_new[2] = key + istride;
    key_new[3] = key + istride + 1;
    
    for (int ic = 0; ic < 4; ic++){
        key = key_new[ic];
        //probe = read_hash(key, h_hash[lev]);
        intintHash_QuerySingle(h_hashTable[lev], key, &probe);
        if (probe >= 0) {
            sum += icells.values[probe];
        } else {
            i = key % istride;
            j = key / istride;
            sum += avg_sub_cells_h_compact(icells, i, j, lev, h_hashTable, ibasesize);
        }
    }

    return sum/4.0;
}

void h_remap (cell_list icells, cell_list ocells) {
    
    int** h_hash = (int **) malloc((icells.levmax+1)*sizeof(uint *));
    
    //initialize 2d array
    //worth checking for an empty level?
    for (uint i = 0; i <= icells.levmax; i++) {
        size_t hash_size = icells.ibasesize*two_to_the(i)*icells.jbasesize*two_to_the(i);
        h_hash[i] = (int *) malloc(hash_size*sizeof(uint));
        //memset(h_hash[i], -2, hash_size*sizeof(uint));
    }
    
    //place the cells and their breadcrumbs 
    for (uint n = 0; n < icells.ncells; n++) {
        uint i = icells.i[n];
        uint j = icells.j[n];
        int lev = icells.level[n];
        uint key = j * icells.ibasesize*two_to_the(lev) + i;
        h_hash[lev][key] = n;
        
        while (i%2 == 0 && j%2 == 0 && lev > 0) {
            i /= 2;
            j /= 2;
            lev--;
            key = j * icells.ibasesize*two_to_the(lev) + i;
            h_hash[lev][key] = -1;
        }
    }
    
    for (uint n = 0; n < ocells.ncells; n++) {
        //printf("cell# %d\n", n);
        //printCell(ocells[n]);
        uint oi = ocells.i[n];
        uint oj = ocells.j[n];
        uint olev = ocells.level[n];
        
        //printf("probe_lev: %d\n", probe_lev);
        //printf("key: %u\n", key);
        //printf("probe: %d\n", probe);
        
        int probe = -1;
        // loop until either we find a valid hash value or the probelev is the output cell level
        for (uint probe_lev = 0; probe < 0 && probe_lev <= olev; probe_lev++){
            //uint key = translate_cell(oi, oj, olev, probe_lev);
            int levdiff = olev - probe_lev;
            uint key = (oj >> levdiff)*icells.ibasesize*two_to_the(probe_lev) + (oi >> levdiff);
            probe = h_hash[probe_lev][key];
        }
        if (probe >= 0) {
            ocells.values[n] = icells.values[probe];
        } else {
            ocells.values[n] = avg_sub_cells_h (icells, oi, oj, olev, h_hash, icells.ibasesize);
        }
        
        //printf ("#%d:\t", n);
        //printCell(ocells[n]);
    }
    
    for (uint i = 0; i <= icells.levmax; i++) {
        free(h_hash[i]);
    }
    free(h_hash);
    
}

//#define HASH_TYPE HASH_ALL_C_HASHES
#define HASH_TYPE (LCG_QUADRATIC_OPEN_COMPACT_HASH_ID | IDENTITY_SENTINEL_PERFECT_HASH_ID)
//#define HASH_TYPE LCG_QUADRATIC_OPEN_COMPACT_HASH_ID
#define HASH_LOAD_FACTOR 0.3333333

void h_remap_compact (cell_list icells, cell_list ocells, intintHash_Factory *factory) {
    
#ifdef DETAILED_TIMING
    struct timeval timer;

    cpu_timer_start(&timer);
#endif

    intintHash_Table** h_hashTable = (intintHash_Table **) malloc((icells.levmax+1)*sizeof(intintHash_Table *));
    
    uint *num_at_level = (uint *)malloc((icells.levmax+1)*sizeof(uint));
    uint *h_hashtype;
    if (DEBUG >= 2) {
       h_hashtype = (uint *)malloc((icells.levmax+1)*sizeof(uint));
    }

    for (uint i = 0; i <= icells.levmax; i++) {
       num_at_level[i] = 0;
    }

    for (uint n = 0; n < icells.ncells; n++) {
       uint lev = icells.level[n];
       num_at_level[lev]++;
    }

    // lev must be int (not uint) to allow -1 for exit
    for (int lev = icells.levmax-1; lev >= 0; lev--) {
       num_at_level[lev] += num_at_level[lev+1]/4;
    }

    //initialize 2d array
    //worth checking for an empty level?
    for (uint i = 0; i <= icells.levmax; i++) {
        size_t hash_size = icells.ibasesize*two_to_the(i)*icells.jbasesize*two_to_the(i);
        h_hashTable[i] = intintHash_CreateTable(factory, HASH_TYPE, hash_size, num_at_level[i], HASH_LOAD_FACTOR);
        if (DEBUG >= 2) {
           h_hashtype[i] = intintHash_GetTableType(h_hashTable[i]);
           if (h_hashtype[i] == IDENTITY_PERFECT_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"IDENTITY_PERFECT_HASH_ID");
           } else if (h_hashtype[i] == IDENTITY_SENTINEL_PERFECT_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"IDENTITY_SENTINEL_PERFECT_HASH_ID");
           } else if (h_hashtype[i] == LCG_LINEAR_OPEN_COMPACT_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"LCG_LINEAR_OPEN_COMPACT_HASH_ID");
           } else if (h_hashtype[i] == LCG_QUADRATIC_OPEN_COMPACT_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"LCG_QUADRATIC_OPEN_COMPACT_HASH_ID");
           }
        }

        //Empty Hash Table
        intintHash_SetupTable(h_hashTable[i]);
        //h_hash[i] = compact_hash_init(icells.ncells, hash_size, 1, 0);
        //memset(h_hash[i], -2, hash_size*sizeof(uint));
    }

    free(num_at_level);
    if (DEBUG >= 2) {
       free(h_hashtype);
    }

#ifdef DETAILED_TIMING
    double setup_time = cpu_timer_stop(timer);
    printf("setup time is %8.4f ms\n",setup_time*1000.0);
    cpu_timer_start(&timer);
#endif

    
    //place the cells and their breadcrumbs 
    for (uint n = 0; n < icells.ncells; n++) {
        uint i = icells.i[n];
        uint j = icells.j[n];
        int lev = icells.level[n];

        uint key = j * icells.ibasesize*two_to_the(lev) + i;
        intintHash_InsertSingle(h_hashTable[lev], key, n);
        //write_hash(n, key, h_hash[lev]);

        while (i%2 == 0 && j%2 == 0 && lev > 0) {
            i /= 2;
            j /= 2;
            lev--;
            key = j * icells.ibasesize*two_to_the(lev) + i;
            intintHash_InsertSingle(h_hashTable[lev], key, -1);
            //write_hash(-1, key, h_hash[lev]);
        }
    }
    
#ifdef DETAILED_TIMING
    double write_time = cpu_timer_stop(timer);
    printf("write time is %8.4f ms\n",write_time*1000.0);
    cpu_timer_start(&timer);
#endif

    for (uint n = 0; n < ocells.ncells; n++) {
        uint oi = ocells.i[n];
        uint oj = ocells.j[n];
        uint olev = ocells.level[n];
        
        int probe = -1;
        for (uint probe_lev = 0; probe < 0 && probe_lev <= olev; probe_lev++){
            int levdiff = olev - probe_lev;
            uint key = (oj >> levdiff)*icells.ibasesize*two_to_the(probe_lev) + (oi >> levdiff);
            //probe = read_hash(key, h_hash[probe_lev]);
            intintHash_QuerySingle(h_hashTable[probe_lev], key, &probe);
        }

        if (probe >= 0) {
            ocells.values[n] = icells.values[probe];
        } else {
            ocells.values[n] = avg_sub_cells_h_compact (icells, oi, oj, olev, h_hashTable, icells.ibasesize);
        }
    }
    
#ifdef DETAILED_TIMING
    double read_time = cpu_timer_stop(timer);
    printf("read time is %8.4f ms\n",read_time*1000.0);
    cpu_timer_start(&timer);
#endif

    for (uint i = 0; i <= icells.levmax; i++) {
        //free(h_hash[i]);
        intintHash_DestroyTable(h_hashTable[i]);
    }
    free(h_hashTable);

#ifdef DETAILED_TIMING
    double cleanup_time = cpu_timer_stop(timer);
    printf("cleanup time is %8.4f ms\n",cleanup_time*1000.0);
#endif
}

#ifdef _OPENMP
void h_remap_openMP (cell_list icells, cell_list ocells) {
    
    int** h_hash = (int **) malloc((icells.levmax+1)*sizeof(uint *));

    //initialize 2d array
    //worth checking for an empty level?
    for (uint i = 0; i <= icells.levmax; i++) {
        size_t hash_size = icells.ibasesize*two_to_the(i)*icells.jbasesize*two_to_the(i);
        h_hash[i] = (int *) malloc(hash_size*sizeof(uint));
        //memset(h_hash[i], -2, hash_size*sizeof(uint));
    }

#pragma omp parallel default(none)  shared (h_hash, icells, ocells)
    {
        uint ilength = icells.ncells;
        uint olength = ocells.ncells;
        uint ibasesize = icells.ibasesize;

    //place the cells and their breadcrumbs
#pragma omp for
    for (uint n = 0; n < ilength; n++) {
        uint i = icells.i[n];
        uint j = icells.j[n];
        int lev = icells.level[n];
        uint key = j * ibasesize*two_to_the(lev) + i;
        h_hash[lev][key] = n;

        while (i%2 == 0 && j%2 == 0 && lev > 0) {
            i /= 2;
            j /= 2;
            lev--;
            key = j * ibasesize*two_to_the(lev) + i;
            h_hash[lev][key] = -1;
        }
    }
    
#pragma omp for
        for (uint n = 0; n < olength; n++) {
            uint oi = ocells.i[n];
            uint oj = ocells.j[n];
            uint olev = ocells.level[n];
        
            int probe = -1;
            // loop until either we find a valid hash value or the probelev is the output cell level
            for (uint probe_lev = 0; probe == -1 && probe_lev <= olev; probe_lev++){
                int levdiff = olev - probe_lev;
                uint key = (oj >> levdiff)*icells.ibasesize*two_to_the(probe_lev) + (oi >> levdiff);
                probe = h_hash[probe_lev][key];
            }

            if (probe >= 0) {
                ocells.values[n] = icells.values[probe];
            } else {
                ocells.values[n] = avg_sub_cells_h (icells, oi, oj, olev, h_hash, icells.ibasesize);
            }
        }
    }
    
    for (uint i = 0; i <= icells.levmax; i++) {
        free(h_hash[i]);
    }
    free(h_hash);
}

//#define HASH_OPENMP_TYPE HASH_ALL_OPENMP_HASHES
#define HASH_OPENMP_TYPE (LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID | IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID)
//#define HASH_OPENMP_TYPE LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID 
void h_remap_compact_openMP (cell_list icells, cell_list ocells, intintHash_Factory *factory) {
    
#ifdef DETAILED_TIMING
    struct timeval timer;

    cpu_timer_start(&timer);
#endif

    intintHash_Table** h_hashTable = (intintHash_Table **) malloc((icells.levmax+1)*sizeof(intintHash_Table *));
    
    uint *num_at_level = (uint *)malloc((icells.levmax+1)*sizeof(uint));
    uint *h_hashtype;
    if (DEBUG >= 2) {
       h_hashtype = (uint *)malloc((icells.levmax+1)*sizeof(uint));
    }

    for (uint i = 0; i <= icells.levmax; i++) {
       num_at_level[i] = 0;
    }

//#pragma omp parallel for reduction(+:num_at_level)
    for (uint n = 0; n < icells.ncells; n++) {
       uint lev = icells.level[n];
       num_at_level[lev]++;
    }

    // lev must be int (not uint) to allow -1 for exit
    for (int lev = icells.levmax-1; lev >= 0; lev--) {
       num_at_level[lev] += num_at_level[lev+1]/4;
    }

    //initialize 2d array
    //worth checking for an empty level?
    for (uint i = 0; i <= icells.levmax; i++) {
        size_t hash_size = icells.ibasesize*two_to_the(i)*icells.jbasesize*two_to_the(i);
        //h_hashTable[i] = intintHash_CreateTable(factory, LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID, hash_size, num_at_level[i], HASH_LOAD_FACTOR);
        h_hashTable[i] = intintHash_CreateTable(factory, HASH_OPENMP_TYPE, hash_size, num_at_level[i], HASH_LOAD_FACTOR);
        if (DEBUG >= 2) {
           h_hashtype[i] = intintHash_GetTableType(h_hashTable[i]);
           if (h_hashtype[i] == IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID");
           } else if (h_hashtype[i] == LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID");
           } else if (h_hashtype[i] == LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID");
           } else {
              printf("Type of hash for lev %d is unknown type %d\n",i,h_hashtype[i]);
           }
        }

        //Empty Hash Table
        intintHash_SetupTable(h_hashTable[i]);
        //h_hash[i] = compact_hash_init(icells.ncells, hash_size, 1, 0);
        //memset(h_hash[i], -2, hash_size*sizeof(uint));
    }

    free(num_at_level);
    if (DEBUG >= 2) {
       free(h_hashtype);
    }
    
#ifdef DETAILED_TIMING
    double setup_time = cpu_timer_stop(timer);
    printf("setup time is %8.4f ms\n",setup_time*1000.0);
    cpu_timer_start(&timer);
#endif

    //place the cells and their breadcrumbs 
#ifdef DETAILED_TIMING
#pragma omp parallel default(none) shared(icells, ocells, h_hashTable) firstprivate(timer)
#else
#pragma omp parallel default(none) shared(icells, ocells, h_hashTable)
#endif
    {
        uint ilength = icells.ncells;
        uint olength = ocells.ncells;

#pragma omp for
         for (uint n = 0; n < ilength; n++) {
             uint i = icells.i[n];
             uint j = icells.j[n];
             int lev = icells.level[n];

             uint key = j * icells.ibasesize*two_to_the(lev) + i;
             intintHash_InsertSingle(h_hashTable[lev], key, n);
             //write_hash(n, key, h_hash[lev]);

             while (i%2 == 0 && j%2 == 0 && lev > 0) {
                 i /= 2;
                 j /= 2;
                 lev--;
                 key = j * icells.ibasesize*two_to_the(lev) + i;
                 intintHash_InsertSingle(h_hashTable[lev], key, -1);
                 //write_hash(-1, key, h_hash[lev]);
             }
         }
    
#ifdef DETAILED_TIMING
    double write_time = cpu_timer_stop(timer);
#pragma omp master
    {
       printf("write time is %8.4f ms\n",write_time*1000.0);
    }
    cpu_timer_start(&timer);
#endif

#pragma omp for
         for (uint n = 0; n < olength; n++) {
             uint oi = ocells.i[n];
             uint oj = ocells.j[n];
             uint olev = ocells.level[n];
        
             int probe = -1;
             for (uint probe_lev = 0; probe < 0 && probe_lev <= olev; probe_lev++){
                 int levdiff = olev - probe_lev;
                 uint key = (oj >> levdiff)*icells.ibasesize*two_to_the(probe_lev) + (oi >> levdiff);
                 //probe = read_hash(key, h_hash[probe_lev]);
                 intintHash_QuerySingle(h_hashTable[probe_lev], key, &probe);
             }

             if (probe >= 0) {
                 ocells.values[n] = icells.values[probe];
             } else {
                 ocells.values[n] = avg_sub_cells_h_compact (icells, oi, oj, olev, h_hashTable, icells.ibasesize);
             }
        }
    } // end omp parallel
    
#ifdef DETAILED_TIMING
    double read_time = cpu_timer_stop(timer);
    printf("read time is %8.4f ms\n",read_time*1000.0);
    cpu_timer_start(&timer);
#endif

    for (uint i = 0; i <= icells.levmax; i++) {
        //free(h_hash[i]);
        intintHash_DestroyTable(h_hashTable[i]);
    }
    free(h_hashTable);
    
#ifdef DETAILED_TIMING
    double cleanup_time = cpu_timer_stop(timer);
    printf("cleanup time is %8.4f ms\n",cleanup_time*1000.0);
#endif
}

#endif

