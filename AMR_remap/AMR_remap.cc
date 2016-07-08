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

#include "AMR_remap.h"
#include "timer.h"

#include <sys/resource.h>

#ifdef USE_ASSERT
#include <assert.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

enum meshgen_type {
   HIERARCHICAL_MESHGEN = 0,
   SPARSE_MESHGEN,
   ADAPT_MESHGEN };

#include "brute_force_remap.h"
#include "kdtree_remap.h"
#include "full_perfect_remap.h"
#include "singlewrite_remap.h"
#include "hierarchical_remap.h"

#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#include "h_remap_gpu.h"
#include "simplehash/simplehash.h"
#endif

#include "HashFactory/HashFactory.h"

#ifndef DONT_CATCH_SIGNALS
#include <signal.h>

void handle_signal(int signal);
#endif

int TILE_SIZE = 128;
intintHash_Factory *factory;
intintHash_Factory *OpenMPfactory;
intintHash_Factory *CLFactory;

struct timeval timer;

#ifndef DONT_CATCH_SIGNALS
volatile sig_atomic_t quit_loop = 0; 
#endif

int main (int argc, char** argv) {
    const rlim_t kStackSize = 16 * 1024 * 1024;   // min stack size = 16 MB
    struct rlimit rl;
    int result;

    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
    {
        if (rl.rlim_cur < kStackSize)
        {
            rl.rlim_cur = kStackSize;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
            {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
            }
        }
    }

    int run_brute = 1;
    int run_tree = 1;
    int run_tests = 1;
    int plot_file = 0;
    int meshgen = HIERARCHICAL_MESHGEN;
    if (argc < 6) {
       printf("Usage -- ./AMR_remap <i_level_diff> <ilength> <o_level_diff> <olength> <num_rep> [-no-brute,-no-test,-no-tree,-plot-file]\n");
       printf("   or\n");
       printf("Usage -- ./AMR_remap_openMP <i_level_diff> <ilength> <o_level_diff> <olength> <num_rep> [-no-brute,-no-test,-no-tree,-plot-file]\n");
       printf("   Alternate mesh generation usage:\n");
       printf("Usage -- ./AMR_remap <size_base_mesh> <levmax> <refine_threshold> 0 <num_rep> -adapt-meshgen [-no-brute,-no-test,-no-tree,-plot-file]\n");
       printf("   or\n");
       printf("Usage -- ./AMR_remap_openMP <size_base_mesh> <levmax> <refine_threshold> 0 <num_rep> -adapt-meshgen [-no-brute,-no-test,-no-tree,-plot-file]\n");
       exit(-1);
    }
    if (argc>6){
        for (int i = 6; i < argc; i++){
            char* arg = argv[i];
            if (strcmp(arg,"-no-brute")==0){
                run_brute = 0;
            } else
            if (strcmp(arg,"-no-test")==0){
                run_tests = 0;
            } else
            if (strcmp(arg,"-no-tree")==0){
                run_tree = 0;
            } else 
            if (strcmp(arg,"-plot-file")==0){
                plot_file = 1;
            } else 
            if (strcmp(arg,"-sparse-meshgen")==0){
                meshgen = SPARSE_MESHGEN;
            } else 
            if (strcmp(arg,"-adapt-meshgen")==0){
                meshgen = ADAPT_MESHGEN;
            } else
            printf ("Invalid Argument: %s\n", arg);
        }
    }

#ifdef _OPENMP
   int nt = 0;
   int tid = 0;

   nt = omp_get_max_threads();
   tid = omp_get_thread_num();
   if (0 == tid) {
        printf("--- max num openmp threads: %d\n", nt);
   }
#pragma omp parallel
   {
      nt = omp_get_num_threads();
      tid = omp_get_thread_num();

#pragma omp master
      printf("--- num openmp threads in parallel region: %d\n", nt);
   }
#endif

    //uint num_divisions1 = atoi (argv[1]);
    //uint num_divisions2 = atoi (argv[2]);
    uint num_rep = atoi (argv[5]);
    //uint num_divisions1 = 3;
    //uint num_divisions2 = 3;
    //int num_rep = 1;
    //int target_level;

    double full_perfect_remap_time = 0.0;
    double singlewrite_remap_time = 0.0;
    double compact_singlewrite_remap_time = 0.0;
    double h_remap_time = 0.0;
    double compact_h_remap_time = 0.0;
#ifdef _OPENMP
    double full_perfect_remap_openMP_time = 0.0;
    double singlewrite_remap_openMP_time = 0.0;
    double compact_singlewrite_remap_openMP_time = 0.0;
    double h_remap_openMP_time = 0.0;
    double compact_h_remap_openMP_time = 0.0;
#endif
    double brute_force_time = 0.0;
    double kd_tree_time = 0.0;

#ifdef HAVE_OPENCL
    double gpu_full_perfect_remap_time = 0.0; 
    double gpu_singlewrite_remap_time = 0.0; 
    double gpu_compact_singlewrite_remap_time = 0.0; 
    double gpu_hierarchical_remap_time = 0.0; 
    double gpu_compact_hierarchical_remap_time = 0.0; 
#endif
    
    uint ilength = atoi (argv[2]);
    uint i_level_diff = atoi (argv[1]);
    uint olength = atoi (argv[4]);
    uint o_level_diff = atoi (argv[3]);
    uint i_max_level;
    uint o_max_level;
    uint i_min_level;
    uint o_min_level;
    //i for in o for out
    cell_list icells;
    cell_list ocells;
    
#ifdef _OPENMP
    cell_list icells_openmp;
    cell_list ocells_openmp;
#endif

    //ilength = (num_divisions1 * 3) + 1;
    //olength = (num_divisions2 * 3) + 1;
    
    
    //uint i [31] = {0, 2, 3, 2, 6, 14, 15, 14, 15, 6, 7, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 1, 2, 6, 7, 6, 7};
    //uint j [31] = {1, 3, 3, 2, 5, 11, 11, 10, 10, 4, 4, 1, 1, 1, 1, 3, 3, 3, 3, 2, 2, 2, 2, 0, 0, 0, 0, 1, 1, 0, 0};
    //uint level [31] = {1, 2, 2, 2, 3, 4, 4, 4, 4, 3, 3, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2, 2, 3, 3, 3, 3};
    //double val [31] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
    
    //icells = new_cell_list (i, j, level, val);
    
    
    
    //uint i2 [16] = {0, 0, 1, 2, 2, 3, 6, 6, 7, 7, 2, 3, 2, 3, 0, 1};
    //uint j2 [16] = {3, 0, 0, 3, 2, 2, 7, 6, 6, 7, 7, 7, 6, 6, 2, 2};
    //uint level2 [16] = {2, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2};
    //double val2 [16] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    
    //ocells1 = new_cell_list (i2, j2, level2, val2);
    //srand (time(NULL));
    srand(0);
    
    int emptyNeighborValue = -5;
    factory = intintHash_CreateFactory(HASH_ALL_C_HASHES, &emptyNeighborValue, 0, NULL, NULL);
    //OpenMPfactory = intintHash_CreateFactory(LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID, &emptyNeighborValue, 0, NULL, NULL);
    OpenMPfactory = intintHash_CreateFactory(HASH_ALL_OPENMP_HASHES, &emptyNeighborValue, 0, NULL, NULL);
    //factory = intintHash_CreateFactory(HASH_ALL_C_HASHES, NULL, 0, NULL, NULL);

#ifdef HAVE_OPENCL
    int ierr = ezcl_devtype_init(CL_DEVICE_TYPE_GPU);
    if (ierr == EZCL_NODEVICE) {
       printf ("\nNo GPU found - trying CPU\n");
       ierr = ezcl_devtype_init(CL_DEVICE_TYPE_CPU);
       if (ierr == EZCL_NODEVICE) {
          printf("Warning -- no OpenCL device available\n");
       } else {
          printf("CPU successful\n");
       }
    }
    setup_cl();
    hash_lib_init();

    cl_context       context = ezcl_get_context();
    cl_command_queue queue   = ezcl_get_command_queue();
    
    uint lws = TILE_SIZE;
    CLFactory = intintHash_CreateFactory(HASH_ALL_CL_HASHES, &emptyNeighborValue, lws, &context, &queue);
#endif

    
    printf("                      Input mesh                       Output mesh\n");
    printf("           --------------------------------  --------------------------------\n");
    printf("run num    sparsity percent compressibility  sparsity percent compressibility\n");
    printf("-------    ---------------- ---------------  ---------------- ---------------\n");

    size_t sum_ncells = 0;
    size_t save_num_fine_cells=0;

    int mesh_size, levmax;
    
    #ifndef DONT_CATCH_SIGNALS
    signal (SIGINT, handle_signal);
    #endif

    for (uint n = 0; n < num_rep; n++) {
    
        #ifndef DONT_CATCH_SIGNALS
        if (quit_loop){
            printf ("Keyboard Interupt\n");
            num_rep = n;
            break;
        }
        #endif
    
        //if (n % 10 == 0) {
        printf("run #%i", n);
        //}
        //printf("Trying icells construction\n");

        if (meshgen == HIERARCHICAL_MESHGEN) {
           ilength = atoi (argv[2]);
           i_level_diff = atoi (argv[1]);
           olength = atoi (argv[4]);
           o_level_diff = atoi (argv[3]);

           icells = mesh_maker_sparsity(icells, i_level_diff, &ilength, &i_max_level, &i_min_level, 1);
        
           uint num_fine_cells = four_to_the(i_max_level) * icells.ibasesize * icells.ibasesize;
           printf("         %f",(float)(num_fine_cells-icells.ncells)/(float)num_fine_cells*100.0);
           printf("         %f",(float)num_fine_cells/(float)icells.ncells);
           //printf("Trying ocells construction\n");
           ocells = mesh_maker_sparsity(ocells, o_level_diff, &olength, &o_max_level, &o_min_level, 1);
           num_fine_cells = four_to_the(o_max_level) * ocells.ibasesize * ocells.ibasesize;
           printf("         %f",(float)(num_fine_cells-ocells.ncells)/(float)num_fine_cells*100.0);
           printf("         %f",(float)num_fine_cells/(float)ocells.ncells);
           printf("\n");

           icells.ncells    = ilength;
           icells.jbasesize = icells.ibasesize;

           ocells.ncells    = olength;
           ocells.jbasesize = ocells.ibasesize;
           
           if (icells.ibasesize != ocells.ibasesize) {
                printf("Meshes of incompatible size. Exiting.\n");
                exit(0);
           }

#ifdef _OPENMP
           icells_openmp.ncells    = ilength;
           icells_openmp.ibasesize = icells.ibasesize;
           icells_openmp.jbasesize = icells.ibasesize;
           icells_openmp.levmax    = icells.levmax;

           ocells_openmp.ncells    = olength;
           ocells_openmp.ibasesize = icells.ibasesize;
           ocells_openmp.jbasesize = icells.ibasesize;
           ocells_openmp.levmax    = icells.levmax;
#endif
        } else if (meshgen == SPARSE_MESHGEN){
        } else if (meshgen == ADAPT_MESHGEN){
           float threshold = 1.0;
           mesh_size = atoi (argv[1]);
           levmax = atoi (argv[2]);
           threshold = atof (argv[3]);
           int target_ncells = atoi (argv[4]);

           //o_level_diff = atoi (argv[3]);
           //uint num_fine_cells = (double)ilength/i_sparsity;
         
           //mesh_size = (int)(sqrt((double)num_fine_cells)/(double)two_to_the(levmax));
           //printf("DEBUG -- ilength %d num_cells %d four to the levels %d size %d\n",ilength,(int)((double)ilength/i_sparsity),two_to_the(levmax),mesh_size);

           icells = adaptiveMeshConstructorWij(icells, mesh_size, levmax, threshold, target_ncells);
           sum_ncells += icells.ncells;

           ocells = adaptiveMeshConstructorWij(ocells, mesh_size, levmax, threshold, target_ncells);
           sum_ncells += ocells.ncells;

           size_t num_fine_cells = (size_t)mesh_size*(size_t)two_to_the(levmax)*(size_t)mesh_size*(size_t)two_to_the(levmax);;
           save_num_fine_cells = num_fine_cells;

           printf("         %f",(float)(num_fine_cells-icells.ncells)/(float)num_fine_cells*100.0);
           printf("         %f",(float)num_fine_cells/(float)icells.ncells);

           printf("         %f",(float)(num_fine_cells-ocells.ncells)/(float)num_fine_cells*100.0);
           printf("         %f",(float)num_fine_cells/(float)ocells.ncells);
           //printf("\n");
           printf(" fine cells %lu ncells in %u ncells out %u\n",num_fine_cells,icells.ncells,ocells.ncells);

#ifdef _OPENMP
           icells_openmp.ncells    = icells.ncells;
           icells_openmp.ibasesize = mesh_size;
           icells_openmp.jbasesize = mesh_size;
           icells_openmp.levmax    = levmax;

           ocells_openmp.ncells    = ocells.ncells;
           ocells_openmp.ibasesize = mesh_size;
           ocells_openmp.jbasesize = mesh_size;
           ocells_openmp.levmax    = levmax;
#endif

        }
        
        double *val_test         = NULL;
        double *val_test_brute   = NULL;
        double *val_test_kdtree  = NULL;
        double *val_test_perfect = NULL;

        double *val_test_answer  = NULL;
    
        for (uint n = 0; n < icells.ncells; n++) {
            icells.values[n] = rand () % 100;
        }
/*        print_cell_list(icells, ilength);*/
/*        printf("\n\n");*/
/*        print_cell_list(ocells1, olength);*/

        memset(ocells.values,  0xFFFFFFFF, ocells.ncells*sizeof(double));

        // Save original val array to restore later
        val_test = ocells.values;

        ilength = icells.ncells;
        olength = ocells.ncells;

// Brute Force Remap
        
        if (run_brute){
            val_test_brute = (double*)malloc(olength*sizeof(double));
            memset(val_test_brute, 0xFFFFFFFF, olength*sizeof(double));
            ocells.values = val_test_brute;
                
            cpu_timer_start(&timer);
            brute_force_remap (icells, ocells);
            brute_force_time += cpu_timer_stop(timer);
            //print_cell_list (ocells, olength);

            val_test_answer = val_test_brute;
        }
        
// KD Tree Remap

        if (run_tree){
            val_test_kdtree = (double*)malloc(olength*sizeof(double));
            memset(val_test_kdtree, 0, olength*sizeof(double));
            ocells.values = val_test_kdtree;

            cpu_timer_start(&timer);
            remap_kDtree2d(icells, ocells);
            kd_tree_time += cpu_timer_stop(timer);

            if (val_test_answer == NULL) {
                val_test_answer = val_test_kdtree;
            } else if (run_tests) {
                check_output("KD Tree Remap", olength, ocells.values, val_test_answer);
            }
        }

// Full Perfect Remap

        val_test_perfect = (double*)malloc(olength*sizeof(double));
        memset(val_test_perfect, 0xFFFFFFFF, olength*sizeof(double));
        ocells.values = val_test_perfect;
        
        cpu_timer_start(&timer);
        full_perfect_remap (icells, ocells);
        full_perfect_remap_time += cpu_timer_stop(timer);

        if (val_test_answer == NULL) {
            val_test_answer = val_test_perfect;
        } else if (run_tests) {
            check_output("Full Perfect Remap", olength, ocells.values, val_test_answer);
        }

        ocells.values = val_test;

// Single-write Remap
        
        memset(ocells.values,  0xFFFFFFFF, olength*sizeof(double));
        
        cpu_timer_start(&timer);
        singlewrite_remap (icells, ocells);
        singlewrite_remap_time += cpu_timer_stop(timer);

        if (run_tests) check_output("Single-write Remap", olength, ocells.values, val_test_answer);
        
// Hierarchical Remap
        
        memset(ocells.values,  0xFFFFFFFF, olength*sizeof(double));
        
        cpu_timer_start(&timer);
        h_remap (icells, ocells);
        h_remap_time += cpu_timer_stop(timer);
        
        if (run_tests) check_output("Hierarchical Remap", olength, ocells.values, val_test_answer);
        

// Compact Single-write Remap
        
        memset(ocells.values,  0xFFFFFFFF, olength*sizeof(double));
        
        cpu_timer_start(&timer);
        singlewrite_remap_compact (icells, ocells);
        compact_singlewrite_remap_time += cpu_timer_stop(timer);

        if (run_tests) check_output("Compact Single-write Remap", olength, ocells.values, val_test_answer);
        
// Compact Hierarchical Remap
        
        memset(ocells.values,  0xFFFFFFFF, olength*sizeof(double));

        cpu_timer_start(&timer);
        h_remap_compact (icells, ocells, factory);
        compact_h_remap_time += cpu_timer_stop(timer);

        if (run_tests) check_output("Compact Hierarchical Remap", olength, ocells.values, val_test_answer);
        

#ifdef _OPENMP

// Setup for OpenMP to get memory aligned with processors

        icells_openmp.i      = (uint *)  malloc(ilength*sizeof(uint));
        icells_openmp.j      = (uint *)  malloc(ilength*sizeof(uint));
        icells_openmp.level  = (uint *)  malloc(ilength*sizeof(uint));
        icells_openmp.values = (double *)malloc(ilength*sizeof(double));

        ocells_openmp.i      = (uint *)  malloc(olength*sizeof(uint));
        ocells_openmp.j      = (uint *)  malloc(olength*sizeof(uint));
        ocells_openmp.level  = (uint *)  malloc(olength*sizeof(uint));
        ocells_openmp.values = (double *)malloc(olength*sizeof(double));

#pragma omp parallel default(none) firstprivate(ilength, olength, i_min_level) shared(icells_openmp, ocells_openmp, icells, ocells)
        {
#pragma omp for 
           for (uint ic=0; ic < ilength; ic++){
              icells_openmp.i[ic]      = icells.i[ic];
              icells_openmp.j[ic]      = icells.j[ic];
              icells_openmp.level[ic]  = icells.level[ic];
              icells_openmp.values[ic] = icells.values[ic];
           }
#pragma omp for 
           for (uint ic=0; ic < olength; ic++){
              ocells_openmp.i[ic]      = ocells.i[ic];
              ocells_openmp.j[ic]      = ocells.j[ic];
              ocells_openmp.level[ic]  = ocells.level[ic];
              ocells_openmp.values[ic] = -1;
           }

        }

// Full Perfect Remap OpenMP

#pragma omp parallel default(none) firstprivate(ilength, olength) shared(icells_openmp, ocells_openmp, icells, ocells)
        {

#pragma omp for
           for (uint ic=0; ic < olength; ic++){
              ocells_openmp.values[ic] = -1;
           }
        }

        cpu_timer_start(&timer);
        full_perfect_remap_openMP (icells_openmp, ocells_openmp);
        full_perfect_remap_openMP_time += cpu_timer_stop(timer);

        if (run_tests) check_output("Full Perfect Remap OpenMP", ocells.ncells, ocells_openmp.values, val_test_answer);
        
// Single-write Remap OpenMP

#pragma omp parallel default(none) firstprivate(ilength, olength) shared(icells_openmp, ocells_openmp, icells, ocells)
        {

#pragma omp for
           for (uint ic=0; ic < olength; ic++){
              ocells_openmp.values[ic] = -1;
           }
        }

        cpu_timer_start(&timer);
        singlewrite_remap_openMP (icells_openmp, ocells_openmp);
        singlewrite_remap_openMP_time += cpu_timer_stop(timer);

        if (run_tests) check_output("Single-write Remap OpenMP", ocells.ncells, ocells_openmp.values, val_test_answer);
        
// Hierarchical Remap OpenMP

#pragma omp parallel default(none) firstprivate(ilength, olength) shared(icells_openmp, ocells_openmp, icells, ocells)
        {

#pragma omp for
           for (uint ic=0; ic < olength; ic++){
              ocells_openmp.values[ic] = -1;
           }
        }

        cpu_timer_start(&timer);
        h_remap_openMP (icells_openmp, ocells_openmp);
        h_remap_openMP_time += cpu_timer_stop(timer);

        if (run_tests) check_output("Hierarchical Remap OpenMP", olength, ocells_openmp.values, val_test_answer);

// Compact Single-write Remap OpenMP

#pragma omp parallel default(none) firstprivate(ilength, olength) shared(icells_openmp, ocells_openmp, icells, ocells)
        {

#pragma omp for
           for (uint ic=0; ic < olength; ic++){
              ocells_openmp.values[ic] = -1;
           }
        }

        cpu_timer_start(&timer);
        singlewrite_remap_compact_openMP (icells_openmp, ocells_openmp);
        compact_singlewrite_remap_openMP_time += cpu_timer_stop(timer);

        if (run_tests) check_output("Compact Single-write Remap OpenMP", ocells.ncells, ocells_openmp.values, val_test_answer);
        
// Compact Hierarchical Remap OpenMP

#pragma omp parallel default(none) firstprivate(ilength, olength) shared(icells_openmp, ocells_openmp, icells, ocells)
        {

#pragma omp for
           for (uint ic=0; ic < olength; ic++){
              ocells_openmp.values[ic] = -1;
           }
        }

        cpu_timer_start(&timer);
        h_remap_compact_openMP (icells_openmp, ocells_openmp, OpenMPfactory);
        compact_h_remap_openMP_time += cpu_timer_stop(timer);

        if (run_tests) check_output("Compact Hierarchical Remap OpenMP", olength, ocells_openmp.values, val_test_answer);

        free(icells_openmp.i);
        free(icells_openmp.j);
        free(icells_openmp.level);
        free(icells_openmp.values);

        free(ocells_openmp.i);
        free(ocells_openmp.j);
        free(ocells_openmp.level);
        free(ocells_openmp.values);
#endif

#ifdef HAVE_OPENCL

// Full Perfect Remap GPU
        memset(ocells.values,  0xFFFFFFFF, olength*sizeof(double));

        gpu_full_perfect_remap_time+=cl_full_perfect_remap(icells, ocells, run_tests);

        if (run_tests) check_output("GPU Full Perfect Remap", ocells.ncells, ocells.values, val_test_answer);
        

// Single-write Remap GPU
        memset(ocells.values,  0xFFFFFFFF, olength*sizeof(double));

        gpu_singlewrite_remap_time+=cl_singlewrite_remap(icells, ocells, run_tests);

        if (run_tests) check_output("GPU Singlewrite Remap", ocells.ncells, ocells.values, val_test_answer);
        
// Hierarchical Remap GPU
        memset(ocells.values,  0xFFFFFFFF, olength*sizeof(double));

        gpu_hierarchical_remap_time+=cl_hierarchical_remap(icells, ocells, run_tests);

        if (run_tests) check_output("GPU Hierarchical Remap", ocells.ncells, ocells.values, val_test_answer);
        
// Compact Single-write Remap GPU
        memset(ocells.values,  0xFFFFFFFF, olength*sizeof(double));

        gpu_compact_singlewrite_remap_time+=cl_compact_singlewrite_remap(icells, ocells, run_tests);

        if (run_tests) check_output("GPU Compact Singlewrite Remap", ocells.ncells, ocells.values, val_test_answer);
        
// Compact Hierarchical Remap GPU
        memset(ocells.values,  0xFFFFFFFF, olength*sizeof(double));

        gpu_compact_hierarchical_remap_time+=cl_compact_hierarchical_remap(icells, ocells, CLFactory, run_tests);

        if (run_tests) check_output("GPU Compact Hierarchical Remap", ocells.ncells, ocells.values, val_test_answer);
        
#endif

        if (run_brute)  free(val_test_brute);
        if (run_tree)   free(val_test_kdtree);
        free(val_test_perfect);

        destroy(icells);  
        destroy(ocells);
    }

    intintHash_DestroyFactory(factory);
    intintHash_DestroyFactory(OpenMPfactory);
    
    printf("~~~~~~~~~~~~~~~~Averages~~ ~~~~~~~~~~~~~~\n");
    size_t average_ncells = sum_ncells/num_rep/2;

    printf("sparsity percent compressibility    cells in fine mesh average ncells\n");
    printf("---------------- ---------------      --------------     ----------  \n");
    printf("     %f",(float)(save_num_fine_cells-average_ncells)/(float)save_num_fine_cells*100.0);
    printf("     %f",(float)save_num_fine_cells/(float)average_ncells);
    printf("            %ld                %ld",save_num_fine_cells,average_ncells);
    printf("\n");
    printf(" --------------------------------------------------------------------\n");

    if (plot_file) {
       char filename[40];
       if (meshgen == ADAPT_MESHGEN){
          sprintf(filename,"rundata%3d.dat", mesh_size);
       } else {
          sprintf(filename,"rundata%3d.dat", i_level_diff);
       }

       FILE *fout = fopen(filename,"a");
       fprintf(fout,"%2d,\t", levmax);
       if (run_brute) {
          fprintf(fout,"%9.3f,\t", brute_force_time/num_rep*1000);
       } else {
          fprintf(fout,"%9.3f,\t", 0.0);
       }
       if (run_tree) {
          fprintf(fout,"%9.3f,\t", kd_tree_time/num_rep*1000);
       } else {
          fprintf(fout,"%9.3f,\t", 0.0);
       }
       fprintf(fout,"%9.3f,\t", full_perfect_remap_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", singlewrite_remap_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", h_remap_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", compact_singlewrite_remap_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", compact_h_remap_time/num_rep*1000);
#ifdef _OPENMP
       fprintf(fout,"%9.3f,\t", full_perfect_remap_openMP_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", singlewrite_remap_openMP_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", h_remap_openMP_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", compact_singlewrite_remap_openMP_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", compact_h_remap_openMP_time/num_rep*1000);
#endif
       fprintf(fout,"%9.3f,\t", gpu_full_perfect_remap_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", gpu_singlewrite_remap_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", gpu_hierarchical_remap_time/num_rep*1000);
       fprintf(fout,"%9.3f,\t", gpu_compact_singlewrite_remap_time/num_rep*1000);

       fprintf(fout,"%8lu,\t",average_ncells);
       fprintf(fout,"%8.2f,\t",(float)save_num_fine_cells/(float)average_ncells);
       fprintf(fout,"%12lu",save_num_fine_cells);
       fprintf(fout,"\n");
    }

    if (run_brute)
       printf("Brute Force:\t\t\t\t%10.4f ms\n", brute_force_time/num_rep*1000);
    if (run_tree)
       printf("KD Tree Remap:\t\t\t\t%10.4f ms\n", kd_tree_time/num_rep*1000);

    printf("Full Perfect Remap:\t\t\t%10.4f ms\n", full_perfect_remap_time/num_rep*1000);
    printf("Singlewrite Remap:\t\t\t%10.4f ms Speedup relative to full hash %8.2lf\n",
           singlewrite_remap_time/num_rep*1000, full_perfect_remap_time/singlewrite_remap_time);
    printf("Hierarchical Remap:\t\t\t%10.4f ms Speedup relative to full hash %8.2lf\n",
           h_remap_time/num_rep*1000, full_perfect_remap_time/h_remap_time);
    printf("Compact Singlewrite Remap:\t\t%10.4f ms\n", compact_singlewrite_remap_time/num_rep*1000);
    printf("Compact Hierarchical Remap:\t\t%10.4f ms\n", compact_h_remap_time/num_rep*1000);
#ifdef _OPENMP
    printf("\nOpenMP Full Perfect Remap:\t\t%10.4f ms speedup \t%8.2lf\n",
           full_perfect_remap_openMP_time/num_rep*1000, full_perfect_remap_time/full_perfect_remap_openMP_time);
    printf("OpenMP Singlewrite Remap:\t\t%10.4f ms speedup \t%8.2lf\n",
           singlewrite_remap_openMP_time/num_rep*1000,singlewrite_remap_time/singlewrite_remap_openMP_time);
    printf("OpenMP Hierarchical Remap:\t\t%10.4f ms speedup \t%8.2lf\n",
           h_remap_openMP_time/num_rep*1000, h_remap_time/h_remap_openMP_time);
    printf("OpenMP Compact Singlewrite Remap:\t%10.4f ms speedup \t%8.2lf\n",
           compact_singlewrite_remap_openMP_time/num_rep*1000,compact_singlewrite_remap_time/compact_singlewrite_remap_openMP_time);
    printf("OpenMP Compact Hierarchical Remap:\t%10.4f ms speedup \t%8.2lf\n",
           compact_h_remap_openMP_time/num_rep*1000, compact_h_remap_time/compact_h_remap_openMP_time);
#endif
#ifdef HAVE_OPENCL
    printf("\nGPU Full Perfect Remap:\t\t\t%10.4f ms\n", gpu_full_perfect_remap_time/num_rep*1000);
    printf("GPU Singlewrite Remap:\t\t\t%10.4f ms\n", gpu_singlewrite_remap_time/num_rep*1000);
    printf("GPU Hierarchical Remap:\t\t\t%10.4f ms\n", gpu_hierarchical_remap_time/num_rep*1000);
    printf("GPU Compact Singlewrite Remap:\t\t%10.4f ms\n", gpu_compact_singlewrite_remap_time/num_rep*1000);
    printf("GPU Compact Hierarchical Remap:\t\t%10.4f ms\n", gpu_compact_hierarchical_remap_time/num_rep*1000);

    hash_lib_terminate();
    cleanup_cl();
    ezcl_terminate();

    ezcl_mem_walk_all();
#endif
    
    return 0;
}

#ifndef DONT_CATCH_SIGNALS
void handle_signal(int signal){
    (void) signal;
    quit_loop = 1;
}
#endif

void check_output(const char *string, uint olength, double *output_val, double *val_test_answer){
    //printf("Checking %s\n",string);
    int icount = 0;
    for (uint m = 0; m < olength; m++) {
        if (output_val[m] == -1 || output_val[m] != val_test_answer[m]) {
            printf ("%s failed at cell %u\nExpected %f, but found %f\n",
                string, m, val_test_answer[m], output_val[m]);
            icount++;
            if (icount > 6) return;
        }
    }
}
