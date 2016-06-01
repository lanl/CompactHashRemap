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

#define TILE_SIZE 256

#ifdef HAVE_OPENCL
#ifdef __APPLE_CC__
#include <OpenCL/OpenCL.h>
#else
#include <CL/cl.h>
#endif
#endif
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include "meshgen/meshgen.h"
#include "h_remap_gpu.h"
#include "timer.h"

#ifdef HAVE_OPENCL
#include "simplehash/simplehash.h"
#include "ezcl/ezcl.h"

#include "h_remap_kern.inc"

#ifndef DEBUG
#define DEBUG 0
#endif

//#define DETAILED_TIMING 1

cl_platform_id platform_id;
cl_device_id device_id;
cl_context context;
cl_command_queue queue;
cl_kernel full_perfect_hash_setup_kernel;
cl_kernel full_perfect_hash_query_kernel;
cl_kernel singlewrite_hash_init_kernel;
cl_kernel singlewrite_hash_setup_kernel;
cl_kernel singlewrite_hash_query_kernel;
cl_kernel compact_singlewrite_hash_setup_kernel;
cl_kernel compact_singlewrite_hash_query_kernel;
cl_kernel hierarchical_insert_kernel;
cl_kernel hierarchical_probe_kernel;
cl_kernel hierarchical_count_levels_stage1of2_kernel;
cl_kernel hierarchical_count_levels_stage2of2_kernel;
cl_kernel hierarchical_compact_insert_kernel;
cl_kernel hierarchical_compact_probe_kernel;

void setup_cl (){
    size_t program_size = strlen(h_remap_kern_source);
    size_t hashlib_size = strlen(get_hash_kernel_source_string());
    size_t HashFactory_size = strlen(Hash_GetKernelSourceString());

    //printf("DEBUG -- sizes are program_size %lu hashlib_size %lu HashFactory_size %lu\n",program_size,hashlib_size,HashFactory_size);
    size_t all_sources_size = program_size + hashlib_size + HashFactory_size;
    char *all_sources = (char *)malloc((all_sources_size+1)*sizeof(char));

    strcpy(all_sources, Hash_GetKernelSourceString());
    strcat(all_sources, get_hash_kernel_source_string());
    strcat(all_sources, h_remap_kern_source);
    all_sources[all_sources_size] = '\0';
    
    cl_context context ezcl_get_context();
    cl_program program = ezcl_create_program_wsource(context, NULL, (const char *)all_sources);

    full_perfect_hash_setup_kernel      = ezcl_create_kernel_wprogram(program, "full_perfect_hash_setup");
    full_perfect_hash_query_kernel      = ezcl_create_kernel_wprogram(program, "full_perfect_hash_query");
    singlewrite_hash_init_kernel          = ezcl_create_kernel_wprogram(program, "singlewrite_hash_init");
    singlewrite_hash_setup_kernel         = ezcl_create_kernel_wprogram(program, "singlewrite_hash_setup");
    singlewrite_hash_query_kernel         = ezcl_create_kernel_wprogram(program, "singlewrite_hash_query");
    compact_singlewrite_hash_setup_kernel = ezcl_create_kernel_wprogram(program, "compact_singlewrite_hash_setup");
    compact_singlewrite_hash_query_kernel = ezcl_create_kernel_wprogram(program, "compact_singlewrite_hash_query");
    hierarchical_insert_kernel          = ezcl_create_kernel_wprogram(program, "hierarchical_cell_insert");
    hierarchical_probe_kernel           = ezcl_create_kernel_wprogram(program, "hierarchical_hash_probe");
    hierarchical_count_levels_stage1of2_kernel = ezcl_create_kernel_wprogram(program, "hierarchical_count_levels_stage1of2");
    hierarchical_count_levels_stage2of2_kernel = ezcl_create_kernel_wprogram(program, "hierarchical_count_levels_stage2of2");
    hierarchical_compact_insert_kernel  = ezcl_create_kernel_wprogram(program, "hierarchical_compact_insert");
    hierarchical_compact_probe_kernel   = ezcl_create_kernel_wprogram(program, "hierarchical_compact_probe");

    free(all_sources);

    ezcl_program_release(program);
}

void cleanup_cl(){
    ezcl_kernel_release(full_perfect_hash_setup_kernel);
    ezcl_kernel_release(full_perfect_hash_query_kernel);
    ezcl_kernel_release(singlewrite_hash_init_kernel);
    ezcl_kernel_release(singlewrite_hash_setup_kernel);
    ezcl_kernel_release(singlewrite_hash_query_kernel);
    ezcl_kernel_release(compact_singlewrite_hash_setup_kernel);
    ezcl_kernel_release(compact_singlewrite_hash_query_kernel);
    ezcl_kernel_release(hierarchical_insert_kernel);
    ezcl_kernel_release(hierarchical_probe_kernel);
    ezcl_kernel_release(hierarchical_count_levels_stage1of2_kernel);
    ezcl_kernel_release(hierarchical_count_levels_stage2of2_kernel);
    ezcl_kernel_release(hierarchical_compact_insert_kernel);
    ezcl_kernel_release(hierarchical_compact_probe_kernel);
 }

// returns the timing and the results in oval if run_tests flag is set
double cl_full_perfect_remap (cell_list icells, cell_list ocells, int run_tests){
                 
    cl_context       context = ezcl_get_context();
    cl_command_queue queue   = ezcl_get_command_queue();

    struct timeval timer;
#ifdef DETAILED_TIMING
    struct timeval timer1;
#endif

    // Input mesh
    cl_mem icelli_buffer = ezcl_device_memory_malloc(context, NULL, "icelli", icells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem icellj_buffer = ezcl_device_memory_malloc(context, NULL, "icellj", icells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ilevel_buffer = ezcl_device_memory_malloc(context, NULL, "ilevel", icells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ival_buffer   = ezcl_device_memory_malloc(context, NULL, "ival",   icells.ncells, sizeof(double), CL_MEM_READ_WRITE, 0);

    ezcl_enqueue_write_buffer(queue, icelli_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.i,      NULL);
    ezcl_enqueue_write_buffer(queue, icellj_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.j,      NULL);
    ezcl_enqueue_write_buffer(queue, ilevel_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.level,  NULL);
    ezcl_enqueue_write_buffer(queue, ival_buffer,   CL_TRUE,  0, icells.ncells*sizeof(double), icells.values, NULL);
    
    // The output mesh
    cl_mem ocelli_buffer = ezcl_device_memory_malloc(context, NULL, "ocelli", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ocellj_buffer = ezcl_device_memory_malloc(context, NULL, "ocellj", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem olevel_buffer = ezcl_device_memory_malloc(context, NULL, "olevel", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem oval_buffer   = ezcl_device_memory_malloc(context, NULL, "oval",   ocells.ncells, sizeof(double), CL_MEM_READ_WRITE, 0);

    ezcl_enqueue_write_buffer(queue, ocelli_buffer, CL_FALSE, 0, ocells.ncells*sizeof(uint), ocells.i,     NULL);
    ezcl_enqueue_write_buffer(queue, ocellj_buffer, CL_FALSE, 0, ocells.ncells*sizeof(uint), ocells.j,     NULL);
    ezcl_enqueue_write_buffer(queue, olevel_buffer, CL_TRUE,  0, ocells.ncells*sizeof(uint), ocells.level, NULL);
    
    //START TIMER
    cpu_timer_start(&timer);
#ifdef DETAILED_TIMING
    cpu_timer_start(&timer1);
#endif
    
    size_t global_work_size[1];
    size_t local_work_size[1];
    
    local_work_size[0] = TILE_SIZE;

    global_work_size[0] = ((local_work_size[0]+icells.ncells-1)/local_work_size[0])*local_work_size[0];

    // The hash is the size of the sum of all the points in the array - so enough to hold all the scanned powers of 4
    uint hash_size = icells.ibasesize*two_to_the(icells.levmax)*icells.jbasesize*two_to_the(icells.levmax);
    cl_mem hash_buffer = ezcl_device_memory_malloc(context, NULL, "hash", hash_size, sizeof(int), CL_MEM_READ_WRITE, 0);

    ezcl_set_kernel_arg(full_perfect_hash_setup_kernel, 0, sizeof(cl_uint), &icells.ncells);
    ezcl_set_kernel_arg(full_perfect_hash_setup_kernel, 1, sizeof(cl_uint), &icells.ibasesize);
    ezcl_set_kernel_arg(full_perfect_hash_setup_kernel, 2, sizeof(cl_uint), &icells.levmax);
    ezcl_set_kernel_arg(full_perfect_hash_setup_kernel, 3, sizeof(cl_mem),  &icelli_buffer);
    ezcl_set_kernel_arg(full_perfect_hash_setup_kernel, 4, sizeof(cl_mem),  &icellj_buffer);
    ezcl_set_kernel_arg(full_perfect_hash_setup_kernel, 5, sizeof(cl_mem),  &ilevel_buffer);
    ezcl_set_kernel_arg(full_perfect_hash_setup_kernel, 6, sizeof(cl_mem),  &hash_buffer);

#ifdef DETAILED_TIMING
    cl_event hash_setup_event;
    ezcl_enqueue_ndrange_kernel(queue, full_perfect_hash_setup_kernel, 1, 0, global_work_size, local_work_size, &hash_setup_event);
    long long gpu_time = ezcl_timer_calc(&hash_setup_event, &hash_setup_event);
    double setup_time = cpu_timer_stop(timer1);
    printf("setup time is %8.4f ms %8.4f ms\n",1.0e-6*(double)gpu_time, setup_time*1000.0);
    cpu_timer_start(&timer1);
    //clReleaseEvent(hash_setup_event);
#else
    ezcl_enqueue_ndrange_kernel(queue, full_perfect_hash_setup_kernel, 1, 0, global_work_size, local_work_size, NULL);
#endif

    global_work_size[0] = ((local_work_size[0]+ocells.ncells-1)/local_work_size[0])*local_work_size[0];
    
    ezcl_set_kernel_arg(full_perfect_hash_query_kernel, 0, sizeof(cl_uint), &ocells.ncells);
    ezcl_set_kernel_arg(full_perfect_hash_query_kernel, 1, sizeof(cl_uint), &ocells.ibasesize);
    ezcl_set_kernel_arg(full_perfect_hash_query_kernel, 2, sizeof(cl_uint), &ocells.levmax);
    ezcl_set_kernel_arg(full_perfect_hash_query_kernel, 3, sizeof(cl_mem),  &hash_buffer);
    ezcl_set_kernel_arg(full_perfect_hash_query_kernel, 4, sizeof(cl_mem),  &ival_buffer);
    ezcl_set_kernel_arg(full_perfect_hash_query_kernel, 5, sizeof(cl_mem),  &ocelli_buffer);
    ezcl_set_kernel_arg(full_perfect_hash_query_kernel, 6, sizeof(cl_mem),  &ocellj_buffer);
    ezcl_set_kernel_arg(full_perfect_hash_query_kernel, 7, sizeof(cl_mem),  &olevel_buffer);
    ezcl_set_kernel_arg(full_perfect_hash_query_kernel, 8, sizeof(cl_mem),  &oval_buffer);

#ifdef DETAILED_TIMING
    cl_event hash_query_event;
    ezcl_enqueue_ndrange_kernel(queue, full_perfect_hash_query_kernel, 1, 0, global_work_size, local_work_size, &hash_query_event);
    gpu_time = ezcl_timer_calc(&hash_query_event, &hash_query_event);
    double query_time = cpu_timer_stop(timer1);
    printf("query time is %8.4f ms %8.4f ms\n",1.0e-6*(double)gpu_time, query_time*1000.0);
    //clReleaseEvent(hash_query_event);
#else
    ezcl_enqueue_ndrange_kernel(queue, full_perfect_hash_query_kernel, 1, 0, global_work_size, local_work_size, NULL);
#endif
    
    ezcl_finish(queue);
    
    //END TIMER
    double time = cpu_timer_stop(timer);
#ifdef DETAILED_TIMING
    printf("perfect hash time is %8.4f\n", time*1000.0);
#endif
    
    if (run_tests) ezcl_enqueue_read_buffer(queue, oval_buffer, CL_TRUE, 0, ocells.ncells*sizeof(double), ocells.values, NULL);

    ezcl_device_memory_delete(icelli_buffer);
    ezcl_device_memory_delete(icellj_buffer);
    ezcl_device_memory_delete(ilevel_buffer);
    ezcl_device_memory_delete(ival_buffer);
    
    ezcl_device_memory_delete(ocelli_buffer);
    ezcl_device_memory_delete(ocellj_buffer);
    ezcl_device_memory_delete(olevel_buffer);
    ezcl_device_memory_delete(oval_buffer);
    
    ezcl_device_memory_delete(hash_buffer);

    return time;
}

// returns the timing and the results in oval if run_tests flag is set
double cl_singlewrite_remap (cell_list icells, cell_list ocells, int run_tests){

    cl_context       context = ezcl_get_context();
    cl_command_queue queue   = ezcl_get_command_queue();

    struct timeval timer;
#ifdef DETAILED_TIMING
    struct timeval timer1;
#endif

    // Input mesh
    cl_mem icelli_buffer = ezcl_device_memory_malloc(context, NULL, "icelli", icells.ncells, sizeof(uint), CL_MEM_READ_WRITE, 0);
    cl_mem icellj_buffer = ezcl_device_memory_malloc(context, NULL, "icellj", icells.ncells, sizeof(uint), CL_MEM_READ_WRITE, 0);
    cl_mem ilevel_buffer = ezcl_device_memory_malloc(context, NULL, "ilevel", icells.ncells, sizeof(uint), CL_MEM_READ_WRITE, 0);
    cl_mem ival_buffer   = ezcl_device_memory_malloc(context, NULL, "ival",   icells.ncells, sizeof(double), CL_MEM_READ_WRITE, 0);

    ezcl_enqueue_write_buffer(queue, icelli_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.i,      NULL);
    ezcl_enqueue_write_buffer(queue, icellj_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.j,      NULL);
    ezcl_enqueue_write_buffer(queue, ilevel_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.level,  NULL);
    ezcl_enqueue_write_buffer(queue, ival_buffer,   CL_TRUE,  0, icells.ncells*sizeof(double), icells.values, NULL);
    
    // The output mesh
    cl_mem ocelli_buffer = ezcl_device_memory_malloc(context, NULL, "ocelli", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ocellj_buffer = ezcl_device_memory_malloc(context, NULL, "ocellj", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem olevel_buffer = ezcl_device_memory_malloc(context, NULL, "olevel", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem oval_buffer   = ezcl_device_memory_malloc(context, NULL, "oval",   ocells.ncells, sizeof(double), CL_MEM_READ_WRITE, 0);

    ezcl_enqueue_write_buffer(queue, ocelli_buffer, CL_FALSE, 0, ocells.ncells*sizeof(uint),   ocells.i,     NULL);
    ezcl_enqueue_write_buffer(queue, ocellj_buffer, CL_FALSE, 0, ocells.ncells*sizeof(uint),   ocells.j,     NULL);
    ezcl_enqueue_write_buffer(queue, olevel_buffer, CL_TRUE,  0, ocells.ncells*sizeof(uint),   ocells.level, NULL);
    
    //START TIMER
    cpu_timer_start(&timer);
#ifdef DETAILED_TIMING
    cpu_timer_start(&timer1);
#endif
    
    size_t global_work_size[1];
    size_t local_work_size[1];
    
    local_work_size[0] = TILE_SIZE;

    // The hash is the size of the sum of all the points in the array - so enough to hold all the scanned powers of 4
    uint hash_size = icells.ibasesize*two_to_the(icells.levmax)*icells.jbasesize*two_to_the(icells.levmax);
    cl_mem hash_buffer = ezcl_device_memory_malloc(context, NULL, "hash", hash_size, sizeof(int), CL_MEM_READ_WRITE, 0);

    global_work_size[0] = ((local_work_size[0]+hash_size-1)/local_work_size[0])*local_work_size[0];

    ezcl_set_kernel_arg(singlewrite_hash_init_kernel, 0, sizeof(cl_uint), &hash_size);
    ezcl_set_kernel_arg(singlewrite_hash_init_kernel, 1, sizeof(cl_mem),  &hash_buffer);

#ifdef DETAILED_TIMING
    cl_event hash_init_event;
    ezcl_enqueue_ndrange_kernel(queue, singlewrite_hash_init_kernel, 1, 0, global_work_size, local_work_size, &hash_init_event);
    long long gpu_time = ezcl_timer_calc(&hash_init_event, &hash_init_event);
    double init_time = cpu_timer_stop(timer1);
    printf("init time is %8.4f ms %8.4f ms\n",1.0e-6*(double)gpu_time, init_time*1000.0);
    cpu_timer_start(&timer1);
    //clReleaseEvent(hash_init_event);
#else
    ezcl_enqueue_ndrange_kernel(queue, singlewrite_hash_init_kernel, 1, 0, global_work_size, local_work_size, NULL);
#endif

    global_work_size[0] = ((local_work_size[0]+icells.ncells-1)/local_work_size[0])*local_work_size[0];

    ezcl_set_kernel_arg(singlewrite_hash_setup_kernel, 0, sizeof(cl_uint), &icells.ncells);
    ezcl_set_kernel_arg(singlewrite_hash_setup_kernel, 1, sizeof(cl_uint), &icells.ibasesize);
    ezcl_set_kernel_arg(singlewrite_hash_setup_kernel, 2, sizeof(cl_uint), &icells.levmax);
    ezcl_set_kernel_arg(singlewrite_hash_setup_kernel, 3, sizeof(cl_mem),  &icelli_buffer);
    ezcl_set_kernel_arg(singlewrite_hash_setup_kernel, 4, sizeof(cl_mem),  &icellj_buffer);
    ezcl_set_kernel_arg(singlewrite_hash_setup_kernel, 5, sizeof(cl_mem),  &ilevel_buffer);
    ezcl_set_kernel_arg(singlewrite_hash_setup_kernel, 6, sizeof(cl_mem),  &hash_buffer);

#ifdef DETAILED_TIMING
    cl_event hash_setup_event;
    ezcl_enqueue_ndrange_kernel(queue, singlewrite_hash_setup_kernel, 1, 0, global_work_size, local_work_size, &hash_setup_event);
    gpu_time = ezcl_timer_calc(&hash_setup_event, &hash_setup_event);
    double setup_time = cpu_timer_stop(timer1);
    printf("setup time is %8.4f ms %8.4f ms\n",1.0e-6*(double)gpu_time, setup_time*1000.0);
    cpu_timer_start(&timer1);
    //clReleaseEvent(hash_setup_event);
#else
    ezcl_enqueue_ndrange_kernel(queue, singlewrite_hash_setup_kernel, 1, 0, global_work_size, local_work_size, NULL);
#endif

    global_work_size[0] = ((local_work_size[0]+ocells.ncells-1)/local_work_size[0])*local_work_size[0];
    
    ezcl_set_kernel_arg(singlewrite_hash_query_kernel, 0, sizeof(cl_uint), &ocells.ncells);
    ezcl_set_kernel_arg(singlewrite_hash_query_kernel, 1, sizeof(cl_uint), &ocells.ibasesize);
    ezcl_set_kernel_arg(singlewrite_hash_query_kernel, 2, sizeof(cl_uint), &ocells.levmax);
    ezcl_set_kernel_arg(singlewrite_hash_query_kernel, 3, sizeof(cl_mem),  &hash_buffer);
    ezcl_set_kernel_arg(singlewrite_hash_query_kernel, 4, sizeof(cl_mem),  &ilevel_buffer);
    ezcl_set_kernel_arg(singlewrite_hash_query_kernel, 5, sizeof(cl_mem),  &ival_buffer);
    ezcl_set_kernel_arg(singlewrite_hash_query_kernel, 6, sizeof(cl_mem),  &ocelli_buffer);
    ezcl_set_kernel_arg(singlewrite_hash_query_kernel, 7, sizeof(cl_mem),  &ocellj_buffer);
    ezcl_set_kernel_arg(singlewrite_hash_query_kernel, 8, sizeof(cl_mem),  &olevel_buffer);
    ezcl_set_kernel_arg(singlewrite_hash_query_kernel, 9, sizeof(cl_mem),  &oval_buffer);

#ifdef DETAILED_TIMING
    cl_event hash_query_event;
    ezcl_enqueue_ndrange_kernel(queue, singlewrite_hash_query_kernel, 1, 0, global_work_size, local_work_size, &hash_query_event);
    gpu_time = ezcl_timer_calc(&hash_query_event, &hash_query_event);
    double query_time = cpu_timer_stop(timer1);
    printf("query time is %8.4f ms %8.4f ms\n",1.0e-6*(double)gpu_time, query_time*1000.0);
    //clReleaseEvent(hash_query_event);
#else
    ezcl_enqueue_ndrange_kernel(queue, singlewrite_hash_query_kernel, 1, 0, global_work_size, local_work_size, NULL);
#endif
    
    ezcl_finish(queue);
    
    //END TIMER
    double time = cpu_timer_stop(timer);
#ifdef DETAILED_TIMING
    printf("singlewrite hash time is %8.4f\n", time*1000.0);
#endif
    
    if (run_tests) ezcl_enqueue_read_buffer(queue, oval_buffer, CL_TRUE, 0, ocells.ncells*sizeof(double), ocells.values, NULL);

    ezcl_device_memory_delete(icelli_buffer);
    ezcl_device_memory_delete(icellj_buffer);
    ezcl_device_memory_delete(ilevel_buffer);
    ezcl_device_memory_delete(ival_buffer);
    
    ezcl_device_memory_delete(ocelli_buffer);
    ezcl_device_memory_delete(ocellj_buffer);
    ezcl_device_memory_delete(olevel_buffer);
    ezcl_device_memory_delete(oval_buffer);
    
    ezcl_device_memory_delete(hash_buffer);

    return time;
}

// returns the timing and the results in oval if run_tests flag is set
double cl_compact_singlewrite_remap (cell_list icells, cell_list ocells, int run_tests){
                 
    cl_context       context = ezcl_get_context();
    cl_command_queue queue   = ezcl_get_command_queue();

    struct timeval timer;
#ifdef DETAILED_TIMING
    struct timeval timer1;
#endif

    // Input mesh
    cl_mem icelli_buffer = ezcl_device_memory_malloc(context, NULL, "icelli", icells.ncells, sizeof(uint), CL_MEM_READ_WRITE, 0);
    cl_mem icellj_buffer = ezcl_device_memory_malloc(context, NULL, "icellj", icells.ncells, sizeof(uint), CL_MEM_READ_WRITE, 0);
    cl_mem ilevel_buffer = ezcl_device_memory_malloc(context, NULL, "ilevel", icells.ncells, sizeof(uint), CL_MEM_READ_WRITE, 0);
    cl_mem ival_buffer   = ezcl_device_memory_malloc(context, NULL, "ival",   icells.ncells, sizeof(double), CL_MEM_READ_WRITE, 0);

    ezcl_enqueue_write_buffer(queue, icelli_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.i,      NULL);
    ezcl_enqueue_write_buffer(queue, icellj_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.j,      NULL);
    ezcl_enqueue_write_buffer(queue, ilevel_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.level,  NULL);
    ezcl_enqueue_write_buffer(queue, ival_buffer,   CL_TRUE,  0, icells.ncells*sizeof(double), icells.values, NULL);
    
    // The output mesh
    cl_mem ocelli_buffer = ezcl_device_memory_malloc(context, NULL, "ocelli", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ocellj_buffer = ezcl_device_memory_malloc(context, NULL, "ocellj", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem olevel_buffer = ezcl_device_memory_malloc(context, NULL, "olevel", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem oval_buffer   = ezcl_device_memory_malloc(context, NULL, "oval",   ocells.ncells, sizeof(double), CL_MEM_READ_WRITE, 0);

    ezcl_enqueue_write_buffer(queue, ocelli_buffer, CL_FALSE, 0, ocells.ncells*sizeof(uint),   ocells.i,   NULL);
    ezcl_enqueue_write_buffer(queue, ocellj_buffer, CL_FALSE, 0, ocells.ncells*sizeof(uint),   ocells.j,   NULL);
    ezcl_enqueue_write_buffer(queue, olevel_buffer, CL_TRUE,  0, ocells.ncells*sizeof(uint),   ocells.level, NULL);
    
    //START TIMER
    cpu_timer_start(&timer);
#ifdef DETAILED_TIMING
    cpu_timer_start(&timer1);
#endif

    size_t global_work_size[1];
    size_t local_work_size[1];

    local_work_size[0] = TILE_SIZE;

    int gpu_hash_method = QUADRATIC;
    uint hash_report_level = 0;
    ulong gpu_hash_table_size = 0;
    size_t hashsize;
    cl_mem dev_hash_header = NULL;
    cl_mem dev_hash = gpu_compact_hash_init(icells.ncells, icells.ibasesize*two_to_the(icells.levmax),
       icells.jbasesize*two_to_the(icells.levmax), 0, gpu_hash_method, hash_report_level,
       &gpu_hash_table_size, &hashsize, &dev_hash_header);
#ifdef DETAILED_TIMING
    double init_time = cpu_timer_stop(timer1);
    printf("init time is %8.4f ms\n",init_time*1000.0);
    cpu_timer_start(&timer1);
#endif

    global_work_size[0] = ((local_work_size[0]+icells.ncells-1)/local_work_size[0])*local_work_size[0];

    ezcl_set_kernel_arg(compact_singlewrite_hash_setup_kernel, 0, sizeof(cl_uint), &icells.ncells);
    ezcl_set_kernel_arg(compact_singlewrite_hash_setup_kernel, 1, sizeof(cl_uint), &icells.ibasesize);
    ezcl_set_kernel_arg(compact_singlewrite_hash_setup_kernel, 2, sizeof(cl_uint), &icells.levmax);
    ezcl_set_kernel_arg(compact_singlewrite_hash_setup_kernel, 3, sizeof(cl_mem),  &icelli_buffer);
    ezcl_set_kernel_arg(compact_singlewrite_hash_setup_kernel, 4, sizeof(cl_mem),  &icellj_buffer);
    ezcl_set_kernel_arg(compact_singlewrite_hash_setup_kernel, 5, sizeof(cl_mem),  &ilevel_buffer);
    ezcl_set_kernel_arg(compact_singlewrite_hash_setup_kernel, 6, sizeof(cl_mem),  &dev_hash_header);
    ezcl_set_kernel_arg(compact_singlewrite_hash_setup_kernel, 7, sizeof(cl_mem),  &dev_hash);

#ifdef DETAILED_TIMING
    cl_event hash_setup_event;
    ezcl_enqueue_ndrange_kernel(queue, compact_singlewrite_hash_setup_kernel, 1, 0, global_work_size, local_work_size, &hash_setup_event);
    long long gpu_time = ezcl_timer_calc(&hash_setup_event, &hash_setup_event);
    double setup_time = cpu_timer_stop(timer1);
    printf("setup time is %8.4f ms %8.4f ms\n",1.0e-6*(double)gpu_time, setup_time*1000.0);
    cpu_timer_start(&timer1);
    //clReleaseEvent(hash_setup_event);
#else
    ezcl_enqueue_ndrange_kernel(queue, compact_singlewrite_hash_setup_kernel, 1, 0, global_work_size, local_work_size, NULL);
#endif

    global_work_size[0] = ((local_work_size[0]+ocells.ncells-1)/local_work_size[0])*local_work_size[0];

    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel,  0, sizeof(cl_uint), &ocells.ncells);
    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel,  1, sizeof(cl_uint), &ocells.ibasesize);
    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel,  2, sizeof(cl_uint), &ocells.levmax);
    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel,  3, sizeof(cl_mem),  &dev_hash_header);
    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel,  4, sizeof(cl_mem),  &dev_hash);
    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel,  5, sizeof(cl_mem),  &ilevel_buffer);
    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel,  6, sizeof(cl_mem),  &ival_buffer);
    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel,  7, sizeof(cl_mem),  &ocelli_buffer);
    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel,  8, sizeof(cl_mem),  &ocellj_buffer);
    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel,  9, sizeof(cl_mem),  &olevel_buffer);
    ezcl_set_kernel_arg(compact_singlewrite_hash_query_kernel, 10, sizeof(cl_mem),  &oval_buffer);

#ifdef DETAILED_TIMING
    cl_event hash_query_event;
    ezcl_enqueue_ndrange_kernel(queue, compact_singlewrite_hash_query_kernel, 1, 0, global_work_size, local_work_size, &hash_query_event);
    gpu_time = ezcl_timer_calc(&hash_query_event, &hash_query_event);
    double query_time = cpu_timer_stop(timer1);
    printf("query time is %8.4f ms %8.4f ms\n",1.0e-6*(double)gpu_time, query_time*1000.0);
    //clReleaseEvent(hash_query_event);
#else
    ezcl_enqueue_ndrange_kernel(queue, compact_singlewrite_hash_query_kernel, 1, 0, global_work_size, local_work_size, NULL);
#endif

    ezcl_finish(queue);
    
    //END TIMER
    double time = cpu_timer_stop(timer);
#ifdef DETAILED_TIMING
    printf("compact singlewrite hash time is %8.4f\n", time*1000.0);
#endif
    
    if (run_tests) ezcl_enqueue_read_buffer(queue, oval_buffer, CL_TRUE, 0, ocells.ncells*sizeof(double), ocells.values, NULL);

    ezcl_device_memory_delete(icelli_buffer);
    ezcl_device_memory_delete(icellj_buffer);
    ezcl_device_memory_delete(ilevel_buffer);
    ezcl_device_memory_delete(ival_buffer);
    
    ezcl_device_memory_delete(ocelli_buffer);
    ezcl_device_memory_delete(ocellj_buffer);
    ezcl_device_memory_delete(olevel_buffer);
    ezcl_device_memory_delete(oval_buffer);
    
    ezcl_device_memory_delete(dev_hash_header);
    ezcl_device_memory_delete(dev_hash);

    return time;
}


// returns the timing and the results in oval if run_tests flag is set
double cl_hierarchical_remap (cell_list icells, cell_list ocells, int run_tests){
                 
    cl_context       context = ezcl_get_context();
    cl_command_queue queue   = ezcl_get_command_queue();

    struct timeval timer;
#ifdef DETAILED_TIMING
    struct timeval timer1;
#endif

    // Input mesh
    cl_mem icelli_buffer = ezcl_device_memory_malloc(context, NULL, "icelli", icells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem icellj_buffer = ezcl_device_memory_malloc(context, NULL, "icellj", icells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ilevel_buffer = ezcl_device_memory_malloc(context, NULL, "ilevel", icells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ival_buffer   = ezcl_device_memory_malloc(context, NULL, "ival",   icells.ncells, sizeof(double), CL_MEM_READ_WRITE, 0);

    ezcl_enqueue_write_buffer(queue, icelli_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.i,      NULL);
    ezcl_enqueue_write_buffer(queue, icellj_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.j,      NULL);
    ezcl_enqueue_write_buffer(queue, ilevel_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.level,  NULL);
    ezcl_enqueue_write_buffer(queue, ival_buffer,   CL_TRUE,  0, icells.ncells*sizeof(double), icells.values, NULL);
    
    // The output mesh
    cl_mem ocelli_buffer = ezcl_device_memory_malloc(context, NULL, "ocelli", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ocellj_buffer = ezcl_device_memory_malloc(context, NULL, "ocellj", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem olevel_buffer = ezcl_device_memory_malloc(context, NULL, "olevel", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem oval_buffer   = ezcl_device_memory_malloc(context, NULL, "oval",   ocells.ncells, sizeof(double), CL_MEM_READ_WRITE, 0);

    ezcl_enqueue_write_buffer(queue, ocelli_buffer, CL_FALSE, 0, ocells.ncells*sizeof(uint), ocells.i,     NULL);
    ezcl_enqueue_write_buffer(queue, ocellj_buffer, CL_FALSE, 0, ocells.ncells*sizeof(uint), ocells.j,     NULL);
    ezcl_enqueue_write_buffer(queue, olevel_buffer, CL_TRUE,  0, ocells.ncells*sizeof(uint), ocells.level, NULL);
    
    //START TIMER
    cpu_timer_start(&timer);
#ifdef DETAILED_TIMING
    cpu_timer_start(&timer1);
#endif
    
    size_t global_work_size[1];
    size_t local_work_size[1];
    
    local_work_size[0] = TILE_SIZE;
    global_work_size[0] = TILE_SIZE;
            
    /*
    * initialize array
    */
    // This is a simple prefix scan - it is done on the cpu because it also 
    // determines the amount of memory that needs to be alloc'd, so the memory
    // would have to be moved cpu/gpu anyway. Should not be many operations.
    uint index_cume = 0;
    uint* hash_memory_indices = (uint*)malloc((icells.levmax+2)*sizeof(uint));
    memset (hash_memory_indices,0,(icells.levmax+2)*sizeof(uint));
    for (uint i = 0; i <= icells.levmax; i++){
        size_t hash_size = icells.ibasesize*two_to_the(i)*icells.jbasesize*two_to_the(i);
        index_cume += hash_size;
        // the last value (at index maxLev + 1) is the size of the hash, so kept
        hash_memory_indices[i+1] = index_cume;
    }
    
    cl_mem hash_memory_indices_buffer = ezcl_device_memory_malloc(context, NULL, "hash_memory_indices_buffer", (icells.levmax+1), sizeof(uint), CL_MEM_READ_WRITE, 0);
    ezcl_enqueue_write_buffer(queue, hash_memory_indices_buffer, CL_TRUE, 0, (icells.levmax+1)*sizeof(uint), hash_memory_indices, NULL);
    
    // The hash is the size of the sum of all the points in the array - so enough to hold all the scanned powers of 4
    cl_mem hhash = ezcl_device_memory_malloc(context, NULL, "hhash", hash_memory_indices[icells.levmax+1], sizeof(int), CL_MEM_READ_WRITE, 0);
    
    global_work_size[0] = ((local_work_size[0]+icells.ncells-1)/local_work_size[0])*local_work_size[0];
    
    /*
    * insert the cells into the hash
    */
    
    ezcl_set_kernel_arg(hierarchical_insert_kernel, 0, sizeof(cl_mem),(void*)&hhash);
    ezcl_set_kernel_arg(hierarchical_insert_kernel, 1, sizeof(cl_mem), (void*)&hash_memory_indices_buffer);
    ezcl_set_kernel_arg(hierarchical_insert_kernel, 2, sizeof(cl_mem), (void*)&icelli_buffer);
    ezcl_set_kernel_arg(hierarchical_insert_kernel, 3, sizeof(cl_mem), (void*)&icellj_buffer);
    ezcl_set_kernel_arg(hierarchical_insert_kernel, 4, sizeof(cl_mem), (void*)&ilevel_buffer);
    ezcl_set_kernel_arg(hierarchical_insert_kernel, 5, sizeof(uint), &icells.ncells);
    ezcl_set_kernel_arg(hierarchical_insert_kernel, 6, sizeof(uint), &icells.ibasesize);
    
#ifdef DETAILED_TIMING
    cl_event insert_event;
    ezcl_enqueue_ndrange_kernel(queue, hierarchical_insert_kernel, 1, 0, global_work_size, local_work_size, &insert_event);
    long long gpu_time = ezcl_timer_calc(&insert_event, &insert_event);
    double insert_time = cpu_timer_stop(timer1);
    printf("insert time is %8.4f ms %8.4f ms\n",1.0e-6*(double)gpu_time, insert_time*1000.0);
#else
    ezcl_enqueue_ndrange_kernel(queue, hierarchical_insert_kernel, 1, 0, global_work_size, local_work_size, NULL);
#endif

    // precompute the powers so that we don't have to recalculate for every cell    
    /*
    uint secondPowers [maxLev+1];
    for (int i = 0; i <= maxLev; i++){
        secondPowers[i] = two_to_the (i);
    }
    ezcl_enqueue_write_buffer(queue, secondPower_buffer, CL_TRUE, 0, (maxLev+1)*sizeof(uint), secondPowers, NULL);
    */
	
    global_work_size[0] = ((local_work_size[0]+ocells.ncells-1)/local_work_size[0])*local_work_size[0];
    
    ezcl_set_kernel_arg(hierarchical_probe_kernel, 0, sizeof(cl_mem),(void*)&hhash);
    ezcl_set_kernel_arg(hierarchical_probe_kernel, 1, sizeof(cl_mem), (void*)&hash_memory_indices_buffer);
    ezcl_set_kernel_arg(hierarchical_probe_kernel, 2, sizeof(cl_mem), (void*)&ocelli_buffer);
    ezcl_set_kernel_arg(hierarchical_probe_kernel, 3, sizeof(cl_mem), (void*)&ocellj_buffer);
    ezcl_set_kernel_arg(hierarchical_probe_kernel, 4, sizeof(cl_mem), (void*)&olevel_buffer);
    ezcl_set_kernel_arg(hierarchical_probe_kernel, 5, sizeof(cl_mem), (void*)&ival_buffer);
    ezcl_set_kernel_arg(hierarchical_probe_kernel, 6, sizeof(cl_mem), (void*)&oval_buffer);
    ezcl_set_kernel_arg(hierarchical_probe_kernel, 7, sizeof(uint), &ocells.ncells);
    ezcl_set_kernel_arg(hierarchical_probe_kernel, 8, sizeof(uint), &ocells.ibasesize);
    
#ifdef DETAILED_TIMING
    cl_event probe_event;
    ezcl_enqueue_ndrange_kernel(queue, hierarchical_probe_kernel, 1, 0, global_work_size, local_work_size, &probe_event);
    gpu_time = ezcl_timer_calc(&probe_event, &probe_event);
    double probe_time = cpu_timer_stop(timer1);
    printf("probe time is %8.4f ms %8.4f ms\n",1.0e-6*(double)gpu_time, probe_time*1000.0);
#else
    ezcl_enqueue_ndrange_kernel(queue, hierarchical_probe_kernel, 1, 0, global_work_size, local_work_size, NULL);
#endif
    
    ezcl_finish(queue);
    
    //END TIMER
    double time = cpu_timer_stop(timer);
    
    if (run_tests) ezcl_enqueue_read_buffer(queue, oval_buffer, CL_TRUE, 0, ocells.ncells*sizeof(double), ocells.values, NULL);

    ezcl_device_memory_delete(icelli_buffer);
    ezcl_device_memory_delete(icellj_buffer);
    ezcl_device_memory_delete(ilevel_buffer);
    ezcl_device_memory_delete(ival_buffer);
    
    ezcl_device_memory_delete(ocelli_buffer);
    ezcl_device_memory_delete(ocellj_buffer);
    ezcl_device_memory_delete(olevel_buffer);
    ezcl_device_memory_delete(oval_buffer);
    
    ezcl_device_memory_delete(hhash);

    ezcl_device_memory_delete(hash_memory_indices_buffer);    
    
    free(hash_memory_indices);
    
    return time;
}

// returns the timing and the results in oval if run_tests flag is set
double cl_compact_hierarchical_remap (cell_list icells, cell_list ocells,
       intintHash_Factory *CLFactory, int run_tests){
                 
    cl_context       context = ezcl_get_context();
    cl_command_queue queue   = ezcl_get_command_queue();

    struct timeval timer;
#ifdef DETAILED_TIMING
    struct timeval timer1;
#endif

    // Input mesh
    cl_mem icelli_buffer = ezcl_device_memory_malloc(context, NULL, "icelli", icells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem icellj_buffer = ezcl_device_memory_malloc(context, NULL, "icellj", icells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ilevel_buffer = ezcl_device_memory_malloc(context, NULL, "ilevel", icells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ival_buffer   = ezcl_device_memory_malloc(context, NULL, "ival",   icells.ncells, sizeof(double), CL_MEM_READ_WRITE, 0);
    // TODO: ierr should be moved to local memory if possible; a buffer of this size is not necessary,
    // and should not be accessed through global memory
    cl_mem ierr_buffer = ezcl_device_memory_malloc(context, NULL, "ierr", icells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);

    ezcl_enqueue_write_buffer(queue, icelli_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.i,      NULL);
    ezcl_enqueue_write_buffer(queue, icellj_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.j,      NULL);
    ezcl_enqueue_write_buffer(queue, ilevel_buffer, CL_FALSE, 0, icells.ncells*sizeof(uint),   icells.level,  NULL);
    ezcl_enqueue_write_buffer(queue, ival_buffer,   CL_TRUE,  0, icells.ncells*sizeof(double), icells.values, NULL);
    
    
    
    // The output mesh
    cl_mem ocelli_buffer = ezcl_device_memory_malloc(context, NULL, "ocelli", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem ocellj_buffer = ezcl_device_memory_malloc(context, NULL, "ocellj", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem olevel_buffer = ezcl_device_memory_malloc(context, NULL, "olevel", ocells.ncells, sizeof(uint),   CL_MEM_READ_WRITE, 0);
    cl_mem oval_buffer   = ezcl_device_memory_malloc(context, NULL, "oval",   ocells.ncells, sizeof(double), CL_MEM_READ_WRITE, 0);

    ezcl_enqueue_write_buffer(queue, ocelli_buffer, CL_FALSE, 0, ocells.ncells*sizeof(uint), ocells.i,     NULL);
    ezcl_enqueue_write_buffer(queue, ocellj_buffer, CL_FALSE, 0, ocells.ncells*sizeof(uint), ocells.j,     NULL);
    ezcl_enqueue_write_buffer(queue, olevel_buffer, CL_TRUE,  0, ocells.ncells*sizeof(uint), ocells.level, NULL);
    
    //START TIMER
    cpu_timer_start(&timer);
#ifdef DETAILED_TIMING
    cpu_timer_start(&timer1);
#endif
    
    size_t global_work_size[1];
    size_t local_work_size[1];
    
    local_work_size[0] = TILE_SIZE;
    global_work_size[0] = ((local_work_size[0]+icells.ncells-1)/local_work_size[0])*local_work_size[0];

    size_t block_size = global_work_size[0]/local_work_size[0];
    cl_mem dev_redscratch = ezcl_device_memory_malloc(context, NULL, "dev_redscratch", block_size*(icells.levmax+1),
       sizeof(cl_uint), CL_MEM_READ_WRITE, 0);

    cl_mem num_at_level_buffer = ezcl_device_memory_malloc(context, NULL, "num_at_level_buffer", icells.levmax+1,
       sizeof(uint), CL_MEM_READ_WRITE, 0);

    size_t tile_size = local_work_size[0]*sizeof(cl_int)*(icells.levmax+1);

    //count_levels_stage1of2
    ezcl_set_kernel_arg(hierarchical_count_levels_stage1of2_kernel, 0, sizeof(uint),   &icells.ncells);
    ezcl_set_kernel_arg(hierarchical_count_levels_stage1of2_kernel, 1, sizeof(uint),   &icells.levmax);
    ezcl_set_kernel_arg(hierarchical_count_levels_stage1of2_kernel, 2, sizeof(cl_mem), &ilevel_buffer);
    ezcl_set_kernel_arg(hierarchical_count_levels_stage1of2_kernel, 3, sizeof(cl_mem), &num_at_level_buffer);
    ezcl_set_kernel_arg(hierarchical_count_levels_stage1of2_kernel, 4, sizeof(cl_mem), &dev_redscratch);
    ezcl_set_kernel_arg(hierarchical_count_levels_stage1of2_kernel, 5, tile_size,      NULL);

    ezcl_enqueue_ndrange_kernel(queue, hierarchical_count_levels_stage1of2_kernel, 1, 0, global_work_size, local_work_size, NULL);
    
    if (block_size > 1) {
       ezcl_set_kernel_arg(hierarchical_count_levels_stage2of2_kernel, 0, sizeof(uint),   &block_size);
       ezcl_set_kernel_arg(hierarchical_count_levels_stage2of2_kernel, 1, sizeof(uint),   &icells.levmax);
       ezcl_set_kernel_arg(hierarchical_count_levels_stage2of2_kernel, 2, sizeof(cl_mem), &num_at_level_buffer);
       ezcl_set_kernel_arg(hierarchical_count_levels_stage2of2_kernel, 3, sizeof(cl_mem), &dev_redscratch);
       ezcl_set_kernel_arg(hierarchical_count_levels_stage2of2_kernel, 4, tile_size,      NULL);

       ezcl_enqueue_ndrange_kernel(queue, hierarchical_count_levels_stage2of2_kernel, 1, 0, local_work_size, local_work_size, NULL);
    }

// This is temporary for checking -- below
    uint *num_at_level = (uint *)malloc((icells.levmax+1)*sizeof(uint));

    for (uint i = 0; i <= icells.levmax; i++) {
       num_at_level[i] = 0;
    }

    for (uint n = 0; n < icells.ncells; n++) {
       uint lev = icells.level[n];
       num_at_level[lev]++;
    }
/*
    for (uint i = 0; i <= icells.levmax; i++) {
       printf("DEBUG -- levsum for lev %u is %u\n",i,num_at_level[i]);
    }
*/    
// This is temporary for checking -- above


// Needed to allocate the hash tables
    ezcl_enqueue_read_buffer(queue, num_at_level_buffer, CL_TRUE, 0, (icells.levmax+1)*sizeof(uint), num_at_level, NULL);

/*
// This is temporary for checking -- below
    for (uint i = 0; i <= icells.levmax; i++) {
       printf("DEBUG -- levsum for lev %u is %u\n",i,num_at_level[i]);
    }
// This is temporary for checking -- above
*/

    // Every fourth cell will write a breadcrumb to the level above
    // lev must be int (not uint) to allow -1 for exit
    for (int lev = icells.levmax-1; lev >= 0; lev--) {
       num_at_level[lev] += num_at_level[lev+1]/4;
    }

#define HASH_CL_TYPE (LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID | IDENTITY_SENTINEL_PERFECT_CL_HASH_ID)
#define HASH_LOAD_FACTOR 0.3333333

    intintHash_Table** h_hashTable = (intintHash_Table **) malloc((icells.levmax+1)*sizeof(intintHash_Table *));

    uint *h_hashtype;
    if (DEBUG >= 2) {
       h_hashtype = (uint *)malloc((icells.levmax+1)*sizeof(uint));
    }

    //initialize 2d array
    for (uint i = 0; i <= icells.levmax; i++) {
        size_t hash_size = icells.ibasesize*two_to_the(i)*icells.jbasesize*two_to_the(i);
        h_hashTable[i] = intintHash_CreateTable(CLFactory, HASH_CL_TYPE, hash_size, num_at_level[i], HASH_LOAD_FACTOR);
        if (DEBUG >= 2) {
           h_hashtype[i] = intintHash_GetTableType(h_hashTable[i]);
           if (h_hashtype[i] == IDENTITY_PERFECT_CL_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"IDENTITY_PERFECT_CL_HASH_ID");
           } else if (h_hashtype[i] == IDENTITY_SENTINEL_PERFECT_CL_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"IDENTITY_SENTINEL_PERFECT_CL_HASH_ID");
           } else if (h_hashtype[i] == LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID");
           } else if (h_hashtype[i] == LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID) {
              printf("Type of hash for lev %d is %s\n",i,"LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID");
           }
        }

        //Empty Hash Table
        intintHash_EmptyTable(h_hashTable[i]);
    }

    ezcl_device_memory_delete(num_at_level_buffer);
    ezcl_device_memory_delete(dev_redscratch);

    free(num_at_level);
    if (DEBUG >=2){
       free(h_hashtype);
    }

    /*
    * Setup hash
    */

    ezcl_set_kernel_arg(hierarchical_compact_insert_kernel, 0, sizeof(uint),   &icells.ncells);
    ezcl_set_kernel_arg(hierarchical_compact_insert_kernel, 1, sizeof(uint),   &icells.ibasesize);
    ezcl_set_kernel_arg(hierarchical_compact_insert_kernel, 2, sizeof(cl_mem), &icelli_buffer);
    ezcl_set_kernel_arg(hierarchical_compact_insert_kernel, 3, sizeof(cl_mem), &icellj_buffer);
    ezcl_set_kernel_arg(hierarchical_compact_insert_kernel, 4, sizeof(cl_mem), &ilevel_buffer);
    for (int ilev = 0; ilev <= (int)icells.levmax; ilev++){
       ezcl_set_kernel_arg(hierarchical_compact_insert_kernel, 5+ilev, sizeof(cl_mem), intintHash_GetTableDataBufferPtr(h_hashTable[ilev]));
    }
    for (int ilev = icells.levmax+1; ilev <= 10; ilev++){
       ezcl_set_kernel_arg(hierarchical_compact_insert_kernel, 5+ilev, sizeof(cl_mem), NULL);
    }
    

    ezcl_enqueue_ndrange_kernel(queue, hierarchical_compact_insert_kernel, 1, 0, global_work_size, local_work_size, NULL);

    /*
    * Query hash
    */

    ezcl_set_kernel_arg(hierarchical_compact_probe_kernel, 0, sizeof(uint),   &ocells.ncells);
    ezcl_set_kernel_arg(hierarchical_compact_probe_kernel, 1, sizeof(uint),   &ocells.ibasesize);
    ezcl_set_kernel_arg(hierarchical_compact_probe_kernel, 2, sizeof(cl_mem), &ocelli_buffer);
    ezcl_set_kernel_arg(hierarchical_compact_probe_kernel, 3, sizeof(cl_mem), &ocellj_buffer);
    ezcl_set_kernel_arg(hierarchical_compact_probe_kernel, 4, sizeof(cl_mem), &olevel_buffer);
    ezcl_set_kernel_arg(hierarchical_compact_probe_kernel, 5, sizeof(cl_mem), &ival_buffer);
    ezcl_set_kernel_arg(hierarchical_compact_probe_kernel, 6, sizeof(cl_mem), &oval_buffer);
    for (int ilev = 0; ilev <= (int)ocells.levmax; ilev++){
       ezcl_set_kernel_arg(hierarchical_compact_probe_kernel, 7+ilev, sizeof(cl_mem), intintHash_GetTableDataBufferPtr(h_hashTable[ilev]));
    }
    for (int ilev = ocells.levmax+1; ilev <= 8; ilev++){
       ezcl_set_kernel_arg(hierarchical_compact_probe_kernel, 7+ilev, sizeof(cl_mem), NULL);
    }
    
    //TODO: See todo above. This should be local memory.
    ezcl_set_kernel_arg(hierarchical_compact_probe_kernel, 16, sizeof(cl_mem), &ierr_buffer);

    ezcl_enqueue_ndrange_kernel(queue, hierarchical_compact_probe_kernel, 1, 0, global_work_size, local_work_size, NULL);

    //END TIMER
    double time = cpu_timer_stop(timer);
    
    if (run_tests) ezcl_enqueue_read_buffer(queue, oval_buffer, CL_TRUE, 0, ocells.ncells*sizeof(double), ocells.values, NULL);

    ezcl_device_memory_delete(icelli_buffer);
    ezcl_device_memory_delete(icellj_buffer);
    ezcl_device_memory_delete(ilevel_buffer);
    ezcl_device_memory_delete(ival_buffer);
    
    ezcl_device_memory_delete(ocelli_buffer);
    ezcl_device_memory_delete(ocellj_buffer);
    ezcl_device_memory_delete(olevel_buffer);
    ezcl_device_memory_delete(oval_buffer);
    
    ezcl_device_memory_delete(ierr_buffer);
    
    return time;
}
#endif
