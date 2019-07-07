/* Copyright 2015-19.  Triad National Security, LLC. This material was produced
 * under U.S. Government contract 89233218CNA000001 for Los Alamos National 
 * Laboratory (LANL), which is operated by Triad National Security, LLC
 * for the U.S. Department of Energy. The U.S. Government has rights to use,
 * reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
 * TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
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
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <float.h>
#include <limits.h>
#include <pthread.h>
#include <omp.h>
#include <CL/cl.h>
//#include <sort.h>

#define MAX_UINT 4294967295

#define TILE_SIZE 128

#define PRINT_UINT_BUFFER(length, buffer) \
do{ \
    uint *arr = (uint*)malloc(((uint)length)*sizeof(uint)); \
    error = clEnqueueReadBuffer(queue, (cl_mem)buffer, CL_TRUE, 0, ((uint)length)*sizeof(uint), arr, 0, NULL, NULL); \
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__); \
    for(uint i = 0; i < ((uint)length); i++){ \
        printf("arr %i: %u\t",i,  arr[i]); \
    } \
    free(arr); \
}while (0)

#define PRINT_DOUBLE_BUFFER(length, buffer) \
do{ \
    double *arr = (double*)malloc(((uint)length)*sizeof(double)); \
    error = clEnqueueReadBuffer(queue, (cl_mem)buffer, CL_TRUE, 0, ((uint)length)*sizeof(double), arr, 0, NULL, NULL); \
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__); \
    for(uint i = 0; i < ((uint)length); i++){ \
        printf("arr: %3.3f\t", arr[i]); \
    } \
    free(arr); \
}while (0)

typedef double real;

typedef unsigned int uint;
struct timeval timer;
double t1, t2;


cl_platform_id platform_id;
cl_device_id device_id;
cl_context context;
cl_command_queue queue;
cl_program program;
cl_kernel scan1_kernel, scan2_kernel, scan3_kernel, mem_init_kernel, reduce_kernel, init_kernel, histogram_kernel;
cl_int error;



//int compare (const void * a, const void * b) { return ((int)(*(double*)a - *(double*)b )); }
uint compare (const void * a, const void * b)
{
  if ( *(double*)a <  *(double*)b ) return -1;
  if ( *(double*)a >  *(double*)b ) return 1;
  return 0;
}

#define ValType double
#define IS_LESS(v1, v2)  (v1 < v2)
 
void siftDown( ValType *a, int start, int count);
 
#define SWAP(r,s)  do{ValType t=r; r=s; s=t; } while(0)
 

cl_mem procrastisort_gpu(uint length, cl_mem xcoor_buffer) {
            //prefix scanning
			//solid section
	    	/*numD = logTwoC(solidLength);
			#pragma omp parallel
			{
			
			for (int i = 0; i < numD; i++) {
				int jInc = solidLength/(1 << (i + 1));
				#pragma omp for
				for (int j = 0; j < jInc; j++ ) {
					int realJ = j * (1 << (i + 1));
					num[realJ + (1 << (i + 1)) - 1] = num[realJ + (1 << i) - 1] + num[realJ + (1 << (i + 1)) - 1]; 
				}
			}
			
			solidTotal = num[solidLength - 1];
			//printf("solidTotal: %i\n", solidTotal);
			if (numD > 0) {
				num[solidLength - 1] = 0;
				
				for (int i = numD - 1; i >= 0; i--) {
					int jInc = solidLength/(1 << (i +1));
					#pragma omp for
					for (int j = 0; j < jInc; j++) {
						int realJ = j * (1 << (i + 1));
						int t = num[realJ + (1 << i) - 1];
						num[realJ + (1 << i) - 1] = num[realJ + (1 << (i + 1)) - 1];
						num[realJ + (1 << (i + 1)) - 1] = t + num[realJ + (1 << (i + 1)) - 1];
					}
				}
			}
			}
			
			//padded section
			if (solidLength == paddedLength){
			
			} else if (solidLength + 1 == paddedLength) {
				num[solidLength] = 0;
			} else {
				numD = logTwoC(paddedLength - solidLength);
				
				for (int i = 0; i < numD; i++) {
					for (int j = solidLength; (j < paddedLength - 1); (j += ((int)pow(2, (i + 1)))) ) {
						num[j + ((int)pow(2,(i + 1))) - 1] = num[j + ((int)pow(2, (i))) - 1] + num[j + ((int)pow(2, (i + 1))) - 1]; 
					}
				}
		
				num[paddedLength - 1] = 0;
	
				for (int i = numD - 1; i >= 0; i--) {
					for (int j = solidLength; j < paddedLength; j += (int)pow(2, (i + 1))) {
						uint t = num[j + ((int)pow(2, i)) - 1];
						num[j + (int)pow(2, i) - 1] = num[j + (int)pow(2, (i + 1)) - 1];
						num[j + (int)pow(2, (i + 1)) - 1] = t + num[j + (int)pow(2, (i + 1)) - 1];
					}
				}
			}
			
			
			if (solidTotal > 0) {
				for (int i = solidLength; i < length; i++) {
					num[i] += solidTotal;
				}
			}*/
	/* 
	* Min/max finder
	*/
	cl_mem num_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, length*sizeof(cl_uint), NULL, &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
	reduce_kernel = clCreateKernel(program, "min_reduce_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    size_t global_work_size[1];
    size_t local_work_size[1];
    
    local_work_size[0] = TILE_SIZE;
    global_work_size[0] = TILE_SIZE;
    
    size_t num_extremes = (global_work_size[0]/local_work_size[0])*2;
    
    cl_mem extremes_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, num_extremes*sizeof(real), NULL, &error);
    if (error != CL_SUCCESS) {
       //printf("Error is %d at line %d\n",error,__LINE__);
       exit(0);
    } 
     
    error = clSetKernelArg(reduce_kernel, 0, sizeof(cl_mem), (void*)&xcoor_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(reduce_kernel, 1, local_work_size[0]*sizeof(real),NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(reduce_kernel, 2, local_work_size[0]*sizeof(real),NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(reduce_kernel, 3, sizeof(cl_uint), &length);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(reduce_kernel, 4, sizeof(cl_mem), (void*)&extremes_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    cl_event reduce_event;
 
    error = clEnqueueNDRangeKernel(queue, reduce_kernel, 1, 0, global_work_size, local_work_size, 0, NULL, &reduce_event);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    real *extremes_arr = (real*)malloc(num_extremes*sizeof(real));
    error = clEnqueueReadBuffer(queue, extremes_buffer, CL_TRUE, 0, num_extremes*sizeof(real), extremes_arr, 0, NULL, &reduce_event);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    real min = extremes_arr[0];
    real max = extremes_arr[1];
    clReleaseMemObject(extremes_buffer);
	free(extremes_arr);
	
	/*
	* Initialize array
	*/
	
	init_kernel = clCreateKernel(program, "init_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    local_work_size[0] = TILE_SIZE;
    global_work_size[0] = TILE_SIZE;
    
    error = clSetKernelArg(init_kernel, 0, sizeof(cl_uint), &length);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(init_kernel, 1, sizeof(cl_mem),(void*)&xcoor_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(init_kernel, 2, sizeof(cl_double), &min);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    
    cl_event init_event;
 
    error = clEnqueueNDRangeKernel(queue, init_kernel, 1, 0, global_work_size, local_work_size, 0, NULL, &init_event);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    clWaitForEvents(1,&init_event);
	
    /*
    * Prefix scan
    */
    
    real range = max - min;
    uint numRec;
    uint binIndex;
    real binSize = (range)/(length-1);
    
    // Find the number of collisions
    
    mem_init_kernel = clCreateKernel(program,"mem_init_kern",&error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    histogram_kernel = clCreateKernel(program, "histogram_kern", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    local_work_size[0] = TILE_SIZE;
    global_work_size[0] = TILE_SIZE;
    
    // first set all histogram values to 0
    
    error = clSetKernelArg(mem_init_kernel, 0, sizeof(cl_uint), &length);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(mem_init_kernel, 1, sizeof(cl_mem),(void*)&num_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    uint histogram_init_value = 0;
    error = clSetKernelArg(mem_init_kernel, 2, sizeof(cl_uint),&histogram_init_value);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clFinish(queue);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    
    
    
    cl_event mem_init_event;
 
    error = clEnqueueNDRangeKernel(queue, mem_init_kernel, 1, 0, global_work_size, local_work_size, 0, NULL, &mem_init_event);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    
    
    // then count the histogram
    
    error = clSetKernelArg(histogram_kernel, 0, sizeof(cl_uint), &length);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(histogram_kernel, 1, sizeof(cl_double), &binSize);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(histogram_kernel, 2, sizeof(cl_mem),(void*)&xcoor_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(histogram_kernel, 3, sizeof(cl_mem),(void*)&num_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    
    cl_event histogram_event;
 
    error = clEnqueueNDRangeKernel(queue, histogram_kernel, 1, 0, global_work_size, local_work_size, 0, NULL, &histogram_event);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clFlush(queue);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);

    // and now the scans
    
    scan1_kernel = clCreateKernel(program, "scan1", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    scan2_kernel = clCreateKernel(program, "scan2", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    scan3_kernel = clCreateKernel(program, "scan3", &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    /* scan 1 */
    
    global_work_size[0] = ((length+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];
 
    uint group_size = (uint)(global_work_size[0]/local_work_size[0]);
    
    cl_mem ioffset_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, group_size*sizeof(uint), NULL, &error);
    if (error != CL_SUCCESS) {
       //printf("Error is %d at line %d\n",error,__LINE__);
       //clReleaseMemObject(hash_buffer);
       return(NULL);
    }
  
    error = clSetKernelArg(scan1_kernel, 0, sizeof(cl_uint), &length);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan1_kernel, 1, sizeof(cl_mem), (void*)&ioffset_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan1_kernel, 2, local_work_size[0]*sizeof(uint), NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan1_kernel, 3, sizeof(cl_mem), (void*)&num_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
  
    cl_event scan1_event;
    
    error = clEnqueueNDRangeKernel(queue, scan1_kernel, 1, 0, global_work_size, local_work_size, 0, NULL, &scan1_event);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clFlush(queue);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    


    /* scan 2 */
    global_work_size[0] = local_work_size[0];

    cl_event scan2_event;

        
    uint elements_per_thread = (group_size+local_work_size[0]-1)/local_work_size[0];
                
    error = clSetKernelArg(scan2_kernel, 0, local_work_size[0]*sizeof(uint), NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan2_kernel, 1, sizeof(cl_mem), (void*)&ioffset_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan2_kernel, 2, sizeof(uint), &group_size);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    error = clEnqueueNDRangeKernel(queue, scan2_kernel, 1, 0, global_work_size, local_work_size, 0, NULL, &scan2_event);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clFinish(queue);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
 
    /* scan 3 */
    
    cl_mem sorted_buffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length*sizeof(real), NULL, &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    cl_mem hash_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, length*sizeof(uint), NULL, &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    global_work_size[0] = ((length+local_work_size[0]-1)/local_work_size[0])*local_work_size[0];
        
    error = clSetKernelArg(scan3_kernel, 0, sizeof(cl_uint), &length);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan3_kernel, 1, sizeof(cl_mem), (void*)&ioffset_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan3_kernel, 2, local_work_size[0]*sizeof(uint), NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan3_kernel, 3, sizeof(cl_mem), (void*)&num_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan3_kernel, 4, sizeof(cl_mem), (void *)&xcoor_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan3_kernel, 5, sizeof(cl_mem), (void *)&sorted_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan3_kernel, 6, sizeof(cl_mem), (void *)&hash_buffer);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clSetKernelArg(scan3_kernel, 7, sizeof(cl_double), &binSize);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    
    cl_event scan3_event;
    
    if (clEnqueueNDRangeKernel(queue, scan3_kernel, 1, 0, global_work_size, local_work_size, 0, NULL, &scan3_event) != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    clFinish(queue);
    PRINT_DOUBLE_BUFFER(length, xcoor_buffer);
    printf("\n\n");
    PRINT_UINT_BUFFER(length, num_buffer);
    printf("\n\n");
    PRINT_UINT_BUFFER(length, hash_buffer);
    printf("\n\n");
    clReleaseMemObject(hash_buffer);
    clReleaseMemObject(ioffset_buffer);
    clReleaseMemObject(num_buffer);
    return sorted_buffer;
}
	
double* hashsort( uint length, double *arr, double min_diff, double min_val, double max_val ) {
    uint hash_size;
    int *hash=NULL;
    double *sorted=NULL;
    
    sorted = (double*)malloc(length*sizeof(double));

    //create hash table with buckets of size min_diff 
    //   -- +2.5 rounds up and adds one space to either side
    hash_size = (uint)((max_val - min_val)/min_diff + 2.5);
    hash = (int*)malloc(hash_size*sizeof(int));

    //set all elements of hash array to -1
    memset(hash, -1, hash_size*sizeof(int));
    
    for(uint i = 0; i < length; i++) {
       //place index of current arr element into hash according to where the arr value
        hash[(int)((arr[i]-min_val)/min_diff)] = i;
    }
    
    int count=0;
    for(uint i = 0; i < hash_size; i++) {
        if(hash[i] >= 0) {
            //sweep through hash and put set values in a sorted array
            sorted[count] = arr[hash[i]];
            count++;
        }
    }
    
    free(hash);
    return sorted;
}

/* generate a randomly mixed up array with size size to be stored in pointer. the elements will have a minimum value min, and
    the difference between elements when sorted will be between mindx and maxdx. the maximum value is recorded in max. */
void generate_array(uint size, double *ptr, double mindx, double maxdx, double min, double *max ) {
    
    double swap;
    int index, front = 0;
    double running_min = maxdx;
        
    ptr[0] = min;        //start the array using the minimum value
    
    /* for each element, add a random value between mindx and maxdx to the previous element's value */
    for(int i = 1; i < size; i++) {
        ptr[i] = ptr[i-1] + mindx + ((double)rand() * (maxdx - mindx) / (double)RAND_MAX);
        if(ptr[i]-ptr[i-1] < running_min) running_min = ptr[i]-ptr[i-1];
    }

    *max = ptr[size-1];                    //set the max value to the last element's value
    //*max = min + (size-1) * maxdx;    //force the range for timings isolating a different variable
    
    /* Mix up the array by selecting elements from shrinking front portion of array and placing them on back end of array */
    for(int i = 0; (i < size) && (size - i != 0) ; i++) {
        index = rand() % (size - i - front) + front;
        swap = ptr[size-i-1];
        ptr[size-i-1] = ptr[index];
        ptr[index] = swap;
    }
}


int main (int argc, char** argv) {
	
	//printf("starting main\n");
	double qsortTime = 0;
	double procrastisortGTime = 0;
	double hashsortTime = 0;
	uint length = atoi (argv[1]);
	uint numRep = atoi (argv[2]);
	double mindx =1, maxdx = 3, min = 100, max=0;
	double *arr = (double*)malloc(length*sizeof(double));
	generate_array(length, arr, mindx, maxdx, min, &max);
	
	for (int i = 0; i < numRep; i++) {
	
	double *sorted=NULL, *sort_test=NULL;
	

	
	//printf("max: %f,\tmin: %f\n", max, min);
	/*for (int i=0; i<length; i++){
		printf("%f, \n", arr[i]);
	}
	printf("\n");*/
	
	//printf("generated array successfully\n");
	
    sorted = (double*)malloc(length*sizeof(double));
    memcpy(sorted, arr, length*sizeof(double));
    gettimeofday(&timer, NULL);
    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
    //heapsort(sorted, length);
    qsort(sorted, length, sizeof(double), compare);
    gettimeofday(&timer, NULL);
    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
    //printf("qsort time: %.6lf,\n", t2 - t1);
    qsortTime += (t2-t1);
	
	 error = clGetPlatformIDs(1,&platform_id,NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    error = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU,1, &device_id, NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    context = clCreateContext(0,1,&device_id,NULL,NULL,&error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    FILE* program_file = fopen ("procrastisort_kern.cl", "r");
    fseek(program_file,0,SEEK_END);
    size_t program_size =ftell (program_file);
    rewind(program_file);
    
    char* program_buffer = (char*)malloc (program_size+1);
    program_buffer[program_size] = '\0';
    fread (program_buffer,sizeof(char),program_size, program_file);
    fclose(program_file);

    queue = clCreateCommandQueue(context,device_id,0,&error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    program = clCreateProgramWithSource(context, 1, (const char**) &program_buffer, &program_size, &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    
    error = clBuildProgram(program,1,&device_id,NULL,NULL,NULL);
    if (error == CL_BUILD_PROGRAM_FAILURE) {
        // Determine the size of the log
        size_t log_size;
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

        // Allocate memory for the log
        char *log = (char *) malloc(log_size);

        // Get the log
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
    
        // Print the log
        printf("%s\n", log);
    }

    cl_mem xcoor_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, length*sizeof(double), NULL, &error);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
    error = clEnqueueWriteBuffer(queue, xcoor_buffer, CL_TRUE, 0, length*sizeof(double), arr, 0, NULL, NULL);
    if (error != CL_SUCCESS) printf("Error is %d at line %d\n",error,__LINE__);
	
	gettimeofday(&timer, NULL);
    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
	cl_mem sort_buffer = procrastisort_gpu(length, xcoor_buffer); 
	gettimeofday(&timer, NULL);
    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
    
    procrastisortGTime += (t2-t1);
    
    PRINT_DOUBLE_BUFFER(length,sort_buffer);
	printf("\n");
	
	/*sort_test=(double*)malloc(length*sizeof(double));
	memcpy(sort_test, arr, length*sizeof(double));
	gettimeofday(&timer, NULL);
    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
	sort_test = procrastisortD(length, arr, min, max); 
	gettimeofday(&timer, NULL);
    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
    
    
    //printf("procrastisort time: %.6lf,\n", t2 - t1);
    procrastisortDTime += (t2-t1);
    
    //icount=0;
    for(uint i = 0; i < length; i++) {
       if (sort_test[i] != sorted[i]) {
          //printf("Check failed for procrastisort CPU index %d procrastisort value %lf gold standard %lf\n",i,sort_test[i],sorted[i]);
		icount++;
       }
    }
    printf("errors: %i\n", icount);
    
	free(sort_test);
	sort_test=NULL;*/
    
    sort_test=(double*)malloc(length*sizeof(double));
	memcpy(sort_test, arr, length*sizeof(double));
	gettimeofday(&timer, NULL);
    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
	//sort_test = hashsort (length, arr, mindx, min, max); 
	gettimeofday(&timer, NULL);
    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
    
    free(sort_test);
    sort_test=NULL;
    
    hashsortTime += (t2 - t1);
    //printf("hashsort time: %.6lf,\n", t2 - t1);
    

    
   
    
	free(sorted);
	}
	free(arr);
	
	printf("gpsort: %f\thashsort: %f\tqsort: %f\n", (procrastisortGTime/numRep), (hashsortTime/numRep), (qsortTime/numRep));
	return 0;
}


