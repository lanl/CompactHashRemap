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

#pragma OPENCL EXTENSION cl_khr_fp64: enable
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable

typedef double real;
typedef unsigned int uint;

#define MAX_UINT 4294967295

/*
* Finds the min *AND* max -- modified from a paper
* by Bryan Catanzaro found here:
* http://developer.amd.com/resources/documentation-articles/articles-whitepapers/opencl-optimization-case-study-simple-reductions/
*/
__kernel void min_reduce_kern(__global real* buffer,
            __local real* min_scratch,
            __local real* max_scratch,
            __const uint length,
            __global real* result) {

  uint global_index = get_global_id(0);
  real min_accumulator;
  real max_accumulator;
  if(global_index<length){
    min_accumulator = buffer[global_index];
    max_accumulator = buffer[global_index];
    global_index += get_global_size(0);
  }
  // Loop sequentially over chunks of input vector
  while (global_index < length) {
    real element = buffer[global_index];
    min_accumulator = (min_accumulator < element) ? min_accumulator : element;
    max_accumulator = (max_accumulator > element) ? max_accumulator : element;
    global_index += get_global_size(0);
  } 
  
  // Perform parallel reduction
  uint local_index = get_local_id(0);
  min_scratch[local_index] = min_accumulator;
  max_scratch[local_index] = max_accumulator;
  global_index = get_global_id(0);
  barrier(CLK_LOCAL_MEM_FENCE);
  for(uint offset = get_local_size(0) / 2;
      offset > 0;
      offset = offset / 2) {
    if (local_index < offset&&global_index+offset<length) {
      real other_min = min_scratch[local_index + offset];
      real my_min = min_scratch[local_index];
      min_scratch[local_index] = (my_min < other_min) ? my_min : other_min;
      
      real other_max = max_scratch[local_index + offset];
      real my_max = max_scratch[local_index];
      max_scratch[local_index] = (my_max > other_max) ? my_max : other_max;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }
  
  if (local_index == 0) {
    result[get_group_id(0)*2] = min_scratch[0];
    result[get_group_id(0)*2+1] = max_scratch[0];
  }
}


// This should use atomic operators to sum up all the collisions into a buffer
__kernel void histogram_kern(__const uint length, __const real bin_size, __global real* buffer, volatile __global uint* collisions){
    uint global_index = get_global_id(0);
    while (global_index < length) {
        uint bin_index = (uint) ((buffer[global_index])/bin_size);
		atom_inc(&collisions[bin_index]);
		//collisions[global_index]+=bin_index;
        global_index += get_global_size(0);
    } 
}

// Initializes the array by subracting the min from all values
__kernel void init_kern(__const uint length, __global real* buffer, __const real min){
    uint global_index = get_global_id(0);
    while (global_index < length) {
        buffer[global_index] -= min;
        global_index += get_global_size(0);
    } 
}

// Sets all buffer values to the specified value
__kernel void mem_init_kern(__const uint length, __global uint* buffer, __const uint value){
    uint global_index = get_global_id(0);
    while (global_index < length) {
        buffer[global_index] = value;
        global_index += get_global_size(0);
    } 
}

/* Prefix scan */

__kernel void scan1(
	const uint isize,
	__global uint *ioffset,
	__local volatile uint *itile,
	__global const uint *temp) {
		
	const uint giX = get_global_id(0);
	const uint tiX = get_local_id(0);
	const uint ntX = get_local_size(0);
	const uint group_id = get_group_id(0);

        uint temp_val = temp[giX];

        itile[tiX] = temp_val;
	barrier(CLK_GLOBAL_MEM_FENCE);
	
	for(uint offset = ntX >> 1; offset > 32; offset >>= 1) {
		if(tiX < offset) {
			itile[tiX] += itile[tiX+offset];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

        if(giX >= isize) return;
    
    //  Unroll the remainder of the loop as 32 threads must proceed in lockstep.
    if (tiX < 32)
    {  itile[tiX] += itile[tiX+32];
       itile[tiX] += itile[tiX+16];
       itile[tiX] += itile[tiX+8];
       itile[tiX] += itile[tiX+4];
       itile[tiX] += itile[tiX+2];
       itile[tiX] += itile[tiX+1]; }
    
    if(tiX == 0) {
        ioffset[group_id] = itile[0];
    }
}

inline uint scan_warp_exclusive(__local volatile uint *input, const uint idx, const uint lane) {
    if (lane > 0 ) input[idx] += input[idx - 1];
    if (lane > 1 ) input[idx] += input[idx - 2];
    if (lane > 3 ) input[idx] += input[idx - 4];
    if (lane > 7 ) input[idx] += input[idx - 8];
    if (lane > 15) input[idx] += input[idx - 16];
    
    return (lane > 0) ? input[idx-1] : 0;
}

inline uint scan_warp_inclusive(__local volatile uint *input, const uint idx, const uint lane) {
    if (1) {
       if (lane > 0 ) input[idx] += input[idx - 1];
       if (lane > 1 ) input[idx] += input[idx - 2];
       if (lane > 3 ) input[idx] += input[idx - 4];
       if (lane > 7 ) input[idx] += input[idx - 8];
       if (lane > 15) input[idx] += input[idx - 16];
       return input[idx];
    }
}

inline uint scan_workgroup_exclusive(
    __local uint* itile,
    const uint tiX,
    const uint lane,
    const uint warpID) {
    
    // Step 1: scan each warp
    uint val = scan_warp_exclusive(itile, tiX, lane);
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Step 2: Collect per-warp sums
    if (lane == 31) itile[warpID] = itile[tiX];
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Step 3: Use 1st warp to scan per-warp sums
    if (warpID == 0) scan_warp_inclusive(itile, tiX, lane);
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Step 4: Accumulate results from Steps 1 and 3
    if (warpID > 0) val += itile[warpID-1];
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Step 6: Write and return the final result
    itile[tiX] = val;
    barrier(CLK_LOCAL_MEM_FENCE);
    
    return val;
}

__kernel void scan2(
    __local uint* itile,
    __global uint* ioffset,
    const uint size) {

    size_t tiX = get_local_id(0);
    const uint gID = get_group_id(0);
    const uint ntX = get_local_size(0);
    
    const uint lane = tiX & 31;
    const uint warpID = tiX >> 5;
    const uint EPT = (size+ntX-1)/ntX; //elements_per_thread;
    
    uint reduceValue = 0;
    
//  #pragma unroll 4
    for(uint i = 0; i < EPT; ++i)
    {
       uint offsetIdx = i * ntX + tiX;

#ifdef IS_NVIDIA
//     if (offsetIdx >= size) return;
#endif
        
       // Step 1: Read ntX elements from global (off-chip) memory to local memory (on-chip)
       uint input = 0;
       if (offsetIdx < size) input = ioffset[offsetIdx];           
       itile[tiX] = input;           
       barrier(CLK_LOCAL_MEM_FENCE);
        
       // Step 2: Perform scan on ntX elements
       uint val = scan_workgroup_exclusive(itile, tiX, lane, warpID);
        
       // Step 3: Propagate reduced result from previous block of ntX elements
       val += reduceValue;
        
       // Step 4: Write out data to global memory
       if (offsetIdx < size) ioffset[offsetIdx] = val;
     
       // Step 5: Choose reduced value for next iteration
       if (tiX == (ntX-1)) itile[tiX] = input + val;
       barrier(CLK_LOCAL_MEM_FENCE);
        
       reduceValue = itile[ntX-1];
       barrier(CLK_LOCAL_MEM_FENCE);
    }
}

inline uint do_element_pass(uint offsetIdx, uint ntX, uint tiX, uint lane, uint warpID,
      uint reduceValue, uint size, __global uint *ioffset, __local uint *itile) {
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 1: Read ntX elements from global (off-chip) memory to local memory (on-chip)
    uint input = 0;
    if (offsetIdx < size) input = ioffset[offsetIdx];           
    itile[tiX] = input;           
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Step 2: Perform scan on ntX elements
    uint val = scan_workgroup_exclusive(itile, tiX, lane, warpID);
   
    // Step 3: Propagate reduced result from previous block of ntX elements
    val += reduceValue;
  
    // Step 4: Write out data to global memory
    if (offsetIdx < size) ioffset[offsetIdx] = val;
  
    // Step 5: Choose reduced value for next iteration
    if (tiX == (ntX-1)) itile[tiX] = input + val;
    barrier(CLK_LOCAL_MEM_FENCE);

    reduceValue = itile[ntX-1];

    return(reduceValue);
}

__kernel void scan_lev(
    __local uint* itile,
    __global uint* ioffset,
    __global uint* workgroup_results,
    const uint size) {

    uint tiX = get_local_id(0);
    uint giX = get_global_id(0);
    const uint gID = get_group_id(0);
    const uint ntX = get_local_size(0);
    
    const uint lane = tiX & 31;
    const uint warpID = tiX >> 5;

    workgroup_results[gID] = 0;

    // Step 1: Read ntX elements from global (off-chip) memory to local memory (on-chip)
    uint input = 0;
    if (giX < size) input = ioffset[giX];
    itile[tiX] = input;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 2: Perform scan on ntX elements
    uint val = scan_workgroup_exclusive(itile, tiX, lane, warpID);

    // Step 3: Collect per-workgroup partial results
    workgroup_results[gID] = itile[tiX];
}

__kernel void scan_workgroup_results(
    __global uint* workgroup_results)
{
    uint tiX = get_local_id(0);
    
    const uint lane = tiX & 31;
    const uint warpID = tiX >> 5;

    // Step 4: Use 1st warp to scan workgroup_results
    //if (warpID == 0) scan_warp_inclusive(workgroup_results, tiX, lane);
}

__kernel void accumulate_workgroup_results()
{
    uint tiX = get_local_id(0);
    
    const uint warpID = tiX >> 5;

    // Step 5: Accumulate results from steps 2 and 4
    //uint val += itile[warpID-1];
}

__kernel void scan3 (
    const uint isize,
    __global const uint *ioffset,
    __local uint *itile,
    __global uint *temp,
    __global const real *arr,
    __global volatile real *sorted,
    __global uint *hash,
    const real bin_size) {
    
    const uint giX = get_global_id(0);
    const uint tiX = get_local_id(0);
    const uint group_id = get_group_id(0);

    const uint lane   = tiX & 31;
    const uint warpid = tiX >> 5;

    // Step 1: load global data into tile
    int temp_val = -1;
    if (giX < isize) temp_val = temp[giX];
    itile[tiX] = 0;
    if (temp_val >= 0) itile[tiX] = temp[tiX];
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 2: scan each warp
    uint val = scan_warp_exclusive(itile, tiX, lane);
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 3: Collect per-warp sums
    if (lane == 31) itile[warpid] = itile[tiX];
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 4: Use 1st warp to scan per-warp sums
    if (warpid == 0) scan_warp_inclusive(itile, tiX, lane);
    barrier(CLK_LOCAL_MEM_FENCE);

    // Step 5: Accumulate results from Steps 2 and 4
    if (warpid > 0) val += itile[warpid-1];
    barrier(CLK_LOCAL_MEM_FENCE);

    if (giX >= isize || temp_val < 0) return;

    // Step 6: Write and return the final result
    //itile[tiX] = val;
    //barrier(CLK_LOCAL_MEM_FENCE);
    
    //temp[giX] = itile[tiX]+ioffset[group_id];

    val += ioffset[group_id];
    
    temp[giX] = val;
    hash[giX] = MAX_UINT;
    barrier(CLK_GLOBAL_MEM_FENCE);
    
    /*
    * code seems to work correctly up until this point.
    *
    * After here, most of the code has been changed for debugging.
    */
    
    uint bin_index = (uint) ((arr[giX])/bin_size);
    uint number_stored = 0;
    uint hash_value = giX;
    
    /*
    * This block of code was intended to write to the required bin (using atomics)
    * unless a value was already there. If one was there, it moves the largest over
    * inside the bin until a free space is found, then inserts it. This step is done
    * with the hash instead of directly into the sorted array because atomics can only
    * be used with integer values.
    *
    * The program throws a kernel error instead.
    *
    * Note that large parts of the current running code were put in as part of 
    * debugging (such as putting strange integers into hash and temp at the bottom).
    */
    
    /*
    while(!number_stored){
        uint old_hash = atom_cmpxchg(&hash[bin_index],MAX_UINT,hash_value);
        if (old_hash == MAX_UINT){
            number_stored = 1;
        }else{
            if(arr[old_hash]>arr[hash_value]){
                // when the hash has already changed, someone else is operating and we need to restart (to avoid double write)
                if(atom_cmpxchg(&hash[bin_index],old_hash,hash_value)!=old_hash){
                    //continue;
                }
                hash_value = old_hash;
                
            }
        }
        
        bin_index++;
    }
    }
   
    
    
    sorted[giX] = ((arr[giX])/bin_size);
    */
    barrier(CLK_GLOBAL_MEM_FENCE);
    hash[temp[bin_index]] = giX;
    barrier(CLK_GLOBAL_MEM_FENCE);
    temp[bin_index] = giX;
    if(hash[giX] == MAX_UINT){
        //sorted[giX] = -1;
    }else{
        //sorted[giX] = arr[hash[giX]];
    }
    
    
}
