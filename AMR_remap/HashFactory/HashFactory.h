
/* Copyright (C) 1991-2012 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* We do support the IEC 559 math functionality, real and complex.  */
/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */
/* We do not support C11 <threads.h>.  */
/* Copyright 2013-14.  Los Alamos National Security, LLC. This material was produced
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
 * Under this license, it is required to include a reference to this work. We
 * request that each derivative work contain a reference to LANL Copyright 
 * Disclosure C14043/LA-CC-14-003 so that this work's impact can be roughly
 * measured. In addition, it is requested that a modifier is included as in
 * the following example:
 *
 * //<Uses | improves on | modified from> LANL Copyright Disclosure C14043/LA-CC-14-003
 *
 * This is LANL Copyright Disclosure C14043/LA-CC-14-003
 */
/**
 * @file   HashFactory.h
 * @author Peter Ahrens
 * @date   Thu Jun 6 2013 
 */
//
#ifndef HASHFACTORY_H
#define HASHFACTORY_H
//
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
//
#ifndef UINT_TYPE
#define UINT_TYPE
typedef unsigned int uint;
#endif
//
#ifdef HAVE_OPENCL
#ifdef __APPLE_CC__
#include <OpenCL/OpenCL.h>
#else
#include <CL/cl.h>
#endif
#else
typedef void *cl_mem;
typedef int cl_context;
typedef int cl_command_queue;
typedef int cl_program;
typedef int cl_device_id;
typedef int cl_int;
int clRetainContext(int context);
int clRetainCommandQueue(int command_queue);
int clGetContextInfo(int context, int param, size_t size, void *value,
		     size_t * size_ret);
int clReleaseContext(int context);
int clReleaseCommandQueue(int command_queue);
int clReleaseProgram(int program);
int clRetainProgram(int program);
int clRetainKernel(int kernel);
cl_mem clCreateBuffer(int context, int flags, size_t size, void *value,
		      int *size_ret);
int clEnqueueWriteBuffer(int command_queue, void *buffer, int blocking_write,
			 size_t offset, size_t cb, const void *ptr,
			 uint nevents, const int *wait_list, int *event);
int clEnqueueReadBuffer(int command_queue, void *buffer, int blocking_write,
			size_t offset, size_t cb, const void *ptr, uint nevents,
			const int *wait_list, int *event);
int clCreateKernel(int program, const char *kernel_name, int *errcode_ret);
int clReleaseKernel(int kernel);
int clReleaseMemObject(void *memobj);
int clSetKernelArg(int kernel, uint arg_index, size_t arg_size,
		   const void *arg_value);
int clGetKernelWorkGroupInfo(int kernel, int device, int param_name,
			     size_t size, void *value, size_t * size_ret);
int clEnqueueNDRangeKernel(int command_queue, int kernel, uint work_dim,
			   const size_t * offset, const size_t * size,
			   const size_t * local_size, uint nevents,
			   const int *wait_list, int *event);
int clFinish(int command_queue);
#endif
//
#include "CLHash_Utilities.h"
//
//
//
#define HASH_REPORT_NEVER /**/ 0
#define HASH_REPORT_CYCLE /**/ 1
#define HASH_REPORT_END /****/ 2
//
#define HASH_EXIT_CODE_NORMAL /****************/ -1
#define HASH_EXIT_CODE_ERROR /*****************/ -2
#define HASH_EXIT_CODE_OVERWRITE /*************/ -3
#define HASH_EXIT_CODE_KEY_DNE /***************/ -4
#define HASH_EXIT_CODE_CYCLE /*****************/ -5
#define HASH_EXIT_CODE_MAX_ENTRIES_EXCEEDED /**/ -6
#define HASH_EXIT_CODE_BUCKET_INDEX_OOB /******/ -7
//
#define HASH_SEARCH_CODE_MATCH /*****/ 0
#define HASH_SEARCH_CODE_MISMATCH /**/ 1
#define HASH_SEARCH_CODE_EMPTY /*****/ 2
//
#define HASH_NUM_HASHES /*********/ 12
#define HASH_NUM_C_HASHES /*******/ 4
#define HASH_NUM_CL_HASHES /******/ 4
#define HASH_NUM_OPENMP_HASHES /**/ 4
//
#define IDENTITY_PERFECT_HASH_ID /*******************/ 1
#define IDENTITY_SENTINEL_PERFECT_HASH_ID /**********/ 2
#define LCG_LINEAR_OPEN_COMPACT_HASH_ID /************/ 4
#define LCG_QUADRATIC_OPEN_COMPACT_HASH_ID /*********/ 8
//
#define IDENTITY_PERFECT_CL_HASH_ID /****************/ 16
#define IDENTITY_SENTINEL_PERFECT_CL_HASH_ID /*******/ 32
#define LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID /*********/ 64
#define LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID /******/ 128
//
#define IDENTITY_PERFECT_OPENMP_HASH_ID /************/ 256
#define IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID /***/ 512
#define LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID /*****/ 1024
#define LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID /**/ 2048
//
/**
 * HASH_ALL_C_HASHES hash table types that run on the host.
 */
#define HASH_ALL_C_HASHES (IDENTITY_PERFECT_HASH_ID | IDENTITY_SENTINEL_PERFECT_HASH_ID | LCG_LINEAR_OPEN_COMPACT_HASH_ID | LCG_QUADRATIC_OPEN_COMPACT_HASH_ID)
/**
 * HASH_ALL_CL_HASHES hash table types that run in an OpenCL context.
 */
#define HASH_ALL_CL_HASHES (IDENTITY_PERFECT_CL_HASH_ID | IDENTITY_SENTINEL_PERFECT_CL_HASH_ID | LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID | LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID)
/**
 * HASH_ALL_OPENMP_HASHES hash table types that run with OpenMP.
 */
#define HASH_ALL_OPENMP_HASHES (IDENTITY_PERFECT_OPENMP_HASH_ID | IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID | LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID | LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID)
/**
 * HASH_ALL_HASHES all hash table types.
 */
#define HASH_ALL_HASHES (HASH_ALL_C_HASHES | HASH_ALL_CL_HASHES | HASH_ALL_OPENMP_HASHES)
/**
 * HASH_NOSENTINEL_PERFECT_HASHES perfect hash table types that do not use a sentinel value to mark empty buckets.
 */
#define HASH_NOSENTINEL_PERFECT_HASHES (IDENTITY_PERFECT_HASH_ID | IDENTITY_PERFECT_CL_HASH_ID | IDENTITY_PERFECT_OPENMP_HASH_ID)
/**
 * HASH_SENTINEL_PERFECT_HASHES perfect hash table types that use a sentinel value to mark empty buckets.
 */
#define HASH_SENTINEL_PERFECT_HASHES (IDENTITY_SENTINEL_PERFECT_HASH_ID | IDENTITY_SENTINEL_PERFECT_CL_HASH_ID | IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID)
/**
 * HASH_PERFECT_HASHES perfect hash table types.
 */
#define HASH_PERFECT_HASHES (HASH_NOSENTINEL_PERFECT_HASHES | HASH_SENTINEL_PERFECT_HASHES)
/**
 * HASH_QUADRATIC_COMPACT_HASHES compact hash table types that use a quadratic probe sequence.
 */
#define HASH_QUADRATIC_COMPACT_HASHES (LCG_QUADRATIC_OPEN_COMPACT_HASH_ID | LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID | LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID)
/**
 * HASH_LINEAR_COMPACT_HASHES compact hash table types that use a linear probe sequence.
 */
#define HASH_LINEAR_COMPACT_HASHES (LCG_LINEAR_OPEN_COMPACT_HASH_ID | LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID | LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID)
/**
 * HASH_COMPACT_HASHES compact hash table types that use a linear probe sequence.
 */
#define HASH_COMPACT_HASHES (HASH_QUADRATIC_COMPACT_HASHES | HASH_LINEAR_COMPACT_HASHES)
//
#define HASH_DEFAULT_LOCAL_WORK_SIZE /*******/ 64
#define HASH_DEFAULT_LOAD_FACTOR /*******/ 0.3
#define HASH_MIN_LOAD_FACTOR /*******/ 0.0009
#define HASH_PERFECT_COMPACT_SWITCH_FACTOR /*******/ 20
#define HASH_LCG_A /*******/ 2147483629
#define HASH_LCG_C /*******/ 2147483587
#define HASH_LCG_M /*******/ 2147483647
//
#define HASH_BUCKET_STATUS_EMPTY /**/ -1
#define HASH_BUCKET_STATUS_FULL /***/ -2
#define HASH_BUCKET_STATUS_LOCK /***/ -3
#ifdef __cplusplus
extern "C" {
#endif
/**
 * Hash_ExitCodeString will return a string representation of the given exit
 * code.
 * 
 * @param exitCode
 * 
 * @return A string representation of that exit code.
 */
	char *Hash_ExitCodeString(int exitCode);
/**
 * Hash_ExitCodeDebug will print a string representation of the given exit code
 * if it is not EXIT_CODE_NORMAL.
 * 
 * @param exitCode
 */
	void Hash_ExitCodeDebug(int exitCode);
/**
 * Hash_SetReportLevel sets a static report level variable in hash.c. It should 
 * be called before hash tables are created. 
 *
 * @param level The level of data collection desired.
 */
	void Hash_SetReportLevel(int level);
/**
 * Hash_GetReportLevel gets this variable.
 *
 * @return The current report level.
 *         Special Values: Meaning
 *         HASH_REPORT_NEVER: The hash table will not collect data. This is the
 *                       default. A call to Hash_Report will return an empty 
 *                       string.
 *         HASH_REPORT_END: The hash table will collect data and calls to
 *                     Hash_Report will return summary information.
 *         HASH_REPORT_CYCLE: The hash table will collect data and calls to 
 *                       Hash_Report will return information related to the
 *                       last important call.
 */
	int Hash_GetReportLevel();
	const char *Hash_GetKernelSourceString();
	int smallestProthPrimeAbove(int N);
	int largestProthPrimeUnder(int N);
	typedef struct uintuintHash_Table_ uintuintHash_Table;
	typedef struct uintuintCLHash_Table_ uintuintCLHash_Table;
	typedef struct uintuintHash_Factory_ uintuintHash_Factory;
	typedef struct uintuintCLHash_Factory_ uintuintCLHash_Factory;
	uintuintHash_Factory *uintuintHash_CreateFactory(int HashTypes,
							 uint * emptyValue,
							 size_t localWorkSize,
							 cl_context * context,
							 cl_command_queue *
							 queue);
	int uintuintHash_DestroyFactory(uintuintHash_Factory * factory);
	uintuintHash_Table *uintuintHash_CreateTable(uintuintHash_Factory *
						     factory, int hashTypes,
						     size_t keyRange,
						     size_t numEntries,
						     float loadFactor);
	int uintuintHash_SetupTable(uintuintHash_Table * table);
	int uintuintHash_EmptyTable(uintuintHash_Table * table);
	int uintuintHash_DestroyTable(uintuintHash_Table * table);
	cl_mem uintuintHash_GetTableDataBuffer(uintuintHash_Table * table);
	cl_mem *uintuintHash_GetTableDataBufferPtr(uintuintHash_Table * table);
	int uintuintHash_GetTableType(uintuintHash_Table * table);
	int uintuintHash_Query(uintuintHash_Table * table, size_t numKeys,
			       uint * keys, uint * valuesOutput);
	int uintuintHash_QuerySingle(uintuintHash_Table * table, uint key,
				     uint * valueOutput);
	int uintuintHash_Insert(uintuintHash_Table * table, size_t numEntries,
				uint * keys, uint * values);
	int uintuintHash_InsertSingle(uintuintHash_Table * table, uint key,
				      uint value);
	int uintuintHash_InsertNoOverwrite(uintuintHash_Table * table,
					   size_t numEntries, uint * keys,
					   uint * values);
	int uintuintHash_InsertSingleNoOverwrite(uintuintHash_Table * table,
						 uint key, uint value);
	int uintuintHash_BufferQuery(uintuintHash_Table * table, size_t numKeys,
				     cl_mem keys, cl_mem valuesOutput);
	int uintuintHash_BufferInsert(uintuintHash_Table * table,
				      size_t numEntries, cl_mem keys,
				      cl_mem values);
	int uintuintHash_BufferInsertNoOverwrite(uintuintHash_Table * table,
						 size_t numEntries, cl_mem keys,
						 cl_mem values);
	static inline unsigned int uintuintHash_CompressIdentity(char data,
								 int hashCode) {
		return hashCode;
	}
	typedef struct uintuintHash_CompressLCGData {
		long unsigned int a;
		long unsigned int c;
		unsigned int m;
		unsigned int n;
	} uintuintHash_CompressLCGData;
	static inline unsigned int
	    uintuintHash_CompressLCG(uintuintHash_CompressLCGData
				     compressLCGData, int hashCode) {
		return ((compressLCGData.a * hashCode +
			 compressLCGData.c) % compressLCGData.m) %
		    compressLCGData.n;
	}
	int uintuintIdentityPerfectHash_CreateFactory(uintuintHash_Factory *
						      factory, int hashIndex);
	int uintuintIdentityPerfectHash_DestroyFactory(uintuintHash_Factory *
						       factory, int hashIndex);
	uintuintHash_Table
	    *uintuintIdentityPerfectHash_CreateTable(uintuintHash_Factory *
						     factory, int hashIndex,
						     size_t keyRange,
						     size_t numEntries,
						     float loadFactor);
	int uintuintIdentityPerfectHash_InitTable(uintuintHash_Table * table,
						  va_list args);
	int uintuintIdentityPerfectHash_DestroyTable(uintuintHash_Table *
						     table);
	char *uintuintIdentityPerfectHash_Report(uintuintHash_Table * table);
	int uintuintIdentityPerfectHash_SetupTable(uintuintHash_Table * table);
	int uintuintIdentityPerfectHash_EmptyTable(uintuintHash_Table * table);
	int uintuintIdentityPerfectHash_Query(uintuintHash_Table * table,
					      size_t numKeys, uint * keys,
					      uint * valuesOutput);
	int uintuintIdentityPerfectHash_QuerySingle(uintuintHash_Table * table,
						    uint key,
						    uint * valueOutput);
	int uintuintIdentityPerfectHash_Insert(uintuintHash_Table * table,
					       size_t numEntries, uint * keys,
					       uint * values);
	int uintuintIdentityPerfectHash_InsertSingle(uintuintHash_Table * table,
						     uint key, uint value);
	int uintuintIdentityPerfectHash_InsertNoOverwrite(uintuintHash_Table *
							  table,
							  size_t numEntries,
							  uint * keys,
							  uint * values);
	int uintuintIdentityPerfectHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintIdentityPerfectCLHash_CreateFactory(uintuintHash_Factory *
							factory, int hashIndex);
	int uintuintIdentityPerfectCLHash_DestroyFactory(uintuintHash_Factory *
							 factory,
							 int hashIndex);
	uintuintHash_Table
	    *uintuintIdentityPerfectCLHash_CreateTable(uintuintHash_Factory *
						       factory, int hashIndex,
						       size_t keyRange,
						       size_t numEntries,
						       float loadFactor);
	int uintuintIdentityPerfectCLHash_InitTable(uintuintHash_Table * table,
						    va_list args);
	int uintuintIdentityPerfectCLHash_DestroyTable(uintuintHash_Table *
						       table);
	char *uintuintIdentityPerfectCLHash_Report(uintuintHash_Table * table);
	int uintuintIdentityPerfectCLHash_SetupTable(uintuintHash_Table *
						     table);
	int uintuintIdentityPerfectCLHash_EmptyTable(uintuintHash_Table *
						     table);
	int uintuintIdentityPerfectCLHash_Query(uintuintHash_Table * table,
						size_t numKeys, uint * keys,
						uint * valuesOutput);
	int uintuintIdentityPerfectCLHash_QuerySingle(uintuintHash_Table *
						      table, uint key,
						      uint * valueOutput);
	int uintuintIdentityPerfectCLHash_Insert(uintuintHash_Table * table,
						 size_t numEntries, uint * keys,
						 uint * values);
	int uintuintIdentityPerfectCLHash_InsertSingle(uintuintHash_Table *
						       table, uint key,
						       uint value);
	int uintuintIdentityPerfectCLHash_InsertNoOverwrite(uintuintHash_Table *
							    table,
							    size_t numEntries,
							    uint * keys,
							    uint * values);
	int uintuintIdentityPerfectCLHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintIdentityPerfectCLHash_BufferQuery(uintuintHash_Table *
						      table, size_t numKeys,
						      cl_mem keysBuffer,
						      cl_mem
						      valuesOutputBuffer);
	int uintuintIdentityPerfectCLHash_BufferInsert(uintuintHash_Table *
						       table, size_t numEntries,
						       cl_mem keysBuffer,
						       cl_mem valuesBuffer);
	int uintuintIdentityPerfectCLHash_BufferInsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintIdentityPerfectOpenMPHash_CreateFactory(uintuintHash_Factory
							    * factory,
							    int hashIndex);
	int uintuintIdentityPerfectOpenMPHash_DestroyFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	uintuintHash_Table
	    *uintuintIdentityPerfectOpenMPHash_CreateTable(uintuintHash_Factory
							   * factory,
							   int hashIndex,
							   size_t keyRange,
							   size_t numEntries,
							   float loadFactor);
	int uintuintIdentityPerfectOpenMPHash_InitTable(uintuintHash_Table *
							table, va_list args);
	int uintuintIdentityPerfectOpenMPHash_DestroyTable(uintuintHash_Table *
							   table);
	char *uintuintIdentityPerfectOpenMPHash_Report(uintuintHash_Table *
						       table);
	int uintuintIdentityPerfectOpenMPHash_SetupTable(uintuintHash_Table *
							 table);
	int uintuintIdentityPerfectOpenMPHash_EmptyTable(uintuintHash_Table *
							 table);
	int uintuintIdentityPerfectOpenMPHash_Query(uintuintHash_Table * table,
						    size_t numKeys, uint * keys,
						    uint * valuesOutput);
	int uintuintIdentityPerfectOpenMPHash_QuerySingle(uintuintHash_Table *
							  table, uint key,
							  uint * valueOutput);
	int uintuintIdentityPerfectOpenMPHash_Insert(uintuintHash_Table * table,
						     size_t numEntries,
						     uint * keys,
						     uint * values);
	int uintuintIdentityPerfectOpenMPHash_InsertSingle(uintuintHash_Table *
							   table, uint key,
							   uint value);
	int uintuintIdentityPerfectOpenMPHash_InsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, uint * keys,
	     uint * values);
	int uintuintIdentityPerfectOpenMPHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintIdentityPerfectOpenMPHash_BufferQuery(uintuintHash_Table *
							  table, size_t numKeys,
							  cl_mem keysBuffer,
							  cl_mem
							  valuesOutputBuffer);
	int uintuintIdentityPerfectOpenMPHash_BufferInsert(uintuintHash_Table *
							   table,
							   size_t numEntries,
							   cl_mem keysBuffer,
							   cl_mem valuesBuffer);
	int uintuintIdentityPerfectOpenMPHash_BufferInsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintIdentitySentinelPerfectHash_CreateFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	int uintuintIdentitySentinelPerfectHash_DestroyFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	uintuintHash_Table
	    *uintuintIdentitySentinelPerfectHash_CreateTable
	    (uintuintHash_Factory * factory, int hashIndex, size_t keyRange,
	     size_t numEntries, float loadFactor);
	int uintuintIdentitySentinelPerfectHash_InitTable(uintuintHash_Table *
							  table, va_list args);
	int uintuintIdentitySentinelPerfectHash_DestroyTable(uintuintHash_Table
							     * table);
	char *uintuintIdentitySentinelPerfectHash_Report(uintuintHash_Table *
							 table);
	int uintuintIdentitySentinelPerfectHash_SetupTable(uintuintHash_Table *
							   table);
	int uintuintIdentitySentinelPerfectHash_EmptyTable(uintuintHash_Table *
							   table);
	int uintuintIdentitySentinelPerfectHash_Query(uintuintHash_Table *
						      table, size_t numKeys,
						      uint * keys,
						      uint * valuesOutput);
	int uintuintIdentitySentinelPerfectHash_QuerySingle(uintuintHash_Table *
							    table, uint key,
							    uint * valueOutput);
	int uintuintIdentitySentinelPerfectHash_Insert(uintuintHash_Table *
						       table, size_t numEntries,
						       uint * keys,
						       uint * values);
	int uintuintIdentitySentinelPerfectHash_InsertSingle(uintuintHash_Table
							     * table, uint key,
							     uint value);
	int uintuintIdentitySentinelPerfectHash_InsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, uint * keys,
	     uint * values);
	int uintuintIdentitySentinelPerfectHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintIdentitySentinelPerfectCLHash_CreateFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	int uintuintIdentitySentinelPerfectCLHash_DestroyFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	uintuintHash_Table
	    *uintuintIdentitySentinelPerfectCLHash_CreateTable
	    (uintuintHash_Factory * factory, int hashIndex, size_t keyRange,
	     size_t numEntries, float loadFactor);
	int uintuintIdentitySentinelPerfectCLHash_InitTable(uintuintHash_Table *
							    table,
							    va_list args);
	int uintuintIdentitySentinelPerfectCLHash_DestroyTable
	    (uintuintHash_Table * table);
	char *uintuintIdentitySentinelPerfectCLHash_Report(uintuintHash_Table *
							   table);
	int uintuintIdentitySentinelPerfectCLHash_SetupTable(uintuintHash_Table
							     * table);
	int uintuintIdentitySentinelPerfectCLHash_EmptyTable(uintuintHash_Table
							     * table);
	int uintuintIdentitySentinelPerfectCLHash_Query(uintuintHash_Table *
							table, size_t numKeys,
							uint * keys,
							uint * valuesOutput);
	int uintuintIdentitySentinelPerfectCLHash_QuerySingle(uintuintHash_Table
							      * table, uint key,
							      uint *
							      valueOutput);
	int uintuintIdentitySentinelPerfectCLHash_Insert(uintuintHash_Table *
							 table,
							 size_t numEntries,
							 uint * keys,
							 uint * values);
	int uintuintIdentitySentinelPerfectCLHash_InsertSingle
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintIdentitySentinelPerfectCLHash_InsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, uint * keys,
	     uint * values);
	int uintuintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintIdentitySentinelPerfectCLHash_BufferQuery(uintuintHash_Table
							      * table,
							      size_t numKeys,
							      cl_mem keysBuffer,
							      cl_mem
							      valuesOutputBuffer);
	int uintuintIdentitySentinelPerfectCLHash_BufferInsert
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintIdentitySentinelPerfectOpenMPHash_CreateFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	int uintuintIdentitySentinelPerfectOpenMPHash_DestroyFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	uintuintHash_Table
	    *uintuintIdentitySentinelPerfectOpenMPHash_CreateTable
	    (uintuintHash_Factory * factory, int hashIndex, size_t keyRange,
	     size_t numEntries, float loadFactor);
	int uintuintIdentitySentinelPerfectOpenMPHash_InitTable
	    (uintuintHash_Table * table, va_list args);
	int uintuintIdentitySentinelPerfectOpenMPHash_DestroyTable
	    (uintuintHash_Table * table);
	char *uintuintIdentitySentinelPerfectOpenMPHash_Report
	    (uintuintHash_Table * table);
	int uintuintIdentitySentinelPerfectOpenMPHash_SetupTable
	    (uintuintHash_Table * table);
	int uintuintIdentitySentinelPerfectOpenMPHash_EmptyTable
	    (uintuintHash_Table * table);
	int uintuintIdentitySentinelPerfectOpenMPHash_Query(uintuintHash_Table *
							    table,
							    size_t numKeys,
							    uint * keys,
							    uint *
							    valuesOutput);
	int uintuintIdentitySentinelPerfectOpenMPHash_QuerySingle
	    (uintuintHash_Table * table, uint key, uint * valueOutput);
	int uintuintIdentitySentinelPerfectOpenMPHash_Insert(uintuintHash_Table
							     * table,
							     size_t numEntries,
							     uint * keys,
							     uint * values);
	int uintuintIdentitySentinelPerfectOpenMPHash_InsertSingle
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintIdentitySentinelPerfectOpenMPHash_InsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, uint * keys,
	     uint * values);
	int uintuintIdentitySentinelPerfectOpenMPHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintIdentitySentinelPerfectOpenMPHash_BufferQuery
	    (uintuintHash_Table * table, size_t numKeys, cl_mem keysBuffer,
	     cl_mem valuesOutputBuffer);
	int uintuintIdentitySentinelPerfectOpenMPHash_BufferInsert
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintIdentitySentinelPerfectOpenMPHash_BufferInsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintLCGLinearOpenCompactHash_CreateFactory(uintuintHash_Factory
							   * factory,
							   int hashIndex);
	int uintuintLCGLinearOpenCompactHash_DestroyFactory(uintuintHash_Factory
							    * factory,
							    int hashIndex);
	uintuintHash_Table
	    *uintuintLCGLinearOpenCompactHash_CreateTable(uintuintHash_Factory *
							  factory,
							  int hashIndex,
							  size_t keyRange,
							  size_t numEntries,
							  float loadFactor);
	int uintuintLCGLinearOpenCompactHash_InitTable(uintuintHash_Table *
						       table, va_list args);
	int uintuintLCGLinearOpenCompactHash_DestroyTable(uintuintHash_Table *
							  table);
	char *uintuintLCGLinearOpenCompactHash_Report(uintuintHash_Table *
						      table);
	int uintuintLCGLinearOpenCompactHash_SetupTable(uintuintHash_Table *
							table);
	int uintuintLCGLinearOpenCompactHash_EmptyTable(uintuintHash_Table *
							table);
	int uintuintLCGLinearOpenCompactHash_Query(uintuintHash_Table * table,
						   size_t numKeys, uint * keys,
						   uint * valuesOutput);
	int uintuintLCGLinearOpenCompactHash_QuerySingle(uintuintHash_Table *
							 table, uint key,
							 uint * valueOutput);
	int uintuintLCGLinearOpenCompactHash_Insert(uintuintHash_Table * table,
						    size_t numEntries,
						    uint * keys, uint * values);
	int uintuintLCGLinearOpenCompactHash_InsertSingle(uintuintHash_Table *
							  table, uint key,
							  uint value);
	int uintuintLCGLinearOpenCompactHash_InsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, uint * keys,
	     uint * values);
	int uintuintLCGLinearOpenCompactHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintLCGLinearOpenCompactCLHash_CreateFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	int uintuintLCGLinearOpenCompactCLHash_DestroyFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	uintuintHash_Table
	    *uintuintLCGLinearOpenCompactCLHash_CreateTable(uintuintHash_Factory
							    * factory,
							    int hashIndex,
							    size_t keyRange,
							    size_t numEntries,
							    float loadFactor);
	int uintuintLCGLinearOpenCompactCLHash_InitTable(uintuintHash_Table *
							 table, va_list args);
	int uintuintLCGLinearOpenCompactCLHash_DestroyTable(uintuintHash_Table *
							    table);
	char *uintuintLCGLinearOpenCompactCLHash_Report(uintuintHash_Table *
							table);
	int uintuintLCGLinearOpenCompactCLHash_SetupTable(uintuintHash_Table *
							  table);
	int uintuintLCGLinearOpenCompactCLHash_EmptyTable(uintuintHash_Table *
							  table);
	int uintuintLCGLinearOpenCompactCLHash_Query(uintuintHash_Table * table,
						     size_t numKeys,
						     uint * keys,
						     uint * valuesOutput);
	int uintuintLCGLinearOpenCompactCLHash_QuerySingle(uintuintHash_Table *
							   table, uint key,
							   uint * valueOutput);
	int uintuintLCGLinearOpenCompactCLHash_Insert(uintuintHash_Table *
						      table, size_t numEntries,
						      uint * keys,
						      uint * values);
	int uintuintLCGLinearOpenCompactCLHash_InsertSingle(uintuintHash_Table *
							    table, uint key,
							    uint value);
	int uintuintLCGLinearOpenCompactCLHash_InsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, uint * keys,
	     uint * values);
	int uintuintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintLCGLinearOpenCompactCLHash_BufferQuery(uintuintHash_Table *
							   table,
							   size_t numKeys,
							   cl_mem keysBuffer,
							   cl_mem
							   valuesOutputBuffer);
	int uintuintLCGLinearOpenCompactCLHash_BufferInsert(uintuintHash_Table *
							    table,
							    size_t numEntries,
							    cl_mem keysBuffer,
							    cl_mem
							    valuesBuffer);
	int uintuintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintLCGLinearOpenCompactOpenMPHash_CreateFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	int uintuintLCGLinearOpenCompactOpenMPHash_DestroyFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	uintuintHash_Table
	    *uintuintLCGLinearOpenCompactOpenMPHash_CreateTable
	    (uintuintHash_Factory * factory, int hashIndex, size_t keyRange,
	     size_t numEntries, float loadFactor);
	int uintuintLCGLinearOpenCompactOpenMPHash_InitTable(uintuintHash_Table
							     * table,
							     va_list args);
	int uintuintLCGLinearOpenCompactOpenMPHash_DestroyTable
	    (uintuintHash_Table * table);
	char *uintuintLCGLinearOpenCompactOpenMPHash_Report(uintuintHash_Table *
							    table);
	int uintuintLCGLinearOpenCompactOpenMPHash_SetupTable(uintuintHash_Table
							      * table);
	int uintuintLCGLinearOpenCompactOpenMPHash_EmptyTable(uintuintHash_Table
							      * table);
	int uintuintLCGLinearOpenCompactOpenMPHash_Query(uintuintHash_Table *
							 table, size_t numKeys,
							 uint * keys,
							 uint * valuesOutput);
	int uintuintLCGLinearOpenCompactOpenMPHash_QuerySingle
	    (uintuintHash_Table * table, uint key, uint * valueOutput);
	int uintuintLCGLinearOpenCompactOpenMPHash_Insert(uintuintHash_Table *
							  table,
							  size_t numEntries,
							  uint * keys,
							  uint * values);
	int uintuintLCGLinearOpenCompactOpenMPHash_InsertSingle
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintLCGLinearOpenCompactOpenMPHash_InsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, uint * keys,
	     uint * values);
	int uintuintLCGLinearOpenCompactOpenMPHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintLCGLinearOpenCompactOpenMPHash_BufferQuery
	    (uintuintHash_Table * table, size_t numKeys, cl_mem keysBuffer,
	     cl_mem valuesOutputBuffer);
	int uintuintLCGLinearOpenCompactOpenMPHash_BufferInsert
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintLCGLinearOpenCompactOpenMPHash_BufferInsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintLCGQuadraticOpenCompactHash_CreateFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	int uintuintLCGQuadraticOpenCompactHash_DestroyFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	uintuintHash_Table
	    *uintuintLCGQuadraticOpenCompactHash_CreateTable
	    (uintuintHash_Factory * factory, int hashIndex, size_t keyRange,
	     size_t numEntries, float loadFactor);
	int uintuintLCGQuadraticOpenCompactHash_InitTable(uintuintHash_Table *
							  table, va_list args);
	int uintuintLCGQuadraticOpenCompactHash_DestroyTable(uintuintHash_Table
							     * table);
	char *uintuintLCGQuadraticOpenCompactHash_Report(uintuintHash_Table *
							 table);
	int uintuintLCGQuadraticOpenCompactHash_SetupTable(uintuintHash_Table *
							   table);
	int uintuintLCGQuadraticOpenCompactHash_EmptyTable(uintuintHash_Table *
							   table);
	int uintuintLCGQuadraticOpenCompactHash_Query(uintuintHash_Table *
						      table, size_t numKeys,
						      uint * keys,
						      uint * valuesOutput);
	int uintuintLCGQuadraticOpenCompactHash_QuerySingle(uintuintHash_Table *
							    table, uint key,
							    uint * valueOutput);
	int uintuintLCGQuadraticOpenCompactHash_Insert(uintuintHash_Table *
						       table, size_t numEntries,
						       uint * keys,
						       uint * values);
	int uintuintLCGQuadraticOpenCompactHash_InsertSingle(uintuintHash_Table
							     * table, uint key,
							     uint value);
	int uintuintLCGQuadraticOpenCompactHash_InsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, uint * keys,
	     uint * values);
	int uintuintLCGQuadraticOpenCompactHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintLCGQuadraticOpenCompactCLHash_CreateFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	int uintuintLCGQuadraticOpenCompactCLHash_DestroyFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	uintuintHash_Table
	    *uintuintLCGQuadraticOpenCompactCLHash_CreateTable
	    (uintuintHash_Factory * factory, int hashIndex, size_t keyRange,
	     size_t numEntries, float loadFactor);
	int uintuintLCGQuadraticOpenCompactCLHash_InitTable(uintuintHash_Table *
							    table,
							    va_list args);
	int uintuintLCGQuadraticOpenCompactCLHash_DestroyTable
	    (uintuintHash_Table * table);
	char *uintuintLCGQuadraticOpenCompactCLHash_Report(uintuintHash_Table *
							   table);
	int uintuintLCGQuadraticOpenCompactCLHash_SetupTable(uintuintHash_Table
							     * table);
	int uintuintLCGQuadraticOpenCompactCLHash_EmptyTable(uintuintHash_Table
							     * table);
	int uintuintLCGQuadraticOpenCompactCLHash_Query(uintuintHash_Table *
							table, size_t numKeys,
							uint * keys,
							uint * valuesOutput);
	int uintuintLCGQuadraticOpenCompactCLHash_QuerySingle(uintuintHash_Table
							      * table, uint key,
							      uint *
							      valueOutput);
	int uintuintLCGQuadraticOpenCompactCLHash_Insert(uintuintHash_Table *
							 table,
							 size_t numEntries,
							 uint * keys,
							 uint * values);
	int uintuintLCGQuadraticOpenCompactCLHash_InsertSingle
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, uint * keys,
	     uint * values);
	int uintuintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintLCGQuadraticOpenCompactCLHash_BufferQuery(uintuintHash_Table
							      * table,
							      size_t numKeys,
							      cl_mem keysBuffer,
							      cl_mem
							      valuesOutputBuffer);
	int uintuintLCGQuadraticOpenCompactCLHash_BufferInsert
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_CreateFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_DestroyFactory
	    (uintuintHash_Factory * factory, int hashIndex);
	uintuintHash_Table
	    *uintuintLCGQuadraticOpenCompactOpenMPHash_CreateTable
	    (uintuintHash_Factory * factory, int hashIndex, size_t keyRange,
	     size_t numEntries, float loadFactor);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_InitTable
	    (uintuintHash_Table * table, va_list args);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_DestroyTable
	    (uintuintHash_Table * table);
	char *uintuintLCGQuadraticOpenCompactOpenMPHash_Report
	    (uintuintHash_Table * table);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_SetupTable
	    (uintuintHash_Table * table);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_EmptyTable
	    (uintuintHash_Table * table);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_Query(uintuintHash_Table *
							    table,
							    size_t numKeys,
							    uint * keys,
							    uint *
							    valuesOutput);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_QuerySingle
	    (uintuintHash_Table * table, uint key, uint * valueOutput);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_Insert(uintuintHash_Table
							     * table,
							     size_t numEntries,
							     uint * keys,
							     uint * values);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_InsertSingle
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_InsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, uint * keys,
	     uint * values);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_InsertSingleNoOverwrite
	    (uintuintHash_Table * table, uint key, uint value);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_BufferQuery
	    (uintuintHash_Table * table, size_t numKeys, cl_mem keysBuffer,
	     cl_mem valuesOutputBuffer);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_BufferInsert
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
	int uintuintLCGQuadraticOpenCompactOpenMPHash_BufferInsertNoOverwrite
	    (uintuintHash_Table * table, size_t numEntries, cl_mem keysBuffer,
	     cl_mem valuesBuffer);
#ifdef __cplusplus
}
#endif				/* __cplusplus */
#endif				/* HASHFACTORY_H */
