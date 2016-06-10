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
int intintIdentityPerfectCLHash_InsertSingle(__global char *tableData,
					     int key, int value);
int intintIdentityPerfectCLHash_InnerInsertSingle(__global char *tableData,
						  int key, int value);
int intintHash_InsertSingle(__global char *tableData, int key, int value);
int intintIdentityPerfectCLHash_InnerQuery(__global char *tableData,
					   unsigned int numKeys,
					   __global int *keys,
					   __global int *valuesOutput);
int intintIdentityPerfectCLHash_InnerQuerySingle(__global char *tableData,
						 int key,
						 __global int *valueOutput);
int intintIdentityPerfectCLHash_InnerInsert(__global char *tableData,
					    unsigned int numEntries,
					    __global int *keys,
					    __global int *values);
int intintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(__global char
							     *tableData,
							     int key,
							     int value);
int intintIdentityPerfectCLHash_InnerInsertNoOverwrite(__global char *tableData,
						       unsigned int numEntries,
						       __global int *keys,
						       __global int *values);
int intintIdentityPerfectCLHash_QuerySingle(__global char *tableData, int key,
					    __global int *valueOutput);
int intintIdentityPerfectCLHash_QuerySingle(__global char *tableData, int key,
					    __global int *valueOutput);
int intintIdentityPerfectCLHash_Query(__global char *tableData, size_t numKeys,
				      __global int *keys,
				      __global int *valuesOutput);
int intintIdentityPerfectCLHash_Insert(__global char *tableData,
				       size_t numEntries, __global int *keys,
				       __global int *values);
int intintIdentityPerfectCLHash_InsertSingleNoOverwrite(__global char
							*tableData, int key,
							int value);
int intintIdentityPerfectCLHash_InsertNoOverwrite(__global char *tableData,
						  size_t numEntries,
						  __global int *keys,
						  __global int *values);
int intintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite(__global char
							       *tableData,
							       unsigned int
							       numEntries,
							       __global int
							       *keys,
							       __global int
							       *values);
int intintIdentitySentinelPerfectCLHash_InnerQuerySingle(__global char
							 *tableData, int key,
							 __global int
							 *valueOutput);
int intintIdentitySentinelPerfectCLHash_InnerQuery(__global char *tableData,
						   unsigned int numKeys,
						   __global int *keys,
						   __global int *valuesOutput);
int intintIdentitySentinelPerfectCLHash_InnerInsertSingle(__global char
							  *tableData, int key,
							  int value);
int intintIdentitySentinelPerfectCLHash_InnerInsert(__global char *tableData,
						    unsigned int numEntries,
						    __global int *keys,
						    __global int *values);
int intintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite(__global
								     char
								     *tableData,
								     int key,
								     int value);
int intintIdentitySentinelPerfectCLHash_QuerySingle(__global char *tableData,
						    int key,
						    __global int *valueOutput);
int intintIdentitySentinelPerfectCLHash_Query(__global char *tableData,
					      size_t numKeys,
					      __global int *keys,
					      __global int *valuesOutput);
int intintIdentitySentinelPerfectCLHash_InsertSingle(__global char *tableData,
						     int key, int value);
int intintIdentitySentinelPerfectCLHash_Insert(__global char *tableData,
					       size_t numEntries,
					       __global int *keys,
					       __global int *values);
int intintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite(__global char
								*tableData,
								int key,
								int value);
int intintIdentitySentinelPerfectCLHash_InsertNoOverwrite(__global char
							  *tableData,
							  size_t numEntries,
							  __global int *keys,
							  __global int *values);
int intintLCGLinearOpenCompactCLHash_InnerQuerySingle(__global char *tableData,
						      int key,
						      __global int
						      *valueOutput);
int intintLCGLinearOpenCompactCLHash_QuerySingle(__global char *tableData,
						 int key,
						 __global int *valueOutput);
int intintLCGLinearOpenCompactCLHash_Query(__global char *tableData,
					   size_t numKeys, __global int *keys,
					   __global int *valuesOutput);
int intintLCGLinearOpenCompactCLHash_InsertSingle(__global char *tableData,
						  int key, int value);
int intintLCGLinearOpenCompactCLHash_Insert(__global char *tableData,
					    size_t numEntries,
					    __global int *keys,
					    __global int *values);
int intintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite(__global char
							     *tableData,
							     int key,
							     int value);
int intintLCGLinearOpenCompactCLHash_InsertNoOverwrite(__global char *tableData,
						       size_t numEntries,
						       __global int *keys,
						       __global int *values);
int intintLCGLinearOpenCompactCLHash_InnerQuery(__global char *tableData,
						unsigned int numKeys,
						__global int *keys,
						__global int *valuesOutput);
int intintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(__global char
							    *tableData,
							    unsigned int
							    numEntries,
							    __global int *keys,
							    __global int
							    *values);
int intintLCGLinearOpenCompactCLHash_InnerInsertSingle(__global char *tableData,
						       int key, int value);
int intintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global char
								  *tableData,
								  int key,
								  int value);
int intintLCGLinearOpenCompactCLHash_InnerInsert(__global char *tableData,
						 unsigned int numEntries,
						 __global int *keys,
						 __global int *values);
int intintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(__global char
							 *tableData, int key,
							 __global int
							 *valueOutput);
int intintLCGQuadraticOpenCompactCLHash_InnerQuery(__global char *tableData,
						   unsigned int numKeys,
						   __global int *keys,
						   __global int *valuesOutput);
int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(__global char
							  *tableData, int key,
							  int value);
int intintLCGQuadraticOpenCompactCLHash_InnerInsert(__global char *tableData,
						    unsigned int numEntries,
						    __global int *keys,
						    __global int *values);
int intintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global
								     char
								     *tableData,
								     int key,
								     int value);
int intintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite(__global char
							       *tableData,
							       unsigned int
							       numEntries,
							       __global int
							       *keys,
							       __global int
							       *values);
int intintLCGQuadraticOpenCompactCLHash_QuerySingle(__global char *tableData,
						    int key,
						    __global int *valueOutput);
int intintLCGQuadraticOpenCompactCLHash_Query(__global char *tableData,
					      size_t numKeys,
					      __global int *keys,
					      __global int *valuesOutput);
int intintLCGQuadraticOpenCompactCLHash_InsertSingle(__global char *tableData,
						     int key, int value);
int intintLCGQuadraticOpenCompactCLHash_Insert(__global char *tableData,
					       size_t numEntries,
					       __global int *keys,
					       __global int *values);
int intintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite(__global char
								*tableData,
								int key,
								int value);
int intintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(__global char
							  *tableData,
							  size_t numEntries,
							  __global int *keys,
							  __global int *values);
int intintHash_Query(__global char *tableData, unsigned int numKeys,
		     __global int *keys, __global int *valuesOutput);
int intintHash_QuerySingle(__global char *tableData, int key,
			   __global int *valueOutput);
int intintHash_Insert(__global char *tableData, unsigned int numEntries,
		      __global int *keys, __global int *values);
int intintHash_InsertNoOverwrite(__global char *tableData,
				 unsigned int numEntries, __global int *keys,
				 __global int *values);
int intintHash_InsertSingleNoOverwrite(__global char *tableData, int key,
				       int value);
#define HASH_REPORT_NEVER /**/ 0
#define HASH_REPORT_CYCLE /**/ 1
#define HASH_REPORT_END /****/ 2
//
#define HASH_EXIT_CODE_NORMAL /****************/ UINT_MAX
#define HASH_EXIT_CODE_ERROR /*****************/ UINT_MAX-1
#define HASH_EXIT_CODE_OVERWRITE /*************/ UINT_MAX-2
#define HASH_EXIT_CODE_KEY_DNE /***************/ UINT_MAX-3
#define HASH_EXIT_CODE_CYCLE /*****************/ UINT_MAX-4
#define HASH_EXIT_CODE_MAX_ENTRIES_EXCEEDED /**/ UINT_MAX-5
#define HASH_EXIT_CODE_BUCKET_INDEX_OOB /******/ UINT_MAX-6
//
#define HASH_SEARCH_CODE_MATCH /*****/ 0
#define HASH_SEARCH_CODE_MISMATCH /**/ 1
#define HASH_SEARCH_CODE_EMPTY /*****/ 2
//
#define IDENTITY_PERFECT_CL_HASH_ID /****************/ 16
#define IDENTITY_SENTINEL_PERFECT_CL_HASH_ID /*******/ 32
#define LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID /*********/ 64
#define LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID /******/ 128
//
#define HASH_BUCKET_STATUS_EMPTY /**/ -1
#define HASH_BUCKET_STATUS_FULL /***/ -2
#define HASH_BUCKET_STATUS_LOCK /***/ -3
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
static inline unsigned int uintuintHash_CompressLCG(uintuintHash_CompressLCGData
						    compressLCGData,
						    int hashCode) {
	return ((compressLCGData.a * hashCode +
		 compressLCGData.c) % compressLCGData.m) % compressLCGData.n;
}

typedef struct uintuintIdentityPerfectCLHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
} uintuintIdentityPerfectCLHash_TableData;
typedef struct uintuintIdentityPerfectCLHash_Bucket {
	uint key;
	uint value;
} uintuintIdentityPerfectCLHash_Bucket;
int uintuintIdentityPerfectCLHash_InnerQuerySingle(__global char *tableData,
						   uint key,
						   __local uint *
						   valueOutput) {
	__global uintuintIdentityPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectCLHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((__global
					    uintuintIdentityPerfectCLHash_TableData
					    *) tableData)->compressFuncData,
					  key);
	if ((buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
		if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_MISMATCH;
		}
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int uintuintIdentityPerfectCLHash_InnerQuery(__global char *tableData,
					     unsigned int numKeys,
					     __global uint * keys,
					     __local uint * valuesOutput) {
	__global uintuintIdentityPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectCLHash_TableData)];
	uint key;
	__local uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    uintuintHash_CompressIdentity(((__global
						    uintuintIdentityPerfectCLHash_TableData
						    *) tableData)->
						  compressFuncData, key);
		if ((buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
			if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
			} else {
				exitCode = HASH_SEARCH_CODE_MISMATCH;
			}
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int uintuintIdentityPerfectCLHash_InnerInsertSingle(__global char *tableData,
						    uint key, uint value) {
	__global uintuintIdentityPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectCLHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((__global
					    uintuintIdentityPerfectCLHash_TableData
					    *) tableData)->compressFuncData,
					  key);
	if (((buckets[index].key ==
	      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
					   key,
					   HASH_BUCKET_STATUS_EMPTY) :
	     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
		if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_MISMATCH;
		}
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int uintuintIdentityPerfectCLHash_InnerInsert(__global char *tableData,
					      unsigned int numEntries,
					      __global uint * keys,
					      __global uint * values) {
	__global uintuintIdentityPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((__global
						    uintuintIdentityPerfectCLHash_TableData
						    *) tableData)->
						  compressFuncData, key);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
			if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
			} else {
				exitCode = HASH_SEARCH_CODE_MISMATCH;
			}
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int uintuintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(__global char
							       *tableData,
							       uint key,
							       uint value) {
	__global uintuintIdentityPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectCLHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((__global
					    uintuintIdentityPerfectCLHash_TableData
					    *) tableData)->compressFuncData,
					  key);
	if (((buckets[index].key ==
	      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
					   key,
					   HASH_BUCKET_STATUS_EMPTY) :
	     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
		if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_MISMATCH;
		}
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int uintuintIdentityPerfectCLHash_InnerInsertNoOverwrite(__global char
							 *tableData,
							 unsigned int
							 numEntries,
							 __global uint * keys,
							 __global uint *
							 values) {
	__global uintuintIdentityPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentityPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((__global
						    uintuintIdentityPerfectCLHash_TableData
						    *) tableData)->
						  compressFuncData, key);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) != HASH_BUCKET_STATUS_EMPTY) {
			if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
			} else {
				exitCode = HASH_SEARCH_CODE_MISMATCH;
			}
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int uintuintIdentityPerfectCLHash_QuerySingle(__global char *tableData,
					      uint key,
					      __local uint * valueOutput) {
	return uintuintIdentityPerfectCLHash_InnerQuerySingle(tableData, key,
							      valueOutput);
}
int uintuintIdentityPerfectCLHash_Query(__global char *tableData,
					size_t numKeys, __global uint * keys,
					__local uint * valuesOutput) {
	return uintuintIdentityPerfectCLHash_InnerQuery(tableData, numKeys,
							keys, valuesOutput);
}
int uintuintIdentityPerfectCLHash_InsertSingle(__global char *tableData,
					       uint key, uint value) {
	return uintuintIdentityPerfectCLHash_InnerInsertSingle(tableData, key,
							       value);
}
int uintuintIdentityPerfectCLHash_Insert(__global char *tableData,
					 size_t numEntries,
					 __global uint * keys,
					 __global uint * values) {
	return uintuintIdentityPerfectCLHash_InnerInsert(tableData, numEntries,
							 keys, values);
}
int uintuintIdentityPerfectCLHash_InsertSingleNoOverwrite(__global char
							  *tableData, uint key,
							  uint value) {
	return
	    uintuintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite
	    (tableData, key, value);
}
int uintuintIdentityPerfectCLHash_InsertNoOverwrite(__global char *tableData,
						    size_t numEntries,
						    __global uint * keys,
						    __global uint * values) {
	return uintuintIdentityPerfectCLHash_InnerInsertNoOverwrite(tableData,
								    numEntries,
								    keys,
								    values);
}
__kernel void uintuintIdentityPerfectCLHash_RangeQuerySingle(__global char
							     *tableData,
							     unsigned int
							     numQueries,
							     __global uint *
							     keys,
							     __local uint *
							     valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	uintuintIdentityPerfectCLHash_InnerQuerySingle(tableData, keys[i],
						       valuesOutput + i);
}
__kernel void uintuintIdentityPerfectCLHash_RangeQuery(__global char *tableData,
						       unsigned int numQueries,
						       unsigned int numKeys,
						       __global uint * keys,
						       __local uint *
						       valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	uintuintIdentityPerfectCLHash_InnerQuery(tableData, numKeys,
						 keys + (i * numKeys),
						 valuesOutput + (i * numKeys));
}
__kernel void uintuintIdentityPerfectCLHash_RangeInsertSingle(__global char
							      *tableData,
							      unsigned int
							      numInsertions,
							      __global uint *
							      keys,
							      __global uint *
							      values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintIdentityPerfectCLHash_InnerInsertSingle(tableData, keys[i],
							values[i]);
}
__kernel void uintuintIdentityPerfectCLHash_RangeInsert(__global char
							*tableData,
							unsigned int
							numInsertions,
							unsigned int numEntries,
							__global uint * keys,
							__global uint *
							values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintIdentityPerfectCLHash_InnerInsert(tableData, numEntries,
						  keys + (i * numEntries),
						  values + (i * numEntries));
}
__kernel void
uintuintIdentityPerfectCLHash_RangeInsertSingleNoOverwrite(__global char
							   *tableData,
							   unsigned int
							   numInsertions,
							   __global uint * keys,
							   __global uint *
							   values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintIdentityPerfectCLHash_InnerInsertSingleNoOverwrite(tableData,
								   keys[i],
								   values[i]);
}
__kernel void uintuintIdentityPerfectCLHash_RangeInsertNoOverwrite(__global char
								   *tableData,
								   unsigned int
								   numInsertions,
								   unsigned int
								   numEntries,
								   __global uint
								   * keys,
								   __global uint
								   * values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintIdentityPerfectCLHash_InnerInsertNoOverwrite(tableData,
							     numEntries,
							     keys +
							     (i * numEntries),
							     values +
							     (i * numEntries));
}

typedef struct uintuintIdentitySentinelPerfectCLHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
	uint emptyValue;
} uintuintIdentitySentinelPerfectCLHash_TableData;
typedef struct uintuintIdentitySentinelPerfectCLHash_Bucket {
	uint value;
} uintuintIdentitySentinelPerfectCLHash_Bucket;
int uintuintIdentitySentinelPerfectCLHash_InnerQuerySingle(__global char
							   *tableData, uint key,
							   __local uint *
							   valueOutput) {
	__global uintuintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectCLHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((__global
					    uintuintIdentitySentinelPerfectCLHash_TableData
					    *) tableData)->compressFuncData,
					  key);
	if (buckets[index].value !=
	    ((__global uintuintIdentitySentinelPerfectCLHash_TableData *)
	     tableData)->emptyValue) {
		exitCode = HASH_SEARCH_CODE_MATCH;
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int uintuintIdentitySentinelPerfectCLHash_InnerQuery(__global char *tableData,
						     unsigned int numKeys,
						     __global uint * keys,
						     __local uint * valuesOutput) {
	__global uintuintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectCLHash_TableData)];
	uint key;
	__local uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    uintuintHash_CompressIdentity(((__global
						    uintuintIdentitySentinelPerfectCLHash_TableData
						    *) tableData)->
						  compressFuncData, key);
		if (buckets[index].value !=
		    ((__global uintuintIdentitySentinelPerfectCLHash_TableData
		      *) tableData)->emptyValue) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int uintuintIdentitySentinelPerfectCLHash_InnerInsertSingle(__global char
							    *tableData,
							    uint key,
							    uint value) {
	__global uintuintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectCLHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((__global
					    uintuintIdentitySentinelPerfectCLHash_TableData
					    *) tableData)->compressFuncData,
					  key);
	if (buckets[index].value !=
	    ((__global uintuintIdentitySentinelPerfectCLHash_TableData *)
	     tableData)->emptyValue) {
		exitCode = HASH_SEARCH_CODE_MATCH;
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int uintuintIdentitySentinelPerfectCLHash_InnerInsert(__global char *tableData,
						      unsigned int numEntries,
						      __global uint * keys,
						      __global uint * values) {
	__global uintuintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((__global
						    uintuintIdentitySentinelPerfectCLHash_TableData
						    *) tableData)->
						  compressFuncData, key);
		if (buckets[index].value !=
		    ((__global uintuintIdentitySentinelPerfectCLHash_TableData
		      *) tableData)->emptyValue) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int uintuintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite(__global
								       char
								       *tableData,
								       uint key,
								       uint
								       value) {
	__global uintuintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectCLHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((__global
					    uintuintIdentitySentinelPerfectCLHash_TableData
					    *) tableData)->compressFuncData,
					  key);
	if (buckets[index].value !=
	    ((__global uintuintIdentitySentinelPerfectCLHash_TableData *)
	     tableData)->emptyValue) {
		exitCode = HASH_SEARCH_CODE_MATCH;
	} else {
		exitCode = HASH_SEARCH_CODE_EMPTY;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int uintuintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite(__global char
								 *tableData,
								 unsigned int
								 numEntries,
								 __global uint *
								 keys,
								 __global uint *
								 values) {
	__global uintuintIdentitySentinelPerfectCLHash_Bucket *buckets =
	    (__global uintuintIdentitySentinelPerfectCLHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((__global
						    uintuintIdentitySentinelPerfectCLHash_TableData
						    *) tableData)->
						  compressFuncData, key);
		if (buckets[index].value !=
		    ((__global uintuintIdentitySentinelPerfectCLHash_TableData
		      *) tableData)->emptyValue) {
			exitCode = HASH_SEARCH_CODE_MATCH;
		} else {
			exitCode = HASH_SEARCH_CODE_EMPTY;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int uintuintIdentitySentinelPerfectCLHash_QuerySingle(__global char *tableData,
						      uint key,
						      __local uint *
						      valueOutput) {
	return uintuintIdentitySentinelPerfectCLHash_InnerQuerySingle(tableData,
								      key,
								      valueOutput);
}
int uintuintIdentitySentinelPerfectCLHash_Query(__global char *tableData,
						size_t numKeys,
						__global uint * keys,
						__local uint * valuesOutput) {
	return uintuintIdentitySentinelPerfectCLHash_InnerQuery(tableData,
								numKeys, keys,
								valuesOutput);
}
int uintuintIdentitySentinelPerfectCLHash_InsertSingle(__global char *tableData,
						       uint key, uint value) {
	return
	    uintuintIdentitySentinelPerfectCLHash_InnerInsertSingle(tableData,
								    key, value);
}
int uintuintIdentitySentinelPerfectCLHash_Insert(__global char *tableData,
						 size_t numEntries,
						 __global uint * keys,
						 __global uint * values) {
	return uintuintIdentitySentinelPerfectCLHash_InnerInsert(tableData,
								 numEntries,
								 keys, values);
}
int uintuintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite(__global char
								  *tableData,
								  uint key,
								  uint value) {
	return
	    uintuintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite
	    (tableData, key, value);
}
int uintuintIdentitySentinelPerfectCLHash_InsertNoOverwrite(__global char
							    *tableData,
							    size_t numEntries,
							    __global uint *
							    keys,
							    __global uint *
							    values) {
	return
	    uintuintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite
	    (tableData, numEntries, keys, values);
}
__kernel void uintuintIdentitySentinelPerfectCLHash_RangeQuerySingle(__global
								     char
								     *tableData,
								     unsigned
								     int
								     numQueries,
								     __global
								     uint *
								     keys,
								     __local
								     uint *
								     valuesOutput) 
{
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	uintuintIdentitySentinelPerfectCLHash_InnerQuerySingle(tableData,
							       keys[i],
							       valuesOutput +
							       i);
}
__kernel void uintuintIdentitySentinelPerfectCLHash_RangeQuery(__global char
							       *tableData,
							       unsigned int
							       numQueries,
							       unsigned int
							       numKeys,
							       __global uint *
							       keys,
							       __local uint *
							       valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	uintuintIdentitySentinelPerfectCLHash_InnerQuery(tableData, numKeys,
							 keys + (i * numKeys),
							 valuesOutput +
							 (i * numKeys));
}
__kernel void uintuintIdentitySentinelPerfectCLHash_RangeInsertSingle(__global
								      char
								      *tableData,
								      unsigned
								      int
								      numInsertions,
								      __global
								      uint *
								      keys,
								      __global
								      uint *
								      values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintIdentitySentinelPerfectCLHash_InnerInsertSingle(tableData,
								keys[i],
								values[i]);
}
__kernel void uintuintIdentitySentinelPerfectCLHash_RangeInsert(__global char
								*tableData,
								unsigned int
								numInsertions,
								unsigned int
								numEntries,
								__global uint *
								keys,
								__global uint *
								values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintIdentitySentinelPerfectCLHash_InnerInsert(tableData, numEntries,
							  keys +
							  (i * numEntries),
							  values +
							  (i * numEntries));
}
__kernel void
uintuintIdentitySentinelPerfectCLHash_RangeInsertSingleNoOverwrite(__global char
								   *tableData,
								   unsigned int
								   numInsertions,
								   __global uint
								   * keys,
								   __global uint
								   * values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintIdentitySentinelPerfectCLHash_InnerInsertSingleNoOverwrite
	    (tableData, keys[i], values[i]);
}
__kernel void
uintuintIdentitySentinelPerfectCLHash_RangeInsertNoOverwrite(__global char
							     *tableData,
							     unsigned int
							     numInsertions,
							     unsigned int
							     numEntries,
							     __global uint *
							     keys,
							     __global uint *
							     values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintIdentitySentinelPerfectCLHash_InnerInsertNoOverwrite(tableData,
								     numEntries,
								     keys +
								     (i *
								      numEntries),
								     values +
								     (i *
								      numEntries));
}

typedef struct uintuintLCGLinearOpenCompactCLHash_TableData {
	int hashID;
	unsigned int numBuckets;
	uintuintHash_CompressLCGData compressFuncData;
} uintuintLCGLinearOpenCompactCLHash_TableData;
typedef struct uintuintLCGLinearOpenCompactCLHash_Bucket {
	uint key;
	uint value;
} uintuintLCGLinearOpenCompactCLHash_Bucket;
int uintuintLCGLinearOpenCompactCLHash_InnerQuerySingle(__global char
							*tableData, uint key,
							__local uint *
							valueOutput) {
	__global uintuintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactCLHash_TableData)];
	uint index;
	int exitCode;
	__global uintuintLCGLinearOpenCompactCLHash_TableData *mytableData =
	    (__global uintuintLCGLinearOpenCompactCLHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((__global uintuintLCGLinearOpenCompactCLHash_TableData *)
		      tableData)->numBuckets);
		if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else if ((index == c && iteration > 0)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int uintuintLCGLinearOpenCompactCLHash_InnerQuery(__global char *tableData,
						  unsigned int numKeys,
						  __global uint * keys,
						  __local uint * valuesOutput) {
	__global uintuintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactCLHash_TableData)];
	uint key;
	__local uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		__global uintuintLCGLinearOpenCompactCLHash_TableData
		    *mytableData =
		    (__global uintuintLCGLinearOpenCompactCLHash_TableData *)
		    tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((__global
			       uintuintLCGLinearOpenCompactCLHash_TableData *)
			      tableData)->numBuckets);
			if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else if ((index == c && iteration > 0)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int uintuintLCGLinearOpenCompactCLHash_InnerInsertSingle(__global char
							 *tableData, uint key,
							 uint value) {
	__global uintuintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactCLHash_TableData)];
	uint index;
	int exitCode;
	__global uintuintLCGLinearOpenCompactCLHash_TableData *mytableData =
	    (__global uintuintLCGLinearOpenCompactCLHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((__global uintuintLCGLinearOpenCompactCLHash_TableData *)
		      tableData)->numBuckets);
		if ((atomic_cmpxchg
		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
		      key)) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else if ((index == c && iteration > 0)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int uintuintLCGLinearOpenCompactCLHash_InnerInsert(__global char *tableData,
						   unsigned int numEntries,
						   __global uint * keys,
						   __global uint * values) {
	__global uintuintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		__global uintuintLCGLinearOpenCompactCLHash_TableData
		    *mytableData =
		    (__global uintuintLCGLinearOpenCompactCLHash_TableData *)
		    tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((__global
			       uintuintLCGLinearOpenCompactCLHash_TableData *)
			      tableData)->numBuckets);
			if ((atomic_cmpxchg
			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
			      key)) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else if ((index == c && iteration > 0)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int uintuintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global
								    char
								    *tableData,
								    uint key,
								    uint value) 
{
	__global uintuintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactCLHash_TableData)];
	uint index;
	int exitCode;
	__global uintuintLCGLinearOpenCompactCLHash_TableData *mytableData =
	    (__global uintuintLCGLinearOpenCompactCLHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((__global uintuintLCGLinearOpenCompactCLHash_TableData *)
		      tableData)->numBuckets);
		if ((atomic_cmpxchg
		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
		      key)) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else if ((index == c && iteration > 0)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int uintuintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(__global char
							      *tableData,
							      unsigned int
							      numEntries,
							      __global uint *
							      keys,
							      __global uint *
							      values) {
	__global uintuintLCGLinearOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGLinearOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		__global uintuintLCGLinearOpenCompactCLHash_TableData
		    *mytableData =
		    (__global uintuintLCGLinearOpenCompactCLHash_TableData *)
		    tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((__global
			       uintuintLCGLinearOpenCompactCLHash_TableData *)
			      tableData)->numBuckets);
			if ((atomic_cmpxchg
			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
			      key)) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else if ((index == c && iteration > 0)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int uintuintLCGLinearOpenCompactCLHash_QuerySingle(__global char *tableData,
						   uint key,
						   __local uint *
						   valueOutput) {
	return uintuintLCGLinearOpenCompactCLHash_InnerQuerySingle(tableData,
								   key,
								   valueOutput);
}
int uintuintLCGLinearOpenCompactCLHash_Query(__global char *tableData,
					     size_t numKeys,
					     __global uint * keys,
					     __local uint * valuesOutput) {
	return uintuintLCGLinearOpenCompactCLHash_InnerQuery(tableData, numKeys,
							     keys,
							     valuesOutput);
}
int uintuintLCGLinearOpenCompactCLHash_InsertSingle(__global char *tableData,
						    uint key, uint value) {
	return uintuintLCGLinearOpenCompactCLHash_InnerInsertSingle(tableData,
								    key, value);
}
int uintuintLCGLinearOpenCompactCLHash_Insert(__global char *tableData,
					      size_t numEntries,
					      __global uint * keys,
					      __global uint * values) {
	return uintuintLCGLinearOpenCompactCLHash_InnerInsert(tableData,
							      numEntries, keys,
							      values);
}
int uintuintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite(__global char
							       *tableData,
							       uint key,
							       uint value) {
	return
	    uintuintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite
	    (tableData, key, value);
}
int uintuintLCGLinearOpenCompactCLHash_InsertNoOverwrite(__global char
							 *tableData,
							 size_t numEntries,
							 __global uint * keys,
							 __global uint *
							 values) {
	return
	    uintuintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(tableData,
								      numEntries,
								      keys,
								      values);
}
__kernel void uintuintLCGLinearOpenCompactCLHash_RangeQuerySingle(__global char
								  *tableData,
								  unsigned int
								  numQueries,
								  __global uint
								  * keys,
								  __local uint
								  *
								  valuesOutput) 
{
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	uintuintLCGLinearOpenCompactCLHash_InnerQuerySingle(tableData, keys[i],
							    valuesOutput + i);
}
__kernel void uintuintLCGLinearOpenCompactCLHash_RangeQuery(__global char
							    *tableData,
							    unsigned int
							    numQueries,
							    unsigned int
							    numKeys,
							    __global uint *
							    keys,
							    __local uint *
							    valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	uintuintLCGLinearOpenCompactCLHash_InnerQuery(tableData, numKeys,
						      keys + (i * numKeys),
						      valuesOutput +
						      (i * numKeys));
}
__kernel void uintuintLCGLinearOpenCompactCLHash_RangeInsertSingle(__global char
								   *tableData,
								   unsigned int
								   numInsertions,
								   __global uint
								   * keys,
								   __global uint
								   * values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintLCGLinearOpenCompactCLHash_InnerInsertSingle(tableData, keys[i],
							     values[i]);
}
__kernel void uintuintLCGLinearOpenCompactCLHash_RangeInsert(__global char
							     *tableData,
							     unsigned int
							     numInsertions,
							     unsigned int
							     numEntries,
							     __global uint *
							     keys,
							     __global uint *
							     values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintLCGLinearOpenCompactCLHash_InnerInsert(tableData, numEntries,
						       keys + (i * numEntries),
						       values +
						       (i * numEntries));
}
__kernel void
uintuintLCGLinearOpenCompactCLHash_RangeInsertSingleNoOverwrite(__global char
								*tableData,
								unsigned int
								numInsertions,
								__global uint *
								keys,
								__global uint *
								values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintLCGLinearOpenCompactCLHash_InnerInsertSingleNoOverwrite
	    (tableData, keys[i], values[i]);
}
__kernel void uintuintLCGLinearOpenCompactCLHash_RangeInsertNoOverwrite(__global
									char
									*tableData,
									unsigned
									int
									numInsertions,
									unsigned
									int
									numEntries,
									__global
									uint *
									keys,
									__global
									uint *
									values) 
{
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintLCGLinearOpenCompactCLHash_InnerInsertNoOverwrite(tableData,
								  numEntries,
								  keys +
								  (i *
								   numEntries),
								  values +
								  (i *
								   numEntries));
}

typedef struct uintuintLCGQuadraticOpenCompactCLHash_TableData {
	int hashID;
	unsigned int numBuckets;
	uintuintHash_CompressLCGData compressFuncData;
} uintuintLCGQuadraticOpenCompactCLHash_TableData;
typedef struct uintuintLCGQuadraticOpenCompactCLHash_Bucket {
	uint key;
	uint value;
} uintuintLCGQuadraticOpenCompactCLHash_Bucket;
int uintuintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(__global char
							   *tableData, uint key,
							   __local uint *
							   valueOutput) {
	__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactCLHash_TableData)];
	uint index;
	int exitCode;
	__global uintuintLCGQuadraticOpenCompactCLHash_TableData *mytableData =
	    (__global uintuintLCGQuadraticOpenCompactCLHash_TableData *)
	    tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((__global uintuintLCGQuadraticOpenCompactCLHash_TableData
		       *) tableData)->numBuckets);
		if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((__global
			   uintuintLCGQuadraticOpenCompactCLHash_TableData *)
			  tableData)->numBuckets)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
		*valueOutput = buckets[index].value;
		return HASH_EXIT_CODE_NORMAL;
	case HASH_SEARCH_CODE_MISMATCH:
	case HASH_SEARCH_CODE_EMPTY:
		return HASH_EXIT_CODE_KEY_DNE;
	default:
		return exitCode;
	}
}
int uintuintLCGQuadraticOpenCompactCLHash_InnerQuery(__global char *tableData,
						     unsigned int numKeys,
						     __global uint * keys,
						     __local uint * valuesOutput) {
	__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactCLHash_TableData)];
	uint key;
	__local uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		__global uintuintLCGQuadraticOpenCompactCLHash_TableData
		    *mytableData =
		    (__global uintuintLCGQuadraticOpenCompactCLHash_TableData *)
		    tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((__global
			       uintuintLCGQuadraticOpenCompactCLHash_TableData
			       *) tableData)->numBuckets);
			if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((__global
				   uintuintLCGQuadraticOpenCompactCLHash_TableData
				   *) tableData)->numBuckets)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
			*valueOutput = buckets[index].value;
			break;
		case HASH_SEARCH_CODE_MISMATCH:
		case HASH_SEARCH_CODE_EMPTY:
			resultExitCode = HASH_EXIT_CODE_KEY_DNE;
			break;
		default:
			return exitCode;
		}
	}
	return resultExitCode;
}
int uintuintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(__global char
							    *tableData,
							    uint key,
							    uint value) {
	__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactCLHash_TableData)];
	uint index;
	int exitCode;
	__global uintuintLCGQuadraticOpenCompactCLHash_TableData *mytableData =
	    (__global uintuintLCGQuadraticOpenCompactCLHash_TableData *)
	    tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((__global uintuintLCGQuadraticOpenCompactCLHash_TableData
		       *) tableData)->numBuckets);
		if ((atomic_cmpxchg
		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
		      key)) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((__global
			   uintuintLCGQuadraticOpenCompactCLHash_TableData *)
			  tableData)->numBuckets)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		buckets[index].value = value;
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int uintuintLCGQuadraticOpenCompactCLHash_InnerInsert(__global char *tableData,
						      unsigned int numEntries,
						      __global uint * keys,
						      __global uint * values) {
	__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		__global uintuintLCGQuadraticOpenCompactCLHash_TableData
		    *mytableData =
		    (__global uintuintLCGQuadraticOpenCompactCLHash_TableData *)
		    tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((__global
			       uintuintLCGQuadraticOpenCompactCLHash_TableData
			       *) tableData)->numBuckets);
			if ((atomic_cmpxchg
			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
			      key)) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((__global
				   uintuintLCGQuadraticOpenCompactCLHash_TableData
				   *) tableData)->numBuckets)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int uintuintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite(__global
								       char
								       *tableData,
								       uint key,
								       uint
								       value) {
	__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactCLHash_TableData)];
	uint index;
	int exitCode;
	__global uintuintLCGQuadraticOpenCompactCLHash_TableData *mytableData =
	    (__global uintuintLCGQuadraticOpenCompactCLHash_TableData *)
	    tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((__global uintuintLCGQuadraticOpenCompactCLHash_TableData
		       *) tableData)->numBuckets);
		if ((atomic_cmpxchg
		     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
		      key)) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((__global
			   uintuintLCGQuadraticOpenCompactCLHash_TableData *)
			  tableData)->numBuckets)) {
			exitCode = HASH_EXIT_CODE_CYCLE;
			break;
		}
		iteration++;
	}
	switch (exitCode) {
	case HASH_SEARCH_CODE_MATCH:
	case HASH_SEARCH_CODE_MISMATCH:
		return HASH_EXIT_CODE_OVERWRITE;
	case HASH_SEARCH_CODE_EMPTY:
		buckets[index].value = value;
		return HASH_EXIT_CODE_NORMAL;
	default:
		return exitCode;
	}
}
int uintuintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite(__global char
								 *tableData,
								 unsigned int
								 numEntries,
								 __global uint *
								 keys,
								 __global uint *
								 values) {
	__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *buckets =
	    (__global uintuintLCGQuadraticOpenCompactCLHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactCLHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		__global uintuintLCGQuadraticOpenCompactCLHash_TableData
		    *mytableData =
		    (__global uintuintLCGQuadraticOpenCompactCLHash_TableData *)
		    tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((__global
			       uintuintLCGQuadraticOpenCompactCLHash_TableData
			       *) tableData)->numBuckets);
			if ((atomic_cmpxchg
			     (&(buckets[index].key), HASH_BUCKET_STATUS_EMPTY,
			      key)) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((__global
				   uintuintLCGQuadraticOpenCompactCLHash_TableData
				   *) tableData)->numBuckets)) {
				exitCode = HASH_EXIT_CODE_CYCLE;
				break;
			}
			iteration++;
		}
		switch (exitCode) {
		case HASH_SEARCH_CODE_MATCH:
		case HASH_SEARCH_CODE_MISMATCH:
			resultExitCode = HASH_EXIT_CODE_OVERWRITE;
			break;
		case HASH_SEARCH_CODE_EMPTY:
			buckets[index].value = values[i];
			break;
		default:
			resultExitCode = exitCode;
		}
	}
	return resultExitCode;
}
int uintuintLCGQuadraticOpenCompactCLHash_QuerySingle(__global char *tableData,
						      uint key,
						      __local uint *
						      valueOutput) {
	return uintuintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(tableData,
								      key,
								      valueOutput);
}
int uintuintLCGQuadraticOpenCompactCLHash_Query(__global char *tableData,
						size_t numKeys,
						__global uint * keys,
						__local uint * valuesOutput) {
	return uintuintLCGQuadraticOpenCompactCLHash_InnerQuery(tableData,
								numKeys, keys,
								valuesOutput);
}
int uintuintLCGQuadraticOpenCompactCLHash_InsertSingle(__global char *tableData,
						       uint key, uint value) {
	return
	    uintuintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(tableData,
								    key, value);
}
int uintuintLCGQuadraticOpenCompactCLHash_Insert(__global char *tableData,
						 size_t numEntries,
						 __global uint * keys,
						 __global uint * values) {
	return uintuintLCGQuadraticOpenCompactCLHash_InnerInsert(tableData,
								 numEntries,
								 keys, values);
}
int uintuintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite(__global char
								  *tableData,
								  uint key,
								  uint value) {
	return
	    uintuintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite
	    (tableData, key, value);
}
int uintuintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(__global char
							    *tableData,
							    size_t numEntries,
							    __global uint *
							    keys,
							    __global uint *
							    values) {
	return
	    uintuintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite
	    (tableData, numEntries, keys, values);
}
__kernel void uintuintLCGQuadraticOpenCompactCLHash_RangeQuerySingle(__global
								     char
								     *tableData,
								     unsigned
								     int
								     numQueries,
								     __global
								     uint *
								     keys,
								     __local
								     uint *
								     valuesOutput) 
{
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	uintuintLCGQuadraticOpenCompactCLHash_InnerQuerySingle(tableData,
							       keys[i],
							       valuesOutput +
							       i);
}
__kernel void uintuintLCGQuadraticOpenCompactCLHash_RangeQuery(__global char
							       *tableData,
							       unsigned int
							       numQueries,
							       unsigned int
							       numKeys,
							       __global uint *
							       keys,
							       __local uint *
							       valuesOutput) {
	uint i = get_global_id(0);
	if (i >= numQueries) {
		return;
	}
	uintuintLCGQuadraticOpenCompactCLHash_InnerQuery(tableData, numKeys,
							 keys + (i * numKeys),
							 valuesOutput +
							 (i * numKeys));
}
__kernel void uintuintLCGQuadraticOpenCompactCLHash_RangeInsertSingle(__global
								      char
								      *tableData,
								      unsigned
								      int
								      numInsertions,
								      __global
								      uint *
								      keys,
								      __global
								      uint *
								      values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintLCGQuadraticOpenCompactCLHash_InnerInsertSingle(tableData,
								keys[i],
								values[i]);
}
__kernel void uintuintLCGQuadraticOpenCompactCLHash_RangeInsert(__global char
								*tableData,
								unsigned int
								numInsertions,
								unsigned int
								numEntries,
								__global uint *
								keys,
								__global uint *
								values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintLCGQuadraticOpenCompactCLHash_InnerInsert(tableData, numEntries,
							  keys +
							  (i * numEntries),
							  values +
							  (i * numEntries));
}
__kernel void
uintuintLCGQuadraticOpenCompactCLHash_RangeInsertSingleNoOverwrite(__global char
								   *tableData,
								   unsigned int
								   numInsertions,
								   __global uint
								   * keys,
								   __global uint
								   * values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintLCGQuadraticOpenCompactCLHash_InnerInsertSingleNoOverwrite
	    (tableData, keys[i], values[i]);
}
__kernel void
uintuintLCGQuadraticOpenCompactCLHash_RangeInsertNoOverwrite(__global char
							     *tableData,
							     unsigned int
							     numInsertions,
							     unsigned int
							     numEntries,
							     __global uint *
							     keys,
							     __global uint *
							     values) {
	uint i = get_global_id(0);
	if (i >= numInsertions) {
		return;
	}
	uintuintLCGQuadraticOpenCompactCLHash_InnerInsertNoOverwrite(tableData,
								     numEntries,
								     keys +
								     (i *
								      numEntries),
								     values +
								     (i *
								      numEntries));
}
__kernel void uintuintHash_RangeQuery(__global char *tableData,
				      unsigned int numQueries,
				      unsigned int numKeys,
				      __global uint * keys,
				      __local uint * valuesOutput) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return uintuintIdentityPerfectCLHash_RangeQuery(tableData,
								numQueries,
								numKeys, keys,
								valuesOutput);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentitySentinelPerfectCLHash_RangeQuery(tableData,
								     numQueries,
								     numKeys,
								     keys,
								     valuesOutput);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return uintuintLCGLinearOpenCompactCLHash_RangeQuery(tableData,
								     numQueries,
								     numKeys,
								     keys,
								     valuesOutput);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGQuadraticOpenCompactCLHash_RangeQuery(tableData,
								     numQueries,
								     numKeys,
								     keys,
								     valuesOutput);
	}
}
__kernel void uintuintHash_RangeQuerySingle(__global char *tableData,
					    unsigned int numQueries,
					    __global uint * keys,
					    __local uint * valueOutput) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return uintuintIdentityPerfectCLHash_RangeQuerySingle(tableData,
								      numQueries,
								      keys,
								      valueOutput);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentitySentinelPerfectCLHash_RangeQuerySingle
		    (tableData, numQueries, keys, valueOutput);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGLinearOpenCompactCLHash_RangeQuerySingle
		    (tableData, numQueries, keys, valueOutput);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGQuadraticOpenCompactCLHash_RangeQuerySingle
		    (tableData, numQueries, keys, valueOutput);
	}
}
__kernel void uintuintHash_RangeInsert(__global char *tableData,
				       unsigned int numInsertions,
				       unsigned int numEntries,
				       __global uint * keys,
				       __global uint * values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return uintuintIdentityPerfectCLHash_RangeInsert(tableData,
								 numInsertions,
								 numEntries,
								 keys, values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentitySentinelPerfectCLHash_RangeInsert(tableData,
								      numInsertions,
								      numEntries,
								      keys,
								      values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return uintuintLCGLinearOpenCompactCLHash_RangeInsert(tableData,
								      numInsertions,
								      numEntries,
								      keys,
								      values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGQuadraticOpenCompactCLHash_RangeInsert(tableData,
								      numInsertions,
								      numEntries,
								      keys,
								      values);
	}
}
__kernel void uintuintHash_RangeInsertSingle(__global char *tableData,
					     unsigned int numInsertions,
					     __global uint * keys,
					     __global uint * values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentityPerfectCLHash_RangeInsertSingle(tableData,
								    numInsertions,
								    keys,
								    values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentitySentinelPerfectCLHash_RangeInsertSingle
		    (tableData, numInsertions, keys, values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGLinearOpenCompactCLHash_RangeInsertSingle
		    (tableData, numInsertions, keys, values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGQuadraticOpenCompactCLHash_RangeInsertSingle
		    (tableData, numInsertions, keys, values);
	}
}
__kernel void uintuintHash_RangeInsertNoOverwrite(__global char *tableData,
						  unsigned int numInsertions,
						  unsigned int numEntries,
						  __global uint * keys,
						  __global uint * values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentityPerfectCLHash_RangeInsertNoOverwrite
		    (tableData, numInsertions, numEntries, keys, values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentitySentinelPerfectCLHash_RangeInsertNoOverwrite
		    (tableData, numInsertions, numEntries, keys, values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGLinearOpenCompactCLHash_RangeInsertNoOverwrite
		    (tableData, numInsertions, numEntries, keys, values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGQuadraticOpenCompactCLHash_RangeInsertNoOverwrite
		    (tableData, numInsertions, numEntries, keys, values);
	}
}
__kernel void uintuintHash_RangeInsertSingleNoOverwrite(__global char
							*tableData,
							unsigned int
							numInsertions,
							__global uint * keys,
							__global uint *
							values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentityPerfectCLHash_RangeInsertSingleNoOverwrite
		    (tableData, numInsertions, keys, values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentitySentinelPerfectCLHash_RangeInsertSingleNoOverwrite
		    (tableData, numInsertions, keys, values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGLinearOpenCompactCLHash_RangeInsertSingleNoOverwrite
		    (tableData, numInsertions, keys, values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGQuadraticOpenCompactCLHash_RangeInsertSingleNoOverwrite
		    (tableData, numInsertions, keys, values);
	}
}
int uintuintHash_Query(__global char *tableData, unsigned int numKeys,
		       __global uint * keys, __local uint * valuesOutput) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return uintuintIdentityPerfectCLHash_Query(tableData, numKeys,
							   keys, valuesOutput);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return uintuintIdentitySentinelPerfectCLHash_Query(tableData,
								   numKeys,
								   keys,
								   valuesOutput);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return uintuintLCGLinearOpenCompactCLHash_Query(tableData,
								numKeys, keys,
								valuesOutput);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return uintuintLCGQuadraticOpenCompactCLHash_Query(tableData,
								   numKeys,
								   keys,
								   valuesOutput);
	}
	return HASH_EXIT_CODE_ERROR;
}
int uintuintHash_QuerySingle(__global char *tableData, uint key,
			     __local uint * valueOutput) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return uintuintIdentityPerfectCLHash_QuerySingle(tableData, key,
								 valueOutput);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentitySentinelPerfectCLHash_QuerySingle(tableData,
								      key,
								      valueOutput);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return uintuintLCGLinearOpenCompactCLHash_QuerySingle(tableData,
								      key,
								      valueOutput);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGQuadraticOpenCompactCLHash_QuerySingle(tableData,
								      key,
								      valueOutput);
	}
	return HASH_EXIT_CODE_ERROR;
}
int uintuintHash_Insert(__global char *tableData, unsigned int numEntries,
			__global uint * keys, __global uint * values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return uintuintIdentityPerfectCLHash_Insert(tableData,
							    numEntries, keys,
							    values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return uintuintIdentitySentinelPerfectCLHash_Insert(tableData,
								    numEntries,
								    keys,
								    values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return uintuintLCGLinearOpenCompactCLHash_Insert(tableData,
								 numEntries,
								 keys, values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return uintuintLCGQuadraticOpenCompactCLHash_Insert(tableData,
								    numEntries,
								    keys,
								    values);
	}
	return HASH_EXIT_CODE_ERROR;
}
int uintuintHash_InsertSingle(__global char *tableData, uint key, uint value) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return uintuintIdentityPerfectCLHash_InsertSingle(tableData,
								  key, value);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentitySentinelPerfectCLHash_InsertSingle
		    (tableData, key, value);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGLinearOpenCompactCLHash_InsertSingle(tableData,
								    key, value);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGQuadraticOpenCompactCLHash_InsertSingle
		    (tableData, key, value);
	}
	return HASH_EXIT_CODE_ERROR;
}
int uintuintHash_InsertNoOverwrite(__global char *tableData,
				   unsigned int numEntries,
				   __global uint * keys,
				   __global uint * values) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentityPerfectCLHash_InsertNoOverwrite(tableData,
								    numEntries,
								    keys,
								    values);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentitySentinelPerfectCLHash_InsertNoOverwrite
		    (tableData, numEntries, keys, values);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGLinearOpenCompactCLHash_InsertNoOverwrite
		    (tableData, numEntries, keys, values);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite
		    (tableData, numEntries, keys, values);
	}
	return HASH_EXIT_CODE_ERROR;
}
int uintuintHash_InsertSingleNoOverwrite(__global char *tableData, uint key,
					 uint value) {
	switch (((__global int *)tableData)[0]) {
	case IDENTITY_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentityPerfectCLHash_InsertSingleNoOverwrite
		    (tableData, key, value);
	case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
		return
		    uintuintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite
		    (tableData, key, value);
	case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite
		    (tableData, key, value);
	case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
		return
		    uintuintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite
		    (tableData, key, value);
	}
	return HASH_EXIT_CODE_ERROR;
}
