
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
 * @file   HashFactory.c
 * @author Peter Ahrens
 * @date   Thu Jun 6 2013 
 */
//
#ifdef _OPENMP
#include <omp.h>
#endif
#include "HashFactory.h"
#ifdef HAVE_OPENCL
#ifdef __APPLE_CC__
#include <OpenCL/OpenCL.h>
#else
#include <CL/cl.h>
#endif
#else
typedef int cl_kernel;
#define CL_CONTEXT_DEVICES 0
#define CL_MEM_READ_WRITE 0
#define CL_TRUE 1
#define CL_MEM_READ_ONLY 0
#define CL_MEM_COPY_HOST_PTR 0
#define CL_MEM_WRITE_ONLY 0
#define CL_KERNEL_WORK_GROUP_SIZE 128

int clRetainContext(int context) {
	return context;
}
int clRetainCommandQueue(int command_queue) {
	return command_queue;
}
int clGetContextInfo(int context, int param, size_t size, void *value,
		     size_t * size_ret) {
	return 0;
}
int clReleaseContext(int context) {
	return context;
}
int clReleaseCommandQueue(int command_queue) {
	return command_queue;
}
int clReleaseProgram(int program) {
	return program;
}
int clRetainKernel(int kernel) {
	return kernel;
}
int clRetainProgram(int program) {
	return program;
}
cl_mem clCreateBuffer(int context, int flags, size_t size, void *value,
		      int *size_ret) {
	return 0;
}
int clEnqueueWriteBuffer(int command_queue, void *buffer, int blocking_write,
			 size_t offset, size_t cb, const void *ptr,
			 uint nevents, const int *wait_list, int *event) {
	return 0;
}
int clEnqueueReadBuffer(int command_queue, void *buffer, int blocking_write,
			size_t offset, size_t cb, const void *ptr, uint nevents,
			const int *wait_list, int *event) {
	return 0;
}
int clCreateKernel(int program, const char *kernel_name, int *errcode_ret) {
	return 0;
}
int clReleaseKernel(int kernel) {
	return kernel;
}
int clReleaseMemObject(void *memobj) {
	return 0;
}
int clSetKernelArg(int kernel, uint arg_index, size_t arg_size,
		   const void *arg_value) {
	return 0;
}
int clGetKernelWorkGroupInfo(int kernel, int device, int param_name,
			     size_t size, void *value, size_t * size_ret) {
	return 0;
}
int clEnqueueNDRangeKernel(int command_queue, int kernel, uint work_dim,
			   const size_t * offset, const size_t * size,
			   const size_t * local_size, uint nevents,
			   const int *wait_list, int *event) {
	return 0;
}
int clFinish(int command_queue) {
	return 0;
}
#endif
#define PRIME_NUM_CHECKS 20
#include <math.h>
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
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic push
static int reportLevel = 0;

#include "HashFactory_source.inc"
const char *Hash_GetKernelSourceString() {
	return HashFactory_source;
}
size_t roundUpToNearest(size_t x, size_t r) {
	return (((x - 1) / r) + 1) * r;
}
int modularPow(int base, int exponent, int modulus) {
	int result = 1;
	while (exponent) {
		if (exponent & 1)
			result = ((long long int)result * base) % modulus;
		exponent >>= 1;
		base = ((long long int)base * base) % modulus;
	}
	return result;
}
int largestProthPrimeUnder(int N) {
	if (N < 4) {
		return N;
	}
	//determine the nearest proth number
	int n;
	frexp((double)N, &n);
	n /= 2;
	int s = 1 << n;
	int p = s * ((N - 1) / s) + 1;
	int i;
	int a;
	srand(p);
	while (p > 3) {
		//check if a proth number is prime
		for (i = 0; i < PRIME_NUM_CHECKS; i++) {
			a = rand();
			if (modularPow(a, (p - 1) / 2, p) == p - 1) {
				return p;
			}
		}
		//determine the next proth number
		if (p - 1 == s * s / 4) {
			s /= 2;
		}
		p -= s;
	}
	return 3;
}
int smallestProthPrimeAbove(int N) {
	if (N < 4) {
		return N;
	}
	//determine the nearest proth number
	int n;
	frexp((double)N, &n);
	n /= 2;
	int s = 1 << n;
	int p = s * ((N - 1) / s) + 1;
	int i;
	int a;
	srand(p);
	while (1) {
		//determine the next proth number
		if (p - 1 == s * s) {
			s *= 2;
		}
		p += s;
		//check if a proth number is prime
		for (i = 0; i < PRIME_NUM_CHECKS; i++) {
			a = rand();
			if (modularPow(a, (p - 1) / 2, p) == p - 1) {
				return p;
			}
		}
	}
	return 3;
}
int intLog2(int n) {
	int result = 0;
	while (n >>= 1) {
		result++;
	}
	return result;
}
void Hash_SetReportLevel(int level) {
	reportLevel = level;
}
int Hash_GetReportLevel() {
	return reportLevel;
}
char *Hash_ExitCodeString(int exitCode) {
	switch (exitCode) {
	case HASH_EXIT_CODE_NORMAL:
		return "Normal";
	case HASH_EXIT_CODE_ERROR:
		return "Error";
	case HASH_EXIT_CODE_OVERWRITE:
		return "Overwrite";
	case HASH_EXIT_CODE_KEY_DNE:
		return "Key Does Not Exist";
	case HASH_EXIT_CODE_CYCLE:
		return "Cycle";
	case HASH_EXIT_CODE_MAX_ENTRIES_EXCEEDED:
		return "Maximum Number Of Entries Exceeded";
	default:
		return "Unknown";
	}
}
void Hash_ExitCodeDebug(int exitCode) {
	if (exitCode != HASH_EXIT_CODE_NORMAL) {
		printf("HashExitCode: %s\n", Hash_ExitCodeString(exitCode));
	}
}

struct uintuintHash_Table_ {
	char *tableData;
	cl_mem tableDataBuffer;
	int (*destroyFunc) (uintuintHash_Table *);
	int (*setupFunc) (uintuintHash_Table *);
	int (*emptyFunc) (uintuintHash_Table *);
	int (*queryFunc) (uintuintHash_Table *, size_t, uint *, uint *);
	int (*querySingleFunc) (uintuintHash_Table *, uint, uint *);
	int (*insertFunc) (uintuintHash_Table *, size_t, uint *, uint *);
	int (*insertSingleFunc) (uintuintHash_Table *, uint, uint);
	int (*insertNoOverwriteFunc) (uintuintHash_Table *, size_t, uint *,
				      uint *);
	int (*insertSingleNoOverwriteFunc) (uintuintHash_Table *, uint, uint);
	int (*bufferQueryFunc) (uintuintHash_Table *, size_t, cl_mem, cl_mem);
	int (*bufferInsertFunc) (uintuintHash_Table *, size_t, cl_mem, cl_mem);
	int (*bufferInsertNoOverwriteFunc) (uintuintHash_Table *, size_t,
					    cl_mem, cl_mem);
	cl_context context;
	cl_command_queue queue;
	cl_program utilProgram;
	cl_kernel emptyKernel;
	size_t emptyKernelLocalWorkSize;
	cl_program program;
	cl_kernel querySingleKernel;
	cl_kernel insertSingleKernel;
	cl_kernel insertSingleNoOverwriteKernel;
	size_t localWorkSize;
};
struct uintuintHash_Factory_ {
	cl_context context;
	cl_program program;
	cl_command_queue queue;
	int hashTypesAvailable;
	cl_program utilProgram[HASH_NUM_CL_HASHES];
	cl_kernel emptyKernel[HASH_NUM_CL_HASHES];
	size_t emptyKernelLocalWorkSize[HASH_NUM_CL_HASHES];
	cl_kernel querySingleKernel[HASH_NUM_CL_HASHES];
	cl_kernel insertSingleKernel[HASH_NUM_CL_HASHES];
	cl_kernel insertSingleNoOverwriteKernel[HASH_NUM_CL_HASHES];
	uint emptyValue;
	size_t localWorkSize;
	uintuintHash_Table *(*createFunc[HASH_NUM_HASHES]) (uintuintHash_Factory
							    *, int hashIndex,
							    size_t keyRange,
							    size_t numEntries,
							    float loadFactor);
	int (*destroyFunc[HASH_NUM_HASHES]) (uintuintHash_Factory *,
					     int hashIndex);
};
uintuintHash_Factory *uintuintHash_CreateFactory(int hashTypes,
						 uint * emptyValue,
						 size_t localWorkSize,
						 cl_context * context,
						 cl_command_queue * queue) {
	if (hashTypes == 0) {
		hashTypes = HASH_ALL_C_HASHES;
	}
	if (!(hashTypes & HASH_ALL_HASHES)) {
		printf("Please specify a valid hash type to create.\n");
		exit(1);
	}
	hashTypes &= HASH_ALL_HASHES;
	if ((hashTypes & HASH_SENTINEL_PERFECT_HASHES) == hashTypes
	    && emptyValue == NULL) {
		printf
		    ("emptyValue must be valid if a sentinel perfect hash is the only option available.\n");
		exit(1);
	}
	uintuintHash_Factory *factory =
	    (uintuintHash_Factory *) malloc(sizeof(uintuintHash_Factory));
	if (emptyValue == NULL) {
		hashTypes &= !HASH_SENTINEL_PERFECT_HASHES;
	} else {
		factory->emptyValue = *emptyValue;
	}
	factory->hashTypesAvailable = hashTypes;
	if (hashTypes & HASH_ALL_CL_HASHES) {
		if (localWorkSize == 0) {
			factory->localWorkSize = HASH_DEFAULT_LOCAL_WORK_SIZE;
		} else {
			factory->localWorkSize = 1 << intLog2(localWorkSize);
		}
		if (context == NULL) {
			CLHash_Utilities_CreateContext(&factory->context,
						       &factory->queue);
		} else {
			factory->context = *context;
			clRetainContext(*context);
			if (queue == NULL) {
				printf
				    ("Please specify a command queue for your context.\n");
				exit(-1);
			}
			factory->queue = *queue;
			clRetainCommandQueue(*queue);
		}
		cl_int error;
		cl_device_id device;
		error =
		    clGetContextInfo(factory->context, CL_CONTEXT_DEVICES,
				     sizeof(device), &device, NULL);
		CLHash_Utilities_HandleError(error,
					     "uintuintHash_CreateFactory",
					     "clGetContextInfo");
		factory->program =
		    CLHash_Utilities_BuildProgramString(factory->context,
							device,
							Hash_GetKernelSourceString
							());
	}
	int hashType = 1;
	for (int hashIndex = 0; hashIndex < HASH_NUM_HASHES; hashIndex++) {
		hashType = 1 << hashIndex;
		switch (hashType & hashTypes) {
		case IDENTITY_PERFECT_HASH_ID:
			uintuintIdentityPerfectHash_CreateFactory(factory,
								  hashIndex);
			break;
		case IDENTITY_PERFECT_CL_HASH_ID:
			uintuintIdentityPerfectCLHash_CreateFactory(factory,
								    hashIndex);
			break;
		case IDENTITY_PERFECT_OPENMP_HASH_ID:
			uintuintIdentityPerfectOpenMPHash_CreateFactory(factory,
									hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_HASH_ID:
			uintuintIdentitySentinelPerfectHash_CreateFactory
			    (factory, hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
			uintuintIdentitySentinelPerfectCLHash_CreateFactory
			    (factory, hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID:
			uintuintIdentitySentinelPerfectOpenMPHash_CreateFactory
			    (factory, hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_HASH_ID:
			uintuintLCGLinearOpenCompactHash_CreateFactory(factory,
								       hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
			uintuintLCGLinearOpenCompactCLHash_CreateFactory
			    (factory, hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID:
			uintuintLCGLinearOpenCompactOpenMPHash_CreateFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_HASH_ID:
			uintuintLCGQuadraticOpenCompactHash_CreateFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
			uintuintLCGQuadraticOpenCompactCLHash_CreateFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID:
			uintuintLCGQuadraticOpenCompactOpenMPHash_CreateFactory
			    (factory, hashIndex);
			break;
		}
	}
	return factory;
}
int uintuintHash_DestroyFactory(uintuintHash_Factory * factory) {
	int hashType = 1;
	for (int hashIndex = 0; hashIndex < HASH_NUM_HASHES; hashIndex++) {
		hashType = 1 << hashIndex;
		switch (hashType & factory->hashTypesAvailable) {
		case IDENTITY_PERFECT_HASH_ID:
			uintuintIdentityPerfectHash_DestroyFactory(factory,
								   hashIndex);
			break;
		case IDENTITY_PERFECT_CL_HASH_ID:
			uintuintIdentityPerfectCLHash_DestroyFactory(factory,
								     hashIndex);
			break;
		case IDENTITY_PERFECT_OPENMP_HASH_ID:
			uintuintIdentityPerfectOpenMPHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_HASH_ID:
			uintuintIdentitySentinelPerfectHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_CL_HASH_ID:
			uintuintIdentitySentinelPerfectCLHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID:
			uintuintIdentitySentinelPerfectOpenMPHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_HASH_ID:
			uintuintLCGLinearOpenCompactHash_DestroyFactory(factory,
									hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID:
			uintuintLCGLinearOpenCompactCLHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID:
			uintuintLCGLinearOpenCompactOpenMPHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_HASH_ID:
			uintuintLCGQuadraticOpenCompactHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID:
			uintuintLCGQuadraticOpenCompactCLHash_DestroyFactory
			    (factory, hashIndex);
			break;
		case LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID:
			uintuintLCGQuadraticOpenCompactOpenMPHash_DestroyFactory
			    (factory, hashIndex);
			break;
		}
		hashIndex++;
	}
	if (factory->hashTypesAvailable & HASH_ALL_CL_HASHES) {
		clReleaseContext(factory->context);
		clReleaseCommandQueue(factory->queue);
		clReleaseProgram(factory->program);
	}
	free(factory);
	return (0);
}
uintuintHash_Table *uintuintHash_CreateTable(uintuintHash_Factory * factory,
					     int hashTypes, size_t keyRange,
					     size_t numEntries,
					     float loadFactor) {
	if (loadFactor > 1.0 || loadFactor < HASH_MIN_LOAD_FACTOR) {
		loadFactor = HASH_DEFAULT_LOAD_FACTOR;
	}
	if (hashTypes == 0) {
		hashTypes = factory->hashTypesAvailable;
		if ((hashTypes & HASH_ALL_CL_HASHES)
		    && (hashTypes & HASH_ALL_OPENMP_HASHES)
		    && (hashTypes & HASH_ALL_C_HASHES)) {
			hashTypes &= HASH_ALL_C_HASHES;
		}
	}
	if (!(hashTypes & factory->hashTypesAvailable)) {
		printf
		    ("None of the selected hash types are supported by this factory.\n");
		exit(1);
	}
	hashTypes &= factory->hashTypesAvailable;
	if ((hashTypes & HASH_ALL_CL_HASHES)
	    && (hashTypes & HASH_ALL_OPENMP_HASHES)
	    && (hashTypes & HASH_ALL_C_HASHES)) {
		printf("Please decide between OpenCL, OpenMP or C hash.\n");
		exit(1);
	}
	if ((hashTypes & HASH_PERFECT_HASHES) == hashTypes && keyRange == 0) {
		printf
		    ("keyRange must be set if a perfect hash is the only option available.\n");
		exit(1);
	}
	if ((hashTypes & HASH_COMPACT_HASHES) == hashTypes && numEntries == 0) {
		printf
		    ("numEntries must be set if a compact hash is the only option available.\n");
		exit(1);
	}
	if (numEntries == 0 && keyRange == 0) {
		printf("either numEntries or keyRange must be set.\n");
		exit(1);
	}
	size_t perfectNumBuckets = keyRange;
	size_t compactNumBuckets = (size_t) (numEntries / loadFactor);
	int hashIndex;
	if ((hashTypes & HASH_SENTINEL_PERFECT_HASHES)
	    && ((hashTypes == (hashTypes & HASH_SENTINEL_PERFECT_HASHES))
		|| (compactNumBuckets == 0
		    || (perfectNumBuckets / compactNumBuckets <
			HASH_PERFECT_COMPACT_SWITCH_FACTOR)))) {
		hashIndex = intLog2(hashTypes & HASH_SENTINEL_PERFECT_HASHES);
	} else if ((hashTypes & HASH_NOSENTINEL_PERFECT_HASHES)
		   && ((hashTypes == (hashTypes & HASH_PERFECT_HASHES))
		       || (compactNumBuckets == 0
			   || (perfectNumBuckets / compactNumBuckets <
			       HASH_PERFECT_COMPACT_SWITCH_FACTOR)))) {
		hashIndex = intLog2(hashTypes & HASH_NOSENTINEL_PERFECT_HASHES);
	} else if ((hashTypes & HASH_LINEAR_COMPACT_HASHES)
		   &&
		   ((hashTypes ==
		     (hashTypes &
		      (HASH_PERFECT_HASHES | HASH_LINEAR_COMPACT_HASHES)))
		    ||
		    ((hashTypes & HASH_COMPACT_HASHES &
		      HASH_LINEAR_COMPACT_HASHES) ==
		     (hashTypes & HASH_COMPACT_HASHES) || loadFactor > 0.5))) {
		hashIndex = intLog2(hashTypes & HASH_LINEAR_COMPACT_HASHES);
	} else {
		hashIndex = intLog2(hashTypes & HASH_QUADRATIC_COMPACT_HASHES);
	}
	uintuintHash_Table *table =
	    factory->createFunc[hashIndex] (factory, hashIndex, keyRange,
					    numEntries, loadFactor);
	return table;
}
int uintuintHash_SetupTable(uintuintHash_Table * table) {
	table->setupFunc(table);
	return (0);
}
int uintuintHash_EmptyTable(uintuintHash_Table * table) {
	table->emptyFunc(table);
	return (0);
}
int uintuintHash_DestroyTable(uintuintHash_Table * table) {
	table->destroyFunc(table);
	return (0);
}
cl_mem uintuintHash_GetTableDataBuffer(uintuintHash_Table * table) {
	return table->tableDataBuffer;
}
cl_mem *uintuintHash_GetTableDataBufferPtr(uintuintHash_Table * table) {
	return &table->tableDataBuffer;
}
int uintuintHash_GetTableType(uintuintHash_Table * table) {
	return ((int *)table->tableData)[0];
} int uintuintHash_Query(uintuintHash_Table * table, size_t numKeys,
			 uint * keys, uint * valuesOutput) {
	table->queryFunc(table, numKeys, keys, valuesOutput);
	return (0);
}
int uintuintHash_QuerySingle(uintuintHash_Table * table, uint key,
			     uint * valueOutput) {
	table->querySingleFunc(table, key, valueOutput);
	return (0);
}
int uintuintHash_Insert(uintuintHash_Table * table, size_t numEntries,
			uint * keys, uint * values) {
	table->insertFunc(table, numEntries, keys, values);
	return (0);
}
int uintuintHash_InsertSingle(uintuintHash_Table * table, uint key, uint value) {
	table->insertSingleFunc(table, key, value);
	return (0);
}
int uintuintHash_InsertNoOverwrite(uintuintHash_Table * table,
				   size_t numEntries, uint * keys,
				   uint * values) {
	table->insertNoOverwriteFunc(table, numEntries, keys, values);
	return (0);
}
int uintuintHash_InsertSingleNoOverwrite(uintuintHash_Table * table, uint key,
					 uint value) {
	table->insertSingleNoOverwriteFunc(table, key, value);
	return (0);
}
int uintuintHash_BufferQuery(uintuintHash_Table * table, size_t numKeys,
			     cl_mem keys, cl_mem valuesOutput) {
	table->bufferQueryFunc(table, numKeys, keys, valuesOutput);
	return (0);
}
int uintuintHash_BufferInsert(uintuintHash_Table * table, size_t numEntries,
			      cl_mem keys, cl_mem values) {
	table->bufferInsertFunc(table, numEntries, keys, values);
	return (0);
}
int uintuintHash_BufferInsertNoOverwrite(uintuintHash_Table * table,
					 size_t numEntries, cl_mem keys,
					 cl_mem values) {
	table->bufferInsertNoOverwriteFunc(table, numEntries, keys, values);
	return (0);
}

typedef struct uintuintIdentityPerfectHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
} uintuintIdentityPerfectHash_TableData;
typedef struct uintuintIdentityPerfectHash_Bucket {
	uint key;
	uint value;
} uintuintIdentityPerfectHash_Bucket;
uintuintHash_Table *uintuintIdentityPerfectHash_CreateTable(uintuintHash_Factory
							    * factory,
							    int hashIndex,
							    size_t keyRange,
							    size_t numEntries,
							    float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc = &uintuintIdentityPerfectHash_DestroyTable;
	table->setupFunc = &uintuintIdentityPerfectHash_SetupTable;
	table->emptyFunc = &uintuintIdentityPerfectHash_EmptyTable;
	table->queryFunc = &uintuintIdentityPerfectHash_Query;
	table->querySingleFunc = &uintuintIdentityPerfectHash_QuerySingle;
	table->insertFunc = &uintuintIdentityPerfectHash_Insert;
	table->insertSingleFunc = &uintuintIdentityPerfectHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintIdentityPerfectHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintIdentityPerfectHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(uintuintIdentityPerfectHash_TableData));
	((uintuintIdentityPerfectHash_TableData *) table->tableData)->hashID =
	    IDENTITY_PERFECT_HASH_ID;
	((uintuintIdentityPerfectHash_TableData *) table->tableData)->
	    numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)malloc(sizeof(uintuintIdentityPerfectHash_TableData) +
			   ((uintuintIdentityPerfectHash_TableData *) table->
			    tableData)->numBuckets *
			   sizeof(uintuintIdentityPerfectHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintIdentityPerfectHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int uintuintIdentityPerfectHash_CreateFactory(uintuintHash_Factory * factory,
					      int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintIdentityPerfectHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintIdentityPerfectHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentityPerfectHash_DestroyFactory(uintuintHash_Factory * factory,
					       int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentityPerfectHash_DestroyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintIdentityPerfectHash_SetupTable(uintuintHash_Table * table) {
	int exitCode = 0;
	uintuintIdentityPerfectHash_Bucket *buckets =
	    (uintuintIdentityPerfectHash_Bucket *) & table->
	    tableData[sizeof(uintuintIdentityPerfectHash_TableData)];
	if (uintuintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
		for (uint index = 0;
		     index <
		     ((uintuintIdentityPerfectHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
		}
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintIdentityPerfectHash_EmptyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	uintuintIdentityPerfectHash_Bucket *buckets =
	    (uintuintIdentityPerfectHash_Bucket *) & table->
	    tableData[sizeof(uintuintIdentityPerfectHash_TableData)];
	for (uint index = 0;
	     index <
	     ((uintuintIdentityPerfectHash_TableData *) table->tableData)->
	     numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintIdentityPerfectHash_InnerQuerySingle(char *tableData, uint key,
						 uint * valueOutput) {
	uintuintIdentityPerfectHash_Bucket *buckets =
	    (uintuintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentityPerfectHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectHash_InnerQuery(char *tableData,
					   unsigned int numKeys, uint * keys,
					   uint * valuesOutput) {
	uintuintIdentityPerfectHash_Bucket *buckets =
	    (uintuintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectHash_TableData)];
	uint key;
	uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentityPerfectHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectHash_InnerInsertSingle(char *tableData, uint key,
						  uint value) {
	uintuintIdentityPerfectHash_Bucket *buckets =
	    (uintuintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentityPerfectHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectHash_InnerInsert(char *tableData,
					    unsigned int numEntries,
					    uint * keys, uint * values) {
	uintuintIdentityPerfectHash_Bucket *buckets =
	    (uintuintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentityPerfectHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectHash_InnerInsertSingleNoOverwrite(char *tableData,
							     uint key,
							     uint value) {
	uintuintIdentityPerfectHash_Bucket *buckets =
	    (uintuintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentityPerfectHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectHash_InnerInsertNoOverwrite(char *tableData,
						       unsigned int numEntries,
						       uint * keys,
						       uint * values) {
	uintuintIdentityPerfectHash_Bucket *buckets =
	    (uintuintIdentityPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentityPerfectHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectHash_QuerySingle(uintuintHash_Table * table,
					    uint key, uint * valueOutput) {
	return uintuintIdentityPerfectHash_InnerQuerySingle(table->tableData,
							    key, valueOutput);
}
int uintuintIdentityPerfectHash_Query(uintuintHash_Table * table,
				      size_t numKeys, uint * keys,
				      uint * valuesOutput) {
	return uintuintIdentityPerfectHash_InnerQuery(table->tableData, numKeys,
						      keys, valuesOutput);
}
int uintuintIdentityPerfectHash_InsertSingle(uintuintHash_Table * table,
					     uint key, uint value) {
	return uintuintIdentityPerfectHash_InnerInsertSingle(table->tableData,
							     key, value);
}
int uintuintIdentityPerfectHash_Insert(uintuintHash_Table * table,
				       size_t numEntries, uint * keys,
				       uint * values) {
	return uintuintIdentityPerfectHash_InnerInsert(table->tableData,
						       numEntries, keys,
						       values);
}
int uintuintIdentityPerfectHash_InsertSingleNoOverwrite(uintuintHash_Table *
							table, uint key,
							uint value) {
	return uintuintIdentityPerfectHash_InnerInsertSingleNoOverwrite(table->
									tableData,
									key,
									value);
}
int uintuintIdentityPerfectHash_InsertNoOverwrite(uintuintHash_Table * table,
						  size_t numEntries,
						  uint * keys, uint * values) {
	return uintuintIdentityPerfectHash_InnerInsertNoOverwrite(table->
								  tableData,
								  numEntries,
								  keys, values);
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
uintuintHash_Table
    *uintuintIdentityPerfectCLHash_CreateTable(uintuintHash_Factory * factory,
					       int hashIndex, size_t keyRange,
					       size_t numEntries,
					       float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc = &uintuintIdentityPerfectCLHash_DestroyTable;
	table->setupFunc = &uintuintIdentityPerfectCLHash_SetupTable;
	table->emptyFunc = &uintuintIdentityPerfectCLHash_EmptyTable;
	table->queryFunc = &uintuintIdentityPerfectCLHash_Query;
	table->querySingleFunc = &uintuintIdentityPerfectCLHash_QuerySingle;
	table->insertFunc = &uintuintIdentityPerfectCLHash_Insert;
	table->insertSingleFunc = &uintuintIdentityPerfectCLHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintIdentityPerfectCLHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintIdentityPerfectCLHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(uintuintIdentityPerfectCLHash_TableData));
	((uintuintIdentityPerfectCLHash_TableData *) table->tableData)->hashID =
	    IDENTITY_PERFECT_CL_HASH_ID;
	table->context = factory->context;
	table->queue = factory->queue;
	table->program = factory->program;
	table->localWorkSize = factory->localWorkSize;
	table->utilProgram = factory->utilProgram[hashIndex];
	table->emptyKernel = factory->emptyKernel[hashIndex];
	table->emptyKernelLocalWorkSize =
	    factory->emptyKernelLocalWorkSize[hashIndex];
	table->querySingleKernel = factory->querySingleKernel[hashIndex];
	table->insertSingleKernel = factory->insertSingleKernel[hashIndex];
	table->insertSingleNoOverwriteKernel =
	    factory->insertSingleNoOverwriteKernel[hashIndex];
	clRetainContext(table->context);
	clRetainCommandQueue(table->queue);
	clRetainProgram(table->program);
	clRetainProgram(table->utilProgram);
	clRetainKernel(table->emptyKernel);
	clRetainKernel(table->querySingleKernel);
	clRetainKernel(table->insertSingleKernel);
	clRetainKernel(table->insertSingleNoOverwriteKernel);;
	((uintuintIdentityPerfectCLHash_TableData *) table->tableData)->
	    numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)malloc(sizeof(uintuintIdentityPerfectCLHash_TableData) +
			   ((uintuintIdentityPerfectCLHash_TableData *) table->
			    tableData)->numBuckets *
			   sizeof(uintuintIdentityPerfectCLHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintIdentityPerfectCLHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	cl_int err;
	table->tableDataBuffer =
	    clCreateBuffer(table->context, CL_MEM_READ_WRITE,
			   sizeof(uintuintIdentityPerfectHash_TableData) +
			   ((uintuintIdentityPerfectHash_TableData *) table->
			    tableData)->numBuckets *
			   sizeof(uintuintIdentityPerfectHash_Bucket), NULL,
			   &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_InitTable",
				     "clCreateBuffer");
	err =
	    clEnqueueWriteBuffer(table->queue, table->tableDataBuffer, CL_TRUE,
				 0,
				 sizeof(uintuintIdentityPerfectHash_TableData),
				 table->tableData, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_InitTable",
				     "clEnqueueWriteBuffer");
	return table;
}
int uintuintIdentityPerfectCLHash_CreateFactory(uintuintHash_Factory * factory,
						int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintIdentityPerfectCLHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintIdentityPerfectCLHash_DestroyFactory;
	cl_int error;
	cl_device_id device;
	error =
	    clGetContextInfo(factory->context, CL_CONTEXT_DEVICES,
			     sizeof(device), &device, NULL);
	CLHash_Utilities_HandleError(error, "uintuintHash_CreateFactory",
				     "clGetContextInfo");
	factory->querySingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintIdentityPerfectCLHash_RangeQuerySingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintIdentityPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintIdentityPerfectCLHash_RangeInsertSingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintIdentityPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleNoOverwriteKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintIdentityPerfectCLHash_RangeInsertSingleNoOverwrite",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintIdentityPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->utilProgram[hashIndex] =
	    CLHash_Utilities_BuildProgramString(factory->context, device,
						"static inline unsigned int uintuintHash_CompressIdentity(char data, int hashCode){ return hashCode; } typedef struct uintuintHash_CompressLCGData{ long unsigned int a; long unsigned int c; unsigned int m; unsigned int n; }uintuintHash_CompressLCGData; static inline unsigned int uintuintHash_CompressLCG(uintuintHash_CompressLCGData compressLCGData, int hashCode){ return ((compressLCGData.a * hashCode + compressLCGData.c) % compressLCGData.m) % compressLCGData.n; } typedef struct uintuintIdentityPerfectCLHash_TableData{ int hashID; unsigned int numBuckets; char compressFuncData; }uintuintIdentityPerfectCLHash_TableData; typedef struct uintuintIdentityPerfectCLHash_Bucket{ uint key; uint value; }uintuintIdentityPerfectCLHash_Bucket; __kernel void uintuintIdentityPerfectCLHash_Empty(__global char *tableData){ int index = get_global_id(0); if(index >= ((__global uintuintIdentityPerfectCLHash_TableData*)tableData)->numBuckets){ return; } __global uintuintIdentityPerfectCLHash_Bucket *buckets = (__global uintuintIdentityPerfectCLHash_Bucket*)&tableData[sizeof(uintuintIdentityPerfectCLHash_TableData)]; buckets[index].key = -1;/*HASH_BUCKET_STATUS_EMPTY*/ }");
	factory->emptyKernel[hashIndex] =
	    clCreateKernel(factory->utilProgram[hashIndex],
			   "uintuintIdentityPerfectCLHash_Empty", &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintIdentityPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	error =
	    clGetKernelWorkGroupInfo(factory->emptyKernel[hashIndex], device,
				     CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
				     &factory->
				     emptyKernelLocalWorkSize[hashIndex], NULL);
	CLHash_Utilities_HandleError(error,
				     "uintuintIdentityPerfectCLHash_CreateFactory",
				     "clGetKernelWorkGroupInfo");;;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentityPerfectCLHash_DestroyFactory(uintuintHash_Factory * factory,
						 int hashIndex) {;
	clReleaseKernel(factory->emptyKernel[hashIndex]);
	clReleaseProgram(factory->utilProgram[hashIndex]);
	clReleaseKernel(factory->querySingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleNoOverwriteKernel[hashIndex]);;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentityPerfectCLHash_DestroyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	clReleaseMemObject(table->tableDataBuffer);
	clReleaseContext(table->context);
	clReleaseCommandQueue(table->queue);
	clReleaseProgram(table->utilProgram);
	clReleaseKernel(table->emptyKernel);
	clReleaseProgram(table->program);
	clReleaseKernel(table->querySingleKernel);
	clReleaseKernel(table->insertSingleKernel);
	clReleaseKernel(table->insertSingleNoOverwriteKernel);
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintIdentityPerfectCLHash_SetupTable(uintuintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((uintuintIdentityPerfectHash_TableData *) table->
			      tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int uintuintIdentityPerfectCLHash_EmptyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((uintuintIdentityPerfectHash_TableData *) table->
			      tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int uintuintIdentityPerfectCLHash_QuerySingle(uintuintHash_Table * table,
					      uint key, uint * valueOutput) {
	return uintuintIdentityPerfectCLHash_Query(table, 1, &key, valueOutput);
}
int uintuintIdentityPerfectCLHash_Query(uintuintHash_Table * table,
					size_t numKeys, uint * keys,
					uint * valuesOutput) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numKeys, keys, &err);
	CLHash_Utilities_HandleError(err, "uintuintIdentityPerfectCLHash_Query",
				     "clCreateBuffer");
	cl_mem valuesOutputBuffer =
	    clCreateBuffer(table->context, CL_MEM_WRITE_ONLY,
			   sizeof(uint) * numKeys, NULL, &err);
	CLHash_Utilities_HandleError(err, "uintuintIdentityPerfectCLHash_Query",
				     "clCreateBuffer");
	uintuintIdentityPerfectCLHash_BufferQuery(table, numKeys, keysBuffer,
						  valuesOutputBuffer);
	err =
	    clEnqueueReadBuffer(table->queue, valuesOutputBuffer, CL_TRUE, 0,
				sizeof(uint) * numKeys, valuesOutput, 0, NULL,
				NULL);
	CLHash_Utilities_HandleError(err, "uintuintIdentityPerfectCLHash_Query",
				     "clEnqueueReadBuffer");
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesOutputBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentityPerfectCLHash_BufferQuery(uintuintHash_Table * table,
					      size_t numKeys, cl_mem keysBuffer,
					      cl_mem valuesOutputBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->querySingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 1, sizeof(unsigned int),
			   &numKeys);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 3, sizeof(cl_mem),
			   &valuesOutputBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numKeys, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->querySingleKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferQuery",
				     "clEnqueueNDRangeKernel");
	clFinish(table->queue);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentityPerfectCLHash_InsertSingle(uintuintHash_Table * table,
					       uint key, uint value) {
	return uintuintIdentityPerfectCLHash_Insert(table, 1, &key, &value);
}
int uintuintIdentityPerfectCLHash_Insert(uintuintHash_Table * table,
					 size_t numEntries, uint * keys,
					 uint * values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_Insert",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_Insert",
				     "clCreateBuffer");
	uintuintIdentityPerfectCLHash_BufferInsert(table, numEntries,
						   keysBuffer, valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentityPerfectCLHash_BufferInsert(uintuintHash_Table * table,
					       size_t numEntries,
					       cl_mem keysBuffer,
					       cl_mem valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 1, sizeof(unsigned int),
			   &numEntries);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 3, sizeof(cl_mem),
			   &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->insertSingleKernel, 1,
				   0, &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err, NULL, "clEnqueueNDRangeKernel");
	return (0);
}
int uintuintIdentityPerfectCLHash_InsertSingleNoOverwrite(uintuintHash_Table *
							  table, uint key,
							  uint value) {
	return uintuintIdentityPerfectCLHash_InsertNoOverwrite(table, 1, &key,
							       &value);
}
int uintuintIdentityPerfectCLHash_InsertNoOverwrite(uintuintHash_Table * table,
						    size_t numEntries,
						    uint * keys,
						    uint * values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	uintuintIdentityPerfectCLHash_BufferInsertNoOverwrite(table, numEntries,
							      keysBuffer,
							      valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentityPerfectCLHash_BufferInsertNoOverwrite(uintuintHash_Table *
							  table,
							  size_t numEntries,
							  cl_mem keysBuffer,
							  cl_mem valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 0,
			   sizeof(cl_mem), &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 1,
			   sizeof(unsigned int), &numEntries);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferInsertNoOverwrite",
				     "ClSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 2,
			   sizeof(cl_mem), &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 3,
			   sizeof(cl_mem), &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue,
				   table->insertSingleNoOverwriteKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentityPerfectCLHash_BufferInsertNoOverwrite",
				     "clEnqueueNDRangeKernel");
	return (0);
}

typedef struct uintuintIdentityPerfectOpenMPHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
} uintuintIdentityPerfectOpenMPHash_TableData;
typedef struct uintuintIdentityPerfectOpenMPHash_Bucket {
	uint key;
	uint value;
} uintuintIdentityPerfectOpenMPHash_Bucket;
uintuintHash_Table
    *uintuintIdentityPerfectOpenMPHash_CreateTable(uintuintHash_Factory *
						   factory, int hashIndex,
						   size_t keyRange,
						   size_t numEntries,
						   float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc = &uintuintIdentityPerfectOpenMPHash_DestroyTable;
	table->setupFunc = &uintuintIdentityPerfectOpenMPHash_SetupTable;
	table->emptyFunc = &uintuintIdentityPerfectOpenMPHash_EmptyTable;
	table->queryFunc = &uintuintIdentityPerfectOpenMPHash_Query;
	table->querySingleFunc = &uintuintIdentityPerfectOpenMPHash_QuerySingle;
	table->insertFunc = &uintuintIdentityPerfectOpenMPHash_Insert;
	table->insertSingleFunc =
	    &uintuintIdentityPerfectOpenMPHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintIdentityPerfectOpenMPHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintIdentityPerfectOpenMPHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(uintuintIdentityPerfectOpenMPHash_TableData));
	((uintuintIdentityPerfectOpenMPHash_TableData *) table->tableData)->
	    hashID = IDENTITY_PERFECT_OPENMP_HASH_ID;
	((uintuintIdentityPerfectOpenMPHash_TableData *) table->tableData)->
	    numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)malloc(sizeof(uintuintIdentityPerfectOpenMPHash_TableData) +
			   ((uintuintIdentityPerfectOpenMPHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(uintuintIdentityPerfectOpenMPHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintIdentityPerfectOpenMPHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int uintuintIdentityPerfectOpenMPHash_CreateFactory(uintuintHash_Factory *
						    factory, int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintIdentityPerfectOpenMPHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintIdentityPerfectOpenMPHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentityPerfectOpenMPHash_DestroyFactory(uintuintHash_Factory *
						     factory, int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentityPerfectOpenMPHash_DestroyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintIdentityPerfectOpenMPHash_SetupTable(uintuintHash_Table * table) {
	int exitCode = 0;
	uintuintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentityPerfectOpenMPHash_Bucket *) & table->
	    tableData[sizeof(uintuintIdentityPerfectOpenMPHash_TableData)];
	if (uintuintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
#pragma omp parallel for
		for (uint index = 0;
		     index <
		     ((uintuintIdentityPerfectOpenMPHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
		}
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintIdentityPerfectOpenMPHash_EmptyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	uintuintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentityPerfectOpenMPHash_Bucket *) & table->
	    tableData[sizeof(uintuintIdentityPerfectOpenMPHash_TableData)];
#pragma omp parallel for
	for (uint index = 0;
	     index <
	     ((uintuintIdentityPerfectOpenMPHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintIdentityPerfectOpenMPHash_InnerQuerySingle(char *tableData,
						       uint key,
						       uint * valueOutput) {
	uintuintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectOpenMPHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectOpenMPHash_InnerQuery(char *tableData,
						 unsigned int numKeys,
						 uint * keys,
						 uint * valuesOutput) {
	uintuintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectOpenMPHash_TableData)];
	uint key;
	uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectOpenMPHash_InnerInsertSingle(char *tableData,
							uint key, uint value) {
	uintuintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectOpenMPHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectOpenMPHash_InnerInsert(char *tableData,
						  unsigned int numEntries,
						  uint * keys, uint * values) {
	uintuintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectOpenMPHash_InnerInsertSingleNoOverwrite(char
								   *tableData,
								   uint key,
								   uint value) {
	uintuintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectOpenMPHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectOpenMPHash_InnerInsertNoOverwrite(char *tableData,
							     unsigned int
							     numEntries,
							     uint * keys,
							     uint * values) {
	uintuintIdentityPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentityPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintIdentityPerfectOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentityPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
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
int uintuintIdentityPerfectOpenMPHash_QuerySingle(uintuintHash_Table * table,
						  uint key,
						  uint * valueOutput) {
	return uintuintIdentityPerfectOpenMPHash_InnerQuerySingle(table->
								  tableData,
								  key,
								  valueOutput);
}
int uintuintIdentityPerfectOpenMPHash_Query(uintuintHash_Table * table,
					    size_t numKeys, uint * keys,
					    uint * valuesOutput) {
	return uintuintIdentityPerfectOpenMPHash_InnerQuery(table->tableData,
							    numKeys, keys,
							    valuesOutput);
}
int uintuintIdentityPerfectOpenMPHash_InsertSingle(uintuintHash_Table * table,
						   uint key, uint value) {
	return uintuintIdentityPerfectOpenMPHash_InnerInsertSingle(table->
								   tableData,
								   key, value);
}
int uintuintIdentityPerfectOpenMPHash_Insert(uintuintHash_Table * table,
					     size_t numEntries, uint * keys,
					     uint * values) {
	return uintuintIdentityPerfectOpenMPHash_InnerInsert(table->tableData,
							     numEntries, keys,
							     values);
}
int uintuintIdentityPerfectOpenMPHash_InsertSingleNoOverwrite(uintuintHash_Table
							      * table, uint key,
							      uint value) {
	return
	    uintuintIdentityPerfectOpenMPHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int uintuintIdentityPerfectOpenMPHash_InsertNoOverwrite(uintuintHash_Table *
							table,
							size_t numEntries,
							uint * keys,
							uint * values) {
	return uintuintIdentityPerfectOpenMPHash_InnerInsertNoOverwrite(table->
									tableData,
									numEntries,
									keys,
									values);
}

typedef struct uintuintIdentitySentinelPerfectHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
	uint emptyValue;
} uintuintIdentitySentinelPerfectHash_TableData;
typedef struct uintuintIdentitySentinelPerfectHash_Bucket {
	uint value;
} uintuintIdentitySentinelPerfectHash_Bucket;
uintuintHash_Table
    *uintuintIdentitySentinelPerfectHash_CreateTable(uintuintHash_Factory *
						     factory, int hashIndex,
						     size_t keyRange,
						     size_t numEntries,
						     float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc = &uintuintIdentitySentinelPerfectHash_DestroyTable;
	table->setupFunc = &uintuintIdentitySentinelPerfectHash_SetupTable;
	table->emptyFunc = &uintuintIdentitySentinelPerfectHash_EmptyTable;
	table->queryFunc = &uintuintIdentitySentinelPerfectHash_Query;
	table->querySingleFunc =
	    &uintuintIdentitySentinelPerfectHash_QuerySingle;
	table->insertFunc = &uintuintIdentitySentinelPerfectHash_Insert;
	table->insertSingleFunc =
	    &uintuintIdentitySentinelPerfectHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintIdentitySentinelPerfectHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintIdentitySentinelPerfectHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(uintuintIdentitySentinelPerfectHash_TableData));
	((uintuintIdentitySentinelPerfectHash_TableData *) table->tableData)->
	    hashID = IDENTITY_SENTINEL_PERFECT_HASH_ID;
	((uintuintIdentitySentinelPerfectHash_TableData *) table->tableData)->
	    emptyValue = factory->emptyValue;
	((uintuintIdentitySentinelPerfectHash_TableData *) table->tableData)->
	    numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)malloc(sizeof(uintuintIdentitySentinelPerfectHash_TableData)
			   +
			   ((uintuintIdentitySentinelPerfectHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(uintuintIdentitySentinelPerfectHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintIdentitySentinelPerfectHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int uintuintIdentitySentinelPerfectHash_CreateFactory(uintuintHash_Factory *
						      factory, int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintIdentitySentinelPerfectHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintIdentitySentinelPerfectHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentitySentinelPerfectHash_DestroyFactory(uintuintHash_Factory *
						       factory,
						       int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentitySentinelPerfectHash_DestroyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintIdentitySentinelPerfectHash_SetupTable(uintuintHash_Table * table) {
	int exitCode = 0;
	uintuintIdentitySentinelPerfectHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectHash_Bucket *) & table->
	    tableData[sizeof(uintuintIdentitySentinelPerfectHash_TableData)];
	if (uintuintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
		for (uint index = 0;
		     index <
		     ((uintuintIdentitySentinelPerfectHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].value =
			    ((uintuintIdentitySentinelPerfectHash_TableData *)
			     table->tableData)->emptyValue;
		}
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintIdentitySentinelPerfectHash_EmptyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	uintuintIdentitySentinelPerfectHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectHash_Bucket *) & table->
	    tableData[sizeof(uintuintIdentitySentinelPerfectHash_TableData)];
	for (uint index = 0;
	     index <
	     ((uintuintIdentitySentinelPerfectHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].value =
		    ((uintuintIdentitySentinelPerfectHash_TableData *) table->
		     tableData)->emptyValue;
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintIdentitySentinelPerfectHash_InnerQuerySingle(char *tableData,
							 uint key,
							 uint * valueOutput) {
	uintuintIdentitySentinelPerfectHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((uintuintIdentitySentinelPerfectHash_TableData *) tableData)->
	    emptyValue) {
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
int uintuintIdentitySentinelPerfectHash_InnerQuery(char *tableData,
						   unsigned int numKeys,
						   uint * keys,
						   uint * valuesOutput) {
	uintuintIdentitySentinelPerfectHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectHash_TableData)];
	uint key;
	uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((uintuintIdentitySentinelPerfectHash_TableData *)
		     tableData)->emptyValue) {
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
int uintuintIdentitySentinelPerfectHash_InnerInsertSingle(char *tableData,
							  uint key,
							  uint value) {
	uintuintIdentitySentinelPerfectHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((uintuintIdentitySentinelPerfectHash_TableData *) tableData)->
	    emptyValue) {
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
int uintuintIdentitySentinelPerfectHash_InnerInsert(char *tableData,
						    unsigned int numEntries,
						    uint * keys,
						    uint * values) {
	uintuintIdentitySentinelPerfectHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((uintuintIdentitySentinelPerfectHash_TableData *)
		     tableData)->emptyValue) {
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
int uintuintIdentitySentinelPerfectHash_InnerInsertSingleNoOverwrite(char
								     *tableData,
								     uint key,
								     uint
								     value) {
	uintuintIdentitySentinelPerfectHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((uintuintIdentitySentinelPerfectHash_TableData *) tableData)->
	    emptyValue) {
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
int uintuintIdentitySentinelPerfectHash_InnerInsertNoOverwrite(char *tableData,
							       unsigned int
							       numEntries,
							       uint * keys,
							       uint * values) {
	uintuintIdentitySentinelPerfectHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectHash_Bucket *) &
	    tableData[sizeof(uintuintIdentitySentinelPerfectHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((uintuintIdentitySentinelPerfectHash_TableData *)
		     tableData)->emptyValue) {
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
int uintuintIdentitySentinelPerfectHash_QuerySingle(uintuintHash_Table * table,
						    uint key,
						    uint * valueOutput) {
	return uintuintIdentitySentinelPerfectHash_InnerQuerySingle(table->
								    tableData,
								    key,
								    valueOutput);
}
int uintuintIdentitySentinelPerfectHash_Query(uintuintHash_Table * table,
					      size_t numKeys, uint * keys,
					      uint * valuesOutput) {
	return uintuintIdentitySentinelPerfectHash_InnerQuery(table->tableData,
							      numKeys, keys,
							      valuesOutput);
}
int uintuintIdentitySentinelPerfectHash_InsertSingle(uintuintHash_Table * table,
						     uint key, uint value) {
	return uintuintIdentitySentinelPerfectHash_InnerInsertSingle(table->
								     tableData,
								     key,
								     value);
}
int uintuintIdentitySentinelPerfectHash_Insert(uintuintHash_Table * table,
					       size_t numEntries, uint * keys,
					       uint * values) {
	return uintuintIdentitySentinelPerfectHash_InnerInsert(table->tableData,
							       numEntries, keys,
							       values);
}
int
uintuintIdentitySentinelPerfectHash_InsertSingleNoOverwrite(uintuintHash_Table *
							    table, uint key,
							    uint value) {
	return
	    uintuintIdentitySentinelPerfectHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int uintuintIdentitySentinelPerfectHash_InsertNoOverwrite(uintuintHash_Table *
							  table,
							  size_t numEntries,
							  uint * keys,
							  uint * values) {
	return
	    uintuintIdentitySentinelPerfectHash_InnerInsertNoOverwrite(table->
								       tableData,
								       numEntries,
								       keys,
								       values);
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
uintuintHash_Table
    *uintuintIdentitySentinelPerfectCLHash_CreateTable(uintuintHash_Factory *
						       factory, int hashIndex,
						       size_t keyRange,
						       size_t numEntries,
						       float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc =
	    &uintuintIdentitySentinelPerfectCLHash_DestroyTable;
	table->setupFunc = &uintuintIdentitySentinelPerfectCLHash_SetupTable;
	table->emptyFunc = &uintuintIdentitySentinelPerfectCLHash_EmptyTable;
	table->queryFunc = &uintuintIdentitySentinelPerfectCLHash_Query;
	table->querySingleFunc =
	    &uintuintIdentitySentinelPerfectCLHash_QuerySingle;
	table->insertFunc = &uintuintIdentitySentinelPerfectCLHash_Insert;
	table->insertSingleFunc =
	    &uintuintIdentitySentinelPerfectCLHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintIdentitySentinelPerfectCLHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(uintuintIdentitySentinelPerfectCLHash_TableData));
	((uintuintIdentitySentinelPerfectCLHash_TableData *) table->tableData)->
	    hashID = IDENTITY_SENTINEL_PERFECT_CL_HASH_ID;
	table->context = factory->context;
	table->queue = factory->queue;
	table->program = factory->program;
	table->localWorkSize = factory->localWorkSize;
	table->utilProgram = factory->utilProgram[hashIndex];
	table->emptyKernel = factory->emptyKernel[hashIndex];
	table->emptyKernelLocalWorkSize =
	    factory->emptyKernelLocalWorkSize[hashIndex];
	table->querySingleKernel = factory->querySingleKernel[hashIndex];
	table->insertSingleKernel = factory->insertSingleKernel[hashIndex];
	table->insertSingleNoOverwriteKernel =
	    factory->insertSingleNoOverwriteKernel[hashIndex];
	clRetainContext(table->context);
	clRetainCommandQueue(table->queue);
	clRetainProgram(table->program);
	clRetainProgram(table->utilProgram);
	clRetainKernel(table->emptyKernel);
	clRetainKernel(table->querySingleKernel);
	clRetainKernel(table->insertSingleKernel);
	clRetainKernel(table->insertSingleNoOverwriteKernel);;
	((uintuintIdentitySentinelPerfectCLHash_TableData *) table->tableData)->
	    emptyValue = factory->emptyValue;
	((uintuintIdentitySentinelPerfectCLHash_TableData *) table->tableData)->
	    numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)
	    malloc(sizeof(uintuintIdentitySentinelPerfectCLHash_TableData) +
		   ((uintuintIdentitySentinelPerfectCLHash_TableData *) table->
		    tableData)->numBuckets *
		   sizeof(uintuintIdentitySentinelPerfectCLHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintIdentitySentinelPerfectCLHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	cl_int err;
	table->tableDataBuffer =
	    clCreateBuffer(table->context, CL_MEM_READ_WRITE,
			   sizeof(uintuintIdentitySentinelPerfectHash_TableData)
			   +
			   ((uintuintIdentitySentinelPerfectHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(uintuintIdentitySentinelPerfectHash_Bucket),
			   NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_InitTable",
				     "clCreateBuffer");
	err =
	    clEnqueueWriteBuffer(table->queue, table->tableDataBuffer, CL_TRUE,
				 0,
				 sizeof
				 (uintuintIdentitySentinelPerfectHash_TableData),
				 table->tableData, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_InitTable",
				     "clEnqueueWriteBuffer");
	return table;
}
int uintuintIdentitySentinelPerfectCLHash_CreateFactory(uintuintHash_Factory *
							factory,
							int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintIdentitySentinelPerfectCLHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintIdentitySentinelPerfectCLHash_DestroyFactory;
	cl_int error;
	cl_device_id device;
	error =
	    clGetContextInfo(factory->context, CL_CONTEXT_DEVICES,
			     sizeof(device), &device, NULL);
	CLHash_Utilities_HandleError(error, "uintuintHash_CreateFactory",
				     "clGetContextInfo");
	factory->querySingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintIdentitySentinelPerfectCLHash_RangeQuerySingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintIdentitySentinelPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintIdentitySentinelPerfectCLHash_RangeInsertSingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintIdentitySentinelPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleNoOverwriteKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintIdentitySentinelPerfectCLHash_RangeInsertSingleNoOverwrite",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintIdentitySentinelPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	factory->utilProgram[hashIndex] =
	    CLHash_Utilities_BuildProgramString(factory->context, device,
						"static inline unsigned int uintuintHash_CompressIdentity(char data, int hashCode){ return hashCode; } typedef struct uintuintHash_CompressLCGData{ long unsigned int a; long unsigned int c; unsigned int m; unsigned int n; }uintuintHash_CompressLCGData; static inline unsigned int uintuintHash_CompressLCG(uintuintHash_CompressLCGData compressLCGData, int hashCode){ return ((compressLCGData.a * hashCode + compressLCGData.c) % compressLCGData.m) % compressLCGData.n; } typedef struct uintuintIdentitySentinelPerfectCLHash_TableData{ int hashID; unsigned int numBuckets; char compressFuncData; uint emptyValue; }uintuintIdentitySentinelPerfectCLHash_TableData; typedef struct uintuintIdentitySentinelPerfectCLHash_Bucket{ uint value; }uintuintIdentitySentinelPerfectCLHash_Bucket; __kernel void uintuintIdentitySentinelPerfectCLHash_Empty(__global char *tableData){ int index = get_global_id(0); if(index >= ((__global uintuintIdentitySentinelPerfectCLHash_TableData*)tableData)->numBuckets){ return; } __global uintuintIdentitySentinelPerfectCLHash_Bucket *buckets = (__global uintuintIdentitySentinelPerfectCLHash_Bucket*)&tableData[sizeof(uintuintIdentitySentinelPerfectCLHash_TableData)]; buckets[index].value = ((__global uintuintIdentitySentinelPerfectCLHash_TableData*)tableData)->emptyValue; }");
	factory->emptyKernel[hashIndex] =
	    clCreateKernel(factory->utilProgram[hashIndex],
			   "uintuintIdentitySentinelPerfectCLHash_Empty",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintIdentitySentinelPerfectCLHash_CreateFactory",
				     "clCreateKernel");
	error =
	    clGetKernelWorkGroupInfo(factory->emptyKernel[hashIndex], device,
				     CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
				     &factory->
				     emptyKernelLocalWorkSize[hashIndex], NULL);
	CLHash_Utilities_HandleError(error,
				     "uintuintIdentitySentinelPerfectCLHash_CreateFactory",
				     "clGetKernelWorkGroupInfo");;;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentitySentinelPerfectCLHash_DestroyFactory(uintuintHash_Factory *
							 factory,
							 int hashIndex) {;
	clReleaseKernel(factory->emptyKernel[hashIndex]);
	clReleaseProgram(factory->utilProgram[hashIndex]);
	clReleaseKernel(factory->querySingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleNoOverwriteKernel[hashIndex]);;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentitySentinelPerfectCLHash_DestroyTable(uintuintHash_Table *
						       table) {
	int exitCode = 0;
	clReleaseMemObject(table->tableDataBuffer);
	clReleaseContext(table->context);
	clReleaseCommandQueue(table->queue);
	clReleaseProgram(table->utilProgram);
	clReleaseKernel(table->emptyKernel);
	clReleaseProgram(table->program);
	clReleaseKernel(table->querySingleKernel);
	clReleaseKernel(table->insertSingleKernel);
	clReleaseKernel(table->insertSingleNoOverwriteKernel);
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintIdentitySentinelPerfectCLHash_SetupTable(uintuintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((uintuintIdentitySentinelPerfectHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int uintuintIdentitySentinelPerfectCLHash_EmptyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((uintuintIdentitySentinelPerfectHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int uintuintIdentitySentinelPerfectCLHash_QuerySingle(uintuintHash_Table *
						      table, uint key,
						      uint * valueOutput) {
	return uintuintIdentitySentinelPerfectCLHash_Query(table, 1, &key,
							   valueOutput);
}
int uintuintIdentitySentinelPerfectCLHash_Query(uintuintHash_Table * table,
						size_t numKeys, uint * keys,
						uint * valuesOutput) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numKeys, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_Query",
				     "clCreateBuffer");
	cl_mem valuesOutputBuffer =
	    clCreateBuffer(table->context, CL_MEM_WRITE_ONLY,
			   sizeof(uint) * numKeys, NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_Query",
				     "clCreateBuffer");
	uintuintIdentitySentinelPerfectCLHash_BufferQuery(table, numKeys,
							  keysBuffer,
							  valuesOutputBuffer);
	err =
	    clEnqueueReadBuffer(table->queue, valuesOutputBuffer, CL_TRUE, 0,
				sizeof(uint) * numKeys, valuesOutput, 0, NULL,
				NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_Query",
				     "clEnqueueReadBuffer");
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesOutputBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentitySentinelPerfectCLHash_BufferQuery(uintuintHash_Table *
						      table, size_t numKeys,
						      cl_mem keysBuffer,
						      cl_mem
						      valuesOutputBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->querySingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 1, sizeof(unsigned int),
			   &numKeys);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 3, sizeof(cl_mem),
			   &valuesOutputBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferQuery",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numKeys, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->querySingleKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferQuery",
				     "clEnqueueNDRangeKernel");
	clFinish(table->queue);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentitySentinelPerfectCLHash_InsertSingle(uintuintHash_Table *
						       table, uint key,
						       uint value) {
	return uintuintIdentitySentinelPerfectCLHash_Insert(table, 1, &key,
							    &value);
}
int uintuintIdentitySentinelPerfectCLHash_Insert(uintuintHash_Table * table,
						 size_t numEntries, uint * keys,
						 uint * values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_Insert",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_Insert",
				     "clCreateBuffer");
	uintuintIdentitySentinelPerfectCLHash_BufferInsert(table, numEntries,
							   keysBuffer,
							   valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentitySentinelPerfectCLHash_BufferInsert(uintuintHash_Table *
						       table, size_t numEntries,
						       cl_mem keysBuffer,
						       cl_mem valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 1, sizeof(unsigned int),
			   &numEntries);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 3, sizeof(cl_mem),
			   &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferInsert",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->insertSingleKernel, 1,
				   0, &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err, NULL, "clEnqueueNDRangeKernel");
	return (0);
}
int
uintuintIdentitySentinelPerfectCLHash_InsertSingleNoOverwrite(uintuintHash_Table
							      * table, uint key,
							      uint value) {
	return uintuintIdentitySentinelPerfectCLHash_InsertNoOverwrite(table, 1,
								       &key,
								       &value);
}
int uintuintIdentitySentinelPerfectCLHash_InsertNoOverwrite(uintuintHash_Table *
							    table,
							    size_t numEntries,
							    uint * keys,
							    uint * values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	uintuintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite(table,
								      numEntries,
								      keysBuffer,
								      valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int
uintuintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite(uintuintHash_Table
							      * table,
							      size_t numEntries,
							      cl_mem keysBuffer,
							      cl_mem
							      valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 0,
			   sizeof(cl_mem), &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 1,
			   sizeof(unsigned int), &numEntries);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite",
				     "ClSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 2,
			   sizeof(cl_mem), &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 3,
			   sizeof(cl_mem), &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue,
				   table->insertSingleNoOverwriteKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintIdentitySentinelPerfectCLHash_BufferInsertNoOverwrite",
				     "clEnqueueNDRangeKernel");
	return (0);
}

typedef struct uintuintIdentitySentinelPerfectOpenMPHash_TableData {
	int hashID;
	unsigned int numBuckets;
	char compressFuncData;
	uint emptyValue;
} uintuintIdentitySentinelPerfectOpenMPHash_TableData;
typedef struct uintuintIdentitySentinelPerfectOpenMPHash_Bucket {
	uint value;
} uintuintIdentitySentinelPerfectOpenMPHash_Bucket;
uintuintHash_Table
    *uintuintIdentitySentinelPerfectOpenMPHash_CreateTable(uintuintHash_Factory
							   * factory,
							   int hashIndex,
							   size_t keyRange,
							   size_t numEntries,
							   float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc =
	    &uintuintIdentitySentinelPerfectOpenMPHash_DestroyTable;
	table->setupFunc =
	    &uintuintIdentitySentinelPerfectOpenMPHash_SetupTable;
	table->emptyFunc =
	    &uintuintIdentitySentinelPerfectOpenMPHash_EmptyTable;
	table->queryFunc = &uintuintIdentitySentinelPerfectOpenMPHash_Query;
	table->querySingleFunc =
	    &uintuintIdentitySentinelPerfectOpenMPHash_QuerySingle;
	table->insertFunc = &uintuintIdentitySentinelPerfectOpenMPHash_Insert;
	table->insertSingleFunc =
	    &uintuintIdentitySentinelPerfectOpenMPHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintIdentitySentinelPerfectOpenMPHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintIdentitySentinelPerfectOpenMPHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(uintuintIdentitySentinelPerfectOpenMPHash_TableData));
	((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) table->
	 tableData)->hashID = IDENTITY_SENTINEL_PERFECT_OPENMP_HASH_ID;
	((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) table->
	 tableData)->emptyValue = factory->emptyValue;
	((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) table->
	 tableData)->numBuckets = keyRange + 1;
	char *tempHashData =
	    (char *)
	    malloc(sizeof(uintuintIdentitySentinelPerfectOpenMPHash_TableData) +
		   ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *)
		    table->tableData)->numBuckets *
		   sizeof(uintuintIdentitySentinelPerfectOpenMPHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintIdentitySentinelPerfectOpenMPHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int uintuintIdentitySentinelPerfectOpenMPHash_CreateFactory(uintuintHash_Factory
							    * factory,
							    int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintIdentitySentinelPerfectOpenMPHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintIdentitySentinelPerfectOpenMPHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int
uintuintIdentitySentinelPerfectOpenMPHash_DestroyFactory(uintuintHash_Factory *
							 factory,
							 int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintIdentitySentinelPerfectOpenMPHash_DestroyTable(uintuintHash_Table *
							   table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintIdentitySentinelPerfectOpenMPHash_SetupTable(uintuintHash_Table *
							 table) {
	int exitCode = 0;
	uintuintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectOpenMPHash_Bucket *) & table->
	    tableData[sizeof
		      (uintuintIdentitySentinelPerfectOpenMPHash_TableData)];
	if (uintuintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
#pragma omp parallel for
		for (uint index = 0;
		     index <
		     ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *)
		      table->tableData)->numBuckets; index++) {
			buckets[index].value =
			    ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) table->tableData)->emptyValue;
		}
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintIdentitySentinelPerfectOpenMPHash_EmptyTable(uintuintHash_Table *
							 table) {
	int exitCode = 0;
	uintuintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectOpenMPHash_Bucket *) & table->
	    tableData[sizeof
		      (uintuintIdentitySentinelPerfectOpenMPHash_TableData)];
#pragma omp parallel for
	for (uint index = 0;
	     index <
	     ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].value =
		    ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *)
		     table->tableData)->emptyValue;
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintIdentitySentinelPerfectOpenMPHash_InnerQuerySingle(char *tableData,
							       uint key,
							       uint *
							       valueOutput) {
	uintuintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintIdentitySentinelPerfectOpenMPHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *)
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
int uintuintIdentitySentinelPerfectOpenMPHash_InnerQuery(char *tableData,
							 unsigned int numKeys,
							 uint * keys,
							 uint * valuesOutput) {
	uintuintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintIdentitySentinelPerfectOpenMPHash_TableData)];
	uint key;
	uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *)
		     tableData)->emptyValue) {
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
int uintuintIdentitySentinelPerfectOpenMPHash_InnerInsertSingle(char *tableData,
								uint key,
								uint value) {
	uintuintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintIdentitySentinelPerfectOpenMPHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *)
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
int uintuintIdentitySentinelPerfectOpenMPHash_InnerInsert(char *tableData,
							  unsigned int
							  numEntries,
							  uint * keys,
							  uint * values) {
	uintuintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintIdentitySentinelPerfectOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *)
		     tableData)->emptyValue) {
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
int uintuintIdentitySentinelPerfectOpenMPHash_InnerInsertSingleNoOverwrite(char
									   *tableData,
									   uint
									   key,
									   uint
									   value) 
{
	uintuintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintIdentitySentinelPerfectOpenMPHash_TableData)];
	uint index;
	int exitCode;
	index =
	    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
	if (buckets[index].value !=
	    ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *)
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
int uintuintIdentitySentinelPerfectOpenMPHash_InnerInsertNoOverwrite(char
								     *tableData,
								     unsigned
								     int
								     numEntries,
								     uint *
								     keys,
								     uint *
								     values) {
	uintuintIdentitySentinelPerfectOpenMPHash_Bucket *buckets =
	    (uintuintIdentitySentinelPerfectOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintIdentitySentinelPerfectOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		index =
		    uintuintHash_CompressIdentity(((uintuintIdentitySentinelPerfectOpenMPHash_TableData *) tableData)->compressFuncData, key);
		if (buckets[index].value !=
		    ((uintuintIdentitySentinelPerfectOpenMPHash_TableData *)
		     tableData)->emptyValue) {
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
int uintuintIdentitySentinelPerfectOpenMPHash_QuerySingle(uintuintHash_Table *
							  table, uint key,
							  uint * valueOutput) {
	return
	    uintuintIdentitySentinelPerfectOpenMPHash_InnerQuerySingle(table->
								       tableData,
								       key,
								       valueOutput);
}
int uintuintIdentitySentinelPerfectOpenMPHash_Query(uintuintHash_Table * table,
						    size_t numKeys, uint * keys,
						    uint * valuesOutput) {
	return uintuintIdentitySentinelPerfectOpenMPHash_InnerQuery(table->
								    tableData,
								    numKeys,
								    keys,
								    valuesOutput);
}
int uintuintIdentitySentinelPerfectOpenMPHash_InsertSingle(uintuintHash_Table *
							   table, uint key,
							   uint value) {
	return
	    uintuintIdentitySentinelPerfectOpenMPHash_InnerInsertSingle(table->
									tableData,
									key,
									value);
}
int uintuintIdentitySentinelPerfectOpenMPHash_Insert(uintuintHash_Table * table,
						     size_t numEntries,
						     uint * keys,
						     uint * values) {
	return uintuintIdentitySentinelPerfectOpenMPHash_InnerInsert(table->
								     tableData,
								     numEntries,
								     keys,
								     values);
}
int
uintuintIdentitySentinelPerfectOpenMPHash_InsertSingleNoOverwrite
(uintuintHash_Table * table, uint key, uint value) {
	return
	    uintuintIdentitySentinelPerfectOpenMPHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int
uintuintIdentitySentinelPerfectOpenMPHash_InsertNoOverwrite(uintuintHash_Table *
							    table,
							    size_t numEntries,
							    uint * keys,
							    uint * values) {
	return
	    uintuintIdentitySentinelPerfectOpenMPHash_InnerInsertNoOverwrite
	    (table->tableData, numEntries, keys, values);
}

typedef struct uintuintLCGLinearOpenCompactHash_TableData {
	int hashID;
	unsigned int numBuckets;
	uintuintHash_CompressLCGData compressFuncData;
} uintuintLCGLinearOpenCompactHash_TableData;
typedef struct uintuintLCGLinearOpenCompactHash_Bucket {
	uint key;
	uint value;
} uintuintLCGLinearOpenCompactHash_Bucket;
uintuintHash_Table
    *uintuintLCGLinearOpenCompactHash_CreateTable(uintuintHash_Factory *
						  factory, int hashIndex,
						  size_t keyRange,
						  size_t numEntries,
						  float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc = &uintuintLCGLinearOpenCompactHash_DestroyTable;
	table->setupFunc = &uintuintLCGLinearOpenCompactHash_SetupTable;
	table->emptyFunc = &uintuintLCGLinearOpenCompactHash_EmptyTable;
	table->queryFunc = &uintuintLCGLinearOpenCompactHash_Query;
	table->querySingleFunc = &uintuintLCGLinearOpenCompactHash_QuerySingle;
	table->insertFunc = &uintuintLCGLinearOpenCompactHash_Insert;
	table->insertSingleFunc =
	    &uintuintLCGLinearOpenCompactHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintLCGLinearOpenCompactHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintLCGLinearOpenCompactHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)malloc(sizeof(uintuintLCGLinearOpenCompactHash_TableData));
	((uintuintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    hashID = LCG_LINEAR_OPEN_COMPACT_HASH_ID;
	((uintuintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    numBuckets = (unsigned int)((double)numEntries / loadFactor);
	((uintuintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.a = HASH_LCG_A;
	((uintuintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.c = HASH_LCG_C;
	((uintuintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.m = HASH_LCG_M;
	((uintuintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.n =
	    ((uintuintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	    numBuckets;
	char *tempHashData =
	    (char *)malloc(sizeof(uintuintLCGLinearOpenCompactHash_TableData) +
			   ((uintuintLCGLinearOpenCompactHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(uintuintLCGLinearOpenCompactHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintLCGLinearOpenCompactHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int uintuintLCGLinearOpenCompactHash_CreateFactory(uintuintHash_Factory *
						   factory, int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintLCGLinearOpenCompactHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintLCGLinearOpenCompactHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGLinearOpenCompactHash_DestroyFactory(uintuintHash_Factory *
						    factory, int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGLinearOpenCompactHash_DestroyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintLCGLinearOpenCompactHash_SetupTable(uintuintHash_Table * table) {
	int exitCode = 0;
	uintuintLCGLinearOpenCompactHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactHash_Bucket *) & table->
	    tableData[sizeof(uintuintLCGLinearOpenCompactHash_TableData)];
	if (uintuintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
		for (uint index = 0;
		     index <
		     ((uintuintLCGLinearOpenCompactHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
		}
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintLCGLinearOpenCompactHash_EmptyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	uintuintLCGLinearOpenCompactHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactHash_Bucket *) & table->
	    tableData[sizeof(uintuintLCGLinearOpenCompactHash_TableData)];
	for (uint index = 0;
	     index <
	     ((uintuintLCGLinearOpenCompactHash_TableData *) table->tableData)->
	     numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintLCGLinearOpenCompactHash_InnerQuerySingle(char *tableData, uint key,
						      uint * valueOutput) {
	uintuintLCGLinearOpenCompactHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGLinearOpenCompactHash_TableData *mytableData =
	    (uintuintLCGLinearOpenCompactHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((uintuintLCGLinearOpenCompactHash_TableData *)
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
int uintuintLCGLinearOpenCompactHash_InnerQuery(char *tableData,
						unsigned int numKeys,
						uint * keys,
						uint * valuesOutput) {
	uintuintLCGLinearOpenCompactHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactHash_TableData)];
	uint key;
	uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		uintuintLCGLinearOpenCompactHash_TableData *mytableData =
		    (uintuintLCGLinearOpenCompactHash_TableData *) tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((uintuintLCGLinearOpenCompactHash_TableData *)
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
int uintuintLCGLinearOpenCompactHash_InnerInsertSingle(char *tableData,
						       uint key, uint value) {
	uintuintLCGLinearOpenCompactHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGLinearOpenCompactHash_TableData *mytableData =
	    (uintuintLCGLinearOpenCompactHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((uintuintLCGLinearOpenCompactHash_TableData *)
		      tableData)->numBuckets);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
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
int uintuintLCGLinearOpenCompactHash_InnerInsert(char *tableData,
						 unsigned int numEntries,
						 uint * keys, uint * values) {
	uintuintLCGLinearOpenCompactHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		uintuintLCGLinearOpenCompactHash_TableData *mytableData =
		    (uintuintLCGLinearOpenCompactHash_TableData *) tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((uintuintLCGLinearOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if (((buckets[index].key ==
			      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
							   key,
							   HASH_BUCKET_STATUS_EMPTY)
			     : buckets[index].key) ==
			    HASH_BUCKET_STATUS_EMPTY) {
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
int uintuintLCGLinearOpenCompactHash_InnerInsertSingleNoOverwrite(char
								  *tableData,
								  uint key,
								  uint value) {
	uintuintLCGLinearOpenCompactHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGLinearOpenCompactHash_TableData *mytableData =
	    (uintuintLCGLinearOpenCompactHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((uintuintLCGLinearOpenCompactHash_TableData *)
		      tableData)->numBuckets);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
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
int uintuintLCGLinearOpenCompactHash_InnerInsertNoOverwrite(char *tableData,
							    unsigned int
							    numEntries,
							    uint * keys,
							    uint * values) {
	uintuintLCGLinearOpenCompactHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		uintuintLCGLinearOpenCompactHash_TableData *mytableData =
		    (uintuintLCGLinearOpenCompactHash_TableData *) tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration +
			      c) %
			     ((uintuintLCGLinearOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if (((buckets[index].key ==
			      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
							   key,
							   HASH_BUCKET_STATUS_EMPTY)
			     : buckets[index].key) ==
			    HASH_BUCKET_STATUS_EMPTY) {
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
int uintuintLCGLinearOpenCompactHash_QuerySingle(uintuintHash_Table * table,
						 uint key, uint * valueOutput) {
	return uintuintLCGLinearOpenCompactHash_InnerQuerySingle(table->
								 tableData, key,
								 valueOutput);
}
int uintuintLCGLinearOpenCompactHash_Query(uintuintHash_Table * table,
					   size_t numKeys, uint * keys,
					   uint * valuesOutput) {
	return uintuintLCGLinearOpenCompactHash_InnerQuery(table->tableData,
							   numKeys, keys,
							   valuesOutput);
}
int uintuintLCGLinearOpenCompactHash_InsertSingle(uintuintHash_Table * table,
						  uint key, uint value) {
	return uintuintLCGLinearOpenCompactHash_InnerInsertSingle(table->
								  tableData,
								  key, value);
}
int uintuintLCGLinearOpenCompactHash_Insert(uintuintHash_Table * table,
					    size_t numEntries, uint * keys,
					    uint * values) {
	return uintuintLCGLinearOpenCompactHash_InnerInsert(table->tableData,
							    numEntries, keys,
							    values);
}
int uintuintLCGLinearOpenCompactHash_InsertSingleNoOverwrite(uintuintHash_Table
							     * table, uint key,
							     uint value) {
	return
	    uintuintLCGLinearOpenCompactHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int uintuintLCGLinearOpenCompactHash_InsertNoOverwrite(uintuintHash_Table *
						       table, size_t numEntries,
						       uint * keys,
						       uint * values) {
	return uintuintLCGLinearOpenCompactHash_InnerInsertNoOverwrite(table->
								       tableData,
								       numEntries,
								       keys,
								       values);
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
uintuintHash_Table
    *uintuintLCGLinearOpenCompactCLHash_CreateTable(uintuintHash_Factory *
						    factory, int hashIndex,
						    size_t keyRange,
						    size_t numEntries,
						    float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc = &uintuintLCGLinearOpenCompactCLHash_DestroyTable;
	table->setupFunc = &uintuintLCGLinearOpenCompactCLHash_SetupTable;
	table->emptyFunc = &uintuintLCGLinearOpenCompactCLHash_EmptyTable;
	table->queryFunc = &uintuintLCGLinearOpenCompactCLHash_Query;
	table->querySingleFunc =
	    &uintuintLCGLinearOpenCompactCLHash_QuerySingle;
	table->insertFunc = &uintuintLCGLinearOpenCompactCLHash_Insert;
	table->insertSingleFunc =
	    &uintuintLCGLinearOpenCompactCLHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintLCGLinearOpenCompactCLHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(uintuintLCGLinearOpenCompactCLHash_TableData));
	((uintuintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    hashID = LCG_LINEAR_OPEN_COMPACT_CL_HASH_ID;
	table->context = factory->context;
	table->queue = factory->queue;
	table->program = factory->program;
	table->localWorkSize = factory->localWorkSize;
	table->utilProgram = factory->utilProgram[hashIndex];
	table->emptyKernel = factory->emptyKernel[hashIndex];
	table->emptyKernelLocalWorkSize =
	    factory->emptyKernelLocalWorkSize[hashIndex];
	table->querySingleKernel = factory->querySingleKernel[hashIndex];
	table->insertSingleKernel = factory->insertSingleKernel[hashIndex];
	table->insertSingleNoOverwriteKernel =
	    factory->insertSingleNoOverwriteKernel[hashIndex];
	clRetainContext(table->context);
	clRetainCommandQueue(table->queue);
	clRetainProgram(table->program);
	clRetainProgram(table->utilProgram);
	clRetainKernel(table->emptyKernel);
	clRetainKernel(table->querySingleKernel);
	clRetainKernel(table->insertSingleKernel);
	clRetainKernel(table->insertSingleNoOverwriteKernel);;
	((uintuintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    numBuckets = (unsigned int)((double)numEntries / loadFactor);
	((uintuintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.a = HASH_LCG_A;
	((uintuintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.c = HASH_LCG_C;
	((uintuintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.m = HASH_LCG_M;
	((uintuintLCGLinearOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.n =
	    ((uintuintLCGLinearOpenCompactCLHash_TableData *) table->
	     tableData)->numBuckets;
	char *tempHashData =
	    (char *)malloc(sizeof(uintuintLCGLinearOpenCompactCLHash_TableData)
			   +
			   ((uintuintLCGLinearOpenCompactCLHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(uintuintLCGLinearOpenCompactCLHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintLCGLinearOpenCompactCLHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	cl_int err;
	table->tableDataBuffer =
	    clCreateBuffer(table->context, CL_MEM_READ_WRITE,
			   sizeof(uintuintLCGLinearOpenCompactHash_TableData) +
			   ((uintuintLCGLinearOpenCompactHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(uintuintLCGLinearOpenCompactHash_Bucket),
			   NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_InitTable",
				     "clCreateBuffer");
	err =
	    clEnqueueWriteBuffer(table->queue, table->tableDataBuffer, CL_TRUE,
				 0,
				 sizeof
				 (uintuintLCGLinearOpenCompactHash_TableData),
				 table->tableData, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_InitTable",
				     "clEnqueueWriteBuffer");
	return table;
}
int uintuintLCGLinearOpenCompactCLHash_CreateFactory(uintuintHash_Factory *
						     factory, int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintLCGLinearOpenCompactCLHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintLCGLinearOpenCompactCLHash_DestroyFactory;
	cl_int error;
	cl_device_id device;
	error =
	    clGetContextInfo(factory->context, CL_CONTEXT_DEVICES,
			     sizeof(device), &device, NULL);
	CLHash_Utilities_HandleError(error, "uintuintHash_CreateFactory",
				     "clGetContextInfo");
	factory->querySingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintLCGLinearOpenCompactCLHash_RangeQuerySingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintLCGLinearOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintLCGLinearOpenCompactCLHash_RangeInsertSingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintLCGLinearOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleNoOverwriteKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintLCGLinearOpenCompactCLHash_RangeInsertSingleNoOverwrite",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintLCGLinearOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->utilProgram[hashIndex] =
	    CLHash_Utilities_BuildProgramString(factory->context, device,
						"static inline unsigned int uintuintHash_CompressIdentity(char data, int hashCode){ return hashCode; } typedef struct uintuintHash_CompressLCGData{ long unsigned int a; long unsigned int c; unsigned int m; unsigned int n; }uintuintHash_CompressLCGData; static inline unsigned int uintuintHash_CompressLCG(uintuintHash_CompressLCGData compressLCGData, int hashCode){ return ((compressLCGData.a * hashCode + compressLCGData.c) % compressLCGData.m) % compressLCGData.n; } typedef struct uintuintLCGLinearOpenCompactCLHash_TableData{ int hashID; unsigned int numBuckets; uintuintHash_CompressLCGData compressFuncData; }uintuintLCGLinearOpenCompactCLHash_TableData; typedef struct uintuintLCGLinearOpenCompactCLHash_Bucket{ uint key; uint value; }uintuintLCGLinearOpenCompactCLHash_Bucket; __kernel void uintuintLCGLinearOpenCompactCLHash_Empty(__global char *tableData){ int index = get_global_id(0); if(index >= ((__global uintuintLCGLinearOpenCompactCLHash_TableData*)tableData)->numBuckets){ return; } __global uintuintLCGLinearOpenCompactCLHash_Bucket *buckets = (__global uintuintLCGLinearOpenCompactCLHash_Bucket*)&tableData[sizeof(uintuintLCGLinearOpenCompactCLHash_TableData)]; buckets[index].key = -1;/*HASH_BUCKET_STATUS_EMPTY*/ }");
	factory->emptyKernel[hashIndex] =
	    clCreateKernel(factory->utilProgram[hashIndex],
			   "uintuintLCGLinearOpenCompactCLHash_Empty", &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintLCGLinearOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	error =
	    clGetKernelWorkGroupInfo(factory->emptyKernel[hashIndex], device,
				     CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
				     &factory->
				     emptyKernelLocalWorkSize[hashIndex], NULL);
	CLHash_Utilities_HandleError(error,
				     "uintuintLCGLinearOpenCompactCLHash_CreateFactory",
				     "clGetKernelWorkGroupInfo");;;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGLinearOpenCompactCLHash_DestroyFactory(uintuintHash_Factory *
						      factory, int hashIndex) {;
	clReleaseKernel(factory->emptyKernel[hashIndex]);
	clReleaseProgram(factory->utilProgram[hashIndex]);
	clReleaseKernel(factory->querySingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleNoOverwriteKernel[hashIndex]);;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGLinearOpenCompactCLHash_DestroyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	clReleaseMemObject(table->tableDataBuffer);
	clReleaseContext(table->context);
	clReleaseCommandQueue(table->queue);
	clReleaseProgram(table->utilProgram);
	clReleaseKernel(table->emptyKernel);
	clReleaseProgram(table->program);
	clReleaseKernel(table->querySingleKernel);
	clReleaseKernel(table->insertSingleKernel);
	clReleaseKernel(table->insertSingleNoOverwriteKernel);
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintLCGLinearOpenCompactCLHash_SetupTable(uintuintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((uintuintLCGLinearOpenCompactHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int uintuintLCGLinearOpenCompactCLHash_EmptyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((uintuintLCGLinearOpenCompactHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int uintuintLCGLinearOpenCompactCLHash_QuerySingle(uintuintHash_Table * table,
						   uint key,
						   uint * valueOutput) {
	return uintuintLCGLinearOpenCompactCLHash_Query(table, 1, &key,
							valueOutput);
}
int uintuintLCGLinearOpenCompactCLHash_Query(uintuintHash_Table * table,
					     size_t numKeys, uint * keys,
					     uint * valuesOutput) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numKeys, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_Query",
				     "clCreateBuffer");
	cl_mem valuesOutputBuffer =
	    clCreateBuffer(table->context, CL_MEM_WRITE_ONLY,
			   sizeof(uint) * numKeys, NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_Query",
				     "clCreateBuffer");
	uintuintLCGLinearOpenCompactCLHash_BufferQuery(table, numKeys,
						       keysBuffer,
						       valuesOutputBuffer);
	err =
	    clEnqueueReadBuffer(table->queue, valuesOutputBuffer, CL_TRUE, 0,
				sizeof(uint) * numKeys, valuesOutput, 0, NULL,
				NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_Query",
				     "clEnqueueReadBuffer");
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesOutputBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGLinearOpenCompactCLHash_BufferQuery(uintuintHash_Table * table,
						   size_t numKeys,
						   cl_mem keysBuffer,
						   cl_mem valuesOutputBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->querySingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 1, sizeof(unsigned int),
			   &numKeys);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 3, sizeof(cl_mem),
			   &valuesOutputBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numKeys, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->querySingleKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferQuery",
				     "clEnqueueNDRangeKernel");
	clFinish(table->queue);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGLinearOpenCompactCLHash_InsertSingle(uintuintHash_Table * table,
						    uint key, uint value) {
	return uintuintLCGLinearOpenCompactCLHash_Insert(table, 1, &key,
							 &value);
}
int uintuintLCGLinearOpenCompactCLHash_Insert(uintuintHash_Table * table,
					      size_t numEntries, uint * keys,
					      uint * values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_Insert",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_Insert",
				     "clCreateBuffer");
	uintuintLCGLinearOpenCompactCLHash_BufferInsert(table, numEntries,
							keysBuffer,
							valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGLinearOpenCompactCLHash_BufferInsert(uintuintHash_Table * table,
						    size_t numEntries,
						    cl_mem keysBuffer,
						    cl_mem valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 1, sizeof(unsigned int),
			   &numEntries);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 3, sizeof(cl_mem),
			   &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->insertSingleKernel, 1,
				   0, &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err, NULL, "clEnqueueNDRangeKernel");
	return (0);
}
int
uintuintLCGLinearOpenCompactCLHash_InsertSingleNoOverwrite(uintuintHash_Table *
							   table, uint key,
							   uint value) {
	return uintuintLCGLinearOpenCompactCLHash_InsertNoOverwrite(table, 1,
								    &key,
								    &value);
}
int uintuintLCGLinearOpenCompactCLHash_InsertNoOverwrite(uintuintHash_Table *
							 table,
							 size_t numEntries,
							 uint * keys,
							 uint * values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	uintuintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite(table,
								   numEntries,
								   keysBuffer,
								   valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int
uintuintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite(uintuintHash_Table *
							   table,
							   size_t numEntries,
							   cl_mem keysBuffer,
							   cl_mem valuesBuffer) 
{
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 0,
			   sizeof(cl_mem), &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 1,
			   sizeof(unsigned int), &numEntries);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite",
				     "ClSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 2,
			   sizeof(cl_mem), &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 3,
			   sizeof(cl_mem), &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue,
				   table->insertSingleNoOverwriteKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGLinearOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clEnqueueNDRangeKernel");
	return (0);
}

typedef struct uintuintLCGLinearOpenCompactOpenMPHash_TableData {
	int hashID;
	unsigned int numBuckets;
	uintuintHash_CompressLCGData compressFuncData;
} uintuintLCGLinearOpenCompactOpenMPHash_TableData;
typedef struct uintuintLCGLinearOpenCompactOpenMPHash_Bucket {
	uint key;
	uint value;
} uintuintLCGLinearOpenCompactOpenMPHash_Bucket;
uintuintHash_Table
    *uintuintLCGLinearOpenCompactOpenMPHash_CreateTable(uintuintHash_Factory *
							factory, int hashIndex,
							size_t keyRange,
							size_t numEntries,
							float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc =
	    &uintuintLCGLinearOpenCompactOpenMPHash_DestroyTable;
	table->setupFunc = &uintuintLCGLinearOpenCompactOpenMPHash_SetupTable;
	table->emptyFunc = &uintuintLCGLinearOpenCompactOpenMPHash_EmptyTable;
	table->queryFunc = &uintuintLCGLinearOpenCompactOpenMPHash_Query;
	table->querySingleFunc =
	    &uintuintLCGLinearOpenCompactOpenMPHash_QuerySingle;
	table->insertFunc = &uintuintLCGLinearOpenCompactOpenMPHash_Insert;
	table->insertSingleFunc =
	    &uintuintLCGLinearOpenCompactOpenMPHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintLCGLinearOpenCompactOpenMPHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintLCGLinearOpenCompactOpenMPHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData));
	((uintuintLCGLinearOpenCompactOpenMPHash_TableData *) table->
	 tableData)->hashID = LCG_LINEAR_OPEN_COMPACT_OPENMP_HASH_ID;
	((uintuintLCGLinearOpenCompactOpenMPHash_TableData *) table->
	 tableData)->numBuckets =
(unsigned int)((double)numEntries / loadFactor);
	((uintuintLCGLinearOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.a = HASH_LCG_A;
	((uintuintLCGLinearOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.c = HASH_LCG_C;
	((uintuintLCGLinearOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.m = HASH_LCG_M;
	((uintuintLCGLinearOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.n =
((uintuintLCGLinearOpenCompactOpenMPHash_TableData *) table->tableData)->numBuckets;
	char *tempHashData =
	    (char *)
	    malloc(sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData) +
		   ((uintuintLCGLinearOpenCompactOpenMPHash_TableData *) table->
		    tableData)->numBuckets *
		   sizeof(uintuintLCGLinearOpenCompactOpenMPHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int uintuintLCGLinearOpenCompactOpenMPHash_CreateFactory(uintuintHash_Factory *
							 factory,
							 int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintLCGLinearOpenCompactOpenMPHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintLCGLinearOpenCompactOpenMPHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGLinearOpenCompactOpenMPHash_DestroyFactory(uintuintHash_Factory *
							  factory,
							  int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGLinearOpenCompactOpenMPHash_DestroyTable(uintuintHash_Table *
							table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintLCGLinearOpenCompactOpenMPHash_SetupTable(uintuintHash_Table *
						      table) {
	int exitCode = 0;
	uintuintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactOpenMPHash_Bucket *) & table->
	    tableData[sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData)];
	if (uintuintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
#pragma omp parallel for
		for (uint index = 0;
		     index <
		     ((uintuintLCGLinearOpenCompactOpenMPHash_TableData *)
		      table->tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
		}
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintLCGLinearOpenCompactOpenMPHash_EmptyTable(uintuintHash_Table *
						      table) {
	int exitCode = 0;
	uintuintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactOpenMPHash_Bucket *) & table->
	    tableData[sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData)];
#pragma omp parallel for
	for (uint index = 0;
	     index <
	     ((uintuintLCGLinearOpenCompactOpenMPHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintLCGLinearOpenCompactOpenMPHash_InnerQuerySingle(char *tableData,
							    uint key,
							    uint *
							    valueOutput) {
	uintuintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
	    (uintuintLCGLinearOpenCompactOpenMPHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((uintuintLCGLinearOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		uint old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
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
int uintuintLCGLinearOpenCompactOpenMPHash_InnerQuery(char *tableData,
						      unsigned int numKeys,
						      uint * keys,
						      uint * valuesOutput) {
	uintuintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData)];
	uint key;
	uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		uintuintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
		    (uintuintLCGLinearOpenCompactOpenMPHash_TableData *)
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
			     ((uintuintLCGLinearOpenCompactOpenMPHash_TableData
			       *) tableData)->numBuckets);
			uint old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
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
int uintuintLCGLinearOpenCompactOpenMPHash_InnerInsertSingle(char *tableData,
							     uint key,
							     uint value) {
	uintuintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
	    (uintuintLCGLinearOpenCompactOpenMPHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((uintuintLCGLinearOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		uint old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
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
int uintuintLCGLinearOpenCompactOpenMPHash_InnerInsert(char *tableData,
						       unsigned int numEntries,
						       uint * keys,
						       uint * values) {
	uintuintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		uintuintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
		    (uintuintLCGLinearOpenCompactOpenMPHash_TableData *)
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
			     ((uintuintLCGLinearOpenCompactOpenMPHash_TableData
			       *) tableData)->numBuckets);
			uint old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
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
int uintuintLCGLinearOpenCompactOpenMPHash_InnerInsertSingleNoOverwrite(char
									*tableData,
									uint
									key,
									uint
									value) {
	uintuintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
	    (uintuintLCGLinearOpenCompactOpenMPHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration +
		      c) %
		     ((uintuintLCGLinearOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		uint old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
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
int uintuintLCGLinearOpenCompactOpenMPHash_InnerInsertNoOverwrite(char
								  *tableData,
								  unsigned int
								  numEntries,
								  uint * keys,
								  uint *
								  values) {
	uintuintLCGLinearOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGLinearOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof(uintuintLCGLinearOpenCompactOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		uintuintLCGLinearOpenCompactOpenMPHash_TableData *mytableData =
		    (uintuintLCGLinearOpenCompactOpenMPHash_TableData *)
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
			     ((uintuintLCGLinearOpenCompactOpenMPHash_TableData
			       *) tableData)->numBuckets);
			uint old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
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
int uintuintLCGLinearOpenCompactOpenMPHash_QuerySingle(uintuintHash_Table *
						       table, uint key,
						       uint * valueOutput) {
	return uintuintLCGLinearOpenCompactOpenMPHash_InnerQuerySingle(table->
								       tableData,
								       key,
								       valueOutput);
}
int uintuintLCGLinearOpenCompactOpenMPHash_Query(uintuintHash_Table * table,
						 size_t numKeys, uint * keys,
						 uint * valuesOutput) {
	return uintuintLCGLinearOpenCompactOpenMPHash_InnerQuery(table->
								 tableData,
								 numKeys, keys,
								 valuesOutput);
}
int uintuintLCGLinearOpenCompactOpenMPHash_InsertSingle(uintuintHash_Table *
							table, uint key,
							uint value) {
	return uintuintLCGLinearOpenCompactOpenMPHash_InnerInsertSingle(table->
									tableData,
									key,
									value);
}
int uintuintLCGLinearOpenCompactOpenMPHash_Insert(uintuintHash_Table * table,
						  size_t numEntries,
						  uint * keys, uint * values) {
	return uintuintLCGLinearOpenCompactOpenMPHash_InnerInsert(table->
								  tableData,
								  numEntries,
								  keys, values);
}
int
uintuintLCGLinearOpenCompactOpenMPHash_InsertSingleNoOverwrite
(uintuintHash_Table * table, uint key, uint value) {
	return
	    uintuintLCGLinearOpenCompactOpenMPHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int uintuintLCGLinearOpenCompactOpenMPHash_InsertNoOverwrite(uintuintHash_Table
							     * table,
							     size_t numEntries,
							     uint * keys,
							     uint * values) {
	return
	    uintuintLCGLinearOpenCompactOpenMPHash_InnerInsertNoOverwrite
	    (table->tableData, numEntries, keys, values);
}

typedef struct uintuintLCGQuadraticOpenCompactHash_TableData {
	int hashID;
	unsigned int numBuckets;
	uintuintHash_CompressLCGData compressFuncData;
} uintuintLCGQuadraticOpenCompactHash_TableData;
typedef struct uintuintLCGQuadraticOpenCompactHash_Bucket {
	uint key;
	uint value;
} uintuintLCGQuadraticOpenCompactHash_Bucket;
uintuintHash_Table
    *uintuintLCGQuadraticOpenCompactHash_CreateTable(uintuintHash_Factory *
						     factory, int hashIndex,
						     size_t keyRange,
						     size_t numEntries,
						     float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc = &uintuintLCGQuadraticOpenCompactHash_DestroyTable;
	table->setupFunc = &uintuintLCGQuadraticOpenCompactHash_SetupTable;
	table->emptyFunc = &uintuintLCGQuadraticOpenCompactHash_EmptyTable;
	table->queryFunc = &uintuintLCGQuadraticOpenCompactHash_Query;
	table->querySingleFunc =
	    &uintuintLCGQuadraticOpenCompactHash_QuerySingle;
	table->insertFunc = &uintuintLCGQuadraticOpenCompactHash_Insert;
	table->insertSingleFunc =
	    &uintuintLCGQuadraticOpenCompactHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintLCGQuadraticOpenCompactHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintLCGQuadraticOpenCompactHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(uintuintLCGQuadraticOpenCompactHash_TableData));
	((uintuintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    hashID = LCG_QUADRATIC_OPEN_COMPACT_HASH_ID;
	((uintuintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    numBuckets = (unsigned int)((double)numEntries / loadFactor);
	((uintuintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.a = HASH_LCG_A;
	((uintuintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.c = HASH_LCG_C;
	((uintuintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.m = HASH_LCG_M;
	((uintuintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    compressFuncData.n =
	    ((uintuintLCGQuadraticOpenCompactHash_TableData *) table->
	     tableData)->numBuckets;
	((uintuintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->
	    numBuckets =
	    largestProthPrimeUnder(((uintuintLCGQuadraticOpenCompactHash_TableData *) table->tableData)->numBuckets);
	char *tempHashData =
	    (char *)malloc(sizeof(uintuintLCGQuadraticOpenCompactHash_TableData)
			   +
			   ((uintuintLCGQuadraticOpenCompactHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(uintuintLCGQuadraticOpenCompactHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintLCGQuadraticOpenCompactHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int uintuintLCGQuadraticOpenCompactHash_CreateFactory(uintuintHash_Factory *
						      factory, int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintLCGQuadraticOpenCompactHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintLCGQuadraticOpenCompactHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGQuadraticOpenCompactHash_DestroyFactory(uintuintHash_Factory *
						       factory,
						       int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGQuadraticOpenCompactHash_DestroyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintLCGQuadraticOpenCompactHash_SetupTable(uintuintHash_Table * table) {
	int exitCode = 0;
	uintuintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactHash_Bucket *) & table->
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactHash_TableData)];
	if (uintuintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
		for (uint index = 0;
		     index <
		     ((uintuintLCGQuadraticOpenCompactHash_TableData *) table->
		      tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
		}
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintLCGQuadraticOpenCompactHash_EmptyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	uintuintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactHash_Bucket *) & table->
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactHash_TableData)];
	for (uint index = 0;
	     index <
	     ((uintuintLCGQuadraticOpenCompactHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintLCGQuadraticOpenCompactHash_InnerQuerySingle(char *tableData,
							 uint key,
							 uint * valueOutput) {
	uintuintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGQuadraticOpenCompactHash_TableData *mytableData =
	    (uintuintLCGQuadraticOpenCompactHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((uintuintLCGQuadraticOpenCompactHash_TableData *)
		      tableData)->numBuckets);
		if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((uintuintLCGQuadraticOpenCompactHash_TableData *)
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
int uintuintLCGQuadraticOpenCompactHash_InnerQuery(char *tableData,
						   unsigned int numKeys,
						   uint * keys,
						   uint * valuesOutput) {
	uintuintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactHash_TableData)];
	uint key;
	uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		uintuintLCGQuadraticOpenCompactHash_TableData *mytableData =
		    (uintuintLCGQuadraticOpenCompactHash_TableData *) tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((uintuintLCGQuadraticOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if ((buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((uintuintLCGQuadraticOpenCompactHash_TableData
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
int uintuintLCGQuadraticOpenCompactHash_InnerInsertSingle(char *tableData,
							  uint key,
							  uint value) {
	uintuintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGQuadraticOpenCompactHash_TableData *mytableData =
	    (uintuintLCGQuadraticOpenCompactHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((uintuintLCGQuadraticOpenCompactHash_TableData *)
		      tableData)->numBuckets);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((uintuintLCGQuadraticOpenCompactHash_TableData *)
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
int uintuintLCGQuadraticOpenCompactHash_InnerInsert(char *tableData,
						    unsigned int numEntries,
						    uint * keys,
						    uint * values) {
	uintuintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		uintuintLCGQuadraticOpenCompactHash_TableData *mytableData =
		    (uintuintLCGQuadraticOpenCompactHash_TableData *) tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((uintuintLCGQuadraticOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if (((buckets[index].key ==
			      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
							   key,
							   HASH_BUCKET_STATUS_EMPTY)
			     : buckets[index].key) ==
			    HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((uintuintLCGQuadraticOpenCompactHash_TableData
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
int uintuintLCGQuadraticOpenCompactHash_InnerInsertSingleNoOverwrite(char
								     *tableData,
								     uint key,
								     uint
								     value) {
	uintuintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGQuadraticOpenCompactHash_TableData *mytableData =
	    (uintuintLCGQuadraticOpenCompactHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((uintuintLCGQuadraticOpenCompactHash_TableData *)
		      tableData)->numBuckets);
		if (((buckets[index].key ==
		      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
						   key,
						   HASH_BUCKET_STATUS_EMPTY) :
		     buckets[index].key) == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (key == buckets[index].key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((uintuintLCGQuadraticOpenCompactHash_TableData *)
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
int uintuintLCGQuadraticOpenCompactHash_InnerInsertNoOverwrite(char *tableData,
							       unsigned int
							       numEntries,
							       uint * keys,
							       uint * values) {
	uintuintLCGQuadraticOpenCompactHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactHash_Bucket *) &
	    tableData[sizeof(uintuintLCGQuadraticOpenCompactHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;;
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		uintuintLCGQuadraticOpenCompactHash_TableData *mytableData =
		    (uintuintLCGQuadraticOpenCompactHash_TableData *) tableData;
		uintuintHash_CompressLCGData compressFuncData =
		    mytableData->compressFuncData;
		unsigned int c =
		    uintuintHash_CompressLCG(compressFuncData, key);
		unsigned long int iteration = 0;
		for (;;) {
			index =
			    ((1 * iteration * iteration + 0 * iteration +
			      c) %
			     ((uintuintLCGQuadraticOpenCompactHash_TableData *)
			      tableData)->numBuckets);
			if (((buckets[index].key ==
			      HASH_BUCKET_STATUS_EMPTY) ? (buckets[index].key =
							   key,
							   HASH_BUCKET_STATUS_EMPTY)
			     : buckets[index].key) ==
			    HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (key == buckets[index].key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((uintuintLCGQuadraticOpenCompactHash_TableData
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
int uintuintLCGQuadraticOpenCompactHash_QuerySingle(uintuintHash_Table * table,
						    uint key,
						    uint * valueOutput) {
	return uintuintLCGQuadraticOpenCompactHash_InnerQuerySingle(table->
								    tableData,
								    key,
								    valueOutput);
}
int uintuintLCGQuadraticOpenCompactHash_Query(uintuintHash_Table * table,
					      size_t numKeys, uint * keys,
					      uint * valuesOutput) {
	return uintuintLCGQuadraticOpenCompactHash_InnerQuery(table->tableData,
							      numKeys, keys,
							      valuesOutput);
}
int uintuintLCGQuadraticOpenCompactHash_InsertSingle(uintuintHash_Table * table,
						     uint key, uint value) {
	return uintuintLCGQuadraticOpenCompactHash_InnerInsertSingle(table->
								     tableData,
								     key,
								     value);
}
int uintuintLCGQuadraticOpenCompactHash_Insert(uintuintHash_Table * table,
					       size_t numEntries, uint * keys,
					       uint * values) {
	return uintuintLCGQuadraticOpenCompactHash_InnerInsert(table->tableData,
							       numEntries, keys,
							       values);
}
int
uintuintLCGQuadraticOpenCompactHash_InsertSingleNoOverwrite(uintuintHash_Table *
							    table, uint key,
							    uint value) {
	return
	    uintuintLCGQuadraticOpenCompactHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int uintuintLCGQuadraticOpenCompactHash_InsertNoOverwrite(uintuintHash_Table *
							  table,
							  size_t numEntries,
							  uint * keys,
							  uint * values) {
	return
	    uintuintLCGQuadraticOpenCompactHash_InnerInsertNoOverwrite(table->
								       tableData,
								       numEntries,
								       keys,
								       values);
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
uintuintHash_Table
    *uintuintLCGQuadraticOpenCompactCLHash_CreateTable(uintuintHash_Factory *
						       factory, int hashIndex,
						       size_t keyRange,
						       size_t numEntries,
						       float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc =
	    &uintuintLCGQuadraticOpenCompactCLHash_DestroyTable;
	table->setupFunc = &uintuintLCGQuadraticOpenCompactCLHash_SetupTable;
	table->emptyFunc = &uintuintLCGQuadraticOpenCompactCLHash_EmptyTable;
	table->queryFunc = &uintuintLCGQuadraticOpenCompactCLHash_Query;
	table->querySingleFunc =
	    &uintuintLCGQuadraticOpenCompactCLHash_QuerySingle;
	table->insertFunc = &uintuintLCGQuadraticOpenCompactCLHash_Insert;
	table->insertSingleFunc =
	    &uintuintLCGQuadraticOpenCompactCLHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(uintuintLCGQuadraticOpenCompactCLHash_TableData));
	((uintuintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    hashID = LCG_QUADRATIC_OPEN_COMPACT_CL_HASH_ID;
	table->context = factory->context;
	table->queue = factory->queue;
	table->program = factory->program;
	table->localWorkSize = factory->localWorkSize;
	table->utilProgram = factory->utilProgram[hashIndex];
	table->emptyKernel = factory->emptyKernel[hashIndex];
	table->emptyKernelLocalWorkSize =
	    factory->emptyKernelLocalWorkSize[hashIndex];
	table->querySingleKernel = factory->querySingleKernel[hashIndex];
	table->insertSingleKernel = factory->insertSingleKernel[hashIndex];
	table->insertSingleNoOverwriteKernel =
	    factory->insertSingleNoOverwriteKernel[hashIndex];
	clRetainContext(table->context);
	clRetainCommandQueue(table->queue);
	clRetainProgram(table->program);
	clRetainProgram(table->utilProgram);
	clRetainKernel(table->emptyKernel);
	clRetainKernel(table->querySingleKernel);
	clRetainKernel(table->insertSingleKernel);
	clRetainKernel(table->insertSingleNoOverwriteKernel);;
	((uintuintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    numBuckets = (unsigned int)((double)numEntries / loadFactor);
	((uintuintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.a = HASH_LCG_A;
	((uintuintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.c = HASH_LCG_C;
	((uintuintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.m = HASH_LCG_M;
	((uintuintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    compressFuncData.n =
	    ((uintuintLCGQuadraticOpenCompactCLHash_TableData *) table->
	     tableData)->numBuckets;
	((uintuintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->
	    numBuckets =
	    largestProthPrimeUnder(((uintuintLCGQuadraticOpenCompactCLHash_TableData *) table->tableData)->numBuckets);
	char *tempHashData =
	    (char *)
	    malloc(sizeof(uintuintLCGQuadraticOpenCompactCLHash_TableData) +
		   ((uintuintLCGQuadraticOpenCompactCLHash_TableData *) table->
		    tableData)->numBuckets *
		   sizeof(uintuintLCGQuadraticOpenCompactCLHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintLCGQuadraticOpenCompactCLHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	cl_int err;
	table->tableDataBuffer =
	    clCreateBuffer(table->context, CL_MEM_READ_WRITE,
			   sizeof(uintuintLCGQuadraticOpenCompactHash_TableData)
			   +
			   ((uintuintLCGQuadraticOpenCompactHash_TableData *)
			    table->tableData)->numBuckets *
			   sizeof(uintuintLCGQuadraticOpenCompactHash_Bucket),
			   NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_InitTable",
				     "clCreateBuffer");
	err =
	    clEnqueueWriteBuffer(table->queue, table->tableDataBuffer, CL_TRUE,
				 0,
				 sizeof
				 (uintuintLCGQuadraticOpenCompactHash_TableData),
				 table->tableData, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_InitTable",
				     "clEnqueueWriteBuffer");
	return table;
}
int uintuintLCGQuadraticOpenCompactCLHash_CreateFactory(uintuintHash_Factory *
							factory,
							int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintLCGQuadraticOpenCompactCLHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintLCGQuadraticOpenCompactCLHash_DestroyFactory;
	cl_int error;
	cl_device_id device;
	error =
	    clGetContextInfo(factory->context, CL_CONTEXT_DEVICES,
			     sizeof(device), &device, NULL);
	CLHash_Utilities_HandleError(error, "uintuintHash_CreateFactory",
				     "clGetContextInfo");
	factory->querySingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintLCGQuadraticOpenCompactCLHash_RangeQuerySingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintLCGQuadraticOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintLCGQuadraticOpenCompactCLHash_RangeInsertSingle",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintLCGQuadraticOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->insertSingleNoOverwriteKernel[hashIndex] =
	    clCreateKernel(factory->program,
			   "uintuintLCGQuadraticOpenCompactCLHash_RangeInsertSingleNoOverwrite",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintLCGQuadraticOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	factory->utilProgram[hashIndex] =
	    CLHash_Utilities_BuildProgramString(factory->context, device,
						"static inline unsigned int uintuintHash_CompressIdentity(char data, int hashCode){ return hashCode; } typedef struct uintuintHash_CompressLCGData{ long unsigned int a; long unsigned int c; unsigned int m; unsigned int n; }uintuintHash_CompressLCGData; static inline unsigned int uintuintHash_CompressLCG(uintuintHash_CompressLCGData compressLCGData, int hashCode){ return ((compressLCGData.a * hashCode + compressLCGData.c) % compressLCGData.m) % compressLCGData.n; } typedef struct uintuintLCGQuadraticOpenCompactCLHash_TableData{ int hashID; unsigned int numBuckets; uintuintHash_CompressLCGData compressFuncData; }uintuintLCGQuadraticOpenCompactCLHash_TableData; typedef struct uintuintLCGQuadraticOpenCompactCLHash_Bucket{ uint key; uint value; }uintuintLCGQuadraticOpenCompactCLHash_Bucket; __kernel void uintuintLCGQuadraticOpenCompactCLHash_Empty(__global char *tableData){ int index = get_global_id(0); if(index >= ((__global uintuintLCGQuadraticOpenCompactCLHash_TableData*)tableData)->numBuckets){ return; } __global uintuintLCGQuadraticOpenCompactCLHash_Bucket *buckets = (__global uintuintLCGQuadraticOpenCompactCLHash_Bucket*)&tableData[sizeof(uintuintLCGQuadraticOpenCompactCLHash_TableData)]; buckets[index].key = -1;/*HASH_BUCKET_STATUS_EMPTY*/ }");
	factory->emptyKernel[hashIndex] =
	    clCreateKernel(factory->utilProgram[hashIndex],
			   "uintuintLCGQuadraticOpenCompactCLHash_Empty",
			   &error);
	CLHash_Utilities_HandleError(error,
				     "uintuintLCGQuadraticOpenCompactCLHash_CreateFactory",
				     "clCreateKernel");
	error =
	    clGetKernelWorkGroupInfo(factory->emptyKernel[hashIndex], device,
				     CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
				     &factory->
				     emptyKernelLocalWorkSize[hashIndex], NULL);
	CLHash_Utilities_HandleError(error,
				     "uintuintLCGQuadraticOpenCompactCLHash_CreateFactory",
				     "clGetKernelWorkGroupInfo");;;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGQuadraticOpenCompactCLHash_DestroyFactory(uintuintHash_Factory *
							 factory,
							 int hashIndex) {;
	clReleaseKernel(factory->emptyKernel[hashIndex]);
	clReleaseProgram(factory->utilProgram[hashIndex]);
	clReleaseKernel(factory->querySingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleKernel[hashIndex]);
	clReleaseKernel(factory->insertSingleNoOverwriteKernel[hashIndex]);;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGQuadraticOpenCompactCLHash_DestroyTable(uintuintHash_Table *
						       table) {
	int exitCode = 0;
	clReleaseMemObject(table->tableDataBuffer);
	clReleaseContext(table->context);
	clReleaseCommandQueue(table->queue);
	clReleaseProgram(table->utilProgram);
	clReleaseKernel(table->emptyKernel);
	clReleaseProgram(table->program);
	clReleaseKernel(table->querySingleKernel);
	clReleaseKernel(table->insertSingleKernel);
	clReleaseKernel(table->insertSingleNoOverwriteKernel);
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintLCGQuadraticOpenCompactCLHash_SetupTable(uintuintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((uintuintLCGQuadraticOpenCompactHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int uintuintLCGQuadraticOpenCompactCLHash_EmptyTable(uintuintHash_Table * table) {
	int exitCode = 0;
	cl_int err;
	err =
	    clSetKernelArg(table->emptyKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_EmptyTable",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(((uintuintLCGQuadraticOpenCompactHash_TableData *)
			      table->tableData)->numBuckets,
			     table->emptyKernelLocalWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->emptyKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->
				   emptyKernelLocalWorkSize, 0, NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_EmptyTable",
				     "clEnqueueNDRangeKernel");
	exitCode = HASH_EXIT_CODE_NORMAL;;
	return exitCode;
}
int uintuintLCGQuadraticOpenCompactCLHash_QuerySingle(uintuintHash_Table *
						      table, uint key,
						      uint * valueOutput) {
	return uintuintLCGQuadraticOpenCompactCLHash_Query(table, 1, &key,
							   valueOutput);
}
int uintuintLCGQuadraticOpenCompactCLHash_Query(uintuintHash_Table * table,
						size_t numKeys, uint * keys,
						uint * valuesOutput) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numKeys, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_Query",
				     "clCreateBuffer");
	cl_mem valuesOutputBuffer =
	    clCreateBuffer(table->context, CL_MEM_WRITE_ONLY,
			   sizeof(uint) * numKeys, NULL, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_Query",
				     "clCreateBuffer");
	uintuintLCGQuadraticOpenCompactCLHash_BufferQuery(table, numKeys,
							  keysBuffer,
							  valuesOutputBuffer);
	err =
	    clEnqueueReadBuffer(table->queue, valuesOutputBuffer, CL_TRUE, 0,
				sizeof(uint) * numKeys, valuesOutput, 0, NULL,
				NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_Query",
				     "clEnqueueReadBuffer");
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesOutputBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGQuadraticOpenCompactCLHash_BufferQuery(uintuintHash_Table *
						      table, size_t numKeys,
						      cl_mem keysBuffer,
						      cl_mem
						      valuesOutputBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->querySingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 1, sizeof(unsigned int),
			   &numKeys);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->querySingleKernel, 3, sizeof(cl_mem),
			   &valuesOutputBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferQuery",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numKeys, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->querySingleKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferQuery",
				     "clEnqueueNDRangeKernel");
	clFinish(table->queue);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGQuadraticOpenCompactCLHash_InsertSingle(uintuintHash_Table *
						       table, uint key,
						       uint value) {
	return uintuintLCGQuadraticOpenCompactCLHash_Insert(table, 1, &key,
							    &value);
}
int uintuintLCGQuadraticOpenCompactCLHash_Insert(uintuintHash_Table * table,
						 size_t numEntries, uint * keys,
						 uint * values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_Insert",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_Insert",
				     "clCreateBuffer");
	uintuintLCGQuadraticOpenCompactCLHash_BufferInsert(table, numEntries,
							   keysBuffer,
							   valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGQuadraticOpenCompactCLHash_BufferInsert(uintuintHash_Table *
						       table, size_t numEntries,
						       cl_mem keysBuffer,
						       cl_mem valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleKernel, 0, sizeof(cl_mem),
			   &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 1, sizeof(unsigned int),
			   &numEntries);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 2, sizeof(cl_mem),
			   &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleKernel, 3, sizeof(cl_mem),
			   &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferInsert",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue, table->insertSingleKernel, 1,
				   0, &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err, NULL, "clEnqueueNDRangeKernel");
	return (0);
}
int
uintuintLCGQuadraticOpenCompactCLHash_InsertSingleNoOverwrite(uintuintHash_Table
							      * table, uint key,
							      uint value) {
	return uintuintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(table, 1,
								       &key,
								       &value);
}
int uintuintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite(uintuintHash_Table *
							    table,
							    size_t numEntries,
							    uint * keys,
							    uint * values) {
	cl_int err;
	cl_mem keysBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, keys, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	cl_mem valuesBuffer =
	    clCreateBuffer(table->context,
			   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			   sizeof(uint) * numEntries, values, &err);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_InsertNoOverwrite",
				     "clCreateBuffer");
	uintuintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite(table,
								      numEntries,
								      keysBuffer,
								      valuesBuffer);
	clReleaseMemObject(keysBuffer);
	clReleaseMemObject(valuesBuffer);
	return HASH_EXIT_CODE_NORMAL;
}
int
uintuintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite(uintuintHash_Table
							      * table,
							      size_t numEntries,
							      cl_mem keysBuffer,
							      cl_mem
							      valuesBuffer) {
	cl_int err;
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 0,
			   sizeof(cl_mem), &table->tableDataBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 1,
			   sizeof(unsigned int), &numEntries);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite",
				     "ClSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 2,
			   sizeof(cl_mem), &keysBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	err =
	    clSetKernelArg(table->insertSingleNoOverwriteKernel, 3,
			   sizeof(cl_mem), &valuesBuffer);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clSetKernelArg");
	const size_t groupWorkSize =
	    roundUpToNearest(numEntries, table->localWorkSize);
	err =
	    clEnqueueNDRangeKernel(table->queue,
				   table->insertSingleNoOverwriteKernel, 1, 0,
				   &groupWorkSize,
				   (const size_t *)&table->localWorkSize, 0,
				   NULL, NULL);
	CLHash_Utilities_HandleError(err,
				     "uintuintLCGQuadraticOpenCompactCLHash_BufferInsertNoOverwrite",
				     "clEnqueueNDRangeKernel");
	return (0);
}

typedef struct uintuintLCGQuadraticOpenCompactOpenMPHash_TableData {
	int hashID;
	unsigned int numBuckets;
	uintuintHash_CompressLCGData compressFuncData;
} uintuintLCGQuadraticOpenCompactOpenMPHash_TableData;
typedef struct uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket {
	uint key;
	uint value;
} uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket;
uintuintHash_Table
    *uintuintLCGQuadraticOpenCompactOpenMPHash_CreateTable(uintuintHash_Factory
							   * factory,
							   int hashIndex,
							   size_t keyRange,
							   size_t numEntries,
							   float loadFactor) {
	uintuintHash_Table *table =
	    (uintuintHash_Table *) malloc(sizeof(uintuintHash_Table));
	table->destroyFunc =
	    &uintuintLCGQuadraticOpenCompactOpenMPHash_DestroyTable;
	table->setupFunc =
	    &uintuintLCGQuadraticOpenCompactOpenMPHash_SetupTable;
	table->emptyFunc =
	    &uintuintLCGQuadraticOpenCompactOpenMPHash_EmptyTable;
	table->queryFunc = &uintuintLCGQuadraticOpenCompactOpenMPHash_Query;
	table->querySingleFunc =
	    &uintuintLCGQuadraticOpenCompactOpenMPHash_QuerySingle;
	table->insertFunc = &uintuintLCGQuadraticOpenCompactOpenMPHash_Insert;
	table->insertSingleFunc =
	    &uintuintLCGQuadraticOpenCompactOpenMPHash_InsertSingle;
	table->insertNoOverwriteFunc =
	    &uintuintLCGQuadraticOpenCompactOpenMPHash_InsertNoOverwrite;
	table->insertSingleNoOverwriteFunc =
	    &uintuintLCGQuadraticOpenCompactOpenMPHash_InsertSingleNoOverwrite;
	table->tableData =
	    (char *)
	    malloc(sizeof(uintuintLCGQuadraticOpenCompactOpenMPHash_TableData));
	((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->hashID = LCG_QUADRATIC_OPEN_COMPACT_OPENMP_HASH_ID;
	((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->numBuckets =
(unsigned int)((double)numEntries / loadFactor);
	((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.a = HASH_LCG_A;
	((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.c = HASH_LCG_C;
	((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.m = HASH_LCG_M;
	((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->compressFuncData.n =
((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->tableData)->numBuckets;
	((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	 tableData)->numBuckets =
largestProthPrimeUnder(((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->tableData)->numBuckets);
	char *tempHashData =
	    (char *)
	    malloc(sizeof(uintuintLCGQuadraticOpenCompactOpenMPHash_TableData) +
		   ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		    table->tableData)->numBuckets *
		   sizeof(uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket));
	memcpy(tempHashData, table->tableData,
	       sizeof(uintuintLCGQuadraticOpenCompactOpenMPHash_TableData));
	free(table->tableData);
	table->tableData = tempHashData;
	return table;
}
int uintuintLCGQuadraticOpenCompactOpenMPHash_CreateFactory(uintuintHash_Factory
							    * factory,
							    int hashIndex) {
	factory->createFunc[hashIndex] =
	    &uintuintLCGQuadraticOpenCompactOpenMPHash_CreateTable;
	factory->destroyFunc[hashIndex] =
	    &uintuintLCGQuadraticOpenCompactOpenMPHash_DestroyFactory;;
	return HASH_EXIT_CODE_NORMAL;
}
int
uintuintLCGQuadraticOpenCompactOpenMPHash_DestroyFactory(uintuintHash_Factory *
							 factory,
							 int hashIndex) {;
	return HASH_EXIT_CODE_NORMAL;
}
int uintuintLCGQuadraticOpenCompactOpenMPHash_DestroyTable(uintuintHash_Table *
							   table) {
	int exitCode = 0;
	free(table->tableData);
	free(table);
	return exitCode;
}
int uintuintLCGQuadraticOpenCompactOpenMPHash_SetupTable(uintuintHash_Table *
							 table) {
	int exitCode = 0;
	uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *) & table->
	    tableData[sizeof
		      (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	if (uintuintHash_GetTableType(table) & ~HASH_SENTINEL_PERFECT_HASHES) {
#pragma omp parallel for
		for (uint index = 0;
		     index <
		     ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		      table->tableData)->numBuckets; index++) {
			buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
		}
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintLCGQuadraticOpenCompactOpenMPHash_EmptyTable(uintuintHash_Table *
							 table) {
	int exitCode = 0;
	uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *) & table->
	    tableData[sizeof
		      (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData)];
#pragma omp parallel for
	for (uint index = 0;
	     index <
	     ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) table->
	      tableData)->numBuckets; index++) {
		buckets[index].key = HASH_BUCKET_STATUS_EMPTY;
	}
	exitCode = HASH_EXIT_CODE_NORMAL;
	return exitCode;
}
int uintuintLCGQuadraticOpenCompactOpenMPHash_InnerQuerySingle(char *tableData,
							       uint key,
							       uint *
							       valueOutput) {
	uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		uint old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData
			   *) tableData)->numBuckets)) {
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
int uintuintLCGQuadraticOpenCompactOpenMPHash_InnerQuery(char *tableData,
							 unsigned int numKeys,
							 uint * keys,
							 uint * valuesOutput) {
	uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	uint key;
	uint *valueOutput;
	uint index;
	int exitCode;
	uint i;
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	for (i = 0; i < numKeys; i++) {
		key = keys[i];
		valueOutput = &valuesOutput[i];
		uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData
		    =
		    (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *)
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
			     ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData)->numBuckets);
			uint old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData)->numBuckets)) {
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
int uintuintLCGQuadraticOpenCompactOpenMPHash_InnerInsertSingle(char *tableData,
								uint key,
								uint value) {
	uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		uint old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData
			   *) tableData)->numBuckets)) {
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
int uintuintLCGQuadraticOpenCompactOpenMPHash_InnerInsert(char *tableData,
							  unsigned int
							  numEntries,
							  uint * keys,
							  uint * values) {
	uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData
		    =
		    (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *)
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
			     ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData)->numBuckets);
			uint old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData)->numBuckets)) {
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
int uintuintLCGQuadraticOpenCompactOpenMPHash_InnerInsertSingleNoOverwrite(char
									   *tableData,
									   uint
									   key,
									   uint
									   value) 
{
	uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	uint index;
	int exitCode;
	uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData;
	uintuintHash_CompressLCGData compressFuncData =
	    mytableData->compressFuncData;
	unsigned int c = uintuintHash_CompressLCG(compressFuncData, key);
	unsigned long int iteration = 0;
	for (;;) {
		index =
		    ((1 * iteration * iteration + 0 * iteration +
		      c) %
		     ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *)
		      tableData)->numBuckets);
		uint old_key =
		    __sync_val_compare_and_swap(&buckets[index].key, -1, key);
		if (old_key == HASH_BUCKET_STATUS_EMPTY) {
			exitCode = HASH_SEARCH_CODE_EMPTY;
			break;
		} else if (old_key == key) {
			exitCode = HASH_SEARCH_CODE_MATCH;
			break;
		} else
		    if ((iteration >
			 ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData
			   *) tableData)->numBuckets)) {
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
int uintuintLCGQuadraticOpenCompactOpenMPHash_InnerInsertNoOverwrite(char
								     *tableData,
								     unsigned
								     int
								     numEntries,
								     uint *
								     keys,
								     uint *
								     values) {
	uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *buckets =
	    (uintuintLCGQuadraticOpenCompactOpenMPHash_Bucket *) &
	    tableData[sizeof
		      (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData)];
	int resultExitCode = HASH_EXIT_CODE_NORMAL;
	uint key;
	uint index;
	int exitCode;
	uint i;
#pragma omp parallel for
	for (i = 0; i < numEntries; i++) {
		key = keys[i];
		uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *mytableData
		    =
		    (uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *)
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
			     ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData)->numBuckets);
			uint old_key =
			    __sync_val_compare_and_swap(&buckets[index].key, -1,
							key);
			if (old_key == HASH_BUCKET_STATUS_EMPTY) {
				exitCode = HASH_SEARCH_CODE_EMPTY;
				break;
			} else if (old_key == key) {
				exitCode = HASH_SEARCH_CODE_MATCH;
				break;
			} else
			    if ((iteration >
				 ((uintuintLCGQuadraticOpenCompactOpenMPHash_TableData *) tableData)->numBuckets)) {
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
int uintuintLCGQuadraticOpenCompactOpenMPHash_QuerySingle(uintuintHash_Table *
							  table, uint key,
							  uint * valueOutput) {
	return
	    uintuintLCGQuadraticOpenCompactOpenMPHash_InnerQuerySingle(table->
								       tableData,
								       key,
								       valueOutput);
}
int uintuintLCGQuadraticOpenCompactOpenMPHash_Query(uintuintHash_Table * table,
						    size_t numKeys, uint * keys,
						    uint * valuesOutput) {
	return uintuintLCGQuadraticOpenCompactOpenMPHash_InnerQuery(table->
								    tableData,
								    numKeys,
								    keys,
								    valuesOutput);
}
int uintuintLCGQuadraticOpenCompactOpenMPHash_InsertSingle(uintuintHash_Table *
							   table, uint key,
							   uint value) {
	return
	    uintuintLCGQuadraticOpenCompactOpenMPHash_InnerInsertSingle(table->
									tableData,
									key,
									value);
}
int uintuintLCGQuadraticOpenCompactOpenMPHash_Insert(uintuintHash_Table * table,
						     size_t numEntries,
						     uint * keys,
						     uint * values) {
	return uintuintLCGQuadraticOpenCompactOpenMPHash_InnerInsert(table->
								     tableData,
								     numEntries,
								     keys,
								     values);
}
int
uintuintLCGQuadraticOpenCompactOpenMPHash_InsertSingleNoOverwrite
(uintuintHash_Table * table, uint key, uint value) {
	return
	    uintuintLCGQuadraticOpenCompactOpenMPHash_InnerInsertSingleNoOverwrite
	    (table->tableData, key, value);
}
int
uintuintLCGQuadraticOpenCompactOpenMPHash_InsertNoOverwrite(uintuintHash_Table *
							    table,
							    size_t numEntries,
							    uint * keys,
							    uint * values) {
	return
	    uintuintLCGQuadraticOpenCompactOpenMPHash_InnerInsertNoOverwrite
	    (table->tableData, numEntries, keys, values);
}
