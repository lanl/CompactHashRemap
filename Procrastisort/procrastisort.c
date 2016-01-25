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
#include <stdbool.h>
//#include <sort.h>



typedef unsigned int uint;
struct timeval timer;
double t1, t2;

void swap_double(double* a, double* b) {
   double c = *a;
   *a = *b;
   *b = c;
}


double* insertionsort(uint length, double *arr) {
	int j;
	for (uint i = 0; i < length; i++) {
		j = i;
		while (j > 0 && arr[j-1] > arr[j]) {
			swap_double(&arr[j], &arr[j-1]);
			j--;
		}
	}
}

double* bubblesort(uint length, double *arr) {
	int n = length;
	bool swapped;
	do {
		swapped = false;
		for (int i = 1; i < n; i++) {
			if (arr[i-1] > arr[i]) {
				swap_double(&arr[i-1], &arr[i]);
				swapped = true;
			}
		}
		n--;
	} while (swapped == true);
}

//int compare (const void * a, const void * b) { return ((int)(*(double*)a - *(double*)b )); }
int compare (const void * a, const void * b)
{
  if ( *(double*)a <  *(double*)b ) return -1;
  if ( *(double*)a >  *(double*)b ) return 1;
  return 0;
}

#define ValType double
#define IS_LESS(v1, v2)  (v1 < v2)
 
void siftDown( ValType *a, int start, int count);
 
#define SWAP(r,s)  do{ValType t=r; r=s; s=t; } while(0)
 
void heapsort( ValType *a, int count)
{
    int start, end;
 
    /* heapify */
    for (start = (count-2)/2; start >=0; start--) {
        siftDown( a, start, count);
    }
 
    for (end=count-1; end > 0; end--) {
        SWAP(a[end],a[0]);
        siftDown(a, 0, end);
    }
}
 
void siftDown( ValType *a, int start, int end)
{
    int root = start;
 
    while ( root*2+1 < end ) {
        int child = 2*root + 1;
        if ((child + 1 < end) && IS_LESS(a[child],a[child+1])) {
            child += 1;
        }
        if (IS_LESS(a[root], a[child])) {
            SWAP( a[child], a[root] );
            root = child;
        }
        else
            return;
    }
}

uint logTwoF (uint numIn) {
	static const char LogTable256[256] = 
	{
	#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
	};

	// 32-bit word to find the log of
	uint r;     // r will be lg(v)
	register uint t, tt; // temporaries
	
	if (tt = numIn >> 16)
	{
		r = (t = tt >> 8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	}
	else 
	{
		r = (t = numIn >> 8) ? 8 + LogTable256[t] : LogTable256[numIn];
	}
	return r;
}

uint logTwoC (uint numIn) {
	static const char LogTable256[256] = 
	{
	#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	    LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	    LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
	};

	// 32-bit word to find the log of
	uint r;     // r will be lg(v)
	register uint t, tt; // temporaries
	
	if (tt = numIn >> 16)
	{
		r = (t = tt >> 8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	}
	else 
	{
		r = (t = numIn >> 8) ? 8 + LogTable256[t] : LogTable256[numIn];
	}
	if (1 << r == numIn) {
		return r;
	} else {
		return r+1;
	}
}

/*This function is the recursive procrastisort. the collisions are not handled in parallel to prevent thread overhead*/
double* procrastisortR(uint length, double *arr) {

    if (length > 1) {
		double min = arr[0];
		double max = arr[0];
		
		for (int i = 1; i < length; i++) {
			if (arr[i] < min) {
				min = arr[i];
			}
			if (arr[i] > max) {
				max = arr[i];
			}
		}
		    
    	double range = max - min;
    
    	if (range > 0) {
    		
    		double *sorted = (double*)malloc(length*sizeof(double));
    	
    		int binIndex;
   		 	int numRec;
    		numRec = 0;
    		double binSize = (range)/(length-1);
    	
    		/*subtract out min*/
    		if (min > 0) {
    			for (uint i = 0; i < length; i++) {
    				arr[i] -= min;
    			}
    		}
    
	    	/*get prefix scan*/
	    	int *num = (int*)malloc(length*sizeof(int));
	    	int *sum = (int*)malloc((length+1)*sizeof(int));
    			
			/*initialise array to 0*/
			memset(num, 0, (length)*sizeof(int));
			
			/*get histogram*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = (int) (arr[i]/binSize);
				num[binIndex]++;
				if (num[binIndex] == 2) {
					numRec++;
				}
	    	}
			
			sum[0] = 0;
	    	for (int i = 1; i < length+1; i++) {
				sum[i] = sum[i-1] + num[i-1];
	    	}
	    	free(num);

	    	/*set all elements of hash array to -1*/
	    	int *hash = (int*)malloc((length)*sizeof(int));
	    	memset(hash, -1, (length)*sizeof(int));

	    	for(uint i = 0; i < length; i++) {
				binIndex = sum[(int) (arr[i]/binSize)];
				while (hash[binIndex] != -1 && binIndex < length) { 
					binIndex++;
				} hash[binIndex] = i;
	    	}
    
	    	if (numRec > 0) {
	    		int colLoc;
	    		int colTally = 0;
		 
		 		//#pragma omp parallel for 
				for (int i = 0; i < numRec; i++) {
					double *colliders;
					int prefixJump = sum[colTally+1]-sum[colTally];
					while (prefixJump <= 1) {
	    					colTally++;
	    					prefixJump = sum[colTally+1]-sum[colTally];
	    			} 
	    			colLoc = colTally;
	    			colTally++; 
	    			
		    		colliders = (double*)malloc((prefixJump)*sizeof(double));
		    		for (int j = 0; j < prefixJump; j++) {
		    			double collider = arr[hash[sum[colLoc] + j]];
		    			colliders[j] = collider;
		    		}
		    		
		    		double* tempCols = procrastisortR(prefixJump, colliders);
		    		
	    			for (int j = 0; j < prefixJump; j++) {
	    				arr[hash[sum[colLoc] + j]] = tempCols[j];
	    			}
	    			
	    			if (tempCols == colliders) {
	    				free(tempCols);
	    				tempCols = NULL;
	    			} else {
	    			free(colliders);
	    			free(tempCols);
	    			}
		    	}
			}	
		
		
			for (int i = 0; i < length; i++) {
	    		sorted[i] = arr[hash[i]] + min;
	    	}
	    	free(sum);
	    	free(hash);
	    	return sorted;
	    } else {
	    	return arr;
    	}
	} else {
		return arr;
	}
}

/*procrastisort that demands min and max for the sake of speed*/
double* procrastisortD(uint length, double *arr, double min, double max) {
   
    if (length > 1) {
    
    	double range = max - min;
    
    	if (range > 0) {
    		
    		double *sorted = (double*)malloc(length*sizeof(double));
    	
    		int binIndex;
   		 	int numRec;
    		numRec = 0;
    		double binSize = (range)/(length-1);
			
    		if (min > 0) {
    			#pragma omp parallel for
    			for (uint i = 0; i < length; i++) {
    				arr[i] -= min;
    			}
    		}
    
	    	int *num = (int*)malloc(length*sizeof(int));
	    	int *sum = (int*)malloc((length+1)*sizeof(int));
	    	
			/*set bin count to zero in each array*/
			memset(num, 0, (length)*sizeof(int));
			
			/*tally number in each bin*/
    		//#pragma omp parallel for
	    	for(uint i = 0; i < length; i++) {
				binIndex = (int) ((arr[i])/binSize);
				num[binIndex]++;
				if (num[binIndex] == 2) {
					numRec++;
				}
	    	}
	    	
			/*prefix scanning*/
			sum[0] = 0;
	    	for (int i = 1; i < length+1; i++) {
				sum[i] = sum[i-1] + num[i-1];
	    	}
	    	free(num);

	    	/*set all elements of hash array to -1*/
	    	int *hash = (int*)malloc((length)*sizeof(int));
	    	memset(hash, -1, (length)*sizeof(int));

			/*hashing using the prefix scan*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = sum[(int) ((arr[i])/binSize)];
				while (hash[binIndex] != -1 && binIndex < length) { 
					binIndex++;
				} hash[binIndex] = i;
	    	}
    		
    		/*handle collisions*/
	    	if (numRec > 0) {
	    		int colTally = 0;
	    		
	    		int *colLocs = (int*)malloc(numRec*sizeof(int));
	    		int *pJumps = (int*)malloc(numRec*sizeof(int));
	    		
	    		/*get collision locations and number of collisions there*/
	    		for (int i = 0; i < numRec; i++) {
					int prefixJump = sum[colTally+1]-sum[colTally];
					while (prefixJump <= 1) {
	    					colTally++;
	    					prefixJump = sum[colTally+1]-sum[colTally];
	    			} 
	    			colLocs[i] = colTally;
	    			pJumps[i] = prefixJump;
	    			colTally++;
	    		}
		 		
		 		/*recursively sort collisions*/
		 		#pragma omp parallel for
				for (int i = 0; i < numRec; i++) {
					double *colliders;
					int prefixJump = pJumps[i];
					int colLoc = colLocs[i];
	    			
		    		colliders = (double*)malloc(prefixJump*sizeof(double));
		    		for (int j = 0; j < prefixJump; j++) {
		    			double collider = arr[hash[sum[colLoc] + j]];
		    			colliders[j] = collider;
		    		}
		    		
		    		double* tempCols = procrastisortR(prefixJump, colliders);
		    		
	    			for (int j = 0; j < prefixJump; j++) {
	    				arr[hash[sum[colLoc] + j]] = tempCols[j];
	    			}
	    			
	    			if (tempCols == colliders) {
	    				free(tempCols);
	    				tempCols = NULL;
	    			} else {
	    			free(colliders);
	    			free(tempCols);
	    			}
		    	}
		    	free(colLocs);
		    	free(pJumps);
			}	
			
			/*finish using hash keys*/
			if (min > 0) {
				#pragma omp parallel for
				for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]] + min;
	    			arr[hash[i]] += min;
	    		
	    		}
	    	} else {
	    		#pragma omp parallel for
	    		for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]];
	    		}
	    	}
	    	free(sum);
	    	free(hash);
	    	return sorted;
	    } else {
	    	return arr;
    	}
	} else {
		return arr;
	}
}

/*procrastisort in serial*/
double* procrastisortS(uint length, double *arr) {
   
    if (length > 1) {
    
		double min = arr[0];
		double max = arr[0];
		
		for (int i = 1; i < length; i++) {
			if (arr[i] < min) {
				min = arr[i];
			}
			if (arr[i] > max) {
				max = arr[i];
			}
		}
		    
    	double range = max - min;
    	
    	if (range > 0) {
    		
    		double *sorted = (double*)malloc(length*sizeof(double));
    	
    		int binIndex;
   		 	int numRec;
    		numRec = 0;
    		double binSize = (range)/(length-1);

			/*subtract min out*/
    		if (min > 0) {
    			for (uint i = 0; i < length; i++) {
    				arr[i] -= min;
    			}
    		}
    
	    	/*generate histogram*/
	    	int *num = (int*)malloc(length*sizeof(int));
	    	int *sum = (int*)malloc((length+1)*sizeof(int));
	    	
	    	/*initialize to 0*/
			memset(num, 0, (length)*sizeof(int));
			
			/*tally number in each bin*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = (int) ((arr[i])/binSize);
				num[binIndex]++;
				if (num[binIndex] == 2) {
					numRec++;
				}
	    	}
	    	
			/*prefix scanning*/
			sum[0] = 0;
	    	for (int i = 1; i < length+1; i++) {
				sum[i] = sum[i-1] + num[i-1];
	    	}
	    	free(num);

	    	/*set all elements of hash array to -1*/
	    	int *hash = (int*)malloc((length)*sizeof(int));
	    	memset(hash, -1, (length)*sizeof(int));

			/*hashing using the prefix scan*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = sum[(int) ((arr[i])/binSize)];
				while (hash[binIndex] != -1 && binIndex < length) { 
					binIndex++;
				} hash[binIndex] = i;
	    	}
    
    		/*sort collisions*/
	    	if (numRec > 0) {
	    		int colTally = 0;
	    		
	    		int *colLocs = (int*)malloc(numRec*sizeof(int));
	    		int *pJumps = (int*)malloc(numRec*sizeof(int));
	    		
	    		/*find collision locations and the number of collisions at that location*/
	    		for (int i = 0; i < numRec; i++) {
					int prefixJump = sum[colTally+1]-sum[colTally];
					while (prefixJump <= 1) {
	    					colTally++;
	    					prefixJump = sum[colTally+1]-sum[colTally];
	    			} 
	    			colLocs[i] = colTally;
	    			pJumps[i] = prefixJump;
	    			colTally++;
	    		}
	    		
	    		/*sort colliders recursively but could be done differently*/
				for (int i = 0; i < numRec; i++) {
					double *colliders;
					int prefixJump = pJumps[i];
					int colLoc = colLocs[i];
	    			
		    		colliders = (double*)malloc(prefixJump*sizeof(double));
		    		for (int j = 0; j < prefixJump; j++) {
		    			double collider = arr[hash[sum[colLoc] + j]];
		    			colliders[j] = collider;
		    		}
		    		
		    		double* tempCols = procrastisortR(prefixJump, colliders);
		    		
	    			for (int j = 0; j < prefixJump; j++) {
	    				arr[hash[sum[colLoc] + j]] = tempCols[j];
	    			}
	    			
	    			if (tempCols == colliders) {
	    				free(tempCols);
	    				tempCols = NULL;
	    			} else {
	    			free(colliders);
	    			free(tempCols);
	    			}
		    	}
		    	free(colLocs);
		    	free(pJumps);
			}	
			
			/*finish by using the keys*/
			if (min > 0) {
				for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]] + min;
	    			arr[hash[i]] += min;
	    		
	    		}
	    	} else {
	    		for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]];
	    		}
	    	}
	    	free(sum);
	    	free(hash);
	    	return sorted;
	    } else {
	    	return arr;
    	}
	} else {
		return arr;
	}
}

/*normal procrastisort runs in parallel using open mp*/
double* procrastisort(uint length, double *arr) {
	/*stops the recursion*/
    if (length > 1) {
		double min = arr[0];
		double max = arr[0];
		
		/*finds the min and max. demanding procrastisort has min and max as an argument*/
		for (int i = 1; i < length; i++) {
			if (arr[i] < min) {
				min = arr[i];
			}
			if (arr[i] > max) {
				max = arr[i];
			}
		}
		    
    	double range = max - min;
    	
		/*if range is 0, the numbers are all the same and don't need to be sorted*/
    	if (range > 0) {
    		
    		double *sorted = (double*)malloc(length*sizeof(double));
    	
    		int binIndex;
			/*number of recursions*/
   		 	int numRec = 0;
    		double binSize = (range)/(length-1);

			/*subtract min out*/
    		if (min > 0) {
    			#pragma omp parallel for
    			for (uint i = 0; i < length; i++) {
    				arr[i] -= min;
    			}
    		}
    
			/*get histogram and prefix scan*/
	    	int *num = (int*)malloc(length*sizeof(int));
	    	int *sum = (int*)malloc((length+1)*sizeof(int));
    		
			/*initialise to zero*/
			memset(num, 0, (length)*sizeof(int));
    
    		//#pragma omp parallel for
	    	for(uint i = 0; i < length; i++) {
				binIndex = (int) ((arr[i])/binSize);
				num[binIndex]++;
				if (num[binIndex] == 2) {
					numRec++;
				}
	    	}
	    	
			/*prefix scanning*/
			sum[0] = 0;
	    	for (int i = 1; i < length+1; i++) {
				sum[i] = sum[i-1] + num[i-1];
	    	}
	    	free(num);

	    	/*set all elements of hash array to -1*/
	    	int *hash = (int*)malloc((length)*sizeof(int));
	    	memset(hash, -1, (length)*sizeof(int));

			/*hashing usig the prefix scan*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = sum[(int) ((arr[i])/binSize)];
				while (hash[binIndex] != -1 && binIndex < length) { 
					binIndex++;
				} hash[binIndex] = i;
	    	}
    		
    		/*handle collisions*/
	    	if (numRec > 0) {
	    		int colTally = 0;
	    		
	    		int *colLocs = (int*)malloc(numRec*sizeof(int));
	    		int *pJumps = (int*)malloc(numRec*sizeof(int));
	    		
	    		/*get collision locations and number of collisions there*/
	    		for (int i = 0; i < numRec; i++) {
					int prefixJump = sum[colTally+1]-sum[colTally];
					while (prefixJump <= 1) {
	    					colTally++;
	    					prefixJump = sum[colTally+1]-sum[colTally];
	    			} 
	    			colLocs[i] = colTally;
	    			pJumps[i] = prefixJump;
	    			colTally++;
	    		}
		 
		 
		 		/*sort colliders recursively in parallel*/
		 		#pragma omp parallel for
				for (int i = 0; i < numRec; i++) {
					double *colliders;
					int prefixJump = pJumps[i];
					int colLoc = colLocs[i];
	    			
		    		colliders = (double*)malloc(prefixJump*sizeof(double));
		    		for (int j = 0; j < prefixJump; j++) {
		    			double collider = arr[hash[sum[colLoc] + j]];
		    			colliders[j] = collider;
		    		}
		    		
		    		double* tempCols = procrastisortR(prefixJump, colliders);
		    		
	    			for (int j = 0; j < prefixJump; j++) {
	    				arr[hash[sum[colLoc] + j]] = tempCols[j];
	    			}
	    			
	    			if (tempCols == colliders) {
	    				free(tempCols);
	    				tempCols = NULL;
	    			} else {
	    			free(colliders);
	    			free(tempCols);
	    			}
		    	}
		    	free(colLocs);
		    	free(pJumps);
			}	
			
			/*finsih by using the keys in the hash*/
			if (min > 0) {
				#pragma omp parallel for
				for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]] + min;
	    			arr[hash[i]] += min;
	    		
	    		}
	    	} else {
	    		#pragma omp parallel for
	    		for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]];
	    		}
	    	}
	    	free(sum);
	    	free(hash);
	    	return sorted;
	    } else {
	    	return arr;
    	}
	} else {
		return arr;
	}
}

double* procrastisortQ(uint length, double *arr) {
	/*stops the recursion*/
    if (length > 1) {
		double min = arr[0];
		double max = arr[0];
		
		/*finds the min and max. demanding procrastisort has min and max as an argument*/
		for (int i = 1; i < length; i++) {
			if (arr[i] < min) {
				min = arr[i];
			}
			if (arr[i] > max) {
				max = arr[i];
			}
		}
		    
    	double range = max - min;
    	
		/*if range is 0, the numbers are all the same and don't need to be sorted*/
    	if (range > 0) {
    		
    		double *sorted = (double*)malloc(length*sizeof(double));
    	
    		int binIndex;
			/*number of recursions*/
   		 	int numRec = 0;
    		double binSize = (range)/(length-1);

			/*subtract min out*/
    		if (min > 0) {
    			#pragma omp parallel for
    			for (uint i = 0; i < length; i++) {
    				arr[i] -= min;
    			}
    		}
    
			/*get histogram and prefix scan*/
	    	int *num = (int*)malloc(length*sizeof(int));
	    	int *sum = (int*)malloc((length+1)*sizeof(int));
    		
			/*initialise to zero*/
			memset(num, 0, (length)*sizeof(int));
    
    		//#pragma omp parallel for
	    	for(uint i = 0; i < length; i++) {
				binIndex = (int) ((arr[i])/binSize);
				num[binIndex]++;
				if (num[binIndex] == 2) {
					numRec++;
				}
	    	}
	    	
			/*prefix scanning*/
			sum[0] = 0;
	    	for (int i = 1; i < length+1; i++) {
				sum[i] = sum[i-1] + num[i-1];
	    	}
	    	free(num);

	    	/*set all elements of hash array to -1*/
	    	int *hash = (int*)malloc((length)*sizeof(int));
	    	memset(hash, -1, (length)*sizeof(int));

			/*hashing usig the prefix scan*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = sum[(int) ((arr[i])/binSize)];
				while (hash[binIndex] != -1 && binIndex < length) { 
					binIndex++;
				} hash[binIndex] = i;
	    	}
    		
    		/*handle collisions*/
	    	if (numRec > 0) {
	    		int colTally = 0;
	    		
	    		int *colLocs = (int*)malloc(numRec*sizeof(int));
	    		int *pJumps = (int*)malloc(numRec*sizeof(int));
	    		
	    		/*get collision locations and number of collisions there*/
	    		for (int i = 0; i < numRec; i++) {
					int prefixJump = sum[colTally+1]-sum[colTally];
					while (prefixJump <= 1) {
	    					colTally++;
	    					prefixJump = sum[colTally+1]-sum[colTally];
	    			} 
	    			colLocs[i] = colTally;
	    			pJumps[i] = prefixJump;
	    			colTally++;
	    		}
		 
		 
		 		/*sort colliders recursively in parallel*/
		 		#pragma omp parallel for
				for (int i = 0; i < numRec; i++) {
					double *colliders;
					int prefixJump = pJumps[i];
					int colLoc = colLocs[i];
	    			
		    		colliders = (double*)malloc(prefixJump*sizeof(double));
		    		for (int j = 0; j < prefixJump; j++) {
		    			double collider = arr[hash[sum[colLoc] + j]];
		    			colliders[j] = collider;
		    		}
		    		
		    		qsort(colliders, prefixJump, sizeof(double), compare);
		    		
	    			for (int j = 0; j < prefixJump; j++) {
	    				arr[hash[sum[colLoc] + j]] = colliders[j];
	    			}
	    			
	    			free(colliders);
	    			
		    	}
		    	free(colLocs);
		    	free(pJumps);
			}	
			
			/*finsih by using the keys in the hash*/
			if (min > 0) {
				for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]] + min;
	    			arr[hash[i]] += min;
	    		
	    		}
	    	} else {
	    		for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]];
	    		}
	    	}
	    	free(sum);
	    	free(hash);
	    	return sorted;
	    } else {
	    	return arr;
    	}
	} else {
		return arr;
	}
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

double* procrastisortI(uint length, double *arr) {
	/*stops the recursion*/
    if (length > 1) {
		double min = arr[0];
		double max = arr[0];
		
		/*finds the min and max. demanding procrastisort has min and max as an argument*/
		for (int i = 1; i < length; i++) {
			if (arr[i] < min) {
				min = arr[i];
			}
			if (arr[i] > max) {
				max = arr[i];
			}
		}
		    
    	double range = max - min;
    	
		/*if range is 0, the numbers are all the same and don't need to be sorted*/
    	if (range > 0) {
    		
    		double *sorted = (double*)malloc(length*sizeof(double));
    	
    		int binIndex;
			/*number of recursions*/
   		 	int numRec = 0;
    		double binSize = (range)/(length-1);

			/*subtract min out*/
    		if (min > 0) {
    			for (uint i = 0; i < length; i++) {
    				arr[i] -= min;
    			}
    		}
    
			/*get histogram and prefix scan*/
	    	int *num = (int*)malloc(length*sizeof(int));
	    	int *sum = (int*)malloc((length+1)*sizeof(int));
    		
			/*initialise to zero*/
			memset(num, 0, (length)*sizeof(int));
    
    		//#pragma omp parallel for
	    	for(uint i = 0; i < length; i++) {
				binIndex = (int) ((arr[i])/binSize);
				num[binIndex]++;
				if (num[binIndex] == 2) {
					numRec++;
				}
	    	}
	    	
			/*prefix scanning*/
			sum[0] = 0;
	    	for (int i = 1; i < length+1; i++) {
				sum[i] = sum[i-1] + num[i-1];
	    	}
	    	free(num);

	    	/*set all elements of hash array to -1*/
	    	int *hash = (int*)malloc((length)*sizeof(int));
	    	memset(hash, -1, (length)*sizeof(int));

			/*hashing usig the prefix scan*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = sum[(int) ((arr[i])/binSize)];
				while (hash[binIndex] != -1 && binIndex < length) { 
					binIndex++;
				} hash[binIndex] = i;
	    	}
	    	
	    	/*finsih by using the keys in the hash*/
			if (min > 0) {
				for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]] + min;
	    			arr[hash[i]] += min;
	    		
	    		}
	    	} else {
	    		for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]];
	    		}
	    	}
    		
    		/*handle collisions*/
	    	insertionsort(length, sorted);	

	    	free(sum);
	    	free(hash);
	    	return sorted;
	    } else {
	    	return arr;
    	}
	} else {
		return arr;
	}
}

double* procrastisortI2(uint length, double *arr) {
	/*stops the recursion*/
    if (length > 1) {
		double min = arr[0];
		double max = arr[0];
		
		/*finds the min and max. demanding procrastisort has min and max as an argument*/
		for (int i = 1; i < length; i++) {
			if (arr[i] < min) {
				min = arr[i];
			}
			if (arr[i] > max) {
				max = arr[i];
			}
		}
		    
    	double range = max - min;
    	
		/*if range is 0, the numbers are all the same and don't need to be sorted*/
    	if (range > 0) {
    		
    		double *sorted = (double*)malloc(length*sizeof(double));
    	
    		int binIndex;
			/*number of recursions*/
   		 	int numRec = 0;
    		double binSize = (range)/(length-1);

			/*subtract min out*/
    		if (min > 0) {
    			//#pragma omp parallel for
    			for (uint i = 0; i < length; i++) {
    				arr[i] -= min;
    			}
    		}
    
			/*get histogram and prefix scan*/
	    	int *num = (int*)malloc(length*sizeof(int));
	    	int *sum = (int*)malloc((length+1)*sizeof(int));
    		
			/*initialise to zero*/
			memset(num, 0, (length)*sizeof(int));
    
    		//#pragma omp parallel for
	    	for(uint i = 0; i < length; i++) {
				binIndex = (int) ((arr[i])/binSize);
				num[binIndex]++;
				if (num[binIndex] == 2) {
					numRec++;
				}
	    	}
	    	
			/*prefix scanning*/
			sum[0] = 0;
	    	for (int i = 1; i < length+1; i++) {
				sum[i] = sum[i-1] + num[i-1];
	    	}
	    	free(num);

	    	/*set all elements of hash array to -1*/
	    	int *hash = (int*)malloc((length)*sizeof(int));
	    	memset(hash, -1, (length)*sizeof(int));

			/*hashing usig the prefix scan*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = sum[(int) ((arr[i])/binSize)];
				while (hash[binIndex] != -1 && binIndex < length) { 
					binIndex++;
				} hash[binIndex] = i;
	    	}
    		
    		/*finsih by using the keys in the hash*/
			if (min > 0) {
				for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]] + min;
	    			arr[hash[i]] += min;
	    		
	    		}
	    	} else {
	    		for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]];
	    		}
	    	}
    		
    		/*handle collisions*/
	    	insertionsort(length, sorted);
			
			
	    	free(sum);
	    	free(hash);
	    	return sorted;
	    } else {
	    	return arr;
    	}
	} else {
		return arr;
	}
}

double* procrastisortB(uint length, double *arr) {
	/*stops the recursion*/
    if (length > 1) {
		double min = arr[0];
		double max = arr[0];
		
		/*finds the min and max. demanding procrastisort has min and max as an argument*/
		for (int i = 1; i < length; i++) {
			if (arr[i] < min) {
				min = arr[i];
			}
			if (arr[i] > max) {
				max = arr[i];
			}
		}
		    
    	double range = max - min;
    	
		/*if range is 0, the numbers are all the same and don't need to be sorted*/
    	if (range > 0) {
    		
    		double *sorted = (double*)malloc(length*sizeof(double));
    	
    		int binIndex;
			/*number of recursions*/
   		 	int numRec = 0;
    		double binSize = (range)/(length-1);

			/*subtract min out*/
    		if (min > 0) {
    			for (uint i = 0; i < length; i++) {
    				arr[i] -= min;
    			}
    		}
    
			/*get histogram and prefix scan*/
	    	int *num = (int*)malloc(length*sizeof(int));
	    	int *sum = (int*)malloc((length+1)*sizeof(int));
    		
			/*initialise to zero*/
			memset(num, 0, (length)*sizeof(int));
    
    		//#pragma omp parallel for
	    	for(uint i = 0; i < length; i++) {
				binIndex = (int) ((arr[i])/binSize);
				num[binIndex]++;
				if (num[binIndex] == 2) {
					numRec++;
				}
	    	}
	    	
			/*prefix scanning*/
			sum[0] = 0;
	    	for (int i = 1; i < length+1; i++) {
				sum[i] = sum[i-1] + num[i-1];
	    	}
	    	free(num);

	    	/*set all elements of hash array to -1*/
	    	int *hash = (int*)malloc((length)*sizeof(int));
	    	memset(hash, -1, (length)*sizeof(int));

			/*hashing using the prefix scan*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = sum[(int) ((arr[i])/binSize)];
				while (hash[binIndex] != -1 && binIndex < length) { 
					binIndex++;
				} hash[binIndex] = i;
	    	}
			
			/*finsih by using the keys in the hash*/
			if (min > 0) {
				for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]] + min;
	    			/*re add the minimum value*/
	    			arr[hash[i]] += min;
	    		
	    		}
	    	} else {
	    		for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]];
	    		}
	    	}
	    	free(sum);
	    	free(hash);
	    	
	    	/*handle collisions*/
	    	bubblesort(length, sorted);
	    	
	    	return sorted;
	    } else {
	    	return arr;
    	}
	} else {
		return arr;
	}
}

double* procrastisortB2(uint length, double *arr) {
	/*stops the recursion*/
    if (length > 1) {
		double min = arr[0];
		double max = arr[0];
		
		/*finds the min and max. demanding procrastisort has min and max as an argument*/
		for (int i = 1; i < length; i++) {
			if (arr[i] < min) {
				min = arr[i];
			}
			if (arr[i] > max) {
				max = arr[i];
			}
		}
		    
    	double range = max - min;
    	
		/*if range is 0, the numbers are all the same and don't need to be sorted*/
    	if (range > 0) {
    		
    		double *sorted = (double*)malloc(length*sizeof(double));
    	
    		int binIndex;
    		double binSize = (range)/(length-1);

			/*subtract min out*/
    		if (min > 0) {
    			//#pragma omp parallel for
    			for (uint i = 0; i < length; i++) {
    				arr[i] -= min;
    			}
    		}
    
			/*get histogram and prefix scan*/
	    	int *num = (int*)malloc(length*sizeof(int));
	    	int *sum = (int*)malloc((length+1)*sizeof(int));
    		
			/*initialise to zero*/
			memset(num, 0, (length)*sizeof(int));
    
    		//#pragma omp parallel for
	    	for(uint i = 0; i < length; i++) {
				binIndex = (int) ((arr[i])/binSize);
				num[binIndex]++;
	    	}
	    	
			/*prefix scanning*/
			sum[0] = 0;
	    	for (int i = 1; i < length+1; i++) {
				sum[i] = sum[i-1] + num[i-1];
	    	}
	    	free(num);

	    	/*set all elements of hash array to -1*/
	    	int *hash = (int*)malloc((length)*sizeof(int));
	        //memset(hash, -1, (length)*sizeof(int));

			/*hashing usig the prefix scan*/
	    	for(uint i = 0; i < length; i++) {
				binIndex = sum[(int) ((arr[i])/binSize)];
				hash[binIndex] = i; 
				sum[(int) ((arr[i])/binSize)]++;
	    	}
    		
    		/*finsih by using the keys in the hash*/
			if (min > 0) {
				for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]] + min;
	    			arr[hash[i]] += min;
	    		
	    		}
	    	} else {
	    		for (int i = 0; i < length; i++) {
	    			sorted[i] = arr[hash[i]];
	    		}
	    	}
    		
    		/*handle collisions*/
	    	bubblesort(length, sorted);
			
			
	    	free(sum);
	    	free(hash);
	    	return sorted;
	    } else {
	    	return arr;
    	}
	} else {
		return arr;
	}
}

/* generate a randomly mixed up array with size size to be stored in pointer. the elements will have a minimum value min, and
    the difference between elements when sorted will be between mindx and maxdx. the maximum value is recorded in max. */
void generate_array(uint size, double *ptr, double mindx, double maxdx, double min, double *max ) {

	srand(1);
    
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

void true_random(uint length, double *arr) {
	srand(time(NULL));
	double min = 0;
	double max = INT_MAX;
	
	for (uint i = 0; i < length; i++) {
		arr[i] = ((double)rand()) * max / ((double)RAND_MAX);
	}
}


int main (int argc, char** argv) {
	
	double procrastisortTime = 0;
	double qsortTime = 0;
	double procrastisortDTime = 0;
	double procrastisortSTime = 0;
	double procrastisortQTime = 0;
	double procrastisortITime = 0;
	double procrastisortI2Time = 0;
	double procrastisortBTime = 0;
	double procrastisortB2Time = 0;
	double hashsortTime = 0;
	uint length;
	uint numRep;
	/*length of the random array is arg 1*/
	length = atoi (argv[1]);
	/*arg 2 is the number of repetitions to ensure accurate data*/
	numRep = atoi (argv[2]);
	double mindx = 0, maxdx = 10, min = 0, max=0;
	double *arr = (double*)malloc(length*sizeof(double));
	double *times = malloc(10*sizeof(double));
	memset(times, 0, 10*sizeof(double));
	//generate_array(length, arr, mindx, maxdx, min, &max);
	true_random(length, arr);
	
	/*for (int i=0; i<length; i++){
			//arr[i] = i - 5;
			printf("%.16lf, \n", arr[i]);
		}
		printf("\n");*/
	
	for (int i = 0; i < numRep; i++) {
	
		double *sorted=NULL, *sort_test=NULL;
		int icount = 0;
	

	
		//printf("max: %f,\tmin: %f\n", max, min);
	
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
    	times[8] += (t2-t1);
    	/*for (int i=0; i<length; i++){
			printf("%f, \n", sorted[i]);
		}
		printf("\n");*/
    
    
		//sort_test = (double*)malloc(length*sizeof(double));
		gettimeofday(&timer, NULL);
   		t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
		sort_test = procrastisort(length, arr); 
		gettimeofday(&timer, NULL);
    	t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
    
    
    	//printf("procrastisort time: %.6lf,\n", t2 - t1);
    	times[0] += (t2-t1);
    
    	//icount=0;
		/*for(uint i = 0; i < length; i++) {
			if (sort_test[i] != sorted[i]) {
			//printf("Check failed for procrastisort CPU index %d procrastisort value %0.20lf gold standard %0.16lf\n",i,sort_test[i],sorted[i]);
			icount++;
			}
		}
		if (icount > 0) {
			printf("errors: %i\n", icount);
		}*/
    
    
	    /*for (int i=0; i<length; i++){
			printf("%f, \n", sort_test[i]);
		}
		printf("\n");*/
		
		if (sort_test != arr) {free(sort_test);}
		
		//sort_test=NULL;
		
		//sort_test = (double*)malloc(length*sizeof(double));
		//if (sort_test != arr) {
			//memcpy(sort_test, arr, length*sizeof(double));
		//}
		gettimeofday(&timer, NULL);
	    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
		//sort_test = procrastisortD(length, arr, min, max); 
		gettimeofday(&timer, NULL);
	    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
	    
    
	    //printf("procrastisort time: %.6lf,\n", t2 - t1);
	    //times[1] += (t2-t1);
	    
		//if (sort_test != arr) {free(sort_test);}
		//sort_test=NULL;
		
		//sort_test = (double*)malloc(length*sizeof(double));
		//memcpy(sort_test, arr, length*sizeof(double));
		gettimeofday(&timer, NULL);
	    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
		sort_test = procrastisortQ(length, arr); 
		gettimeofday(&timer, NULL);
	    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
	    
    
	    //printf("procrastisort time: %.6lf,\n", t2 - t1);
	    times[2] += (t2-t1);
	    
		if (sort_test != arr) {free(sort_test);}
		//sort_test=NULL;
		
		//sort_test=(double*)malloc(length*sizeof(double));
		//memcpy(sort_test, arr, length*sizeof(double));
		gettimeofday(&timer, NULL);
	    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
		sort_test = procrastisortS(length, arr); 
		gettimeofday(&timer, NULL);
	    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
	    
	    
	    //printf("procrastisort time: %.6lf,\n", t2 - t1);
		times[7] += (t2-t1);
	    
	    if (sort_test != arr) {free(sort_test);}
		//sort_test=NULL;
	
		//sort_test=(double*)malloc(length*sizeof(double));
		//memcpy(sort_test, arr, length*sizeof(double));
		gettimeofday(&timer, NULL);
	    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
		sort_test = procrastisortI(length, arr); 
		gettimeofday(&timer, NULL);
	    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
    
    
	    //printf("procrastisort time: %.6lf,\n", t2 - t1);
	    times[3] += (t2-t1);
	    
	    if (sort_test != arr) {free(sort_test);}
		//sort_test=NULL;
	
		//sort_test=(double*)malloc(length*sizeof(double));
		//memcpy(sort_test, arr, length*sizeof(double));
		gettimeofday(&timer, NULL);
	    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
		sort_test = procrastisortI2(length, arr); 
		gettimeofday(&timer, NULL);
	    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
	    
    
	    //printf("procrastisort time: %.6lf,\n", t2 - t1);
	    times[4] += (t2-t1);
	    
	    if (sort_test != arr) {free(sort_test);}
		//sort_test=NULL;
	
		//sort_test=(double*)malloc(length*sizeof(double));
		//memcpy(sort_test, arr, length*sizeof(double));
		gettimeofday(&timer, NULL);
	    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
		sort_test = procrastisortB(length, arr); 
		gettimeofday(&timer, NULL);
	    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
    
    
	    //printf("procrastisort time: %.6lf,\n", t2 - t1);
	    times[5] += (t2-t1);
	    
	    if (sort_test != arr) {free(sort_test);}
		//sort_test=NULL;
	
		//sort_test=(double*)malloc(length*sizeof(double));
		//memcpy(sort_test, arr, length*sizeof(double));
		gettimeofday(&timer, NULL);
	    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
		sort_test = procrastisortB2(length, arr); 
		gettimeofday(&timer, NULL);
	    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
	    
    
	    //printf("procrastisort time: %.6lf,\n", t2 - t1);
	    times[6] += (t2-t1);
	    
	    /*icount=0;
	    for(uint i = 0; i < length; i++) {
	       if (sort_test[i] != sorted[i]) {
	          //printf("Check failed for procrastisort CPU index %d procrastisort value %0.20lf gold standard %0.16lf\n",i,sort_test[i],sorted[i]);
			icount++;
	       }
	    }
	    if (icount > 0) {
	    	printf("errors: %i\n", icount);
	    }*/
	    
		if (sort_test != arr) {free(sort_test);}
		//sort_test=NULL;
	    
	    //sort_test=(double*)malloc(length*sizeof(double));
		//memcpy(sort_test, arr, length*sizeof(double));
		gettimeofday(&timer, NULL);
	    t1 = timer.tv_sec+(timer.tv_usec/1000000.0);
		//sort_test = hashsort (length, arr, mindx, min, max); 
		gettimeofday(&timer, NULL);
	    t2 = timer.tv_sec+(timer.tv_usec/1000000.0);
    
	    //if (sort_test != arr) {free(sort_test);}
	    //sort_test=NULL;
	    
	    times[9] += (t2 - t1);
	    //printf("hashsort time: %.6lf,\n", t2 - t1);
		
		free(sorted);
	}
	free(arr);
	
	for (int i = 0; i < 10; i++) {
		times[i] /= numRep;
	}
	
	printf("procrastisort:\t\t\t%.8lf\ndemanding procrastisort:\t%.8lf\nquicksort procrastisort:\t%.8lf\ninsertion1 procrastisort:\t%.8lf\ninsertion2 procrastisort:\t%.8lf\nbubble1 procrastisort:\t\t%.8lf\nbubble2 procrastisort:\t\t%.8lf\nserial procrastisort:\t\t%.8lf\nquicksort:\t\t\t%.8lf\nhashsort:\t\t\t%.8lf\n", (times[0]), (times[1]), (times[2]), (times[3]), (times[4]), (times[5]), (times[6]),  (times[7]), (times[8]), (times[9]));
	
	double minTime = DBL_MAX;
	int minKey;
	for (int i = 0; i < 10; i++) {
		if (times[i] < minTime) {
		    if (i != 1 && i != 9) {
			minKey = i;
			minTime = times[i];
			}
		}
	}
	
	switch (minKey) {
	case 0:
		printf("procrastisort: %.8lf\n", times[0]);
		break;
	case 1:
		printf("demanding procrastisort: %.8lf\n", times[1]);
		break;
	case 2:
		printf("quicksort procrastisort: %.8lf\n", times[2]);
		break;
	case 3:
		printf("insetion 1 procrastisort: %.8lf\n", times[3]);
		break;
	case 4:
		printf("insertion 2 procrastisort: %.8lf\n", times[4]);
		break;
	case 5:
		printf("bubble 1 procrastisort: %.8lf\n", times[5]);
		break;
	case 6:
		printf("bubble 2 procrastisort: %.8lf\n", times[6]);
		break;
	case 7:
		printf("serial procrastisort: %.8lf\n", times[7]);
		break;
	case 8:
		printf("quicksort: %.8lf\n", times[8]);
		break;
	case 9:
		printf("hashsort: %.8lf\n", times[9]);
		break;
	}
	return 0;
}


