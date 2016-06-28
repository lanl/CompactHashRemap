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
#include <math.h>
#include "meshgen.h"

static bool randomize = true;

#define SQ(x) (( (x)*(x) ))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

void swap_double(double** a, double** b) {
  double* c = *a;
  *a = *b;
  *b = c;
}

void swap_int(int** a, int** b) {
  int* c = *a;
  *a = *b;
  *b = c;
}

void swap_uint(uint** a, uint** b) {
  uint* c = *a;
  *a = *b;
  *b = c;
}

int powerOfFour(int n) {
  int result = 1;
  int i;
  for(i = 0; i < n; i++) {
    result *= 4;
  }
  return result;
}

cell_list new_cell_list(uint *x, uint *y, uint *lev, double *values) {
  cell_list a;
  a.i = x;
  a.j = y;
  a.level = lev;
  a.values = values;
  return a;
}

cell_list create_cell_list(cell_list a, uint length) {
    a.ncells = length;
    a.i      = (uint *)   malloc(length * sizeof(uint));
    a.j      = (uint *)   malloc(length * sizeof(uint));
    a.level  = (uint *)   malloc(length * sizeof(uint));
    a.values = (double *) malloc(length * sizeof(double));
    return(a);
}

void destroy(cell_list a) {
    free(a.i);
    free(a.j);
    free(a.level);
    free(a.values);
}

cell_list mesh_maker_level (cell_list clist, uint levels_diff, uint *length, uint *max_level, 
    uint *min_level) {
    *max_level = levels_diff;
    uint cell_count;
    
    //use a local variable because pointers are hard
    uint num_cells = *length;
    
    if ((num_cells - 1) % 3 != 0) {
        num_cells--;
        num_cells /= 3;
        num_cells *= 3;
        num_cells ++;
        printf("\nImpossible number of cells, using %u instead\n", num_cells);
    }
    
    if (num_cells < levels_diff * 3 + 1) {
        num_cells = levels_diff * 3 + 1;
        printf("Impossible number of cells, using %u instead\n", num_cells);
    }
    
    uint max_level_temp = *max_level;
    while (four_to_the(max_level_temp) < (int)*length) {
        max_level_temp++;
    }
    
    *max_level = max_level_temp;
    *length = num_cells;
    
    *min_level = *max_level - levels_diff;
    
    cell_count = 1;
    uint cell_target, current_max_lev = 0, lev;
    
    clist = create_cell_list (clist, *length);
    clist.i[0]      =  0;
    clist.j[0]      =  0;
    clist.level[0]  =  0;
    clist.values[0] = -1;
    
    for (uint n = 0; n < *min_level; n++) {
        for (int m = 0; m < four_to_the(n); m++) {
            divide_cell (clist.i[m], clist.j[m],
                clist.level[m], clist, cell_count, m);
            cell_count += 3;
        }
        current_max_lev ++;
    }
    
    while (current_max_lev < *max_level || cell_count < *length) {
        cell_target = (uint) rand() % (cell_count);
        //print_cell (clist, cell_target);
        lev = clist.level[cell_target];
        
        if (lev < *max_level && (current_max_lev == *max_level || lev == current_max_lev)) {
            divide_cell (clist.i[cell_target], clist.j[cell_target],
                clist.level[cell_target], clist, cell_count, cell_target);
            
            if (lev + 1 > current_max_lev) {
                current_max_lev = lev + 1;
            }
            cell_count += 3;
        }
    }
    clist.ncells = *length;
    clist.ibasesize = 1;//four_to_the(*min_level-1);
    clist.levmax = *max_level;
    return clist;
}

cell_list mesh_maker_sparsity (cell_list clist, uint levels_diff, uint *length, uint *max_level, 
    uint *min_level, double sparsity) {
    *max_level = levels_diff;
    uint cell_count;
    
    //use a local variable because pointers are hard
    uint num_cells = *length;
    
    if ((num_cells - 1) % 3 != 0) {
        num_cells--;
        num_cells /= 3;
        num_cells *= 3;
        num_cells ++;
        printf("Impossible number of cells, using %u instead\n", num_cells);
    }
    
    if (num_cells < levels_diff * 3 + 1) {
        num_cells = levels_diff * 3 + 1;
        printf("Impossible number of cells, using %u instead\n", num_cells);
    }
    
    uint max_level_temp = *max_level;
    while (four_to_the(max_level_temp) < (int)*length / sparsity) {
        max_level_temp++;
    }
    
    *max_level = max_level_temp;
    *length = num_cells;
    
    *min_level = *max_level - levels_diff;
    
    cell_count = 1;
    uint cell_target, current_max_lev = 0, lev;
    
    clist = create_cell_list (clist, *length);
    clist.i[0]      =  0;
    clist.j[0]      =  0;
    clist.level[0]  =  0;
    clist.values[0] = -1;
    
    for (uint n = 0; n < *min_level; n++) {
        for (int m = 0; m < four_to_the(n); m++) {
            divide_cell (clist.i[m], clist.j[m],
                clist.level[m], clist, cell_count, m);
            cell_count += 3;
        }
        current_max_lev ++;
    }
    
    while (current_max_lev < *max_level || cell_count < *length) {
        cell_target = (uint) rand() % (cell_count);
        //print_cell (clist, cell_target);
        lev = clist.level[cell_target];
        
        if (lev < *max_level && (current_max_lev == *max_level || lev == current_max_lev)) {
            divide_cell (clist.i[cell_target], clist.j[cell_target],
                clist.level[cell_target], clist, cell_count, cell_target);
            
            if (lev + 1 > current_max_lev) {
                current_max_lev = lev + 1;
            }
            cell_count += 3;
        }
    }
    return clist;
}

void divide_cell (uint super_i, uint super_j, uint super_level, cell_list cells, 
    uint cell_count, uint cell_id) {
    
    uint key, ibase, jbase;

    key = translate_cell (super_i, super_j, super_level, super_level + 1);
    ibase = key % two_to_the(super_level + 1);
    jbase = key / two_to_the(super_level + 1);
    
    //New cell 1:
    cells.i[cell_id] = ibase;
    cells.j[cell_id] = jbase;
    cells.level[cell_id] = super_level + 1;
    cells.values[cell_id] =  0xFFFFFFFF;
    //New cell 2:
    cells.i[cell_count] = ibase + 1;
    cells.j[cell_count] = jbase;
    cells.level[cell_count] = super_level + 1;
    cells.values[cell_count] =  0xFFFFFFFF;
    //New cell 3:
    cells.i[cell_count + 1] = ibase;
    cells.j[cell_count + 1] = jbase + 1;
    cells.level[cell_count + 1] = super_level + 1;
    cells.values[cell_count + 1] =  0xFFFFFFFF;
    //New cell 4:
    cells.i[cell_count + 2] = ibase + 1;
    cells.j[cell_count + 2] = jbase + 1;
    cells.level[cell_count + 2] = super_level + 1;
    cells.values[cell_count + 2] =  0xFFFFFFFF;
}

uint translate_cell (uint i, uint j, uint lev, uint new_lev) {
    uint j_comp, i_comp;
    if (new_lev < lev) {
        //j_comp = j / two_to_the(lev - new_lev);
        //i_comp = i / two_to_the(lev - new_lev);
        j_comp = j >> (lev - new_lev);
        i_comp = i >> (lev - new_lev);
    } else {
        //j_comp = j * two_to_the(new_lev - lev);
        //i_comp = i * two_to_the(new_lev - lev);
        j_comp = j << (new_lev - lev);
        i_comp = i << (new_lev - lev);
    }

    //printf("j_comp: %u\tx_comp: %u\n", j_comp, i_comp);
    uint key = (j_comp * two_to_the(new_lev)) + i_comp;
    //printf("newKey: %d\n", newKey);
    return key;
}

uint translate_cell (uint i, uint j, uint lev, uint new_lev, int ibasesize) {
    uint j_comp, i_comp;
    if (new_lev < lev) {
        //j_comp = j / two_to_the(lev - new_lev);
        //i_comp = i / two_to_the(lev - new_lev);
        j_comp = j >> (lev - new_lev);
        i_comp = i >> (lev - new_lev);
    } else {
        //j_comp = j * two_to_the(new_lev - lev);
        //i_comp = i * two_to_the(new_lev - lev);
        j_comp = j << (new_lev - lev);
        i_comp = i << (new_lev - lev);
    }

    //printf("j_comp: %u\tx_comp: %u\n", j_comp, i_comp);
    uint key = (j_comp * ibasesize*two_to_the(new_lev)) + i_comp;
    //printf("newKey: %d\n", newKey);
    return key;
}

void print_cell_list (cell_list cells, uint length) {
    printf("#\ti\tj\tlevel\tvalue\n");
    for (uint n = 0; n < length; n++) {
        printf("%d\t%d\t%d\t%d\t%f\n", n, cells.i[n], cells.j[n],         
            cells.level[n], cells.values[n]);
    }
}

#ifndef USE_MACROS
uint key_to_i (uint key, uint lev) {
    return (key % two_to_the(lev)); 
}

uint key_to_j (uint key, uint lev) {
    return (key / two_to_the(lev));
}

uint truncate_base (uint val, uint val2) {
    return (val/val2) * val2;
}

#ifdef USE_ASSERT
uint four_to_the (int val) {
    assert(val >=0);
    return (1 << (val*2));
}
#else
uint four_to_the (int val) {
    if (val >= 0) {
        return (1 << (val*2));
    } else {
        perror("val is negative");
        exit(-1);
    }
}
#endif

#ifdef USE_ASSERT
int two_to_the (int val) {
    assert(val >=0);
    return (1 << val);
}
#else
int two_to_the (int val) {
    if (val >= 0) {
        return (1 << val);
    } else {
        perror("val is negative");
        exit(-1);
        return -1;
    }
}
#endif

#endif

// adaptiveMeshConstructor()
// Inputs: n (width/height of the square mesh), l (maximum level of refinement),
//         pointers for the level, x, and y arrays (should be NULL for all three)
// Output: number of cells in the adaptive mesh
//
cell_list adaptiveMeshConstructorWij(cell_list icells, const int n, const int levmax, float threshold, int target_ncells) {
  int ncells = SQ(n);

  // ints used for for() loops later
  int ic, xc, yc, xlc, ylc, nlc;

  //printf("\nBuilding the mesh...\n");

  // Initialize Coarse Mesh
  uint*  level = (uint*)  malloc(sizeof(uint)*ncells);
  uint*  i     = (uint*)  malloc(sizeof(uint)*ncells);
  uint*  j     = (uint*)  malloc(sizeof(uint)*ncells);
  for(yc = 0; yc < n; yc++) {
    for(xc = 0; xc < n; xc++) {
      level[n*yc+xc] = 0;
      i[n*yc+xc]     = xc;
      j[n*yc+xc]     = yc;
    }
  }
  //printf("Coarse mesh initialized.\n");

  // Randomly Set Level of Refinement
  //unsigned int iseed = (unsigned int)time(NULL);
  //srand (iseed);
  //srand (0);
  for(int ii = levmax; ii >= 0; ii--) {
    float lev_threshold = threshold*(float)ii/(float)levmax;
    for(ic = 0; ic < ncells; ic++) {
      float jj = (100.0*(float)rand() / ((float)RAND_MAX));
      if(jj<lev_threshold && level[ic] == 0) level[ic] = ii;
    }
  }

  //printf("Levels of refinement randomly set.\n");

  // Smooth the Refinement
  int newcount = -1;
  while(newcount != 0) {
    newcount = 0;
    uint lev = 0;
    for(ic = 0; ic < ncells; ic++) {
      lev = level[ic];
      lev++;
      // Check bottom neighbor
      if(ic - n >= 0) {
        if(level[ic-n] > lev) {
          level[ic] = lev;
          newcount++;
          continue;
        }
      }
      // Check top neighbor
      if(ic + n < ncells) {
        if(level[ic+n] > lev) {
          level[ic] = lev;
          newcount++;
          continue;
        }
      }
      // Check left neighbor
      if((ic%n)-1 >= 0) {
        if(level[ic-1] > lev) {
          level[ic] = lev;
          newcount++;
          continue;
        }
      }
      // Check right neighbor
      if((ic%n)+1 < n) {
        if(level[ic+1] > lev) {
          level[ic] = lev;
          newcount++;
          continue;
        }
      }
    }
  }

  //printf("\nDEBUG -- ncells %d target_ncells %ld fine mesh size %ld\n",ncells,target_ncells,n*two_to_the(levmax)*n*two_to_the(levmax));
  if (target_ncells > ncells && target_ncells < n*two_to_the(levmax)*n*two_to_the(levmax)) {
    int icount = 0;
    int newcount = 0;
    for(ic = 0; ic < ncells; ic++) {newcount += (powerOfFour(level[ic]) - 1);}

    while ( abs((ncells+newcount) - target_ncells) > MAX(5,target_ncells/10000) && icount < 40) {
      icount++;
      //printf("DEBUG -- Adjusting cell count %ld target %ld diff %ld\n",ncells+newcount, target_ncells, abs(ncells+newcount - target_ncells));

      if (ncells+newcount > target_ncells){
        int reduce_count = ((ncells+newcount) - target_ncells);
        //printf("DEBUG -- Too many cells -- need to reduce by %d\n",reduce_count);
        int jcount = 0;
        while (reduce_count > 0 && jcount < ncells) {
          int jj = 1 + (int)((float)ncells*rand() / (RAND_MAX+1.0));
          if(jj>0 && jj<ncells && level[jj] > 0) {
             reduce_count-=4;
          //   printf("DEBUG reducing level for ic %d level %d reduce_count %d jj %d\n",jj,level[jj],reduce_count,jj);
             level[jj]--;
          }
          jcount++;
        }
      } else {
        int increase_count = (target_ncells - (ncells+newcount));
        increase_count /= (levmax*4);
        //printf("DEBUG -- Too few cells -- need to increase by %d\n",increase_count);
        int jcount = 0;
        while (increase_count > 0 && jcount < ncells) {
          int jj = 1 + (int)((float)ncells*rand() / (RAND_MAX+1.0));
          if(jj>0 && jj<ncells && level[jj] < levmax) {
            increase_count-=4;
            level[jj]++;
          }
          jcount++;
        }
      }

      // Smooth the Refinement
      newcount = -1;
      while(newcount != 0) {
        newcount = 0;
        uint lev = 0;
        for(ic = 0; ic < ncells; ic++) {
          lev = level[ic];
          lev++;
          // Check bottom neighbor
          if(ic - n >= 0) {
            if(level[ic-n] > lev) {
              level[ic] = lev;
              newcount++;
              continue;
            }
          }
          // Check top neighbor
          if(ic + n < ncells) {
            if(level[ic+n] > lev) {
              level[ic] = lev;
              newcount++;
              continue;
            }
          }
          // Check left neighbor
          if((ic%n)-1 >= 0) {
            if(level[ic-1] > lev) {
              level[ic] = lev;
              newcount++;
              continue;
            }
          }
          // Check right neighbor
          if((ic%n)+1 < n) {
            if(level[ic+1] > lev) {
              level[ic] = lev;
              newcount++;
              continue;
            }
          }
        }
      } // while(newcount != 0) {
      newcount = 0;
      for(ic = 0; ic < ncells; ic++) {newcount += (powerOfFour(level[ic]) - 1);}
    } // while ( abs(ncells+newcount - target_ncells) > 10 && icount < 10) {

  } //if (target_ncells > 0) {

  //printf("Refinement smoothed.\n");
  int small_cells = 0;
  for(ic = 0; ic < ncells; ic++) {
    if (level[ic] == (uint)levmax) {
      small_cells++;
    }
  }
  //printf("%8d small cells, ", small_cells);
  
  // Allocate Space for the Adaptive Mesh
  newcount = 0;
  for(ic = 0; ic < ncells; ic++) {newcount += (powerOfFour(level[ic]) - 1);}
  //printf("DEBUG -- Exiting cell adjustment with ncells %ld target %ld diff %ld\n\n",ncells+newcount, target_ncells, abs(ncells+newcount - target_ncells));

  uint*  level_temp = (uint*)  malloc(sizeof(uint)*(ncells+newcount));
  uint*  i_temp     = (uint*)  malloc(sizeof(uint)*(ncells+newcount));
  uint*  j_temp     = (uint*)  malloc(sizeof(uint)*(ncells+newcount));

  // Set the Adaptive Mesh
  int offset = 0;
  for(yc = 0; yc < n; yc++) {
    for(xc = 0; xc < n; xc++) {
      ic = n*yc + xc;
      nlc = (int) sqrt( (double) powerOfFour(level[ic]) );
      for(ylc = 0; ylc < nlc; ylc++) {
        for(xlc = 0; xlc < nlc; xlc++) {
          level_temp[ic + offset + (nlc*ylc + xlc)] = level[ic];
          i_temp[ic + offset + (nlc*ylc + xlc)] = i[ic]*pow(2,level[ic]) + xlc;
          j_temp[ic + offset + (nlc*ylc + xlc)] = j[ic]*pow(2,level[ic]) + ylc;
        }         
      }
      offset += powerOfFour(level[ic])-1;
    }
  }
  //printf("Adaptive mesh built.\n");

  // Swap pointers and free memory used by Coarse Mesh
  swap_uint(&level, &level_temp);
  swap_uint(&i, &i_temp);
  swap_uint(&j, &j_temp);
  free(level_temp);
  free(i_temp);
  free(j_temp);

  //printf("Old ncells: %d", ncells);
  // Update ncells
  ncells += newcount;
  //printf("\tNew ncells: %d\n", ncells);

  if (randomize) {
    // Randomize the order of the arrays

    int* random = (int*) malloc(sizeof(int)*ncells);
    uint* temp1 = (uint*) malloc(sizeof(uint)*ncells);
    uint* temp2 = (uint*) malloc(sizeof(uint)*ncells*2);
    // XXX Want better randomization? XXX
    // XXX Why is the time between printf() statements the longest part? XXX
    //printf("Shuffling");
    //fflush(stdout);
    for(ic = 0; ic < ncells; ic++) {random[ic] = ic;}
    //iseed = (unsigned int)time(NULL);
    //srand (iseed);
    srand(0);
    nlc = 0;
    for(int ii = 0; ii < 7; ii++) {
      for(ic = 0; ic < ncells; ic++) {
    
        int jj = (int)( (double)ncells*((double)rand() / (double)(RAND_MAX+1.0) ) );
        // occasionally jj will be ncells and random ratio is 1.0
        if (jj >= ncells) jj=ncells-1;
        nlc = random[jj];
        random[jj] = random[ic];
        random[ic] = nlc;
         if (random[ic] >= ncells) {
           exit(0);
         }
      }
      //printf(".");
      //fflush(stdout);
    }
    //printf("\n");

    for(ic = 0; ic < ncells; ic++) {
      temp1[ic] = level[random[ic]];
      temp2[2*ic] = i[random[ic]];
      temp2[2*ic+1] = j[random[ic]];
    }
    for(ic = 0; ic < ncells; ic++) {
      level[ic] = temp1[ic];
      i[ic]     = temp2[2*ic];
      j[ic]     = temp2[2*ic+1];
    }

    free(temp1);
    free(temp2);
    free(random);
    //printf("Adaptive mesh randomized.\n");
  } // End of if randomize

  icells.ncells = ncells;
  icells.ibasesize = n;
  icells.jbasesize = n;
  icells.levmax = levmax;
  icells.i = i;
  icells.j = j;
  icells.level = level;

  double *values = (double *)malloc(sizeof(double) * ncells);
  icells.values = values;

  //printf("Adaptive mesh construction complete.\n");
  return icells;
}
