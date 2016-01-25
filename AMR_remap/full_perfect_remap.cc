/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */

#include "meshgen/meshgen.h"
#include "genmalloc/genmalloc.h"
#include "full_perfect_remap.h"

void full_perfect_remap (cell_list icells, cell_list ocells) {

    // Allocate a hash table the size of the finest level of the grid
    uint i_max = icells.ibasesize*two_to_the(icells.levmax);
    uint j_max = icells.jbasesize*two_to_the(icells.levmax);
    int **hash = (int **)genmatrix(j_max, i_max, sizeof(int));

    // Fill Hash Table from Input mesh
    for (uint ic = 0; ic < icells.ncells; ic++){
        uint lev = icells.level[ic];
        uint i = icells.i[ic];
        uint j = icells.j[ic];
        // If at the maximum level just set the one cell
        if (lev == icells.levmax) {
            hash[j][i] = ic;
        } else {
            // Set the square block of cells at the finest level
            // to the index number
            uint lev_mod = two_to_the(icells.levmax - lev);
            for (uint jj = j*lev_mod; jj < (j+1)*lev_mod; jj++) {
                for (uint ii = i*lev_mod; ii < (i+1)*lev_mod; ii++) {
                    hash[jj][ii] = ic;
                }
            }
        }
    }

    // Use Hash Table to Perform Remap
    for (uint ic = 0; ic < ocells.ncells; ic++){
        uint lev = ocells.level[ic];
        uint i = ocells.i[ic];
        uint j = ocells.j[ic];

        // If at the finest level, get the index number and
        // get the value of the input mesh at that index
        if (lev == ocells.levmax) {
            ocells.values[ic] = icells.values[hash[j][i]];
        } else {
            // Sum up the values in the underlying block of
            // cells at the finest level and average
            uint lev_mod = two_to_the(ocells.levmax - lev);
            ocells.values[ic] = 0.0;
            for (uint jj = j*lev_mod; jj < (j+1)*lev_mod; jj++) {
                for (uint ii = i*lev_mod; ii < (i+1)*lev_mod; ii++) {
                    ocells.values[ic] += icells.values[hash[jj][ii]];
                }
            }
            // Get average by dividing by number of cells
            ocells.values[ic] /= (double)(lev_mod*lev_mod);
        }
    }

    // Deallocate hash table
    genmatrixfree((void **)hash);
}

#ifdef _OPENMP
void full_perfect_remap_openMP (cell_list icells, cell_list ocells) {

    // Allocate a hash table the size of the finest level of the grid
    uint i_max = icells.ibasesize*two_to_the(icells.levmax);
    uint j_max = icells.jbasesize*two_to_the(icells.levmax);
    int **hash = (int **)genmatrix(j_max, i_max, sizeof(int));

#pragma omp parallel default(none) shared(icells, ocells, hash)
    {
        uint ilength = icells.ncells;
        uint olength = ocells.ncells;
        uint max_lev = icells.levmax;

        // Fill Hash Table from Input mesh
#pragma omp for
        for (uint ic = 0; ic < ilength; ic++){
            uint lev = icells.level[ic];
            uint i = icells.i[ic];
            uint j = icells.j[ic];
            // If at the maximum level just set the one cell
            if (lev == max_lev) {
                hash[j][i] = ic;
            } else {
                // Set the square block of cells at the finest level
                // to the index number
                uint lev_mod = two_to_the(max_lev - lev);
                for (uint jj = j*lev_mod; jj < (j+1)*lev_mod; jj++) {
                    for (uint ii = i*lev_mod; ii < (i+1)*lev_mod; ii++) {
                        hash[jj][ii] = ic;
                    }
                }
            }
        }

    // Use Hash Table to Perform Remap
#pragma omp for
        for (uint ic = 0; ic < olength; ic++){
            uint lev = ocells.level[ic];
            uint i = ocells.i[ic];
            uint j = ocells.j[ic];

            // If at the finest level, get the index number and
            // get the value of the input mesh at that index
            if (lev == max_lev) {
                ocells.values[ic] = icells.values[hash[j][i]];
            } else {
                // Sum up the values in the underlying block of
                // cells at the finest level and average
                uint lev_mod = two_to_the(max_lev - lev);
                ocells.values[ic] = 0.0;
                for (uint jj = j*lev_mod; jj < (j+1)*lev_mod; jj++) {
                    for (uint ii = i*lev_mod; ii < (i+1)*lev_mod; ii++) {
                        ocells.values[ic] += icells.values[hash[jj][ii]];
                    }
                }
                // Get average by dividing by number of cells
                ocells.values[ic] /= (double)(lev_mod*lev_mod);
            }
        }
    }

    // Deallocate hash table
    genmatrixfree((void **)hash);
}
#endif
