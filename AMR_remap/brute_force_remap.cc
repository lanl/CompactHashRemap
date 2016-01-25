/* LANL Copyright 2015
 *
 *
 *  Authors: 
 *        Gerald Collom        gcollom@lanl.gov
 *        Colin Redman        credman@lanl.gov
 */

#include "meshgen/meshgen.h"
#include "brute_force_remap.h"

void brute_force_remap (cell_list icells, cell_list ocells) {
    
    for(uint o = 0; o<ocells.ncells; o++){
        ocells.values[o] = 0;
    }
    
    for (uint o = 0; o < ocells.ncells; o++){
        for (uint i = 0; i < icells.ncells; i++){
            uint lev_out = ocells.level[o];
            uint lev_in = icells.level[i];
            uint i_out = ocells.i[o];
            uint j_out = ocells.j[o];
            uint i_in = icells.i[i];
            uint j_in = icells.j[i];
            if (lev_out > lev_in) {
                uint key_in = j_in * icells.ibasesize*two_to_the(lev_in) + i_in;
                uint key_out = translate_cell(i_out, j_out, lev_out, lev_in, ocells.ibasesize);
                if (key_in == key_out){
                    ocells.values[o] = icells.values[i];
                    continue;
                }
            }else{
                uint key_out = j_out * ocells.ibasesize*two_to_the(lev_out) + i_out;
                uint key_in = translate_cell(i_in, j_in, lev_in, lev_out, icells.ibasesize);
                if (key_out == key_in){
                    ocells.values[o] += (icells.values[i] / four_to_the(lev_in - lev_out));
                }
            }
        }
    }
}

