#include "distribution_plot_gen.h"


typedef unsigned int uint;

const uint basesize = 64;
const uint levmax = 6;
const float adapt_threshhold = 20.0f;
uint numcells = 118;
 

int main (int n_args, char** args){
    int output_mode = 0;
    int adapt_meshgen = 0;
    if (n_args == 3){
        if  (strcmp (args[2],"-output-mode")==0){
            output_mode = 1;
        }
    }
    cell_list ocells;
    ocells = adaptiveMeshConstructorWij(ocells, basesize, levmax, adapt_threshhold, numcells);
    uint levmin = 0;
    uint* olev_count;
    /*levmin = ocells.level[0];
    for (uint i = 1; i < ocells.ncells; i++){
        if (ocells.level[i]<levmin){
            levmin = ocells.level[i];
        }
    }*/
    for (uint run_num = 0; run_num<100; run_num++){
        if (adapt_meshgen){
            printf ("Adapt-meshgen: %u cells.\n", ocells.ncells);
            olev_count = (uint*)malloc(sizeof(uint)*(ocells.levmax+1));
            for (uint i = levmin; i <= ocells.levmax; i++){
                olev_count[i-levmin] = 0;
            }
            for (uint i = levmin; i < ocells.ncells; i++){
                olev_count[ocells.level[i-levmin]]++;
            }
            for (uint i = levmin; i <= ocells.levmax; i++){
                printf ("lev %u: %u\n", i-levmin, olev_count[i]);
            }
            numcells = ocells.ncells;
            PrintMesh(ocells);
            destroy(ocells);
            free (olev_count);
        }else{
        
            uint ilength = numcells;
            uint i_max_level, i_min_level;
            ocells = mesh_maker_level(ocells, levmax, &ilength, &i_max_level, &i_min_level);
            levmin = 0;
            /*levmin = ocells.level[0];
            for (uint i = 1; i < ocells.ncells; i++){
                if (ocells.level[i]<levmin){
                    levmin = ocells.level[i];
                }
            }*/
            if (!output_mode)
                printf ("Levelbased-meshgen: %u cells.\n", ocells.ncells);
            olev_count = (uint*)malloc(sizeof(uint)*(ocells.levmax+1));
            for (uint i = levmin; i <= ocells.levmax; i++){
                olev_count[i-levmin] = 0;
            }
            for (uint i = levmin; i < ocells.ncells; i++){
                olev_count[ocells.level[i]]++;
            }
            for (uint i = levmin; i <= ocells.levmax; i++){
                if (!output_mode){
                    printf ("lev %u: %u\n", i-levmin, olev_count[i]);
                }else{
                    printf ("%u ",olev_count[i]);
                }
            }
            if (output_mode)
                printf ("\n");
            if (!output_mode)
                PrintMesh(ocells);
            destroy(ocells);
            free (olev_count);
        }
    numcells += 3;
    srand (0xDEADBEEF);
    }
    return 0;
}

// Prints the mesh as described by the cell list as if it were in a perfect hash.
void PrintMesh (cell_list icells){
    if (two_to_the(icells.levmax)*icells.ibasesize>64){
        printf("cell list too large! %u level %u base\n", icells.levmax, icells.ibasesize);
        return;
    }
    
    uint* hash = (uint*)malloc (sizeof(uint)*icells.ibasesize*icells.ibasesize*four_to_the(icells.levmax));
    // initialize the hash to aid in checking for errors
    for (uint i = 0; i < icells.ibasesize*four_to_the(icells.levmax); i++){
        hash[i] = 0;
    }
    
    
    // the size across of the finest level of the mesh
    uint fine_size = two_to_the(icells.levmax)*icells.ibasesize;
    
    printf ("finesize: %u from levmax %u and basesize %u\n", fine_size, icells.levmax, icells.ibasesize);
    
    // Add items to the hash
    for (uint ic = 0; ic < icells.ncells; ic++){
        uint lev = icells.level[ic];
        uint i = icells.i[ic];
        uint j = icells.j[ic];
        
        uint lev_mod = two_to_the(icells.levmax - lev);
        // Calculate the position of the cell on the lowest level of the mesh
        uint base_key = i*lev_mod + j*fine_size*lev_mod;
        // 2 to the lev_mod power is four to the difference
        for (uint jc = 0; jc < four_to_the(icells.levmax - lev); jc++){
            uint ii = jc%lev_mod;
            uint jj = jc/lev_mod;
            uint key = base_key + ii + jj*fine_size;
            hash [key] = ic;
        }
    }
    
    for (uint ii = 0; ii < fine_size; ii++){
        for (uint jj = 0; jj < fine_size; jj++){    
            uint key = ii + jj*fine_size;
            
            printf ("%u ", hash[key]);
            // if the number has too few digits, add space
            if(hash[key]/100 == 0){
                printf (" ");
            }
            if(hash[key]/10 == 0){
                printf (" ");
            }
        }
        printf ("\n");
    }
    free(hash);
    
    
}
