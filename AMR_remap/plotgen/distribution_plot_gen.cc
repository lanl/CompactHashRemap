#include "distribution_plot_gen.h"
#include <time.h>


typedef unsigned int uint;

uint basesize = 2;
uint levmax = 6;
float adapt_threshhold = 20.0f;
uint numcells = 16399;
uint num_runs = 1;
uint cell_inc = 3;
float threshhold_inc = 0.0;

long seed = 0xDEADBEEF;

uint output_mode = 0;
uint print_mode = 0;
uint data_mode = 0;
uint adapt_meshgen = 0;
uint force_seed = 1;

double sparsity = 0.1;
 

int main (int argc, char** argv){
    if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)){
        printf ("Usage: %s [--help | -h] [-output, -print-levels, -adapt-meshgen, -ncells <ncells>, -lev <level difference>, adapt-threshhold <adaptive threshhold>, -base-size <base size>, -no-force-seed, -seed <seed>, -cell-inc <cell increment size>, -num-runs <number of runs>, -sparsity <sparsity>]\n", argv[0]);
        printf ("Any variables not set will be run with reasonable defaults\n");
        return 0;
    }
    if (argc>1){
        for (int i = 1; i < argc; i++){
            char* arg = argv[i];
            if (strcmp(arg,"-output")==0){
                output_mode = 1;
            } else
            if (strcmp(arg,"-print-levels")==0){
                print_mode = 1;
            } else
            if (strcmp(arg,"-data-format")==0){
                data_mode = 1;
            } else
            if (strcmp(arg,"-adapt-meshgen")==0){
                adapt_meshgen = 1;
            } else
            if (strcmp(arg,"-threshhold-inc")==0){
                i++;
                threshhold_inc = atof(argv[i]);
            } else
            if (strcmp(arg,"-ncells")==0){
                i++;
                numcells = atoi(argv[i]);
                if ((numcells - 1) % 3 != 0) {
                    numcells--;
                    numcells /= 3;
                    numcells *= 3;
                    numcells ++;
                    //printf("\nImpossible number of cells, using %u instead\n", numcells);
                }
            } else
            if (strcmp(arg,"-lev")==0){
                i++;
                levmax = atoi(argv[i]);
            } else
            if (strcmp(arg,"-adapt-threshhold")==0){
                i++;
                adapt_threshhold = atof(argv[i]);
            } else
            if (strcmp(arg,"-base-size")==0){
                i++;
                basesize = atoi(argv[i]);
            } else
            if (strcmp(arg,"-no-force-seed")==0){
                force_seed = 0;
            } else
            if (strcmp(arg,"-seed")==0){
                i++;
                seed = atol(argv[i]);
            } else
            if (strcmp(arg,"-cell-inc")==0){
                i++;
                cell_inc = atoi(argv[i]);
                if ((cell_inc) % 3 != 0) {
                    cell_inc /= 3;
                    cell_inc *= 3;
                    printf("\nImpossible increment, using %u instead\n", cell_inc);
                }
            } else
            if (strcmp(arg,"-num-runs")==0){
                i++;
                num_runs = atoi(argv[i]);
            } else
            if (strcmp(arg,"-sparsity")==0){
                i++;
                sparsity = atof(argv[i]);
            } else
            printf ("Invalid Argument: %s\n", arg);
        }
    }
    cell_list ocells;
    uint levmin = 0;
    uint* olev_count;
    if (!force_seed){
        srand ((unsigned)time(NULL));
    }
    for (uint run_num = 0; run_num<num_runs; run_num++){
        if (force_seed) {
            srand(seed);
        }
        if (adapt_meshgen){
            ocells = adaptiveMeshConstructorWij(ocells, basesize, levmax, adapt_threshhold, numcells);
            if (print_mode){
                printf ("Adapt-meshgen: %u cells.\n", ocells.ncells);
                if (threshhold_inc!=0){
                    printf ("Threshhold = %f\n", adapt_threshhold);
                }
                olev_count = (uint*)malloc(sizeof(uint)*(ocells.levmax));
                for (uint i = 0; i < ocells.levmax; i++){
                    olev_count[i-levmin] = 0;
                }
                for (uint i = 0; i < ocells.ncells; i++){
                    olev_count[ocells.level[i]]++;
                }
                for (uint i = 0; i < ocells.levmax; i++){
                    printf ("lev %u: %u\n", i-levmin, olev_count[i]);
                }
                free (olev_count);
            }
            if (output_mode){
                PrintMesh(ocells);
            }
            destroy(ocells);
            
        }else{
        
            uint ilength = numcells;
            uint i_max_level;
            ocells = mesh_maker(ocells, levmax, &ilength, &i_max_level, sparsity, basesize);
            numcells = ocells.ncells;
            ocells = shuffle_cell_list(ocells, ilength*0xFF);
            //printf("Max lev: %u\n", i_max_level);
            levmin = 0;
            //levmin = ocells.level[0];
            /*for (uint i = 1; i < ocells.ncells; i++){
                if (ocells.level[i]<levmin){
                    levmin = ocells.level[i];
                }
            }*/
            if (print_mode) {
                printf ("Levelbased-meshgen: %u cells.\n", ocells.ncells);
                olev_count = (uint*)malloc(sizeof(uint)*(ocells.levmax+1));
                
                for (uint i = levmin; i <= ocells.levmax; i++){
                    olev_count[i] = 0;
                }
                
                for (uint i = 0; i < ocells.ncells; i++){
                    olev_count[ocells.level[i]]++;
                }
                for (uint i = levmin; i <= ocells.levmax; i++){
                    printf ("lev %u: %u\n", i-levmin, olev_count[i]);
                }
            }
            
            if (data_mode) {
                olev_count = (uint*)malloc(sizeof(uint)*(ocells.levmax+1));
                
                for (uint i = levmin; i <= ocells.levmax; i++){
                    olev_count[i] = 0;
                }
                
                for (uint i = 0; i < ocells.ncells; i++){
                    olev_count[ocells.level[i]]++;
                }
                printf("%u ", ocells.ibasesize);
                for (uint i = levmin; i <= ocells.levmax; i++){
                    printf ("%u ", olev_count[i]);
                }
            }
          
            if (!output_mode)
                printf ("\n");
            if (output_mode)
                PrintMesh(ocells);
            destroy(ocells);
            if (print_mode)
                free (olev_count);
        }
    numcells += cell_inc;
    adapt_threshhold+=threshhold_inc;
    }
    return 0;
}

// Prints the mesh as described by the cell list as if it were in a perfect hash.
void PrintMesh (cell_list icells){
    //if (two_to_the(icells.levmax)*icells.ibasesize>32){
    //    printf("cell list too large! %u level %u base\n", icells.levmax, icells.ibasesize);
    //    return;
    //}
    
    uint* hash = (uint*)malloc (sizeof(uint)*icells.ibasesize*icells.ibasesize*four_to_the(icells.levmax));
    // initialize the hash to aid in checking for errors
    for (uint i = 0; i < icells.ibasesize*four_to_the(icells.levmax); i++){
        hash[i] = 0;
    }
    
    
    // the size across of the finest level of the mesh
    uint fine_size = two_to_the(icells.levmax)*icells.ibasesize;
    
    //printf ("finesize: %u from levmax %u and basesize %u\n", fine_size, icells.levmax, icells.ibasesize);
    
    // Add items to the hash
    for (uint ic = 0; ic < icells.ncells; ic++){
        uint lev = icells.level[ic];
        uint i = icells.i[ic];
        uint j = icells.j[ic];
        
        uint lev_mod = two_to_the(icells.levmax - lev);
        // Calculate the position of the cell on the lowest level of the mesh
        uint base_key = i*lev_mod + j*fine_size*lev_mod;
        // 2 to the lev_mod power is four to the difference
        for (uint jc = 0; jc < (uint)four_to_the(icells.levmax - lev); jc++){
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
