# CompactHashRemap
Fast mesh remapping algorithm based on hashing techniques

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

To checkout the code:

   git clone git@github.com:losalamos/CompactHashRemap

Common commands for code developement

   git add <file>  // Adding a new file or change to local repository
   git commit -a   // Committing all changes to local repository
   git pull        // Pulling in changes from remote respository
   git push        // Pushing changes to remote repository

To update from remote repository when local branch has modifications
   git stash
   git pull
   git stash pop

To setup the environment:

   module commands -- need list
   module load cmake/3.0.2
   //find MPI command

The Nvidia Kepler K40 with a 32-cpu host node is the Darwin node for timing.

salloc -p SAND_16_2:MEM64.0:Accelerator-1:TeslaK40c:10GbE::

Options to cmake for the configure step is
   cmake .                         // in-tree build
   cmake <path-to-src>             // out-of-tree build
   cmake -DCMAKE_BUILD_TYPE=debug <path-to-src>
   cmake -DCMAKE_BUILD_TYPE=release <path-to-src>

   cmake . -DCMAKE_C_COMPILER=bf-clang -DCMAKE_CXX_COMPILER=bf-clang++ // for building with byfl

Build the executables:

   make

Clean

   make clean

Running the code

   cd into the AMR_remap directory
   
   ./AMR_remap <i_level_diff> <ilength> <o_level_diff> <olength> <num_rep>
   
   i_level_diff and o_level_diff are the range of levels of refinement for input(i) and output(o)
   ilength and olength are the number of cells within each mesh. If the number of cells is impossible to obtain, a
   possible number of cells is set instead of the argument.
   num_rep is the number of times the code is repeated to get more rounded averages.
   sample: 2 10 2 10 10
   sample: 8 50000 9 200000 1

   ./AMR_remap_openMP 6 1000000 6 1000000 10 -no-brute
   
   An alternate mesh generation is available that can better explore more realistic
   AMR meshes. To run it add -adapt-meshgen and the input options are the size of the base mesh, the number of
   levels of refinement and the refinement probability threshold in percent. A blank field shown as 0 is also
   necessary. The number of repetitions remains in the fifth argument.

   Usage -- ./AMR_remap_openMP <size_base_mesh> <levmax> <refine_threshold> 0 <num_rep> -adapt-meshgen [-no-brute,-no-test,-no-tree]

   A nice example for a test case would be:

   ./AMR_remap_openMP 128 6 20 0 10 -adapt-meshgen -no-brute

   cd into the Unstruct_remap directory
   
   ./parse_test
   ./read_test
   ./write_test
   ./search_test
   
   each of the above will test the part of the unstructured search that is specified. no arguments are required.


