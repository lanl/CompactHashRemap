cmake_minimum_required(VERSION 2.6)

project (CompactHashLib)
# The version number.
set (COMPACTHASHLIB_VERSION_MAJOR 2)
set (COMPACTHASHLIB_VERSION_MINOR 0.7)
set (COMPACTHASHLIB_VERSION "${COMPACTHASHLIB_VERSION_MAJOR}.${COMPACTHASHLIB_VERSION_MINOR}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

if (NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE RelWithDebInfo)
   set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

#find_program(BF_CLANG bf-clang)

#message("CMAKE_C_COMPILER FOUND ${CMAKE_C_COMPILER}")

#if(BF_CLANG AND "${CMAKE_C_COMPILER}" MATCHES "bf-clang")
#   #message("BF_CLANG FOUND ${BF_CLANG}")
#   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I$LLVM/include")
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I$LLVM/include")
#   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bf-include=_Z14remap_kDtree2d9cell_listS_")
#   set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -bf-include=_Z17singlewrite_remap9cell_listS_")
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bf-include=_Z17singlewrite_remap9cell_listS_")
#endif()

# uncomment to build debug
#set(CMAKE_BUILD_TYPE Debug)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -std=c99")
   #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C
#  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -vec-report3 -std=c99")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
endif()

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELDEBINFO} -g -O3")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   # using Clang
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   # using GCC
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize")
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=5")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
   # using Intel C++
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -vec-report3")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
   # using Visual Studio C++
endif()

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELDEBINFO} -g -O3")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0 -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELDEBINFO} -g -O3")

if (MIC_NATIVE)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static-intel -mmic")
   link_directories("/projects/opt/intel/compilers/parallel_studio_xe_2013/composer_xe_2013_sp1.2.144/compiler/lib/mic")
endif (MIC_NATIVE)

#message("CMAKE_C_FLAGS ${CMAKE_C_FLAGS}")

#if (NOT MIC_NATIVE)
   find_package(OpenCL)
   #  OpenCL_FOUND          - True if OpenCL was found
   #  OpenCL_INCLUDE_DIRS   - include directories for OpenCL
   #  OpenCL_LIBRARIES      - link against this library to use OpenCL
   #  OpenCL_VERSION_STRING - Highest supported OpenCL version (eg. 1.2)
   #  OpenCL_VERSION_MAJOR  - The major version of the OpenCL implementation
   #  OpenCL_VERSION_MINOR  - The minor version of the OpenCL implementation
   if (OpenCL_FOUND)
      include_directories(${OpenCL_INCLUDE_DIRS})
   endif (OpenCL_FOUND)
#endif (NOT MIC_NATIVE)

#find_package(MPI)
#if (MPI_FOUND)
#   include_directories(${MPI_INCLUDE_PATH})
#endif (MPI_FOUND)

#if (BF_CLANG-NOTFOUND)
   if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.6.0)
     cmake_policy(PUSH)
     # Suppress warnings in FindOpenMP
     cmake_policy(SET CMP0012 NEW)
   endif (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.6.0)
   find_package( OpenMP)
   if (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.6.0)
     cmake_policy(POP)
   endif (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} VERSION_GREATER 3.6.0)
#endif()

add_subdirectory(AMR_remap)
add_subdirectory(Unstruct_remap)

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "")
