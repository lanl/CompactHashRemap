cmake_minimum_required(VERSION 2.6)

add_subdirectory(genmalloc)
add_subdirectory(meshgen)
if (OPENCL_FOUND)
   add_subdirectory(ezcl)
endif (OPENCL_FOUND)
add_subdirectory(simplehash)
add_subdirectory(HashFactory)
add_subdirectory(KDTree)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -bf-include=_Z17brute_force_remap9cell_listS_")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bf-include=_Z17brute_force_remap9cell_listS_")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -bf-include=_Z14remap_kDtree2d9cell_listS_")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bf-include=_Z14remap_kDtree2d9cell_listS_")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -bf-include=_Z18full_perfect_remap9cell_listS_")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bf-include=_Z18full_perfect_remap9cell_listS_")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -bf-include=_Z17singlewrite_remap9cell_listS_")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bf-include=_Z17singlewrite_remap9cell_listS_")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -bf-include=_Z25singlewrite_remap_compact9cell_listS_")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bf-include=_Z25singlewrite_remap_compact9cell_listS_")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -bf-include=_Z7h_remap9cell_listS_")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bf-include=_Z7h_remap9cell_listS_")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -bf-include=_Z15h_remap_compact9cell_listS_P19intintHash_Factory_")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -bf-include=_Z15h_remap_compact9cell_listS_P19intintHash_Factory_")

########### global settings ###############

set(AMR_REMAP_SRCS AMR_remap.cc h_remap_gpu.cc hierarchical_remap.cc kdtree_remap.cc
  singlewrite_remap.cc full_perfect_remap.cc brute_force_remap.cc timer.cc)
set(AMR_REMAP_HDRS AMR_remap.h h_remap_gpu.h hierarchical_remap.h kdtree_remap.h
  singlewrite_remap.h full_perfect_remap.h brute_force_remap.h timer.h)

include_directories(.)
########### embed source target ############

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/h_remap_kern.inc
                   COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl 
                      ${CMAKE_CURRENT_SOURCE_DIR}/h_remap_kern.cl > 
                      ${CMAKE_CURRENT_BINARY_DIR}/h_remap_kern.inc
                   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/h_remap_kern.cl ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl)
add_custom_target(h_remap_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/h_remap_kern.inc)

########### AMR_remap target ###############

add_executable(AMR_remap ${AMR_REMAP_SRCS} ${AMR_REMAP_HDRS})
target_link_libraries(AMR_remap genmalloc meshgen simplehash HashFactory kdtree)
if (OPENCL_FOUND)
   set_target_properties(AMR_remap PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
   target_link_libraries(AMR_remap ezcl)
   target_link_libraries(AMR_remap ${OPENCL_LIBRARIES})
   add_dependencies(AMR_remap h_remap_kernel_source)
endif(OPENCL_FOUND)

########### AMR_remap_openMP target ###############
if (OPENMP_FOUND)
   add_executable(AMR_remap_openMP ${AMR_REMAP_SRCS} ${AMR_REMAP_HDRS})

   set_target_properties(AMR_remap_openMP PROPERTIES COMPILE_FLAGS "-I.")
   set_target_properties(AMR_remap_openMP PROPERTIES COMPILE_FLAGS "${OpenMP_C_FLAGS}")
   set_target_properties(AMR_remap_openMP PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")
   target_link_libraries(AMR_remap_openMP genmalloc meshgen simplehash_openmp HashFactory_openmp kdtree)
   if (OPENCL_FOUND)
      set_target_properties(AMR_remap_openMP PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
      target_link_libraries(AMR_remap_openMP ezcl)
      target_link_libraries(AMR_remap_openMP ${OPENCL_LIBRARIES})
      add_dependencies(AMR_remap h_remap_kernel_source)
   endif(OPENCL_FOUND)
endif (OPENMP_FOUND)

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "")

