set(libsimplehash_LIB_SRCS simplehash.cc simplehash.h)

########### embed source target ###############
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/simplehashlib_kern.inc
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl
                          ${CMAKE_CURRENT_SOURCE_DIR}/simplehashlib_kern.cl > simplehashlib_kern.inc
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/simplehashlib_kern.cl
                          ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl)
add_custom_target(simplehashlib_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/simplehashlib_kern.inc)

add_custom_command(OUTPUT
${CMAKE_CURRENT_BINARY_DIR}/simplehash_kern.inc
                  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl
                          ${CMAKE_CURRENT_SOURCE_DIR}/simplehash_kern.cl > simplehash_kern.inc
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/simplehash_kern.cl
                          ${CMAKE_CURRENT_SOURCE_DIR}/embed_source.pl)
add_custom_target(simplehash_kernel_source ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/simplehash_kern.inc)

add_library(simplehash STATIC ${libsimplehash_LIB_SRCS})
set_target_properties(simplehash PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/.. -I.")
target_link_libraries(simplehash genmalloc)
if (OpenCL_FOUND)
    target_link_libraries(simplehash ezcl)
    set_target_properties(simplehash PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
endif (OpenCL_FOUND)
add_dependencies(simplehash simplehash_kernel_source)
add_dependencies(simplehash simplehashlib_kernel_source)

if (OPENMP_FOUND)
    add_library(simplehash_openmp STATIC ${libsimplehash_LIB_SRCS})
    set_target_properties(simplehash_openmp PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/.. -I. ${OpenMP_C_FLAGS}")
#   set_target_properties(simplehash_openmp PROPERTIES COMPILE_FLAGS "${OpenMP_C_FLAGS}")
#   set_target_properties(simplehash_openmp PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")

    target_link_libraries(simplehash_openmp genmalloc)
    if (OpenCL_FOUND)
        target_link_libraries(simplehash_openmp ezcl)
        set_target_properties(simplehash_openmp PROPERTIES COMPILE_DEFINITIONS "HAVE_OPENCL")
    endif (OpenCL_FOUND)
    add_dependencies(simplehash_openmp simplehash_kernel_source)
    add_dependencies(simplehash_openmp simplehashlib_kernel_source)
endif (OPENMP_FOUND)


