cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
if (NOT TARGET distclean)
   INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)
endif (NOT TARGET distclean)

# enable @rpath in the install name for any shared library being built
set(CMAKE_MACOSX_RPATH 1)

project (EZCL)
# The version number.
set (EZCL_VERSION_MAJOR 2)
set (EZCL_VERSION_MINOR 0.7)
set (EZCL_VERSION "${EZCL_VERSION_MAJOR}.${EZCL_VERSION_MINOR}")

include_directories(${PROJECT_SOURCE_DIR})

########### global settings ###############
if (OpenCL_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DOPENCL_VERSION_MAJOR=${OpenCL_VERSION_MAJOR} -DOPENCL_VERSION_MINOR=${OpenCL_VERSION_MINOR} -I${OpenCL_INCLUDE_DIRS}")
else (OpenCL_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I${OpenCL_INCLUDE_DIRS}")
endif (OpenCL_FOUND)

set(H_SRCS ezcl.h)
set(C_SRCS ezcl.c)

########### ezcl target ###############
set(ezcl_LIB_SRCS ${C_SRCS} ${H_SRCS})

add_library(ezcl SHARED ${ezcl_LIB_SRCS})

set_target_properties(ezcl PROPERTIES VERSION ${EZCL_VERSION} SOVERSION ${EZCL_VERSION_MAJOR})
target_link_libraries(ezcl ${OpenCL_LIBRARIES})
install(TARGETS ezcl DESTINATION lib)

########### install files ###############

install(FILES  ezcl.h DESTINATION include/ezcl)

########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "")
