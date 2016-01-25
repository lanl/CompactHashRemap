########### global settings ###############

set(UNSTRUCTURED_REMAP_SRCS unstructured_search.cc unstructured_parse.cc unstructured_write.cc unstructured_read.cc
	local_hash.cc)
set(UNSTRUCTURED_REMAP_HDRS unstructured_types.h unstructured_search.h unstructured_parse.h unstructured_write.h
	unstructured_read.h local_hash.h)

########### unstructured_remap target ###############

add_library(unstructured_remap ${UNSTRUCTURED_REMAP_SRCS} ${UNSTRUCTURED_REMAP_HDRS})

add_executable(read_test unstructured_read_main.cc)
target_link_libraries(read_test unstructured_remap)

add_executable(write_test unstructured_write_main.cc)
target_link_libraries(write_test unstructured_remap)

add_executable(parse_test unstructured_parse_main.cc)
target_link_libraries(parse_test unstructured_remap)

add_executable(search_test unstructured_search_main.cc)
target_link_libraries(search_test unstructured_remap)


########### clean files ################
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
   "")

