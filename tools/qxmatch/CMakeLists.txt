cmake_minimum_required(VERSION 2.6)
project(qxmatch)

include_directories(${CFITSIO_INCLUDES})
include_directories(${CCFITS_INCLUDES})

add_executable(qxmatch qxmatch.cpp)

target_link_libraries(qxmatch ${CCFITS_LIBRARIES})
target_link_libraries(qxmatch ${CFITSIO_LIBRARIES})
if (TCMALLOC_LIBRARY)
    target_link_libraries(qxmatch ${TCMALLOC_LIBRARIES})
endif()
