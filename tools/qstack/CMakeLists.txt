cmake_minimum_required(VERSION 2.6)
project(qstack)

include_directories(${CFITSIO_INCLUDES})
include_directories(${CCFITS_INCLUDES})
include_directories(${WCSLIB_INCLUDES})

add_executable(qstack qstack.cpp)

target_link_libraries(qstack ${CCFITS_LIBRARIES})
target_link_libraries(qstack ${CFITSIO_LIBRARIES})
target_link_libraries(qstack ${WCSLIB_LIBRARIES})
