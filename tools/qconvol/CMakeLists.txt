cmake_minimum_required(VERSION 2.6)
project(qconvol)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(qconvol qconvol.cpp)
target_link_libraries(qconvol ${PHYPP_LIBRARIES})
install(TARGETS qconvol DESTINATION bin)
