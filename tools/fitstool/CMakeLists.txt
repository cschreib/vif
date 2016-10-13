cmake_minimum_required(VERSION 2.6)
project(fitstool)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(fitstool fitstool.cpp)
target_link_libraries(fitstool ${PHYPP_LIBRARIES})
install(TARGETS fitstool DESTINATION bin)
