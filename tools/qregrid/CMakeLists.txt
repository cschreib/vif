cmake_minimum_required(VERSION 2.6)
project(qregrid)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(qregrid qregrid.cpp)
target_link_libraries(qregrid ${PHYPP_LIBRARIES})
install(TARGETS qregrid DESTINATION bin)
