cmake_minimum_required(VERSION 2.6)
project(qstack2)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(qstack2 qstack2.cpp)
target_link_libraries(qstack2 ${PHYPP_LIBRARIES})
install(TARGETS qstack2 DESTINATION bin)
