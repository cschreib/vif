cmake_minimum_required(VERSION 2.6)
project(qstack2)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(qstack2 qstack2.cpp)
target_link_libraries(qstack2 ${VIF_LIBRARIES})
install(TARGETS qstack2 DESTINATION bin)
