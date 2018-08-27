cmake_minimum_required(VERSION 2.6)
project(photinfo)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(photinfo photinfo.cpp)
target_link_libraries(photinfo ${VIF_LIBRARIES})
install(TARGETS photinfo DESTINATION bin)
