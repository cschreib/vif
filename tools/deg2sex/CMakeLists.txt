cmake_minimum_required(VERSION 2.6)
project(deg2sex)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(deg2sex deg2sex.cpp)
target_link_libraries(deg2sex ${VIF_LIBRARIES})
install(TARGETS deg2sex DESTINATION bin)
