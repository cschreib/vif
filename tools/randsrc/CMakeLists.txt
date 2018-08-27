cmake_minimum_required(VERSION 2.6)
project(randsrc)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(randsrc randsrc.cpp)
target_link_libraries(randsrc ${VIF_LIBRARIES})
install(TARGETS randsrc DESTINATION bin)
