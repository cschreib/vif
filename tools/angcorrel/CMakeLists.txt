cmake_minimum_required(VERSION 2.6)
project(angcorrel)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(angcorrel angcorrel.cpp)
target_link_libraries(angcorrel ${VIF_LIBRARIES})
install(TARGETS angcorrel DESTINATION bin)
