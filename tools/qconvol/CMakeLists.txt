cmake_minimum_required(VERSION 2.6)
project(qconvol)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(qconvol qconvol.cpp)
target_link_libraries(qconvol ${VIF_LIBRARIES})
install(TARGETS qconvol DESTINATION bin)
