cmake_minimum_required(VERSION 2.6)
project(getgal)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(getgal getgal.cpp)
target_link_libraries(getgal ${VIF_LIBRARIES})
install(TARGETS getgal DESTINATION bin)
