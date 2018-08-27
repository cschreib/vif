cmake_minimum_required(VERSION 2.6)
project(catinfo)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(catinfo catinfo.cpp)
target_link_libraries(catinfo ${VIF_LIBRARIES})
install(TARGETS catinfo DESTINATION bin)
