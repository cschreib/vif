cmake_minimum_required(VERSION 2.6)
project(psffit)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(psffit psffit.cpp)
target_link_libraries(psffit ${VIF_LIBRARIES})
install(TARGETS psffit DESTINATION bin)
