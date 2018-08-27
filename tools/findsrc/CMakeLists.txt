cmake_minimum_required(VERSION 2.6)
project(findsrc)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(findsrc findsrc.cpp)
target_link_libraries(findsrc ${VIF_LIBRARIES})
install(TARGETS findsrc DESTINATION bin)
