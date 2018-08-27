cmake_minimum_required(VERSION 2.6)
project(subsrc)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(subsrc subsrc.cpp)
target_link_libraries(subsrc ${VIF_LIBRARIES})
install(TARGETS subsrc DESTINATION bin)
