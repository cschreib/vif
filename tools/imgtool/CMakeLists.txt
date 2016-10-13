cmake_minimum_required(VERSION 2.6)
project(imgtool)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(imgtool imgtool.cpp)
target_link_libraries(imgtool ${PHYPP_LIBRARIES})
install(TARGETS imgtool DESTINATION bin)
