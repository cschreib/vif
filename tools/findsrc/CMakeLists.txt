cmake_minimum_required(VERSION 2.6)
project(findsrc)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(findsrc findsrc.cpp)
target_link_libraries(findsrc ${PHYPP_LIBRARIES})
install(TARGETS findsrc DESTINATION bin)
