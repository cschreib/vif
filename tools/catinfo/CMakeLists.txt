cmake_minimum_required(VERSION 2.6)
project(catinfo)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(catinfo catinfo.cpp)
target_link_libraries(catinfo ${PHYPP_LIBRARIES})
install(TARGETS catinfo DESTINATION bin)
