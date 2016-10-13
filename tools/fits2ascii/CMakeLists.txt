cmake_minimum_required(VERSION 2.6)
project(fits2ascii)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(fits2ascii fits2ascii.cpp)
target_link_libraries(fits2ascii ${PHYPP_LIBRARIES})
install(TARGETS fits2ascii DESTINATION bin)

