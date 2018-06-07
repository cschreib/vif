cmake_minimum_required(VERSION 2.6)
project(qaper)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(qaper qaper.cpp)
target_link_libraries(qaper ${PHYPP_LIBRARIES})
install(TARGETS qaper DESTINATION bin)
