cmake_minimum_required(VERSION 2.6)
project(qmask)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(qmask qmask.cpp)
target_link_libraries(qmask ${PHYPP_LIBRARIES})
install(TARGETS qmask DESTINATION bin)
