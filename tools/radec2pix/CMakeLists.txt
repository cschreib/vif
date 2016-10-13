cmake_minimum_required(VERSION 2.6)
project(radec2pix)

include_directories(${PHYPP_INCLUDE_DIRS})
add_executable(radec2pix radec2pix.cpp)
target_link_libraries(radec2pix ${PHYPP_LIBRARIES})
install(TARGETS radec2pix DESTINATION bin)
