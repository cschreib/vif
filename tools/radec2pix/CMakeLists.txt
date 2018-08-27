cmake_minimum_required(VERSION 2.6)
project(radec2pix)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(radec2pix radec2pix.cpp)
target_link_libraries(radec2pix ${VIF_LIBRARIES})
install(TARGETS radec2pix DESTINATION bin)
