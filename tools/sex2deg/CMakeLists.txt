cmake_minimum_required(VERSION 2.6)
project(sex2deg)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(sex2deg sex2deg.cpp)
target_link_libraries(sex2deg ${VIF_LIBRARIES})
install(TARGETS sex2deg DESTINATION bin)
