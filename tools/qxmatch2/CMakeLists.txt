cmake_minimum_required(VERSION 2.6)
project(qxmatch2)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(qxmatch2 qxmatch2.cpp)
target_link_libraries(qxmatch2 ${VIF_LIBRARIES})
install(TARGETS qxmatch2 DESTINATION bin)
