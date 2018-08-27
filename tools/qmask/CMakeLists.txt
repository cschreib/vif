cmake_minimum_required(VERSION 2.6)
project(qmask)

include_directories(${VIF_INCLUDE_DIRS})
add_executable(qmask qmask.cpp)
target_link_libraries(qmask ${VIF_LIBRARIES})
install(TARGETS qmask DESTINATION bin)
