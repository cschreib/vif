cmake_minimum_required(VERSION 2.6)
project(refgen)

include_directories(${CLANG_INCLUDE})

add_executable(refgen refgen.cpp)

target_link_libraries(refgen ${CLANG_LIBS})
