cmake_minimum_required(VERSION 2.6)
project(refgen)

include_directories(${CLANG_INCLUDE})

add_definitions(-D__STDC_CONSTANT_MACROS)
add_definitions(-D__STDC_FORMAT_MACROS)
add_definitions(-D__STDC_LIMIT_MACROS)

add_executable(vif-refgen refgen.cpp)

target_link_libraries(vif-refgen ${CLANG_LIBS})
if (TCMALLOC_LIBRARY)
    target_link_libraries(vif-refgen ${TCMALLOC_LIBRARIES})
endif()

install(TARGETS vif-refgen RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
