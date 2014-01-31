cmake_minimum_required(VERSION 2.6)
project(phy++)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 3.2)
        message(STATUS "clang version >= 3.3 (${CMAKE_CXX_COMPILER_VERSION})")
    else()
        message(ERROR ": phy++ requires advanced features from the C++14 norm that are only available with clang 3.3 or higher (your version: ${CMAKE_CXX_COMPILER_VERSION}). Please upgrade your compiler.")
    endif()

    add_definitions(-Weverything)
    add_definitions(-Wno-c++98-compat-pedantic)
    add_definitions(-Wno-c++98-compat)
    add_definitions(-Wno-unused-parameter)
    add_definitions(-Wno-sign-conversion)
    add_definitions(-Wno-conversion)
    add_definitions(-Wno-missing-variable-declarations)
    add_definitions(-Wno-missing-prototypes)
    add_definitions(-Wno-padded)
    add_definitions(-Wno-float-equal)
    add_definitions(-Wno-unused-variable)
    add_definitions(-Wno-global-constructors)
    add_definitions(-Wno-exit-time-destructors)
    add_definitions(-Wno-weak-vtables)
    add_definitions(-Wno-covered-switch-default)
    add_definitions(-Wno-documentation-unknown-command)
    add_definitions(-Wno-unneeded-internal-declaration)
    add_definitions(-Wno-unused-function)
    add_definitions(-Wno-unused-macros)
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_definitions(-g)
    else()
        add_definitions(-O3)
    endif()
    add_definitions(-std=c++1y)
    add_definitions(-ftemplate-backtrace-limit=0)
    add_definitions(-ferror-limit=5)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.7)
        message(STATUS "gcc version >= 4.8 (${CMAKE_CXX_COMPILER_VERSION})")
    else()
        message(ERROR ": phy++ requires advanced features from the C++14 norm that are only available with gcc 4.8 or higher (your version: ${CMAKE_CXX_COMPILER_VERSION}). Please upgrade your compiler.")
    endif()

    add_definitions(-Wall)
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_definitions(-g)
    else()
        add_definitions(-O3)
    endif()
    add_definitions(-std=c++1y)
    add_definitions(-fmax-errors=5)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  message(ERROR "Intel C++ compiler is not supported")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  message(ERROR "Microsoft Visual C++ compiler is not supported")
endif()

set(PHYPP_INCLUDES "${PROJECT_SOURCE_DIR}/include")

file(MAKE_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

find_package(LibClang 3)
find_package(LLVM 3)
find_package(CFITSIO)
find_package(CCFITS)
find_package(WCSLib)
find_package(LAPACK)
find_package(GSL)

# build the refgen tool
add_subdirectory(tools/refgen)

# build pre-phy++ tools (on demand)
add_subdirectory(tools/qstack)
add_subdirectory(tools/qxmatch)

# build configuration
set(PHYPP_OUTPUT_CONFIG_FILE ".phypprc")
set(PHYPP_OUTPUT_CONFIG "${CMAKE_BINARY_DIR}/${PHYPP_OUTPUT_CONFIG_FILE}")
file(WRITE ${PHYPP_OUTPUT_CONFIG} "# phy++ configuration file\n\n")

file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_PATH=${PROJECT_SOURCE_DIR}\n")
file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_DATA_DIR=\${PHYPP_PATH}/data\n\n")

file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_REFGEN=${CMAKE_INSTALL_PREFIX}/bin/phy++-refgen\n")
file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_TEMP_DIR=/tmp/\n")
file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_REFGEN_FLAGS=\"-DREFLECTION_STAGE -ferror-limit=5 -std=c++1y -Weverything -Wno-c++98-compat-pedantic -Wno-c++98-compat -Wno-unused-parameter -Wno-sign-conversion -Wno-conversion -Wno-missing-prototypes -Wno-padded -Wno-float-equal -Wno-unused-variable -Wno-global-constructors -Wno-exit-time-destructors\"\n\n")

file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER=${CMAKE_CXX_COMPILER}\n")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_FLAGS=\"-ferror-limit=5 -std=c++1y -lstdc++ -lcfitsio -lCCfits -lwcs -llapack -lgsl -Weverything -Wno-c++98-compat-pedantic -Wno-c++98-compat -Wno-unused-parameter -Wno-sign-conversion -Wno-conversion -Wno-missing-prototypes -Wno-padded -Wno-float-equal -Wno-unused-variable -Wno-global-constructors -Wno-exit-time-destructors\"\n")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_DEBUG_FLAGS=-g\n")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_OPT_FLAGS=\"-O3 -ltcmalloc\"\n")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_PROF_FLAGS=-lprofiler\n")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_FLAGS=\"-fmax-errors=5 -std=c++1y -lstdc++ -lcfitsio -lCCfits -lwcs -llapack -lgsl -Wall\"\n")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_DEBUG_FLAGS=-g\n")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_OPT_FLAGS=\"-O3 -ltcmalloc\"\n")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_PROF_FLAGS=-lprofiler\n")
endif()

# install headers
file(GLOB files ${PROJECT_SOURCE_DIR}/include/phypp/*.hpp)
install(FILES ${PROJECT_SOURCE_DIR}/include/phypp.hpp DESTINATION include)
install(FILES ${files} DESTINATION include/phypp)

# install build scripts
install(PROGRAMS bin/phy++ bin/ophy++ bin/gphy++ bin/gophy++ bin/prophy++
    DESTINATION bin COMPONENT scripts)

# install config
install(FILES ${PHYPP_OUTPUT_CONFIG} DESTINATION $ENV{HOME} COMPONENT config)
install(CODE "message(\"NOTE: to finalize the installation, please source $ENV{HOME}/${PHYPP_OUTPUT_CONFIG_FILE} in your startup script (~/.bashrc).\")")
