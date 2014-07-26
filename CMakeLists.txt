cmake_minimum_required(VERSION 2.6)
project(phy++)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.3))
        message(STATUS "clang version >= 3.3 (${CMAKE_CXX_COMPILER_VERSION})")
    else()
        message(FATAL_ERROR "phy++ requires advanced features from the C++14 norm that are only available with clang 3.3 or higher (your version: ${CMAKE_CXX_COMPILER_VERSION}). Please upgrade your compiler.")
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
    add_definitions(-Wno-weak-vtables)
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_definitions(-g)
    else()
        add_definitions(-O3)
    endif()
    add_definitions(-std=c++11)
    add_definitions(-ftemplate-backtrace-limit=0)
    add_definitions(-ferror-limit=5)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if(NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8))
        message(STATUS "gcc version >= 4.8 (${CMAKE_CXX_COMPILER_VERSION})")
    else()
        message(FATAL_ERROR "phy++ requires advanced features from the C++14 norm that are only available with gcc 4.8 or higher (your version: ${CMAKE_CXX_COMPILER_VERSION}). Please upgrade your compiler.")
    endif()

    add_definitions(-Wall)
    if(CMAKE_BUILD_TYPE MATCHES Debug)
        add_definitions(-g)
    else()
        add_definitions(-O3)
    endif()
    add_definitions(-std=c++11)
    add_definitions(-fmax-errors=5)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  message(ERROR "Intel C++ compiler is not supported")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  message(ERROR "Microsoft Visual C++ compiler is not supported")
endif()

set(PHYPP_INCLUDES "${PROJECT_SOURCE_DIR}/include")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

# find required libraries
find_package(LibClang 3)
find_package(CFITSIO)
find_package(CCFITS)
find_package(WCSLib)
find_package(LAPACK)

# find optional libraries
find_package(GSL)
find_package(FFTW 3)
find_package(GooglePerfTools)

set(PHYPP_ADD_COMPILER_FLAGS "")
set(PHYPP_ADD_OPT_COMPILER_FLAGS "")

# handle conditional GSL support
if (NOT GSL_FOUND)
    message("note: the GSL library could not be found: certain mathematical functions will not be available, but apart from that the library will function properly")
endif()
if (NO_GSL)
    message("note: the GSL library has been disabled: certain mathematical functions will not be available, but apart from that the library will function properly")
endif()
if (NOT GSL_FOUND OR NO_GSL)
    add_definitions(-DNO_GSL)
    set(PHYPP_ADD_COMPILER_FLAGS ${PHYPP_ADD_COMPILER_FLAGS} " -DNO_GSL")
endif()

# handle conditional FFTW support

if (NOT FFTW_FOUND)
    message("note: the FFTW library cound not be found: certain mathematical functions will be slow, but apart from that the library will function properly")
endif()
if (NO_FFTW)
    message("note: the FFTW library has been disabled: certain mathematical functions will be slow, but apart from that the library will function properly")
endif()
if (NOT FFTW_FOUND OR NO_FFTW)
    add_definitions(-DNO_FFTW)
    set(PHYPP_ADD_COMPILER_FLAGS ${PHYPP_ADD_COMPILER_FLAGS} " -DNO_FFTW")
else()
    set(PHYPP_ADD_COMPILER_FLAGS ${PHYPP_ADD_COMPILER_FLAGS} " -lfftw3")
endif()

# handle conditional Google perftools support
if (NOT TCMALLOC_LIBRARY)
    message("note: could not find Google's tcmalloc library: installing this library will make phy++ programs slightly faster")
else()
    set(PHYPP_ADD_OPT_COMPILER_FLAGS ${PHYPP_ADD_OPT_COMPILER_FLAGS} " -ltcmalloc")
endif()

if (NOT PROFILER_LIBRARY)
    message("note: could not find Google's profiler library: will not be able to profile phy++ programs, but apart from that the library will function properly")
endif()

# build the refgen tool
add_subdirectory(tools/refgen)

# build configuration
set(PHYPP_OUTPUT_CONFIG_FILE ".phypprc")
set(PHYPP_OUTPUT_CONFIG "${CMAKE_BINARY_DIR}/${PHYPP_OUTPUT_CONFIG_FILE}")
file(WRITE ${PHYPP_OUTPUT_CONFIG} "# phy++ configuration file\n\n")

file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_INCLUDE_DIR=-I\"${CMAKE_INSTALL_PREFIX}/include\"\n")

file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_PATH=${PROJECT_SOURCE_DIR}/\n")
file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_DATA_DIR=\${PHYPP_PATH}data/\n\n")

file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_REFGEN=${CMAKE_INSTALL_PREFIX}/bin/phy++-refgen\n")
file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_TEMP_DIR=/tmp/\n")
file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_REFGEN_FLAGS=\"-DREFLECTION_STAGE -ferror-limit=5 -std=c++11 -Weverything -Wno-c++98-compat-pedantic -Wno-c++98-compat -Wno-unused-parameter -Wno-sign-conversion -Wno-conversion -Wno-missing-prototypes -Wno-padded -Wno-float-equal -Wno-unused-variable -Wno-global-constructors -Wno-exit-time-destructors -Wno-missing-variable-declarations -Wno-weak-vtables\"\n\n")

file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER=${CMAKE_CXX_COMPILER}\n")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_FLAGS=\"-ferror-limit=5 -std=c++11 -lstdc++ -lcfitsio -lCCfits -lwcs -llapack -lgsl -Weverything -Wno-c++98-compat-pedantic -Wno-c++98-compat -Wno-unused-parameter -Wno-sign-conversion -Wno-conversion -Wno-missing-prototypes -Wno-padded -Wno-float-equal -Wno-unused-variable -Wno-global-constructors -Wno-exit-time-destructors -Wno-missing-variable-declarations -Wno-weak-vtables ${PHYPP_ADD_COMPILER_FLAGS}\"\n")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_DEBUG_FLAGS=-g\n")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_OPT_FLAGS=\"-O3 ${PHYPP_ADD_OPT_COMPILER_FLAGS}\"\n")
    if (PROFILER_LIBRARY)
        file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_PROF_FLAGS=-lprofiler\n")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_FLAGS=\"-fmax-errors=5 -std=c++11 -lstdc++ -lcfitsio -lCCfits -lwcs -llapack -lgsl -Wall ${PHYPP_ADD_COMPILER_FLAGS}\"\n")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_DEBUG_FLAGS=-g\n")
    file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_OPT_FLAGS=\"-O3 ${PHYPP_ADD_OPT_COMPILER_FLAGS}\"\n")
    if (PROFILER_LIBRARY)
        file(APPEND ${PHYPP_OUTPUT_CONFIG} "export PHYPP_COMPILER_PROF_FLAGS=-lprofiler\n")
    endif()
endif()

# install headers
file(GLOB PHYPP_INCLUDES ${PROJECT_SOURCE_DIR}/include/phypp/*.hpp)
install(FILES ${PHYPP_INCLUDES} DESTINATION include/phypp)
install(FILES ${PROJECT_SOURCE_DIR}/include/phypp.hpp DESTINATION include)

# install build scripts
install(PROGRAMS ${PROJECT_SOURCE_DIR}/bin/phy++ ${PROJECT_SOURCE_DIR}/bin/ophy++
    ${PROJECT_SOURCE_DIR}/bin/gphy++ ${PROJECT_SOURCE_DIR}/bin/gophy++
    ${PROJECT_SOURCE_DIR}/bin/cphy++
    DESTINATION bin COMPONENT scripts)
if (PROFILER_LIBRARIES)
    install(PROGRAMS ${PROJECT_SOURCE_DIR}/bin/prophy++ DESTINATION bin COMPONENT scripts)
endif()

# install config
install(FILES ${PHYPP_OUTPUT_CONFIG} DESTINATION $ENV{HOME} COMPONENT config)
install(CODE "message(\"note: to finalize the installation, please source $ENV{HOME}/${PHYPP_OUTPUT_CONFIG_FILE} in your startup script (~/.bashrc).\")")
