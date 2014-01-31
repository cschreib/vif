# Find the native LLVM includes and library
#
# LLVM_INCLUDE_DIR - where to find llvm include files
# LLVM_LIBRARY_DIR - where to find llvm libs
# LLVM_MODULE_LIBS - list of llvm libs for working with modules.
# LLVM_FOUND - True if llvm found.
#
# Note: code from https://github.com/cloudera/Impala
# Note: modified by C. Schreiber

find_program(LLVM_CONFIG_EXECUTABLE
  NAMES llvm-config llvm-config-3.0 llvm-config-3.1 llvm-config-3.2
        llvm-config-3.3 llvm-config-3.4 llvm-config-3.5
  PATHS
  $ENV{LLVM_HOME}
  /usr/bin
  NO_DEFAULT_PATH
)
find_program(LLVM_CONFIG_EXECUTABLE llvm-config)

if (NOT LLVM_CONFIG_EXECUTABLE)
  message(FATAL_ERROR "Could not find llvm-config")
endif (NOT LLVM_CONFIG_EXECUTABLE)

message(STATUS "LLVM llvm-config found at: ${LLVM_CONFIG_EXECUTABLE}")

execute_process(
  COMMAND ${LLVM_CONFIG_EXECUTABLE} --version
  OUTPUT_VARIABLE LLVM_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND ${LLVM_CONFIG_EXECUTABLE} --includedir
  OUTPUT_VARIABLE LLVM_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND ${LLVM_CONFIG_EXECUTABLE} --libdir
  OUTPUT_VARIABLE LLVM_LIBRARY_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the link libs we need. llvm has many and we don't want to link all of the libs
# if we don't need them.
execute_process(
  COMMAND ${LLVM_CONFIG_EXECUTABLE} --libnames core
  OUTPUT_VARIABLE LLVM_MODULE_LIBS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# CMake really doesn't like adding link directories and wants absolute paths
# Reconstruct it with LLVM_MODULE_LIBS and LLVM_LIBRARY_DIR
string(REPLACE " " ";" LIBS_LIST ${LLVM_MODULE_LIBS})
set (LLVM_MODULE_LIBS "-ldl")
foreach (LIB ${LIBS_LIST})
  set(LLVM_MODULE_LIBS ${LLVM_MODULE_LIBS} "${LLVM_LIBRARY_DIR}/${LIB}")
endforeach(LIB)

message(STATUS "LLVM include dir: ${LLVM_INCLUDE_DIR}")
message(STATUS "LLVM lib dir: ${LLVM_LIBRARY_DIR}")
message(STATUS "LLVM libs: ${LLVM_MODULE_LIBS}")
