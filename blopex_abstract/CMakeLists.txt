# TOP LEVEL CMakeLists.txt
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "DO NOT BUILD in-tree.")
endif()

# CMake 2.8.8 introduced an "OBJECT" library type that can be created
# with add_library. An object library compiles source files but does not
# archive or link their object files into a library.
# So we require minimum version 2.8.8.
cmake_minimum_required (VERSION 2.8.8 FATAL_ERROR)

project(BLOPEX_ABSTRACT C) #'C' indicates that this is a C project

# Prepare libraries.
if(SHARED)
  set(BLOPEX_ABSTRACT_LIBRARY_TYPE SHARED)
else()
  set(BLOPEX_ABSTRACT_LIBRARY_TYPE STATIC)
endif()

# The POSITION_INDEPENDENT_CODE property determines whether position
# independent executables or shared libraries will be created.
# This property is True by default for SHARED and MODULE library targets and
# False otherwise.
# So in the C compiler flags we do not explicitly initialize it.
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang") # using Clang
  set(compileFlags "-Wall -pedantic -g -march=native")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") # using GCC
  set(compileFlags "-Wall -pedantic -g -march=native")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC") # using Visual Studio C/C++
  set(compileFlags "/Za /W4 /EHsc /arch:AVX2")
# elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") # using Intel C++
#   set(compileFlags " ")
endif()

if (NOT CONFIGURED_ONCE)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${compileFlags}"
       CACHE STRING "Flags used by the compiler during all build types."
       FORCE
  )
  set(CMAKE_BUILD_TYPE Release CACHE STRING "default Release." FORCE)
endif()

# All DIRs added via add_subdirectory() must have a CMakeLists.txt
add_subdirectory(multivector) # make obj libs: multivector, temp_multivector
add_subdirectory(krylov)      # make obj lib: lobpcg
add_subdirectory(utilities)   # make obj lib: fortran_matrix

# Make the BLOPEX library from the corresponding object libraries.
add_library(BLOPEX
  $<TARGET_OBJECTS:lobpcg>
  $<TARGET_OBJECTS:multivector>
  $<TARGET_OBJECTS:temp_multivector>
  $<TARGET_OBJECTS:fortran_matrix>
)

# Install rule.
# ARCHIVE: Static libraries, DLL import libraries (DLL platforms)
# LIBRARY: Module libraries, shared libraries (non-DLL platforms)
# RUNTIME: Executables, DLL part of a shared library (for DLL platforms)
install(TARGETS BLOPEX
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)
# Install rules for header files is written in sub-directory CMakeLists.txt

set (CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing that CMake has configured at least once."
)
