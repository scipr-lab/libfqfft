cmake_minimum_required(VERSION 2.8)

project(libfqfft)

option(
  MULTICORE
  "Enable parallelized execution, using OpenMP"
  OFF
)

option(
  WITH_PROCPS
  "Use procps for memory profiling"
  ON
)

set(
  OPT_FLAGS
  ""
  CACHE
  STRING
  "Override C++ compiler optimization flags"
)

option(
  PROF_DOUBLE
  "Use Double for profiling"
  ON
)

if(CMAKE_COMPILER_IS_GNUCXX OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # Common compilation flags and warning configuration
  set(
    CMAKE_CXX_FLAGS

    "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Wfatal-errors -pthread"
  )

  if("${MULTICORE}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  endif()

   # Default optimizations flags (to override, use -DOPT_FLAGS=...)
  if("${OPT_FLAGS}" STREQUAL "")
    set(
      OPT_FLAGS
      "-ggdb3 -O2 -march=native -mtune=native"
    )
  endif()
endif()

if("${MULTICORE}")
  add_definitions(-DMULTICORE=1)
endif()

set(
  CMAKE_CXX_FLAGS
  "${CMAKE_CXX_FLAGS} ${OPT_FLAGS}"
)

include(FindPkgConfig)
if("${WITH_PROCPS}")
  pkg_check_modules(
    PROCPS
    REQUIRED

    libprocps
  )
else()
  add_definitions(
    -DNO_PROCPS
  )
endif()

include_directories(.)

# GMP
find_path(GMP_INCLUDE_DIR NAMES gmp.h)
find_library(GMP_LIBRARIES NAMES gmp libgmp)
find_library(GMPXX_LIBRARIES NAMES gmpxx libgmpxx)

# libff
find_path(LIBFF_INCLUDE_DIR NAMES libff)
include_directories(${LIBFF_INCLUDE_DIR})
find_library(LIBFF_LIBRARIES NAMES ff libff)

# Google Test
enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

add_custom_target(
  check

  COMMAND
  ${CMAKE_CTEST_COMMAND}
)
add_custom_target(
  tutorials

  COMMAND lagrange_polynomial_evaluation
  COMMAND polynomial_evaluation
  COMMAND polynomial_multiplication
)

add_subdirectory(libfqfft)
add_subdirectory(tutorials)
