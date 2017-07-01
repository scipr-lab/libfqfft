install(
  DIRECTORY "" DESTINATION "include/libfqfft"
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.tcc"
  PATTERN "logs" EXCLUDE
  PATTERN "tests" EXCLUDE
  PATTERN "tutorials" EXCLUDE
)

# GTests
add_executable(
  gtests
  EXCLUDE_FROM_ALL

  tests/init_test
  tests/evaluation_domain_test
  tests/polynomial_arithmetic_test
  tests/kronecker_substitution_test
)
target_link_libraries(
  gtests

  ${GTEST_BOTH_LIBRARIES}
  ${GMP_LIBRARIES}
  ${GMPXX_LIBRARIES}
  ${PROCPS_LIBRARIES}
  ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/libff/lib/libff.a
)
target_include_directories(
  gtests

  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/libff/include
)

# Profiling with Double by default
add_executable(
  profiling_menu
  EXCLUDE_FROM_ALL

  profiling/profiling_menu.cpp
)
set_target_properties(
  profiling_menu

  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY
  ${CMAKE_HOME_DIRECTORY}
)

add_executable(
  profiler
  EXCLUDE_FROM_ALL

  profiling/profile/profile.cpp
)
target_compile_definitions(
  profiler

  PUBLIC
  -DPROFILE_OP_COUNTS
  -DMULTICORE
)
set_target_properties(
  profiler

  PROPERTIES
  COMPILE_FLAGS "-fopenmp"
  LINK_FLAGS "-fopenmp"
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}
)
target_link_libraries(
  profiler

  ${GMP_LIBRARIES}
  ${GMPXX_LIBRARIES}
  ${PROCPS_LIBRARIES}
  ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/libff/lib/libff.a
)
target_include_directories(
  profiler

  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../third_party/libff/include
)

if(${PROF_DOUBLE})
  target_compile_definitions(
    profiling_menu
    PUBLIC -DPROF_DOUBLE
  )
  target_compile_definitions(
    profiler
    PUBLIC -DPROF_DOUBLE
  )
endif()

include(CTest)
add_test(
  NAME gtests
  COMMAND gtests
)

add_dependencies(check gtests)
