include_directories(.)

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

  ${CMAKE_HOME_DIRECTORY}/src/tests/init_test
  ${CMAKE_HOME_DIRECTORY}/src/tests/evaluation_domain_test
  ${CMAKE_HOME_DIRECTORY}/src/tests/polynomial_arithmetic_test
  ${CMAKE_HOME_DIRECTORY}/src/tests/kronecker_substitution_test
)
target_link_libraries(
  gtests

  ${GTEST_BOTH_LIBRARIES}
  ${LIBFF_LIBRARIES}
  ${GMP_LIBRARIES}
  ${GMPXX_LIBRARIES}
  ${PROCPS_LIBRARIES}
)

# Profiling with Double by default
add_executable(
  profiling_menu

  ${CMAKE_HOME_DIRECTORY}/src/profiling/profiling_menu.cpp
)
set_target_properties(
  profiling_menu

  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY
  ${CMAKE_HOME_DIRECTORY}
)

add_executable(
  profiler

  ${CMAKE_HOME_DIRECTORY}/src/profiling/profile/profile.cpp
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

  ${LIBFF_LIBRARIES}
  ${GMP_LIBRARIES}
  ${GMPXX_LIBRARIES}
  ${PROCPS_LIBRARIES}
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