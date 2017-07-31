add_subdirectory(gtest EXCLUDE_FROM_ALL)

# TODO(scripr-lab): upstream these buildsystems to the respective
# projects.


if(${CURVE} STREQUAL "BN128")
  include_directories(ate-pairing/include)
  include_directories(xbyak)
  add_library(
    zm
    STATIC

    ate-pairing/src/zm.cpp
    ate-pairing/src/zm2.cpp
  )
endif()

include(ExternalProject)
ExternalProject_Add(
  ff

  SOURCE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/libff

  CMAKE_ARGS
  -DDEBUG=${DEBUG}
  -DLOWMEM=${LOWMEM}
  -DMULTICORE=${MULTICORE}
  -DBINARY_OUTPUT=${BINARY_OUTPUT}
  -DMONTGOMERY_OUTPUT=${MONTGOMERY_OUTPUT}
  -DUSE_PT_COMPRESSION=${USE_PT_COMPRESSION}
  -DPROFILE_OP_COUNTS=${PROFILE_OP_COUNTS}
  -DUSE_MIXED_ADDITION=${USE_MIXED_ADDITION}
  -DWITH_PROCPS=${WITH_PROCPS}
  -DCPPDEBUG=${CPPDEBUG}
  -DPERFORMANCE=${PERFORMANCE}
  -DUSE_ASM=${USE_ASM}
  -DOPT_FLAGS=${OPT_FLAGS}
  -DCMAKE_INSTALL_PREFIX=${DEPENDS_DIR}/libff
  -DIS_PARENT_CMAKE=OFF

  INSTALL_DIR
  ${DEPENDS_DIR}/libff
)