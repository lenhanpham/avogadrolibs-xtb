# Tests for the ported xtb module. Unlike the calc suite these need no
# external data files: reference values come from the Fortran sources
# (see the parity harness documented in coordinationtest.cpp).
set(tests
  Coordination
  Eeq
  Elements
  GffAngle
  GffBond
  GffEgBond
  GffEgNonbond
  GffData
  GffDisp
  GffHb
  GffHbSet
  GffSetup
  GffDriver
  GffEeq
  GffNb
  GffGraph
  GffHuckel
  GffHyb
  GffOop
  GffTorsion
  GffVbond
  GfnTopo
)

# Build up the source file names.
set(testSrcs "")
foreach(TestName ${tests})
  message(STATUS "Adding ${TestName} test.")
  string(TOLOWER ${TestName} testname)
  list(APPEND testSrcs ${testname}test.cpp)
endforeach()

add_executable(AvogadroXtbTests ${testSrcs})
target_link_libraries(AvogadroXtbTests Avogadro::Xtb
  ${GTEST_BOTH_LIBRARIES} ${EXTRA_LINK_LIB})

foreach(TestName ${tests})
  add_test(NAME "Xtb-${TestName}"
    COMMAND AvogadroXtbTests "--gtest_filter=${TestName}Test.*")
endforeach()
