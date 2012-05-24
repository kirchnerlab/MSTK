#############################################################################
# memory testing
#############################################################################
# check if we have valgrind
SET(HAVE_VALGRIND TRUE)
FIND_PROGRAM(VALGRIND_EXECUTABLE valgrind)
IF(VALGRIND_EXECUTABLE STREQUAL "VALGRIND_EXECUTABLE-NOTFOUND")
    MESSAGE(STATUS "WARNING: Could not find valgrind. Will NOT build memory tests.")
    SET(HAVE_VALGRIND FALSE)
ENDIF(VALGRIND_EXECUTABLE STREQUAL "VALGRIND_EXECUTABLE-NOTFOUND")

#############################################################################
# code coverage analysis
#############################################################################
IF (ENABLE_COVERAGE)
    # adjust the compile flags
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    # make the script available
    CONFIGURE_FILE(${MSTK_SOURCE_DIR}/cmake/scripts/coverage.sh.in
        ${MSTK_BINARY_DIR}/cmake/scripts/coverage.sh
        @ONLY IMMEDIATE
    )
    MESSAGE(STATUS "MSTK components available for coverage analysis: ${MSTK_COMPONENTS}")
    FOREACH(mstkComponent ${MSTK_COMPONENTS})
        # add a custom coverage target for every build target 
        SET(mstkComponentName "${mstkComponent}_coverage")
        ADD_CUSTOM_TARGET(${mstkComponentName} 
            COMMAND /bin/bash ${MSTK_BINARY_DIR}/cmake/scripts/coverage.sh ${mstkComponent}
        )
        MESSAGE(STATUS "Adding coverage target: ${mstkComponentName}")
        # keep a list of coverage targets
        SET(COVERAGE_TARGETS ${COVERAGE_TARGETS} ${mstkComponentName})
    ENDFOREACH(mstkComponent)
    # add a global coverage target that will determine coverage for every
    # active build target
    MESSAGE(STATUS "Global coverage now includes: ${COVERAGE_TARGETS}")
    ADD_CUSTOM_TARGET(coverage)
    ADD_DEPENDENCIES(coverage ${COVERAGE_TARGETS})
ENDIF(ENABLE_COVERAGE)

#############################################################################
# test-specific includes
#############################################################################
INCLUDE_DIRECTORIES(
    ${MSTK_SOURCE_DIR}/lib/vigra/include
    ${MSTK_SOURCE_DIR}/lib/vigra/test/include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
)

LINK_DIRECTORIES(
    ${MSTK_BINARY_DIR}/src
    ${MSTK_BINARY_DIR}/src/common
)

# subdirs
FOREACH(mstkComponent ${MSTK_COMPONENTS})
    ADD_SUBDIRECTORY(${mstkComponent})
ENDFOREACH(mstkComponent)

#### Sources

SET(TEST_LIBS mstk 
    mstk-common
)
IF (ENABLE_COVERAGE)
    SET(TEST_LIBS ${TEST_LIBS} gcov)
ENDIF(ENABLE_COVERAGE)
#### Unit tests
