# test data
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/testdata.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/testdata.hpp
    @ONLY IMMEDIATE
)
INCLUDE_DIRECTORIES(
    ${VIGRA_INCLUDE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
)

# Configure libs for tests
SET(TEST_LIBS mstk-psf mstk-common)
#########  List of tests
ADD_MSTK_TEST("psf" "PeakParameter" PeakParameter-test.cpp)
ADD_MSTK_TEST("psf" "PeakShapeFunction" PeakShapeFunction-test.cpp)
ADD_MSTK_TEST("psf" "PeakShape" PeakShape-test.cpp )
ADD_MSTK_TEST("psf" "SpectrumAlgorithm" SpectrumAlgorithm-test.cpp)

MESSAGE(STATUS "Tests for 'psf': ${MSTK_psf_TEST_NAMES}")
MESSAGE(STATUS "Memory tests for 'psf': ${MSTK_psf_MEMTEST_NAMES}")

ADD_CUSTOM_TARGET(psf_test 
    DEPENDS ${MSTK_psf_TEST_NAMES} ${MSTK_psf_MEMTEST_NAMES}
)

