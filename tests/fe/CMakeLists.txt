# External libs
FIND_PACKAGE(LIBFBI REQUIRED) 
INCLUDE(${LIBFBI_USE_FILE})
# Test utility library
SET(SRCS_MSTK_fe_TEST_UTILS utilities.cpp)
ADD_LIBRARY(mstk-fe-test-utils ${SRCS_MSTK_fe_TEST_UTILS})
ADD_DEPENDENCIES(mstk-fe-test-utils mstk-fe)
TARGET_LINK_LIBRARIES(mstk-fe-test-utils mstk-fe)

# Configure libs for tests
SET(TEST_LIBS mstk-fe mstk-common mstk-fe-test-utils)
#########  List of tests
ADD_MSTK_TEST("fe" "Centroid" Centroid-test.cpp)
ADD_MSTK_TEST("fe" "Centroider" Centroider-test.cpp)
ADD_MSTK_TEST("fe" "CentroidWeightedMeanDisambiguator" CentroidWeightedMeanDisambiguator-test.cpp )
ADD_MSTK_TEST("fe" "GaussianMeanAccumulator" GaussianMeanAccumulator-test.cpp)
ADD_MSTK_TEST("fe" "IsotopePattern" IsotopePattern-test.cpp)
ADD_MSTK_TEST("fe" "IsotopePatternExtractor" IsotopePatternExtractor-test.cpp)
ADD_MSTK_TEST("fe" "QuickCharge" QuickCharge-test.cpp)
ADD_MSTK_TEST("fe" "RunningMeanSmoother" RunningMeanSmoother-test.cpp)
ADD_MSTK_TEST("fe" "SimpleBumpFinder" SimpleBumpFinder-test.cpp)
ADD_MSTK_TEST("fe" "Spectrum" Spectrum-test.cpp)
ADD_MSTK_TEST("fe" "Splitter" Splitter-test.cpp)
ADD_MSTK_TEST("fe" "SumAbundanceAccumulator" SumAbundanceAccumulator-test.cpp)
ADD_MSTK_TEST("fe" "UncenteredCorrelation" UncenteredCorrelation-test.cpp)
ADD_MSTK_TEST("fe" "Xic" Xic-test.cpp)
ADD_MSTK_TEST("fe" "XicExtractor" XicExtractor-test.cpp)
ADD_MSTK_TEST("fe" "XicLocalMinSplitter" XicLocalMinSplitter-test.cpp)

MESSAGE(STATUS "Tests for 'fe': ${MSTK_fe_TEST_NAMES}")
MESSAGE(STATUS "Memory tests for 'fe': ${MSTK_fe_MEMTEST_NAMES}")

ADD_CUSTOM_TARGET(fe_test 
    DEPENDS ${MSTK_fe_TEST_NAMES} ${MSTK_fe_MEMTEST_NAMES}
)

