
SET(TEST_LIBS mstk-common)
#########  List of tests
ADD_MSTK_TEST("common" "Error" Error-test.cpp)
ADD_MSTK_TEST("common" "Log" Log-test.cpp)

MESSAGE(STATUS "Tests for 'common': ${MSTK_common_TEST_NAMES}")
MESSAGE(STATUS "Memory tests for 'common': ${MSTK_common_MEMTEST_NAMES}")

ADD_CUSTOM_TARGET(common_test 
    DEPENDS ${MSTK_common_TEST_NAMES} ${MSTK_common_MEMTEST_NAMES}
)

