# Configure libs for tests
SET(TEST_LIBS mstk-ipaca mstk-common)
#########  List of tests
ADD_MSTK_TEST("ipaca" "Mercury7" Mercury7-test.cpp)
ADD_MSTK_TEST("ipaca" "Mercury7Impl" Mercury7Impl-test.cpp)
ADD_MSTK_TEST("ipaca" "Stoichiometry" Stoichiometry-test.cpp)

MESSAGE(STATUS "Tests for 'ipaca': ${MSTK_ipaca_TEST_NAMES}")
MESSAGE(STATUS "Memory tests for 'ipaca': ${MSTK_ipaca_MEMTEST_NAMES}")

ADD_CUSTOM_TARGET(ipaca_test 
    DEPENDS ${MSTK_ipaca_TEST_NAMES} ${MSTK_ipaca_MEMTEST_NAMES}
)

