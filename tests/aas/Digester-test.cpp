/*
 * Digester-test.cpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
 *
 */

#include <MSTK/aas/tools/Digester.hpp>

#include "unittest.hxx"

#include <iostream>

using namespace mstk;
using namespace aas;
using namespace aas::tools;

/** Short description.
 * Long description.
 */
struct DigesterTestSuite : vigra::test_suite
{
    /** Constructor.
     * The DigesterTestSuite constructor adds all Digester tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    DigesterTestSuite() :
            vigra::test_suite("Digester")
    {
        add(testCase(&DigesterTestSuite::testDigester));
    }

    void testDigester()
    {
        failTest("not implemented yet");
    }

};

/** The main function that runs the tests for class Digester.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    DigesterTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

