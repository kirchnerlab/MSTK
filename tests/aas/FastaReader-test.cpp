/*
 * FastaReader-test.cpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
 *
 */

#include "unittest.hxx"

#include <iostream>

/** Short description.
 * Long description.
 */
struct FastaReaderTestSuite : vigra::test_suite
{
    /** Constructor.
     * The FastaReaderTestSuite constructor adds all FastaReader tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    FastaReaderTestSuite() :
            vigra::test_suite("FastaReader")
    {
        add(testCase(&FastaReaderTestSuite::testFastaReader));
    }

    void testFastaReader()
    {
        failTest("not implemented yet");
    }
};

/** The main function that runs the tests for class FastaReader.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    FastaReaderTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

