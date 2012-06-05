/*
 * FastaReader-test.cpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
 *
 * This file is part of the Mass Spectrometry Toolkit (MSTK).
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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

