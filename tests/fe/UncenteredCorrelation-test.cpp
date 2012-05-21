/*
 * UncenteredCorrelation-test.hpp
 *
 * Copyright (C) 2012 Marc Kirchner
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
#include "utilities.hpp"
#include <MSTK/fe/UncenteredCorrelation.hpp>
#include <MSTK/fe/types/Xic.hpp>
#include <iostream>

using namespace mstk::fe;

namespace { // local namespace

// Concrete class for testing.
class Correlator : public UncenteredCorrelation
{
public:
    template<class IT>
    UncenteredCorrelation::ThresholdType run(IT lhsFirst, IT lhsLast,
        IT rhsFirst, IT rhsLast)
    {
        return correlate(lhsFirst, lhsLast, rhsFirst, rhsLast);
    }
};

}

struct UncenteredCorrelationTestSuite : vigra::test_suite
{
    UncenteredCorrelationTestSuite() :
            vigra::test_suite("UncenteredCorrelation")
    {
        add(testCase(&UncenteredCorrelationTestSuite::testCorrelate));
    }

    void testCorrelate()
    {
        double mz[] = { 99.9, 100.0, 100.1 };
        double rt[] = { 350.0, 352.0, 354.0 };
        unsigned int sn[] = { 42, 43, 44 };
        double ab[] = { 1.0, 1.0, 1.0 };
        Xic x1 = makeXic(3, mz, rt, sn, ab);
        unsigned int sn2[] = { 42, 43, 44 };
        double ab2[] = { 1.0, 1.0, 1.0 };
        Xic x2 = makeXic(3, mz, rt, sn2, ab2);
        Correlator cor;
        shouldEqualTolerance(
            cor.run(x1.begin(), x1.end(), x2.begin(), x2.end()), 1.0, 1e-12);
        shouldEqualTolerance(
            cor.run(x2.begin(), x2.end(), x1.begin(), x1.end()), 1.0, 1e-12);
        unsigned int sn3[] = { 45, 46, 47 };
        double rt3[] = { 356.0, 358.0, 360.0 };
        Xic x3 = makeXic(3, mz, rt3, sn3, ab2);
        shouldEqualTolerance(
            cor.run(x1.begin(), x1.end(), x3.begin(), x3.end()), 0.0, 1e-12);
        shouldEqualTolerance(
            cor.run(x3.begin(), x3.end(), x1.begin(), x1.end()), 0.0, 1e-12);
    }
};

int main()
{
    UncenteredCorrelationTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

