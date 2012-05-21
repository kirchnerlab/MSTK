/*
 * RunningMeanSmoother-test.hpp
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
#include <MSTK/fe/RunningMeanSmoother.hpp>
#include <MSTK/fe/types/Xic.hpp>
#include <iostream>

using namespace mstk::fe;

namespace { // local namespace

// Concrete class for testing.
class Smoother : public RunningMeanSmoother
{
public:
    template <class IT>
    void run(IT first, IT last) {
        smooth(first, last);
    }
};

}

struct RunningMeanSmootherTestSuite : vigra::test_suite {
    RunningMeanSmootherTestSuite() : vigra::test_suite("RunningMeanSmoother") {
        add(testCase(&RunningMeanSmootherTestSuite::testGetSmoothedXic));
    }

    void testGetSmoothedXic()
    {
        double mzs1[] =
                { 100.001, 100.003, 100.002, 100.005, 100.001, 100.003 };
        double rts1[] = { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0 };
        unsigned int sns1[] = { 1, 2, 3, 4, 5, 6 };
        double abs1[] = { 1.0, 2.0, 3.0, 2.0, 1.0, 0.5 };
        Xic xic = makeXic(6, mzs1, rts1, sns1, abs1);
        Xic sxic(xic);
        Smoother smoother;
        smoother.run(sxic.begin(), sxic.end());
        double expectedAbs[] = { 1.0, 2.0, 2.333333, 2, 1.166666, 0.5 };
        for (size_t i = 0; i < 6; ++i) {
            shouldEqualTolerance(sxic[i].getAbundance(), expectedAbs[i], 1e-6);
        }
    }
};

int main()
{
    RunningMeanSmootherTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}


