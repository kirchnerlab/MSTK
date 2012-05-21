/*
 * QuickCharge-test.cpp
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
#include <MSTK/fe/QuickCharge.hpp>
#include <MSTK/fe/types/Xic.hpp>
#include <MSTK/fe/XicTraits.hpp>
#include <MSTK/common/Types.hpp>
#include <iostream>
#include <set>

using namespace mstk::fe;
using namespace mstk;

struct QuickChargeTestSuite : vigra::test_suite {
    QuickChargeTestSuite() : vigra::test_suite("QuickCharge") {
        add(testCase(&QuickChargeTestSuite::test));
    }

    void test() {
        double mzs[] =
                { 100.001, 100.2502, 100.33, 100.501, 101.001 };
        int expectedCharges[] = { 4, 3, 2, 13, 6 };
        std::vector<double> masses(&mzs[0], &mzs[4]);
        IsotopePattern ip = makeIsotopePattern(masses, 3324.0, 1);
        QuickCharge qc;
        std::vector<int> charges;
        qc(ip.begin(), ip.end(), std::back_inserter(charges));
        shouldEqual(charges.size(), static_cast<Size>(5));
        for (Size k = 0; k < charges.size(); ++k) {
            shouldEqual(charges[k], expectedCharges[k]);
        }
    }
};

int main()
{
    QuickChargeTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}


