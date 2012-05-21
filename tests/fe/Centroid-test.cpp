/*
 * Centroid-test.cpp
 *
 * Copyright (C) 2011 Marc Kirchner
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
#include <iostream>
#include <MSTK/common/Types.hpp>
#include <MSTK/fe/types/Centroid.hpp>
#include <MSTK/fe/types/Spectrum.hpp>
#include "unittest.hxx"

using namespace mstk::fe;
using namespace mstk;

struct CentroidTestSuite : vigra::test_suite {
    CentroidTestSuite() : vigra::test_suite("Centroid") {
        add(testCase(&CentroidTestSuite::test));
    }

    void test() {
        // standard tests
        Centroid c;
        // setget
        c.setRetentionTime(c.getRetentionTime());
        c.setMz(c.getMz());
        c.setScanNumber(c.getScanNumber());
        c.setAbundance(c.getAbundance());
        c.setRawData(c.getRawData());
        //
        Spectrum s;
        Centroid k(10.0, 400.0, 5, 1e6, s.begin(), s.end());
        shouldEqual(k.getRetentionTime(), 10.0);
        shouldEqual(k.getMz(), 400.0);
        shouldEqual(k.getScanNumber(), static_cast<UnsignedInt>(5));
        shouldEqual(k.getAbundance(), 1e6);
        shouldEqual(k.getRawData().size(), static_cast<Size>(0));
    }
};

int main()
{
    CentroidTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}


