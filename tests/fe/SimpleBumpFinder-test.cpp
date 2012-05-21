/*
 * SimpleBumpFinder-test.cpp
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
#include <MSTK/config.hpp>
#include "unittest.hxx"

#include <MSTK/fe/SimpleBumpFinder.hpp>
#include <MSTK/fe/types/Spectrum.hpp>
#include <MSTK/common/Log.hpp>
#include <iostream>
#include <boost/assign.hpp>

using namespace mstk;
using namespace mstk::fe;
using namespace boost::assign;

class MyBumpFinder : public fe::SimpleBumpFinder
{
public:
    typedef Spectrum::const_iterator SSCI;
    MyBumpFinder() :
        fe::SimpleBumpFinder()
    {
    }

    std::pair<SSCI, SSCI> dance(SSCI first, SSCI last)
    {
        return findBump(first, last);
    }
};

struct SimpleBumpFinderTestSuite : vigra::test_suite
{
    // typedefs
    typedef Spectrum::const_iterator SSCI;
    typedef std::pair<SSCI, SSCI> Range;

    SimpleBumpFinderTestSuite() :
        vigra::test_suite("SimpleBumpFinder")
    {
        add(testCase(&SimpleBumpFinderTestSuite::testSinglet));
        add(testCase(&SimpleBumpFinderTestSuite::testDoublets));
        add(testCase(&SimpleBumpFinderTestSuite::testTriplets));
        add(testCase(&SimpleBumpFinderTestSuite::testUseCase01));
    }

    void testSinglet()
    {
        MSTK_LOG(logINFO)
                << "testSinglet()";
        MyBumpFinder bf;
        Spectrum ss;
        ss.push_back(SpectrumElement(400.0, 1.0));
        Range r = bf.dance(ss.begin(), ss.end());
        shouldEqual(r.first == ss.begin(), true);
        shouldEqual(r.second == ss.end(), true);
    }

    void testDoublets()
    {
        MSTK_LOG(logINFO)
                << "testDoublets(): up";
        MyBumpFinder bf;
        Spectrum ss;
        // up
        ss.push_back(SpectrumElement(400.0, 1.0));
        ss.push_back(SpectrumElement(400.001, 2.0));
        Range r = bf.dance(ss.begin(), ss.end());
        shouldEqual(r.first == ss.begin(), true);
        shouldEqual(r.second == ss.end(), true);
        // down
        MSTK_LOG(logINFO)
                << "testDoublets(): down";
        ss.clear();
        ss.push_back(SpectrumElement(400.0, 2.0));
        ss.push_back(SpectrumElement(400.001, 1.0));
        r = bf.dance(ss.begin(), ss.end());
        shouldEqual(r.first == ss.begin(), true);
        shouldEqual(r.second == ss.end(), true);
    }

    void testTriplets()
    {
        MSTK_LOG(logINFO)
                << "testTriplets()";
        MyBumpFinder bf;
        Spectrum ss;
        // up
        ss.push_back(SpectrumElement(400.0, 1.0));
        ss.push_back(SpectrumElement(400.001, 2.0));
        ss.push_back(SpectrumElement(400.002, 3.0));
        Range r = bf.dance(ss.begin(), ss.end());
        shouldEqual(r.first == ss.begin(), true);
        shouldEqual(r.second == ss.end(), true);
        // down
        ss.clear();
        ss.push_back(SpectrumElement(400.0, 3.0));
        ss.push_back(SpectrumElement(400.001, 2.0));
        ss.push_back(SpectrumElement(400.002, 1.0));
        r = bf.dance(ss.begin(), ss.end());
        shouldEqual(r.first == ss.begin(), true);
        shouldEqual(r.second == ss.end(), true);
        // local max
        ss.clear();
        ss.push_back(SpectrumElement(400.0, 1.0));
        ss.push_back(SpectrumElement(400.001, 2.0));
        ss.push_back(SpectrumElement(400.002, 1.0));
        r = bf.dance(ss.begin(), ss.end());
        shouldEqual(r.first == ss.begin(), true);
        shouldEqual(r.second == ss.end(), true);
        // local min
        ss.clear();
        ss.push_back(SpectrumElement(400.0, 2.0));
        ss.push_back(SpectrumElement(400.001, 1.0));
        ss.push_back(SpectrumElement(400.002, 2.0));
        r = bf.dance(ss.begin(), ss.end());
        shouldEqual(r.first == ss.begin(), true);
        shouldEqual(r.second == ss.begin()+2, true);
    }

    void testUseCase01()
    {
        std::vector<double> mz;
        mz += 564.777283669, 564.77981869, 564.782353728, 564.784888783, 564.787423855, 564.789958944, 564.792494051;
        std::vector<double> ab;
        ab += 29282.6054688, 84599.90625, 141463.609375, 161800.828125, 131417.28125, 73301.59375, 22801.6269531;
        Spectrum ss(mz, ab);
        MyBumpFinder bf;
        Range r = bf.dance(ss.begin(), ss.end());
        shouldEqual(r.first == ss.begin(), true);
        shouldEqual(r.second == ss.end(), true);
    }
};

int main()
{
    SimpleBumpFinderTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

