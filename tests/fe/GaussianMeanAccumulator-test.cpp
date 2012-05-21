/*
 * GaussianMeanAccumulator-test.cpp
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
#include <MSTK/fe/types/Spectrum.hpp>
#include <MSTK/fe/GaussianMeanAccumulator.hpp>
#include <MSTK/common/Log.hpp>
#include <iostream>
#include <iomanip>

using namespace mstk;
using namespace mstk::fe;

//
// Accessor class for our spectrum type
//
struct SpectrumAccessor
{
    double getAbundance(const Spectrum::Element& e) const
    {
        return e.abundance;
    }
    double getMz(const Spectrum::Element& e) const
    {
        return e.mz;
    }
};

//
// derive our own accumulator (essentially exposing the otherwise
// protected interface).
//
class MyAccumulator : public fe::GaussianMeanAccumulator
{
public:
    typedef Spectrum::const_iterator SSCI;

    MyAccumulator() :
        fe::GaussianMeanAccumulator()
    {
    }

    double myMean(SSCI first, SSCI last)
    {
        return mean(first, last);
    }

    void myTrimAndMax(SSCI first, SSCI last, SSCI& start, SSCI& stop, SSCI& m)
    {
        return trimAndMax(first, last, start, stop, m);
    }
};

struct GaussianMeanAccumulatorTestSuite : vigra::test_suite
{
    GaussianMeanAccumulatorTestSuite() :
        vigra::test_suite("GaussianMeanAccumulator")
    {
        add(testCase(&GaussianMeanAccumulatorTestSuite::testTrimAndMax));
        add(testCase(&GaussianMeanAccumulatorTestSuite::testSingle));
        add(testCase(&GaussianMeanAccumulatorTestSuite::testAverageOfTwo));
        add(testCase(&GaussianMeanAccumulatorTestSuite::testRamp));
        add(testCase(&GaussianMeanAccumulatorTestSuite::testGaussianFit));
    }

    void testTrimAndMax()
    {
        MSTK_LOG(logINFO)
                << "testTrimAndMax()";
        MyAccumulator a;
        typedef Spectrum::const_iterator SSCI;
        {
            // empty spectrum should leave all iterators pointing to last
            Spectrum ss;
            SSCI first, last, m;
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first==last, true);
            shouldEqual(m==last, true);
        }
        {
            // same with zero-only abundances
            Spectrum ss;
            // one entry
            ss.push_back(SpectrumElement(400.0, 0.0));
            SSCI first, last, m;
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first==last, true);
            shouldEqual(m==last, true);
            // two entries
            ss.push_back(SpectrumElement(401.0, 0.0));
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first==last, true);
            shouldEqual(m==last, true);
            // many entries
            ss.push_back(SpectrumElement(402.0, 0.0));
            ss.push_back(SpectrumElement(403.0, 0.0));
            ss.push_back(SpectrumElement(404.0, 0.0));
            ss.push_back(SpectrumElement(405.0, 0.0));
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first==last, true);
            shouldEqual(m==last, true);
        }
        {
            // for non-zero abundances, iterators should be left unaltered
            Spectrum ss;
            // one entry
            ss.push_back(SpectrumElement(400.0, 1.0));
            SSCI first, last, m;
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first == ss.begin(), true);
            shouldEqual(m == ss.begin(), true);
            shouldEqual(last == ss.end(), true);
            // two entries
            ss.push_back(SpectrumElement(401.0, 1.0));
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first == ss.begin(), true);
            shouldEqual(m == ss.begin(), true);
            shouldEqual(last == ss.end(), true);
            // many entries
            ss.push_back(SpectrumElement(402.0, 1.0));
            ss.push_back(SpectrumElement(403.0, 1.0));
            ss.push_back(SpectrumElement(404.0, 1.0));
            ss.push_back(SpectrumElement(405.0, 1.0));
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first == ss.begin(), true);
            shouldEqual(m == ss.begin(), true);
            shouldEqual(last == ss.end(), true);
        }
        {
            // a real use case: remove leading zeros
            // one entry
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 0.0));
            ss.push_back(SpectrumElement(401.0, 1.0));
            SSCI first, last, m;
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first == ++(ss.begin()), true);
            shouldEqual(m == ++(ss.begin()), true);
            shouldEqual(last == ss.end(), true);
            // two entries
            ss.clear();
            ss.push_back(SpectrumElement(400.0, 0.0));
            ss.push_back(SpectrumElement(401.0, 0.0));
            ss.push_back(SpectrumElement(402.0, 1.0));
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first == ss.begin()+2, true);
            shouldEqual(m == ss.begin()+2, true);
            shouldEqual(last == ss.end(), true);
        }
        {
            // a real use case: remove trailing zeros
            // one entry
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 1.0));
            ss.push_back(SpectrumElement(401.0, 0.0));
            SSCI first, last, m;
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first == ss.begin(), true);
            shouldEqual(m == ss.begin(), true);
            shouldEqual(last == ss.begin()+1, true);
            // two entries
            ss.clear();
            ss.push_back(SpectrumElement(400.0, 1.0));
            ss.push_back(SpectrumElement(401.0, 0.0));
            ss.push_back(SpectrumElement(402.0, 0.0));
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first == ss.begin(), true);
            shouldEqual(m == ss.begin(), true);
            shouldEqual(last == ss.begin()+1, true);
        }
        {
            // a real use case: remove leading and trailing zeros
            // one entry
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 0.0));
            ss.push_back(SpectrumElement(401.0, 1.0));
            ss.push_back(SpectrumElement(402.0, 0.0));
            SSCI first, last, m;
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first == ss.begin()+1, true);
            shouldEqual(m == ss.begin()+1, true);
            shouldEqual(last == ss.begin()+2, true);
            // two entries
            ss.clear();
            ss.push_back(SpectrumElement(400.0, 0.0));
            ss.push_back(SpectrumElement(401.0, 0.0));
            ss.push_back(SpectrumElement(402.0, 1.0));
            ss.push_back(SpectrumElement(403.0, 1.0));
            ss.push_back(SpectrumElement(404.0, 0.0));
            ss.push_back(SpectrumElement(405.0, 0.0));
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(first == ss.begin()+2, true);
            shouldEqual(m == ss.begin()+2, true);
            shouldEqual(last == ss.begin()+4, true);
        }
        {
            // test the extreme cases for the max functionality
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 1.0));
            SSCI first, last, m;
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(m == ss.begin(), true);
            // two entries: '\'
            ss.clear();
            ss.push_back(SpectrumElement(401.0, 1.0));
            ss.push_back(SpectrumElement(402.0, 0.5));
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(m == ss.begin(), true);
            // two entries: '/'
            ss.clear();
            ss.push_back(SpectrumElement(401.0, 0.5));
            ss.push_back(SpectrumElement(402.0, 1.0));
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(m == ss.begin()+1, true);
            // two entries: '--'
            ss.clear();
            ss.push_back(SpectrumElement(401.0, 1.0));
            ss.push_back(SpectrumElement(402.0, 1.0));
            a.myTrimAndMax(ss.begin(), ss.end(), first, last, m);
            shouldEqual(m == ss.begin(), true);
        }
    }

    void testSingle()
    {
        MSTK_LOG(logINFO)
                << "testSingle()";
        MyAccumulator a;
        Spectrum ss;
        ss.push_back(SpectrumElement(400.0, 1.0));
        double m = a.myMean(ss.begin(), ss.end());
        shouldEqual(m, 400.0);
    }

    void testAverageOfTwo()
    {
        MSTK_LOG(logINFO)
                << "testAverageofTwo()";
        MyAccumulator a;
        {
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 1.0));
            ss.push_back(SpectrumElement(401.0, 2.0));
            double expectedMass = (400.0 + 2.0 * 401.0) / 3.0;
            double m = a.myMean(ss.begin(), ss.end());
            shouldEqual(m, expectedMass);
        }
        {
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 2.0));
            ss.push_back(SpectrumElement(401.0, 1.0));
            double expectedMass = (2.0 * 400.0 + 401.0) / 3.0;
            double m = a.myMean(ss.begin(), ss.end());
            shouldEqual(m, expectedMass);
        }
        {
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 1.0));
            ss.push_back(SpectrumElement(401.0, 1.0));
            double expectedMass = (400.0 + 401.0) / 2.0;
            double m = a.myMean(ss.begin(), ss.end());
            shouldEqual(m, expectedMass);
        }
    }

    void testRamp()
    {
        MSTK_LOG(logINFO)
                << "testRamp()";
        MyAccumulator a;
        {
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 1.0));
            ss.push_back(SpectrumElement(401.0, 2.0));
            ss.push_back(SpectrumElement(402.0, 3.0));
            double expectedMass = (2.0 * 401.0 + 3.0 * 402.0) / 5.0;
            double m = a.myMean(ss.begin(), ss.end());
            shouldEqualTolerance(m, expectedMass, 1e-12);
        }
        {
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 3.0));
            ss.push_back(SpectrumElement(401.0, 2.0));
            ss.push_back(SpectrumElement(402.0, 1.0));
            double expectedMass = (3.0 * 400.0 + 2.0 * 401.0) / 5.0;
            double m = a.myMean(ss.begin(), ss.end());
            shouldEqualTolerance(m, expectedMass, 1e-12);
        }
    }

    void testGaussianFit()
    {
        MSTK_LOG(logINFO)
                << "testGaussianFit()";
        MyAccumulator a;
        {
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 1.0));
            ss.push_back(SpectrumElement(401.0, 2.0));
            ss.push_back(SpectrumElement(402.0, 1.0));
            double expectedMass = 401.0;
            double m = a.myMean(ss.begin(), ss.end());
            shouldEqualTolerance(m, expectedMass, 1e-12);
        }
        {
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 1.0));
            ss.push_back(SpectrumElement(401.0, 3.0));
            ss.push_back(SpectrumElement(402.0, 2.0));
            double expectedMass = 401.23042271;
            double m = a.myMean(ss.begin(), ss.end());
            shouldEqualTolerance(m, expectedMass, 1e-12);
        }
    }
};

int main()
{
    GaussianMeanAccumulatorTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

