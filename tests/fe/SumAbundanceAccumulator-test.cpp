/*
 * SumAbundanceAccumulator-test.cpp
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
#include <MSTK/fe/SumAbundanceAccumulator.hpp>
#include <MSTK/common/Log.hpp>
#include <iostream>
#include <iomanip>

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
class MyAccumulator : public SumAbundanceAccumulator
{
public:
    typedef Spectrum::const_iterator SSCI;

    MyAccumulator() :
        SumAbundanceAccumulator()
    {
    }

    double myMean(SSCI first, SSCI last)
    {
        return abundance(first, last);
    }
};

struct SumAbundanceAccumulatorTestSuite : vigra::test_suite
{
    SumAbundanceAccumulatorTestSuite() :
        vigra::test_suite("SumAbundanceAccumulator")
    {
        add(testCase(&SumAbundanceAccumulatorTestSuite::test));
    }

    void test()
    {
        MyAccumulator a;
        {
            Spectrum ss;
            ss.push_back(SpectrumElement(400.0, 1.0));
            ss.push_back(SpectrumElement(401.0, 2.0));
            double ab = a.myMean(ss.begin(), ss.end());
            shouldEqual(ab, 3.0);
        }
        {
            Spectrum ss;
            ss.push_back(SpectrumElement(401.0, 1.0));
            ss.push_back(SpectrumElement(400.0, 2.0));
            double ab = a.myMean(ss.begin(), ss.end());
            shouldEqual(ab, 3.0);
        }
    }
};

int main()
{
    SumAbundanceAccumulatorTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

