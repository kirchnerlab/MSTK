/*
 * Splitter-test.cpp
 *
 * Copyright (c) 2011 Marc Kirchner
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
#include <MSTK/fe/Splitter.hpp>
#include <MSTK/common/Types.hpp>
#include <MSTK/common/Log.hpp>
#include <MSTK/fe/types/Spectrum.hpp>
#include <iostream>
#include <vector>
#include <boost/assign.hpp>
#include "unittest.hxx"

// get vector += in scope
using namespace boost::assign;
using namespace mstk;
using namespace mstk::fe;

struct PpmPeakShapeFunction
{
    PpmPeakShapeFunction(const mstk::Double ppm) : ppm_(ppm) {}
    mstk::Double getSupportThreshold(const mstk::Double mz) const
    {
        return ppm_*mz*1e-6;
    }
    mstk::Double ppm_;
};

struct SplitterTestSuite : vigra::test_suite
{
    SplitterTestSuite() :
        vigra::test_suite("Splitter")
    {
        add(testCase(&SplitterTestSuite::test));
        add(testCase(&SplitterTestSuite::testUseCase01));
    }

    void test()
    {
        PpmPeakShapeFunction psf(5.0);
        Splitter<Spectrum, PpmPeakShapeFunction> splitter(psf);
        Spectrum ss;

        // test: empty
        splitter.assign(ss.begin(), ss.end());
        shouldEqual(splitter.size(), static_cast<Size>(0));

        // test: one group w/ two entries
        ss.push_back(SpectrumElement(546.0, 1.0));
        ss.push_back(SpectrumElement(546.0001, 1.0));
        splitter.assign(ss.begin(), ss.end());
        shouldEqual(splitter.size(), static_cast<Size>(1));

        // test: two groups
        typedef Splitter<Spectrum, PpmPeakShapeFunction
            >::value_type Range;
        typedef Splitter<Spectrum, PpmPeakShapeFunction
            >::const_iterator CI;
        ss.push_back(SpectrumElement(547.0, 1.0));
        splitter.assign(ss.begin(), ss.end());
        shouldEqual(splitter.size(), static_cast<Size>(2));
        std::vector<UnsignedInt> sizes;
        sizes += 2, 1;
        Size n = 0;
        for (CI i = splitter.begin(); i != splitter.end(); ++i, ++n) {
            shouldEqual(std::distance(i->first, i->second), sizes[n]);
        }

        // test: more than two
        ss.push_back(SpectrumElement(548.0, 1.0));
        ss.push_back(SpectrumElement(548.0001, 1.0));
        ss.push_back(SpectrumElement(548.0002, 1.0));
        splitter.assign(ss.begin(), ss.end());
        shouldEqual(splitter.size(), static_cast<Size>(3));
        std::vector<UnsignedInt> sizes2;
        sizes2 += 2, 1, 3;
        n = 0;
        for (CI i = splitter.begin(); i != splitter.end(); ++i, ++n) {
            shouldEqual(std::distance(i->first, i->second), sizes2[n]);
        }
    }

    void testUseCase01()
    {
        MSTK_LOG(logINFO) << "testUseCase01:";
        std::vector<Double> mz;
        mz += 546.762451064, 546.764865768, 546.767280487, 546.769695223,
           546.772109974, 546.774524741, 546.776939525, 546.779354324,
           546.78176914, 546.789013682, 546.791428561, 546.793843457,
           546.796258368, 546.798673295, 546.801088239, 546.803503198,
           546.805918173, 546.808333165, 546.810748172, 546.813163196;
        std::vector<Double> ab;
        ab += 8832.40136719, 54699.3867188, 143617.78125, 236407.625,
           272995.125, 227420.328125, 127339.40625, 30389.2695312,
           4973.97851562, 14187.8740234, 71039.6640625, 358838.6875, 835839.375,
           1241756, 1284578.125, 930160.375, 440593.1875, 112644.351562,
           33855.8945312, 16687.5585938;
        Spectrum ss(mz, ab);
        PpmPeakShapeFunction psf(5.0);
        Splitter<Spectrum, PpmPeakShapeFunction
            > splitter(psf);
        splitter.assign(ss.begin(), ss.end());
        typedef Splitter<Spectrum, PpmPeakShapeFunction
            >::const_iterator CI;
        shouldEqual(splitter.size(), static_cast<Size>(2));
    }
};

int main()
{
    SplitterTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

