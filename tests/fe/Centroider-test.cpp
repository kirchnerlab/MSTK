/*
 * Centroider-test.cpp
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
#include <MSTK/fe/Centroider.hpp>
#include <MSTK/fe/SimpleBumpFinder.hpp>
#include <MSTK/fe/GaussianMeanAccumulator.hpp>
#include <MSTK/fe/SumAbundanceAccumulator.hpp>
#include <MSTK/common/Types.hpp>
#include <MSTK/common/Log.hpp>
#include <MSTK/fe/types/Spectrum.hpp>
#include <MSTK/fe/types/Centroid.hpp>
#include <iostream>
#include <boost/assign.hpp>
#include "unittest.hxx"

using namespace boost::assign;

namespace {

// Mock base classes for testing. The mock bump finder simply
// forwards the range.

class MockBumpFinder
{
protected:
    template<typename InputIterator>
    std::pair<InputIterator, InputIterator> findBump(InputIterator first,
        InputIterator last)
    {
        return std::make_pair(first, last);

    }
};

}

struct CentroiderTestSuite : vigra::test_suite
{
    CentroiderTestSuite() :
        vigra::test_suite("Centroider")
    {
        add(testCase(&CentroiderTestSuite::test));
        add(testCase(&CentroiderTestSuite::testUseCase01));
        add(testCase(&CentroiderTestSuite::testUseCase02));
        add(testCase(&CentroiderTestSuite::testUseCase03));
    }

    void test()
    {
        mstk::fe::Spectrum ss;
        ss.push_back(std::make_pair(400.0, 1.0));
        std::vector<mstk::fe::Centroid> centroids;
        mstk::Double rt = 2435.0;
        mstk::UnsignedInt sn = 432;
        typedef mstk::fe::Centroider<mstk::fe::Centroid,
                mstk::fe::SimpleBumpFinder,
                mstk::fe::GaussianMeanAccumulator,
                mstk::fe::SumAbundanceAccumulator> MyCentroider;
        MyCentroider c;
        c(ss.begin(), ss.end(), rt, sn, std::back_inserter(centroids));
        shouldEqual(centroids.size(), static_cast<mstk::Size>(1));
        shouldEqual(centroids[0].getAbundance(), 1.0);
        shouldEqual(centroids[0].getMz(), 400.0);
        mstk::fe::Spectrum raw = centroids[0].getRawData();
        shouldEqual(raw.size(), static_cast<mstk::Size>(1));
        shouldEqual(raw == ss, true);
    }

    void centroid(const std::vector<mstk::Double>& mz,
        const std::vector<mstk::Double>& ab, std::vector<mstk::fe::Centroid>& centroids)
    {
        mstk::fe::Spectrum ss(mz, ab);
        centroids.clear();
        mstk::Double rt = 0.0;
        mstk::UnsignedInt sn = 0;
        typedef mstk::fe::Centroider<mstk::fe::Centroid,
                mstk::fe::SimpleBumpFinder,
                mstk::fe::GaussianMeanAccumulator,
                mstk::fe::SumAbundanceAccumulator> MyCentroider;
        MyCentroider c;
        c(ss.begin(), ss.end(), rt, sn, std::back_inserter(centroids));
    }

    void testUseCase01()
    {
        std::vector<mstk::Double> mz;
        mz += 564.777283669, 564.77981869, 564.782353728, 564.784888783, 564.787423855, 564.789958944, 564.792494051;
        std::vector<mstk::Double> ab;
        ab += 29282.6054688, 84599.90625, 141463.609375, 161800.828125, 131417.28125, 73301.59375, 22801.6269531;
        std::vector<mstk::fe::Centroid> centroids;
        centroid(mz, ab, centroids);
        shouldEqual(centroids.size(), static_cast<mstk::Size>(1));
    }

    void testUseCase02()
    {
        std::vector<mstk::Double> mz;
        mz += 559.786248544, 559.788750036, 559.791251544, 559.793753069, 559.796254611, 559.79875617, 559.801257745, 559.803759337, 560.286880736, 560.289385584, 560.291890449, 560.294395331, 560.296900229, 560.299405145, 560.301910077, 560.304415026;
        std::vector<mstk::Double> ab;
        ab += 2413.484375, 20380.2675781, 49355.7421875, 74698.1875, 77874.4375, 56932.3984375, 29339.2324219, 10459.6621094, 16606.1914062, 28976.7167969, 40484.2109375, 59069.78125, 76599.578125, 71999.65625, 43608.4375, 12770.8847656;
        std::vector<mstk::fe::Centroid> centroids;
        centroid(mz, ab, centroids);
        shouldEqual(centroids.size(), static_cast<mstk::Size>(2));
    }

    void testUseCase03()
    {
        std::vector<mstk::Double> mz;
        mz += 546.762451064, 546.764865768, 546.767280487, 546.769695223, 546.772109974, 546.774524741, 546.776939525, 546.779354324, 546.78176914, 546.789013682, 546.791428561, 546.793843457, 546.796258368, 546.798673295, 546.801088239, 546.803503198, 546.805918173, 546.808333165, 546.810748172, 546.813163196;
        std::vector<mstk::Double> ab;
        ab += 8832.40136719, 54699.3867188, 143617.78125, 236407.625, 272995.125, 227420.328125, 127339.40625, 30389.2695312, 4973.97851562, 14187.8740234, 71039.6640625, 358838.6875, 835839.375, 1241756, 1284578.125, 930160.375, 440593.1875, 112644.351562, 33855.8945312, 16687.5585938;
        std::vector<mstk::fe::Centroid> centroids;
        centroid(mz, ab, centroids);
        shouldEqual(centroids.size(), static_cast<mstk::Size>(2));
    }
};

int main()
{
    CentroiderTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

