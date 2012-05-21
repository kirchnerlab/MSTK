/*
 * XicExtractor-test.cpp
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
#include "unittest.hxx"
#include "utilities.hpp"
#include <iostream>
#include <vector>
#include <MSTK/common/Types.hpp>
#include <MSTK/fe/XicExtractor.hpp>
#include <MSTK/fe/types/Centroid.hpp>
#include <MSTK/fe/CentroidTraits.hpp>
#include <MSTK/fe/types/CentroidFbiTraits.hpp>
#include <MSTK/fe/types/Xic.hpp>
#include <MSTK/fe/CentroidWeightedMeanDisambiguator.hpp>
#include <MSTK/fe/RunningMeanSmoother.hpp>
#include <MSTK/fe/XicLocalMinSplitter.hpp>
#include <fbi/fbi.h>

using namespace mstk::fe;
using namespace mstk;

struct XicExtractorTestSuite : vigra::test_suite
{
    typedef std::vector<Xic> Xics;
    typedef std::vector<Centroid> Centroids;

    XicExtractorTestSuite() :
            vigra::test_suite("XicExtractor")
    {
        add(testCase(&XicExtractorTestSuite::testNormalMax));
        add(testCase(&XicExtractorTestSuite::testRampDown));
        add(testCase(&XicExtractorTestSuite::testRampUp));
        add(testCase(&XicExtractorTestSuite::testSplit));
        add(testCase(&XicExtractorTestSuite::testSplit2));
    }

    void testNormalMax()
    {
        double mz[] = { 100.0, 100.0, 100.0, 100.0, 100.0 };
        double rt[] = { 350.0, 352.0, 354.0, 356.0, 358.0 };
        unsigned int sn[] = { 42, 43, 44, 45, 46 };
        double ab[] = { 0.2, 0.6, 1.0, 0.6, 0.2 };
        Centroids cs = makeCentroids(5, mz, rt, sn, ab);
        Xics xs;

        typedef XicExtractor<CentroidWeightedMeanDisambiguator,
                RunningMeanSmoother, XicLocalMinSplitter<Xic> > MyXicExtractor;
        MyXicExtractor xe;
        CentroidBoxGenerator bg(3, 10.0);
        Size n = xe(cs, bg, 3, 0.76, xs);
        shouldEqual(xs.size(), n);
        shouldEqual(xs.size(), static_cast<size_t>(1));
    }

    void testRampDown()
    {
        double mz[] = { 100.0, 100.0, 100.0, 100.0, 100.0 };
        double rt[] = { 350.0, 352.0, 354.0, 356.0, 358.0 };
        unsigned int sn[] = { 42, 43, 44, 45, 46 };
        double ab[] = { 1.0, 0.8, 0.6, 0.4, 0.2 };
        Centroids cs = makeCentroids(5, mz, rt, sn, ab);
        Xics xs;

        typedef XicExtractor<CentroidWeightedMeanDisambiguator,
                RunningMeanSmoother, XicLocalMinSplitter<Xic> > MyXicExtractor;
        MyXicExtractor xe;
        CentroidBoxGenerator bg(3, 10.0);

        // require too many centroids
        Size n = xe(cs, bg, 6, 0.76, xs);
        shouldEqual(xs.size(), static_cast<size_t>(0));
        shouldEqual(xs.size(), n);
        // require the right number of centroids
        n = xe(cs, bg, 5, 0.76, xs);
        shouldEqual(xs.size(), n);
        shouldEqual(xs.size(), static_cast<size_t>(1));
        //std::cerr << "Got: " << xs[0].size() << std::endl;
    }

    void testRampUp()
    {
        double mz[] = { 100.0, 100.0, 100.0, 100.0, 100.0 };
        double rt[] = { 350.0, 352.0, 354.0, 356.0, 358.0 };
        unsigned int sn[] = { 42, 43, 44, 45, 46 };
        double ab[] = { 0.2, 0.4, 0.6, 0.8, 1.0 };
        Centroids cs = makeCentroids(5, mz, rt, sn, ab);
        Xics xs;

        typedef XicExtractor<CentroidWeightedMeanDisambiguator,
                RunningMeanSmoother, XicLocalMinSplitter<Xic> > MyXicExtractor;
        MyXicExtractor xe;
        CentroidBoxGenerator bg(3, 10.0);

        // require too many centroids
        Size n = xe(cs, bg, 6, 0.76, xs);
        shouldEqual(xs.size(), static_cast<size_t>(0));
        shouldEqual(xs.size(), n);
        // require the right number of centroids
        n = xe(cs, bg, 5, 0.76, xs);
        shouldEqual(xs.size(), n);
        shouldEqual(xs.size(), static_cast<size_t>(1));
        //std::cerr << "Got: " << xs[0].size() << std::endl;
    }

    void testSplit()
    {
        double mz[] = { 100.001, 100.003, 100.002, 100.005, 100.001, 100.003,
                        100.001, 100.004, 100.0, 100.01 };
        double rt[] = { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
                        19.0 };
        unsigned int sn[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        double ab[] = { 1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 0.5, 0.1 };

        Centroids cs = makeCentroids(10, mz, rt, sn, ab);
        Xics xs;

        typedef XicExtractor<CentroidWeightedMeanDisambiguator,
                RunningMeanSmoother, XicLocalMinSplitter<Xic> > MyXicExtractor;
        MyXicExtractor xe;
        CentroidBoxGenerator bg(2, 100.0);

        Size n = xe(cs, bg, 3, 0.76, xs);
        shouldEqual(xs.size(), n);
        shouldEqual(xs.size(), static_cast<size_t>(2));
    }

    void testSplit2()
    {
        // some real-world data that should lead to a single XIC
        double mz[] = { 548.814, 548.813, 548.814, 548.813, 548.813, 548.813,
                        548.813, 548.813, 548.813, 548.813, 548.813, 548.813,
                        548.813, 548.813, 548.813, 548.813, 548.813, 548.813,
                        548.813, 548.813, 548.814 };
        double rt[] = { 2111.24, 2112.09, 2113.31, 2114.72, 2115.56, 2116.55,
                        2118.23, 2119.06, 2120.6, 2121.75, 2122.92, 2124.33,
                        2125.19, 2126.43, 2127.59, 2128.45, 2129.41, 2130.38,
                        2131.8, 2133.23, 2134.37 };
        unsigned int sn[] = { 4000, 4002, 4005, 4009, 4011, 4013, 4017, 4019,
                              4023, 4026, 4029, 4033, 4035, 4038, 4041, 4043,
                              4045, 4047, 4048, 4052, 4055 };
        double ab[] = { 472009.0, 905473.0, 1291190.0, 1828580.0, 1817710.0,
                        2244620.0, 3388290.0, 4188030.0, 5001190.0, 4322550.0,
                        4969260.0, 4725260.0, 4004990.0, 3754450.0, 2327370.0,
                        2761370.0, 1909930.0, 1926350.0, 756566.0, 400389.0,
                        242239 };
        Centroids cs = makeCentroids(21, mz, rt, sn, ab);
        Xics xs;

        typedef XicExtractor<CentroidWeightedMeanDisambiguator,
                RunningMeanSmoother, XicLocalMinSplitter<Xic> > MyXicExtractor;
        MyXicExtractor xe;
        CentroidBoxGenerator bg(2, 7.0);

        Size n = xe(cs, bg, 3, 0.76, xs);
        shouldEqual(xs.size(), n);
        shouldEqual(xs.size(), static_cast<size_t>(1));
    }
};

int main()
{
    XicExtractorTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

