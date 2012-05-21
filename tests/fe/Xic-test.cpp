/*
 * Xic-test.cpp
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

// Expose complete class for testing.
#define protected public
#define private   public
#include <MSTK/fe/types/Xic.hpp>
#undef protected
#undef private

#include <iostream>

#include "unittest.hxx"
#include "utilities.hpp"

#include <MSTK/fe/types/Centroid.hpp>
#include <numeric>

using namespace mstk::fe;

struct XicTestSuite : vigra::test_suite
{
    XicTestSuite() :
        vigra::test_suite("Xic")
    {
        add(testCase(&XicTestSuite::testConstructor));
        add(testCase(&XicTestSuite::testOps));
        add(testCase(&XicTestSuite::testCorrelate));
        add(testCase(&XicTestSuite::testGetter));
        add(testCase(&XicTestSuite::testRecalculate));
        add(testCase(&XicTestSuite::testMergeDuplicates));
        add(testCase(&XicTestSuite::testSplitRt1));
        add(testCase(&XicTestSuite::testSplitRt2));
        add(testCase(&XicTestSuite::testSplitRt3));
        add(testCase(&XicTestSuite::testSplitRt4));
        add(testCase(&XicTestSuite::testSplitRt5));
        add(testCase(&XicTestSuite::testGetSmoothedXic));
    }

    void testConstructor()
    {
        Xic xic;
        shouldEqual(xic.getAbundance(), 0.0);
        shouldEqual(xic.getMz(), 0.0);
        shouldEqual(xic.getRetentionTime(), 0.0);
        shouldEqual(xic.getMzTolerance(), 0.0);
        shouldEqual(xic.getRetentionTimeTolerance(), 0.0);
    }

    void testOps()
    {
        double mz[] = { 99.9, 100.0, 100.1 };
        double rt[] = { 350.0, 352.0, 354.0 };
        unsigned int sn[] = { 42, 43, 44 };
        double ab[] = { 1.0, 1.0, 1.0 };
        Xic x1 = makeXic(3, mz, rt, sn, ab);
        double mz2[] = { 99.9 + 1.0, 100.0 + 1.0, 100.1 + 1.0 };
        double rt2[] = { 350.0 + 1.0, 352.0 + 1.0, 354.0 + 1.0 };
        unsigned int sn2[] = { 42, 43, 44 };
        double ab2[] = { 2.0, 2.0, 2.0 };
        Xic x2 = makeXic(3, mz2, rt2, sn2, ab2);
        // operator==
        shouldEqual(x1 == x1, true);
        shouldEqual(x2 == x2, true);
        Xic x3 = x1;
        shouldEqual(x1 == x3, true);
        shouldEqual(x1 == x2, false);
        // LessThan...
        Xic::LessThanAbundance lta;
        shouldEqual(lta(x1, x2), true);
        shouldEqual(lta(x1, x1), false);
        shouldEqual(lta(x2, x2), false);
        shouldEqual(lta(x2, x1), false);
        Xic::LessThanRt ltr;
        shouldEqual(ltr(x1, x2), true);
        shouldEqual(ltr(x1, x1), false);
        shouldEqual(ltr(x2, x2), false);
        shouldEqual(ltr(x2, x1), false);
        Xic::LessThanMz ltm;
        shouldEqual(ltm(x1, x2), true);
        shouldEqual(ltm(x1, x1), false);
        shouldEqual(ltm(x2, x2), false);
        shouldEqual(ltm(x2, x1), false);
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
        shouldEqualTolerance(x1.correlate(x2), 1.0, 1e-12);
        shouldEqualTolerance(x2.correlate(x1), 1.0, 1e-12);
        unsigned int sn3[] = { 45, 46, 47 };
        Xic x3 = makeXic(3, mz, rt, sn3, ab2);
        shouldEqualTolerance(x1.correlate(x3), 0.0, 1e-12);
        shouldEqualTolerance(x3.correlate(x1), 0.0, 1e-12);
    }

    void testGetter()
    {
        double mz[] = { 100.0 };
        double rt[] = { 350.0 };
        unsigned int sn[] = { 42 };
        double ab[] = { 1.0 };
        std::vector<Centroid> cs = makeCentroids(1, mz, rt, sn, ab);
        Xic xic;
        xic.insert(xic.end(), cs.begin(), cs.end());
        xic.recalculate();
        shouldEqual(xic.getAbundance(), 1.0);
        shouldEqual(xic.getMz(), 100.0);
        shouldEqual(xic.getRetentionTime(), 350.0);
        shouldEqual(xic.getMzTolerance(), 0.0);
        shouldEqual(xic.getRetentionTimeTolerance(), 0.0);
    }

    void testRecalculate()
    {
        double mz[] = { 99.9, 100.0, 100.1 };
        double rt[] = { 350.0, 352.0, 354.0 };
        unsigned int sn[] = { 42, 43, 44 };
        double ab[] = { 1.0, 1.0, 1.0 };
        std::vector<Centroid> cs = makeCentroids(3, mz, rt, sn, ab);
        Xic xic;
        xic.insert(xic.end(), cs.begin(), cs.end());
        // std::cerr << xic << std::endl;
        xic.recalculate();
        // std::cerr << xic << std::endl;
        shouldEqual(xic.getAbundance(), 4.0);
        shouldEqual(xic.getMz(), 100.0);
        shouldEqual(xic.getRetentionTime(), 352.0);
        shouldEqualTolerance(xic.getMzTolerance(), 0.1, 1e-10);
        shouldEqualTolerance(xic.getRetentionTimeTolerance(), 2.0, 1e-10);
    }

    void testGetSmoothedXic()
    {
        double mzs1[] =
                { 100.001, 100.003, 100.002, 100.005, 100.001, 100.003 };
        double rts1[] = { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0 };
        unsigned int sns1[] = { 1, 2, 3, 4, 5, 6 };
        double abs1[] = { 1.0, 2.0, 3.0, 2.0, 1.0, 0.5 };
        Xic xic = makeXic(6, mzs1, rts1, sns1, abs1);
        Xic sxic;
        xic.getSmoothedXic(sxic);
        double expectedAbs[] = { 1.0, 2.0, 2.333333, 2, 1.166666, 0.5 };
        for (size_t i = 0; i < 6; ++i) {
            shouldEqualTolerance(sxic[i].getAbundance(), expectedAbs[i], 1e-6);
        }
    }

    void testMergeDuplicates()
    {
        {
            // no dupes; should not do anything
            double mzs1[] =
                    { 100.001, 100.003, 100.002, 100.005, 100.001, 100.003 };
            double rts1[] = { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0 };
            unsigned int sns1[] = { 1, 2, 3, 4, 5, 6 };
            double abs1[] = { 1.0, 2.0, 3.0, 2.0, 1.0, 0.5 };
            Xic xic = makeXic(6, mzs1, rts1, sns1, abs1);
            xic.mergeDuplicates();
            shouldEqual(xic.size(), static_cast<size_t>(6));
            for (size_t i = 0; i < 6; ++i) {
                shouldEqualTolerance(xic[i].getAbundance(), abs1[i], 1e-12);
                shouldEqualTolerance(xic[i].getMz(), mzs1[i], 1e-12);
                shouldEqualTolerance(xic[i].getRetentionTime(), rts1[i], 1e-12);
                shouldEqual(xic[i].getScanNumber(), sns1[i]);
            }
        }
        {
            // with dupes
            double mzs1[] =
                    { 100.001, 100.004, 100.002, 100.005, 100.001, 100.003 };
            double rts1[] = { 10.0, 11.0, 11.0, 13.0, 14.0, 15.0 };
            unsigned int sns1[] = { 1, 2, 2, 4, 5, 6 };
            double abs1[] = { 1.0, 2.0, 2.0, 2.0, 1.0, 0.5 };
            Xic xic = makeXic(6, mzs1, rts1, sns1, abs1);
            xic.mergeDuplicates();
            shouldEqual(xic.size(), static_cast<size_t>(5));
            double emzs1[] =
                    { 100.001, 100.003, 100.005, 100.001, 100.003 };
            double erts1[] = { 10.0, 11.0, 13.0, 14.0, 15.0 };
            unsigned int esns1[] = { 1, 2, 4, 5, 6 };
            double eabs1[] = { 1.0, 4.0, 2.0, 1.0, 0.5 };

            for (size_t i = 0; i < 5; ++i) {
                shouldEqualTolerance(xic[i].getAbundance(), eabs1[i], 1e-12);
                shouldEqualTolerance(xic[i].getMz(), emzs1[i], 1e-12);
                shouldEqualTolerance(xic[i].getRetentionTime(), erts1[i], 1e-12);
                shouldEqual(xic[i].getScanNumber(), esns1[i]);
            }
        }
        {
            // with zero abundance dupes
            double mzs1[] =
                    { 100.001, 100.004, 100.002, 100.005, 100.001, 100.003 };
            double rts1[] = { 10.0, 11.0, 11.0, 13.0, 14.0, 15.0 };
            unsigned int sns1[] = { 1, 2, 2, 4, 5, 6 };
            double abs1[] = { 1.0, 0.0, 0.0, 2.0, 1.0, 0.5 };
            Xic xic = makeXic(6, mzs1, rts1, sns1, abs1);
            xic.mergeDuplicates();
            shouldEqual(xic.size(), static_cast<size_t>(5));
            double emzs1[] =
                    { 100.001, 100.003, 100.005, 100.001, 100.003 };
            double erts1[] = { 10.0, 11.0, 13.0, 14.0, 15.0 };
            unsigned int esns1[] = { 1, 2, 4, 5, 6 };
            double eabs1[] = { 1.0, 0.0, 2.0, 1.0, 0.5 };

            for (size_t i = 0; i < 5; ++i) {
                shouldEqualTolerance(xic[i].getAbundance(), eabs1[i], 1e-12);
                shouldEqualTolerance(xic[i].getMz(), emzs1[i], 1e-12);
                shouldEqualTolerance(xic[i].getRetentionTime(), erts1[i], 1e-12);
                shouldEqual(xic[i].getScanNumber(), esns1[i]);
            }
        }
        {
            // with dupes, one of which is zero abundance
            double mzs1[] =
                    { 100.001, 100.004, 100.002, 100.005, 100.001, 100.003 };
            double rts1[] = { 10.0, 11.0, 11.0, 13.0, 14.0, 15.0 };
            unsigned int sns1[] = { 1, 2, 2, 4, 5, 6 };
            double abs1[] = { 1.0, 0.0, 1.0, 2.0, 1.0, 0.5 };
            Xic xic = makeXic(6, mzs1, rts1, sns1, abs1);
            xic.mergeDuplicates();
            shouldEqual(xic.size(), static_cast<size_t>(5));
            double emzs1[] =
                    { 100.001, 100.002, 100.005, 100.001, 100.003 };
            double erts1[] = { 10.0, 11.0, 13.0, 14.0, 15.0 };
            unsigned int esns1[] = { 1, 2, 4, 5, 6 };
            double eabs1[] = { 1.0, 1.0, 2.0, 1.0, 0.5 };

            for (size_t i = 0; i < 5; ++i) {
                shouldEqualTolerance(xic[i].getAbundance(), eabs1[i], 1e-12);
                shouldEqualTolerance(xic[i].getMz(), emzs1[i], 1e-12);
                shouldEqualTolerance(xic[i].getRetentionTime(), erts1[i], 1e-12);
                shouldEqual(xic[i].getScanNumber(), esns1[i]);
            }
        }

    }

    void testSplitRt1()
    {
        double mzs1[] =
                { 100.001, 100.003, 100.002, 100.005, 100.001, 100.003 };
        double rts1[] = { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0 };
        unsigned int sns1[] = { 1, 2, 3, 4, 5, 6 };
        double abs1[] = { 1.0, 2.0, 3.0, 2.0, 1.0, 0.5 };
        Xic xic = makeXic(6, mzs1, rts1, sns1, abs1);
        std::vector<Xic> tmp;
        xic.split(tmp);
        shouldEqual(tmp.size(), static_cast<size_t>(1)); // expect one XIC
        // FIXME: compare all elements
    }

    void testSplitRt2()
    {
        double mzs1[] = { 100.001, 100.003, 100.002, 100.005, 100.001,
                          100.003, 100.001 };
        double rts1[] = { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };
        unsigned int sns1[] = { 1, 2, 3, 4, 5, 6 };
        double abs1[] = { 1.0, 2.0, 3.0, 2.0, 1.0, 0.1, 1.0 };
        Xic xic = makeXic(6, mzs1, rts1, sns1, abs1);
        std::vector<Xic> tmp;
        xic.split(tmp);
        shouldEqual(tmp.size(), static_cast<size_t>(1)); // expect one XIC
        // FIXME: compare all elements
    }

    void testSplitRt3()
    {
        double mzs1[] = { 100.001, 100.003, 100.002, 100.005, 100.001,
                          100.003, 100.001, 100.004, 100.0, 100.01 };
        double rts1[] = { 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
                          18.0, 19.0 };
        unsigned int sns1[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
        double abs1[] = { 1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 0.5, 0.1 };
        Xic xic = makeXic(10, mzs1, rts1, sns1, abs1);
        std::vector<Xic> tmp;
        xic.split(tmp);
        shouldEqual(tmp.size(), static_cast<size_t>(2));
        shouldEqual(tmp[0].size(), static_cast<size_t>(5));
        shouldEqual(tmp[1].size(), static_cast<size_t>(5));
        for (size_t i = 0; i < 5; ++i) {
            shouldEqual(tmp[0][i].getAbundance(), abs1[i]);
        }
        for (size_t i = 0; i < 5; ++i) {
            shouldEqual(tmp[1][i].getAbundance(), abs1[i+5]);
        }

    }

    void testSplitRt4()
    {
        // some real world data that should not cause a split
        double rt[] =
                { 35.054338, 35.082955, 35.124493, 35.148760, 35.187350,
                  35.201423, 35.221842, 35.245268, 35.259378, 35.275875,
                  35.303758, 35.317743, 35.343372, 35.362485, 35.382055,
                  35.405510, 35.419758, 35.440430, 35.459777, 35.474085,
                  35.490142, 35.506402, 35.530067, 35.553845, 35.572822 };
        double ab[] = { 34898.0, 0.0, 40727.0, 59495.0, 135552.0, 225115.0,
                        333659.0, 469826.0, 468061.0, 565953.0, 855597.0,
                        1064007.0, 1252753.0, 1094078.0, 1286880.0, 1220093.0,
                        1003690.0, 968112.0, 589395.0, 704491.0, 480898.0,
                        485505.0, 196695.0, 112505.0, 68079.0 };
        std::vector<double> mz(500.0, 25);
        std::vector<unsigned int> sn(25, 1);
        std::partial_sum(sn.begin(), sn.end(), sn.begin());
        Xic xic = makeXic(25, &mz[0], rt, &sn[0], ab);
        std::vector<Xic> tmp;
        xic.split(tmp);
        //for (size_t i = 0; i < tmp.size(); ++i) {
        //    std::cerr << "--- " << i << " ---\n" << tmp[i];
        //}
        shouldEqual(tmp.size(), static_cast<size_t>(1));
    }

    void testSplitRt5()
    {
        // some real world data that should not cause a split
        double rt[] = { 2111.24, 2112.09, 2113.31, 2114.72, 2115.56, 2116.55,
                        2118.23, 2119.06, 2120.6, 2121.75, 2122.92, 2124.33,
                        2125.19, 2126.43, 2127.59, 2128.45, 2129.41, 2130.38,
                        2131.8, 2133.23, 2134.37 };
        unsigned int sn[] = { 4000, 4002, 4005, 4009, 4011, 4013, 4017, 4019,
                              4023, 4026, 4029, 4033, 4035, 4038, 4041, 4043,
                              4045, 4047, 4048, 4052, 4055 };
        double mz[] = { 548.814, 548.813, 548.814, 548.813, 548.813, 548.813,
                        548.813, 548.813, 548.813, 548.813, 548.813, 548.813,
                        548.813, 548.813, 548.813, 548.813, 548.813, 548.813,
                        548.813, 548.813, 548.814 };
        double ab[] = { 472009.0, 905473.0, 1291190.0, 1828580.0, 1817710.0,
                        2244620.0, 3388290.0, 4188030.0, 5001190.0, 4322550.0,
                        4969260.0, 4725260.0, 4004990.0, 3754450.0, 2327370.0,
                        2761370.0, 1909930.0, 1926350.0, 756566.0, 400389.0,
                        242239.0 };
        Xic xic = makeXic(21, mz, rt, sn, ab);
        std::vector<Xic> tmp;
        xic.split(tmp);
        //for (size_t i = 0; i < tmp.size(); ++i) {
        //    std::cerr << "--- " << i << " ---\n" << tmp[i];
        //}
        shouldEqual(tmp.size(), static_cast<size_t>(1));
    }
};

int main()
{
    XicTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

