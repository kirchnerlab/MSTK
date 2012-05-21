/*
 * XicLocalMinSplitter-test.cpp
 *
 * Copyright (c) 2012 Marc Kirchner
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
#include "utilities.hpp"
#include <MSTK/fe/XicLocalMinSplitter.hpp>
#include <MSTK/fe/RunningMeanSmoother.hpp>
#include <MSTK/common/Types.hpp>
#include <MSTK/common/Log.hpp>
#include <MSTK/fe/types/Xic.hpp>
#include <iostream>
#include <vector>

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

struct XicLocalMinSplitterTestSuite : vigra::test_suite
{
    XicLocalMinSplitterTestSuite() :
        vigra::test_suite("XicLocalMinSplitter")
    {
        add(testCase(&XicLocalMinSplitterTestSuite::testSplitRt1));
        add(testCase(&XicLocalMinSplitterTestSuite::testSplitRt2));
        add(testCase(&XicLocalMinSplitterTestSuite::testSplitRt3));
        add(testCase(&XicLocalMinSplitterTestSuite::testSplitRt4));
        add(testCase(&XicLocalMinSplitterTestSuite::testSplitRt5));
    }

    void split(const Xic& xic, std::vector<Xic>& xics) {
        // get a copy
        Xic smoothXic(xic);
        // smooth the copy
        Smoother smoother;
        //std::cerr << "Xic size (pre-smoothing): " << xic.size() << std::endl;
        smoother.run(smoothXic.begin(), smoothXic.end());
        //std::cerr << "Xic size (post-smoothing): " << xic.size() << std::endl;
        // now split
        XicLocalMinSplitter<Xic> splitter;
        splitter.split(xic.begin(), xic.end(), smoothXic.begin(), smoothXic.end());
        // and store
        xics.clear();
        typedef XicLocalMinSplitter<Xic>::const_iterator IT;
        for (IT i = splitter.begin(); i != splitter.end(); ++i) {
            //std::cerr << "subXic size: " << std::distance(i->first, i->second) << std::endl;
            xics.push_back(Xic(i->first, i->second));
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
            split(xic, tmp);
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
            split(xic, tmp);
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
            split(xic, tmp);
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
            split(xic, tmp);
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
            split(xic, tmp);
            //for (size_t i = 0; i < tmp.size(); ++i) {
            //    std::cerr << "--- " << i << " ---\n" << tmp[i];
            //}
            shouldEqual(tmp.size(), static_cast<size_t>(1));
        }
};

int main()
{
    XicLocalMinSplitterTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

