/*
 * IsotopePattern-test.cpp
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
#include <MSTK/fe/types/IsotopePattern.hpp>
#include <MSTK/fe/types/Xic.hpp>
#include <MSTK/fe/types/Spectrum.hpp>
#include <MSTK/common/Error.hpp>
#include <iostream>
#include <boost/assign.hpp>

#include "unittest.hxx"
#include "utilities.hpp"

using namespace mstk;
using namespace mstk::fe;
using namespace boost::assign;

struct IsotopePatternTestSuite : vigra::test_suite
{
    IsotopePatternTestSuite() :
        vigra::test_suite("IsotopePattern")
    {
        add(testCase(&IsotopePatternTestSuite::testCharges));
        add(testCase(&IsotopePatternTestSuite::testAbundance));
        add(testCase(&IsotopePatternTestSuite::testAsSpectrum));
        add(testCase(&IsotopePatternTestSuite::testSplit));
    }

    void test()
    {
        failTest("Unit test not yet implemented for class IsotopePattern");
    }

    void testCharges()
    {
        IsotopePattern ip;
        // setget
        ip.setCharges(ip.getCharges());
        // setter and getter
        std::set<int> r = ip.getCharges();
        shouldEqual(r.empty(), true);
        std::set<int> c = list_of(1)(2)(3)(4).to_container(c);
        ip.setCharges(c);
        r = ip.getCharges();
        shouldEqual(r.size(), static_cast<size_t>(4));
        for (size_t k = 1; k <= r.size(); ++k) {
            shouldEqual(r.count(k), static_cast<size_t>(1));
        }
    }

    void testAbundance()
    {
        double mz[] = { 99.999, 100.0, 100.001 };
        double rt[] = { 350.0, 352.0, 354.0 };
        unsigned int sn[] = { 42, 43, 44 };
        double ab[] = { 0.5, 1.0, 0.5 };
        Xic x1 = makeXic(3, mz, rt, sn, ab);
        // std::cerr << x1.getAbundance() << std::endl;
        IsotopePattern ip;
        shouldEqual(ip.getAbundance(), 0.0);
        ip.push_back(x1);
        shouldEqual(ip.getAbundance(), x1.getAbundance());
        double mz2[] = { 98.999, 99.0, 99.001 };
        Xic x2 = makeXic(3, mz2, rt, sn, ab);
        ip.push_back(x2);
        shouldEqual(ip.getAbundance(), x1.getAbundance() + x2.getAbundance());

    }

    void testAsSpectrum()
    {
        {
            // empty
            IsotopePattern ip;
            Spectrum ss;
            ip.asSpectrum(ss);
            shouldEqual(ss.empty(), true);
        }
        {
            // single entry
            double mz[] = { 99.999, 100.0, 100.001 };
            double rt[] = { 350.0, 352.0, 354.0 };
            unsigned int sn[] = { 42, 43, 44 };
            double ab[] = { 0.5, 1.0, 0.5 };
            Xic x1 = makeXic(3, mz, rt, sn, ab);
            IsotopePattern ip;
            ip.push_back(x1);
            Spectrum ss;
            ip.asSpectrum(ss);
            shouldEqual(ss.size(), static_cast<size_t>(1));
            shouldEqual(ss[0].mz, x1.getMz());
            shouldEqual(ss[0].abundance, x1.getAbundance());
        }
        {
            // multiple entries
            double mz[] = { 99.999, 100.0, 100.001 };
            double rt[] = { 350.0, 352.0, 354.0 };
            unsigned int sn[] = { 42, 43, 44 };
            double ab[] = { 0.5, 1.0, 0.5 };
            Xic x1 = makeXic(3, mz, rt, sn, ab);
            double mz2[] = { 100.999, 101.0, 101.001 };
            Xic x2 = makeXic(3, mz2, rt, sn, ab);
            double mz3[] = { 101.999, 102.0, 102.001 };
            Xic x3 = makeXic(3, mz3, rt, sn, ab);
            double mz4[] = { 102.999, 103.0, 103.001 };
            Xic x4 = makeXic(3, mz4, rt, sn, ab);
            IsotopePattern ip;
            ip.push_back(x1);
            ip.push_back(x2);
            ip.push_back(x3);
            ip.push_back(x4);
            Spectrum ss;
            ip.asSpectrum(ss);
            shouldEqual(ss.size(), static_cast<size_t>(4));
            shouldEqual(ss[0].mz, x1.getMz());
            shouldEqual(ss[0].abundance, x1.getAbundance());
            shouldEqual(ss[1].mz, x2.getMz());
            shouldEqual(ss[1].abundance, x2.getAbundance());
            shouldEqual(ss[2].mz, x3.getMz());
            shouldEqual(ss[2].abundance, x3.getAbundance());
            shouldEqual(ss[3].mz, x4.getMz());
            shouldEqual(ss[3].abundance, x4.getAbundance());
        }
    }

    void testSplit()
    {
        // this is supposed to fail unless the split function
        // gets implemented.
        bool thrown = false;
        try {
            IsotopePattern ip;
            std::vector<IsotopePattern> ips;
            ip.split(ips);
        } catch(RuntimeError& e) {
            thrown = true;
        }
        shouldEqual(thrown, true);
 /*
        {
            // generate a few patterns that cannot be split
            IsotopePattern ip;
            std::vector<IsotopePattern> ips;
            ip.split(ips);
            shouldEqual(ips.size(), static_cast<size_t>(1));
            shouldEqual(ips[0], ip);
            std::vector<double> mz;
            mz += 100.0;
            ip = makeIsotopePattern(mz, 1300.4, 2);
            ips.clear();
            ip.split(ips);
            shouldEqual(ips.size(), static_cast<size_t>(1));
            shouldEqual(ips[0], ip);
        }
        {

            std::vector<double> mz;
            // 13C distances
            mz += 100.0, 101.0033, 102.0066;
            IsotopePattern ip = makeIsotopePattern(mz, 1300.4, 2);
        }
*/
    }
};

int main()
{
    IsotopePatternTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

