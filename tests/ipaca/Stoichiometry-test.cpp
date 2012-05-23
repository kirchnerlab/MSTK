/*
 * Stoichiometry-test.cpp
 *
 * Copyright (C) 2012 Marc Kirchner
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
#include <MSTK/ipaca/Stoichiometry.hpp>
#include <iostream>
#include "unittest.hxx"

using namespace mstk::ipaca;

/** Tests for all methods and free functions in Stoichiometry.cpp.
 */
struct StoichiometryTestSuite : vigra::test_suite
{
    /** Constructor.
     * The StoichiometryTestSuite constructor adds all Stoichiometry tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    StoichiometryTestSuite() :
        vigra::test_suite("Stoichiometry")
    {
        add(testCase(&StoichiometryTestSuite::testIsPlausibleStoichiometry));
        add(testCase(&StoichiometryTestSuite::testSplitStoichiometry));
    }

    detail::Stoichiometry createH2O()
    {
        detail::Stoichiometry h2o;
        detail::Isotope i;
        detail::Element h;
        double massesH[] = { 1.0, 2.0 };
        double freqsH[] = { 0.99, 0.01 };
        for (size_t k = 0; k < 2; ++k) {
            i.mz = massesH[k];
            i.ab = freqsH[k];
            h.isotopes.push_back(i);
        }
        h.count = 2.0;
        detail::Element o;
        double massesO[] = { 16.0, 17.0, 18.0 };
        double freqsO[] = { 0.97, 0.01, 0.02 };
        for (size_t k = 0; k < 3; ++k) {
            i.mz = massesO[k];
            i.ab = freqsO[k];
            o.isotopes.push_back(i);
        }
        o.count = 1.0;
        h2o.push_back(h);
        h2o.push_back(o);
        return h2o;
    }

    void testIsPlausibleStoichiometry()
    {
        // normal H2O
        detail::Stoichiometry h2o = createH2O();
        shouldEqual(detail::isPlausibleStoichiometry(h2o), true);
        // at least one non-zero entry
        h2o[0].count = 0.0;
        shouldEqual(detail::isPlausibleStoichiometry(h2o), true);
        // one negative entry
        h2o[0].count = -1.0;
        shouldEqual(detail::isPlausibleStoichiometry(h2o), false);
        // all zero entries
        h2o[0].count = 0.0;
        h2o[1].count = 0.0;
        shouldEqual(detail::isPlausibleStoichiometry(h2o), false);
    }

    void testSplitStoichiometry()
    {
        detail::Stoichiometry s = createH2O();
        // make it 2.4 H and 1.3 O atoms
        s[0].count = 2.4;
        s[1].count = 1.3;
        detail::Stoichiometry i, f;
        detail::splitStoichiometry(s, i, f);
        shouldEqual(i.size(), static_cast<size_t>(2));
        shouldEqual(f.size(), static_cast<size_t>(2));
        shouldEqual(i[0].count, 2.0);
        shouldEqualTolerance(f[0].count, 0.4, 1e-12);
        shouldEqual(i[1].count, 1.0);
        shouldEqualTolerance(f[1].count, 0.3, 1e-12);
        // make sure that all isotope info is the same
        for (size_t k = 0; k < 2; ++k) {
            shouldEqual(s[0].isotopes[k].mz, i[0].isotopes[k].mz);
            shouldEqual(s[0].isotopes[k].ab, i[0].isotopes[k].ab);
            shouldEqual(s[0].isotopes[k].mz, f[0].isotopes[k].mz);
            shouldEqual(s[0].isotopes[k].ab, f[0].isotopes[k].ab);
        }
        for (size_t k = 0; k < 3; ++k) {
            shouldEqual(s[1].isotopes[k].mz, i[1].isotopes[k].mz);
            shouldEqual(s[1].isotopes[k].ab, i[1].isotopes[k].ab);
            shouldEqual(s[1].isotopes[k].mz, f[1].isotopes[k].mz);
            shouldEqual(s[1].isotopes[k].ab, f[1].isotopes[k].ab);
        }
    }
};

/** The main function that runs the tests for class Stoichiometry.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    StoichiometryTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

