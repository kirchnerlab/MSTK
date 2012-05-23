/*
 * Mercury7Impl-test.cpp
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
#include <MSTK/common/Types.hpp>
// expose the class
#define private public
#define protected public
#include <MSTK/ipaca/Mercury7Impl.hpp>
#undef private
#undef protected
#include <iostream>
#include "unittest.hxx"

using namespace mstk::ipaca;
using namespace mstk;

/** Tests for the Mercury7Impl algorithm.
 */
struct Mercury7TestSuite : vigra::test_suite
{
    /** Constructor.
     * The Mercury7TestSuite constructor adds all Mercury7Impl tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    Mercury7TestSuite() :
        vigra::test_suite("Mercury7")
    {
        add(testCase(&Mercury7TestSuite::testPrune));
        add(testCase(&Mercury7TestSuite::testConvolve));
        add(testCase(&Mercury7TestSuite::testOperator));
    }

    void testPrune()
    {

        // create a mercury functor
        detail::Mercury7Impl m;

        {
            // prune an empty spectrum
            detail::Spectrum empty;
            m.prune(empty, 1.0);
            shouldEqual(empty.size(), static_cast<Size>(0));
        }
        {
            // prune a spectrum with a single entry
            detail::Spectrum single;
            detail::SpectrumElement e;
            e.mz = 1.0;
            e.ab = 0.1;
            single.push_back(e); // ab == 0.1
            m.prune(single, 0.01);
            shouldEqual(single.size(), static_cast<Size>(1));
            shouldEqualTolerance(single[0].mz, 1.0, 1e-12);
            shouldEqualTolerance(single[0].ab, 0.1, 1e-12);
            single.clear();
            single.push_back(e);
            m.prune(single, 1.0);
            shouldEqual(single.size(), static_cast<Size>(0));
        }
        {
            // prune a spectrum w/ three entries

            // generate a spectrum [ 0.1 1.0 0.1 ]
            detail::Spectrum s;
            detail::SpectrumElement e;
            e.mz = 1.0;
            e.ab = 0.1;
            s.push_back(e);
            e.mz = 2.0;
            e.ab = 1.0;
            s.push_back(e);
            e.mz = 3.0;
            e.ab = 0.1;
            s.push_back(e);
            detail::Spectrum t(s);
            // keep everything
            m.prune(t, 0.01);
            shouldEqual(t.size(), static_cast<Size>(3));
            for (Size k = 0; k < 3; ++k) {
                shouldEqualTolerance(s[k].mz, t[k].mz, 1e-12);
                shouldEqualTolerance(s[k].ab, t[k].ab, 1e-12);
            }
            // prune everything except for the middle peak
            t = s;
            m.prune(t, 0.5);
            shouldEqual(t.size(), static_cast<Size>(1));
            shouldEqualTolerance(s[1].mz, t[0].mz, 1e-12);
            shouldEqualTolerance(s[1].ab, t[0].ab, 1e-12);
            // prune everything
            t = s;
            m.prune(t, 1.0); // also tests < vs <=
            shouldEqual(t.size(), static_cast<Size>(0));
        }
    }

    void testConvolve()
    {
        detail::Mercury7Impl m;
        {
            // convolve two empty spectra (make sure r is clean)
            detail::Spectrum s1, s2, r;
            r.resize(100);
            m.convolve(s1, s2, r);
            shouldEqual(r.size(), static_cast<Size>(0));
        }
        {
            // convolve one empty spectrum (left and right)
            detail::Spectrum s1, s2, r;
            detail::SpectrumElement e;
            e.mz = 1.0;
            e.ab = 0.1;
            s1.push_back(e);
            m.convolve(s1, s2, r);
            shouldEqual(r.size(), static_cast<Size>(1));
            shouldEqualTolerance(r[0].mz, 1.0, 1e-12);
            shouldEqualTolerance(r[0].ab, 0.0, 1e-12);
            r.clear();
            m.convolve(s2, s1, r);
            shouldEqual(r.size(), static_cast<Size>(1));
            shouldEqualTolerance(r[0].mz, 1.0, 1e-12);
            shouldEqualTolerance(r[0].ab, 0.0, 1e-12);
        }
        {
            // convolution results
            detail::Spectrum s1, s2, r;
            detail::SpectrumElement e;
            Double masses[] = { 1.0, 2.0 };
            Double freqs[] = { 0.5, 0.5 };
            for (size_t k = 0; k < 2; ++k) {
                e.mz = masses[k];
                e.ab = freqs[k];
                s1.push_back(e);
            }
            s2 = s1;
            m.convolve(s1, s2, r);
            shouldEqual(r.size(), static_cast<Size>(3));
            Double expectedMasses[] = { 2.0, 3.0, 4.0 };
            Double expectedAbundances[] = { 0.25, 0.5, 0.25 };
            for (Size k = 0; k < 3; ++k) {
                shouldEqualTolerance(r[k].mz, expectedMasses[k], 1e-12);
                shouldEqualTolerance(r[k].ab, expectedAbundances[k], 1e-12);
            }
        }
    }

    detail::Stoichiometry createIntegerH2O()
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

    void testOperator()
    {
        detail::Stoichiometry s = createIntegerH2O();
        detail::Mercury7Impl m;
        detail::Spectrum spectrum = m(s);
        Double expectedMasses[] = { 18.0, 19.0, 20.0, 21.0, 22.0 };
        shouldEqual(spectrum.size(), static_cast<Size>(5));
        for (size_t k = 0; k < 5; ++k) {
            shouldEqualTolerance(spectrum[k].mz, expectedMasses[k], 1e-12);
        }
        detail::Stoichiometry s2;
        // get the H
        s2.push_back(s.front());
        // edit its istope distribution
        s2[0].isotopes.pop_back();
        s2[0].isotopes[0].ab = 1.0;

        // single count
        s2[0].count = 1.0;
        spectrum = m(s2);
        shouldEqual(spectrum.size(), static_cast<Size>(1));
        shouldEqualTolerance(spectrum[0].mz, 1.0, 1e-12);
        shouldEqualTolerance(spectrum[0].ab, 1.0, 1e-12);

        // convolve the single isotope element 1000 times with itself
        s2[0].count = 1000.0;
        spectrum = m(s2);
        shouldEqual(spectrum.size(), static_cast<Size>(1));
        shouldEqualTolerance(spectrum[0].mz, 1000.0, 1e-12);
        shouldEqualTolerance(spectrum[0].ab, 1.0, 1e-12);

        // convolve two elements with equal isotope distributions over two masses
        {
            detail::Element e;
            detail::Isotope i;
            Double masses[] = { 1.0, 2.0 };
            Double freqs[] = { 0.5, 0.5 };
            for (size_t k = 0; k < 2; ++k) {
                i.mz = masses[k];
                i.ab = freqs[k];
                e.isotopes.push_back(i);
            }
            // two copies
            e.count = 2.0;
            detail::Stoichiometry s3;
            s3.push_back(e);
            spectrum = m(s3);
            shouldEqual(spectrum.size(), static_cast<Size>(3));
            Double expectedMasses[] = { 2.0, 3.0, 4.0 };
            Double expectedAbundances[] = { 0.25, 0.5, 0.25 };
            for (Size k = 0; k < 3; ++k) {
                shouldEqualTolerance(spectrum[k].mz, expectedMasses[k], 1e-12);
                shouldEqualTolerance(spectrum[k].ab, expectedAbundances[k], 1e-12);
            }
            // now do the same thing with two entries, each of which
            // has only a single copy
            e.count = 1.0;
            s3.clear();
            s3.push_back(e);
            s3.push_back(e);
            spectrum = m(s3);
            shouldEqual(spectrum.size(), static_cast<Size>(3));
            for (Size k = 0; k < 3; ++k) {
                shouldEqualTolerance(spectrum[k].mz, expectedMasses[k], 1e-12);
                shouldEqualTolerance(spectrum[k].ab, expectedAbundances[k], 1e-12);
            }
        }

        // same idea, with 4 copies
        {
            detail::Element e;
            detail::Isotope i;
            Double masses[] = { 1.0, 2.0 };
            Double freqs[] = { 0.5, 0.5 };
            for (size_t k = 0; k < 2; ++k) {
                i.mz = masses[k];
                i.ab = freqs[k];
                e.isotopes.push_back(i);
            }
            detail::Stoichiometry s;
            e.count = 4.0;
            s.push_back(e);
            spectrum = m(s);
            shouldEqual(spectrum.size(), static_cast<Size>(5));
            Double expectedMasses[] = { 4.0, 5.0, 6.0, 7.0, 8.0 };
            Double expectedAbundances[] =
                    { 0.0625, 0.25, 0.375, 0.25, 0.0625 };
            for (Size k = 0; k < 5; ++k) {
                shouldEqualTolerance(spectrum[k].mz, expectedMasses[k], 1e-12);
                shouldEqualTolerance(spectrum[k].ab, expectedAbundances[k], 1e-12);
            }
            // repeat, but with 4 distinct entries in the stoi vector
            // instead of a count of 4 (the difference is the use of msa
            // and esa registers in the mercury7 code).
            s.clear();
            e.count = 1.0;
            s.push_back(e);
            s.push_back(e);
            s.push_back(e);
            s.push_back(e);
            for (Size k = 0; k < 5; ++k) {
                shouldEqualTolerance(spectrum[k].mz, expectedMasses[k], 1e-12);
                shouldEqualTolerance(spectrum[k].ab, expectedAbundances[k], 1e-12);
            }
            // repeat again, using a mixture of entries.
            s.clear();
            e.count = 1.0;
            s.push_back(e);
            s.push_back(e);
            e.count = 2.0;
            s.push_back(e);
            for (Size k = 0; k < 5; ++k) {
                shouldEqualTolerance(spectrum[k].mz, expectedMasses[k], 1e-12);
                shouldEqualTolerance(spectrum[k].ab, expectedAbundances[k], 1e-12);
            }
        }
    }
};

/** The main function that runs the tests for class Mercury7Impl.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    Mercury7TestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

