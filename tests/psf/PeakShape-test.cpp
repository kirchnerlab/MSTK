/*
 * PeakShape-test.cpp
 *
 * Copyright (C) 2011 Bernhard Kausler
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
#include <cmath>
#include <iostream>
#include <MSTK/common/Error.hpp>
#include <MSTK/psf/Error.hpp>
#include <MSTK/psf/PeakShape.hpp>

#include "unittest.hxx"

using namespace mstk;

struct peakshapeTestSuite : vigra::test_suite {
    peakshapeTestSuite() : vigra::test_suite("peakshape") {
        add( testCase(&peakshapeTestSuite::testGaussianPeakShapeConstruction));
        add( testCase(&peakshapeTestSuite::testGaussianPeakShapeGetterSetter));
        add( testCase(&peakshapeTestSuite::testGaussianPeakShapeSigmaFwhmConversion));
        add( testCase(&peakshapeTestSuite::testGaussianPeakShapeAt));
        add( testCase(&peakshapeTestSuite::testGaussianPeakShapeGetSupportThreshold));
    }

    void testGaussianPeakShapeConstruction() {
        // default constructor
        psf::GaussianPeakShape gps;
        shouldEqual(gps.getSigma(), 0.1);

        // constructor init sigma
        psf::GaussianPeakShape gps_sigma(0.79);
        shouldEqual(gps_sigma.getSigma(), 0.79);

        // illegal sigma
        bool thrown = false;
        try {
            psf::GaussianPeakShape(0);
        } catch (const PreconditionViolation& e){
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        try {
            psf::GaussianPeakShape(-0.34);
        } catch (const PreconditionViolation& e){
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;
    }

    void testGaussianPeakShapeGetterSetter() {
        psf::GaussianPeakShape gps;
        bool thrown = false;

        // sigma
        gps.setSigma(0.5);
        shouldEqual(gps.getSigma(), 0.5);
                
        gps.setSigma(0.7);
        shouldEqual(gps.getSigma(), 0.7);

        try {
            gps.setSigma(0.);
        } catch (const PreconditionViolation& e){
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        try {
            gps.setSigma(-1.7);
        } catch (const PreconditionViolation& e){
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        // fwhm
        gps.setFwhm(0.5);
        shouldEqual(gps.getFwhm(), 0.5);

        gps.setFwhm(0.7);
        shouldEqual(gps.getFwhm(), 0.7);

        try {
            gps.setFwhm(0.);
        } catch (const PreconditionViolation& e){
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        try {
            gps.setFwhm(-1.7);
        } catch (const PreconditionViolation& e){
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        // sigmaFactorForSupportThreshold
        gps.setSigmaFactorForSupportThreshold(0.5);
        shouldEqual(gps.getSigmaFactorForSupportThreshold(), 0.5);

        gps.setSigmaFactorForSupportThreshold(0.7);
        shouldEqual(gps.getSigmaFactorForSupportThreshold(), 0.7);

        try {
            gps.setSigmaFactorForSupportThreshold(0.);
        } catch (const PreconditionViolation& e){
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        try {
            gps.setSigmaFactorForSupportThreshold(-1.7);
        } catch (const PreconditionViolation& e){
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;
    }

    void testGaussianPeakShapeSigmaFwhmConversion() {
        psf::GaussianPeakShape gps;

        // conversion factor
        double conversionFactor = gps.sigmaToFwhmConversionFactor();
        shouldEqualTolerance(conversionFactor, 2.35482, 0.00001);

        // sigma / fwhm conversion
        gps.setSigma(0.5);
        shouldEqual(gps.getFwhm(), conversionFactor * 0.5);

        gps.setFwhm(0.5);
        shouldEqual(gps.getSigma(), 0.5 / conversionFactor);
    }

    void testGaussianPeakShapeAt() {
        struct Gauss {
            static double at(const double x, const double sigma) {
                return std::exp(-(x * x) / (2 * sigma * sigma));
            }
        };

        psf::GaussianPeakShape gps;

        // has to be 1 at xCoordinate=0 
        shouldEqual(gps.at(0), 1);

        // test some values
        gps.setSigma(0.5);
        shouldEqual(gps.at(0.1), Gauss::at(0.1, 0.5));
        shouldEqual(gps.at(3.5), Gauss::at(3.5, 0.5));
        shouldEqual(gps.at(-0.34), Gauss::at(-0.34, 0.5));
        shouldEqual(gps.at(-2.73), Gauss::at(-2.73, 0.5));

        gps.setSigma(0.9);
        shouldEqual(gps.at(0.1), Gauss::at(0.1, 0.9));
        shouldEqual(gps.at(3.5), Gauss::at(3.5, 0.9));
        shouldEqual(gps.at(-0.34), Gauss::at(-0.34, 0.9));
        shouldEqual(gps.at(-2.73), Gauss::at(-2.73, 0.9));
    }

    void testGaussianPeakShapeGetSupportThreshold() {
        psf::GaussianPeakShape gps;

        shouldEqual(gps.getSigmaFactorForSupportThreshold(), 3.0);

        gps.setSigma(1.5);
        shouldEqual(gps.getSupportThreshold(), 4.5);

        gps.setSigma(0.7);
        shouldEqual(gps.getSupportThreshold(), 0.7 * 3.0);
    }
};

int main()
{
    peakshapeTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


