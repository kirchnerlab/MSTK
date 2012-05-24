
/*
 * PeakShapeFunction-test.cpp
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

#include "unittest.hxx"
#include <MSTK/common/Log.hpp>
#include <MSTK/psf/PeakShapeFunction.hpp>
#include <MSTK/psf/types/Spectrum.hpp>
#include "testdata.hpp"

using namespace mstk;

struct PsfTestSuite : vigra::test_suite {
    PsfTestSuite() : vigra::test_suite("PeakShapeFunction") {
        add( testCase(&PsfTestSuite::testPsfTypeEnum) );
        add( testCase(&PsfTestSuite::testPeakShapeFunctionType) );
        add( testCase(&PsfTestSuite::testGetType) );
        add( testCase(&PsfTestSuite::testOrbitrapPeakShapeFunction) );
        add( testCase(&PsfTestSuite::testGaussianPeakShapeFunction) );
        add( testCase(&PsfTestSuite::testOperator));
        add( testCase(&PsfTestSuite::testGetSupportThreshold));
        add( testCase(&PsfTestSuite::testSet_GetMinimalPeakHeightForCalibration));
        add( testCase(&PsfTestSuite::testOrbiFwhmLinearSqrtPeakShape));
    }


    void testPsfTypeEnum() {
        // there should be six different enum values: box, gaussian, orbi, orbiBox, orbiConst, tof
        MSTK_LOG(logINFO) << "Testing the PeakShapeFunctionType enum.";
        psf::PeakShapeFunctionTypes t;
        t = psf::box;
        t = psf::gaussian;
        t = psf::orbi;
        t = psf::orbiBox;
        t = psf::tof;
    }


    void testPeakShapeFunctionType() {
        MSTK_LOG(logINFO) << "Testing class PeakShapeFunctionType";

        psf::PeakShapeFunctionType psfBox = psf::box;
        shouldEqual(psfBox.toEnum(), psf::box);
        should(psfBox.toString() == "box");

        psf::PeakShapeFunctionType psfGaussian = psf::gaussian;
        shouldEqual(psfGaussian.toEnum(), psf::gaussian);
        should(psfGaussian.toString() == "gaussian");

        psf::PeakShapeFunctionType psfOrbi = psf::orbi;
        shouldEqual(psfOrbi.toEnum(), psf::orbi);
        should(psfOrbi.toString() == "orbi");

        psf::PeakShapeFunctionType psfOrbiBox = psf::orbiBox;
        shouldEqual(psfOrbiBox.toEnum(), psf::orbiBox);
        should(psfOrbiBox.toString() == "orbiBox");

        psf::PeakShapeFunctionType psfTof = psf::tof;
        shouldEqual(psfTof.toEnum(), psf::tof);
        should(psfTof.toString() == "time-of-flight");

        // test illegal enum (choose the integer high enough...)
        psf::PeakShapeFunctionType psfIllegal = static_cast<psf::PeakShapeFunctionTypes>(200);
        shouldEqual(psfIllegal.toEnum(), static_cast<psf::PeakShapeFunctionTypes>(200));
        should(psfIllegal.toString() == "unknown");

    }
    
    // We want to test this function for every implementation of the abstract
    // 'PeakShapeFunction' interface.
    // The 'box' type is only used for unit testing and not included in the core
    // library. So, we don't test it.
    void testGetType() {
        // Testing getType(). We use sound values for the constructors.
        MSTK_LOG(logINFO) << "Testing the getType() functions.";
        shouldEqual(psf::GaussianPeakShapeFunction().getType().toEnum(), psf::gaussian);
        shouldEqual(psf::OrbitrapPeakShapeFunction().getType().toEnum(), psf::orbi);
    }



    void testOrbitrapPeakShapeFunction() {
        psf::OrbitrapPeakShapeFunction orbi_psf;
        shouldEqual(orbi_psf.getType().toEnum(), psf::orbi);

        // construction
        psf::OrbitrapPeakShapeFunction psf1(0.1214);
        shouldEqual(psf1.getA(), 0.1214);

        // operator()
        orbi_psf.setA(0.0123);
        shouldEqualTolerance(orbi_psf(400., 402.), 0.998856, 0.000001); 
        shouldEqualTolerance(orbi_psf(400., 397.64), 0.998407, 0.000001);
        shouldEqualTolerance(orbi_psf(600., 602.), 0.999661, 0.000001);
        shouldEqualTolerance(orbi_psf(600., 597.64), 0.999528, 0.000001);

        // test half maximum at full width
        orbi_psf.setA(1.); // corresponds to const. fwhm of 2.0
        double fullMaximum = orbi_psf(1., 1.);
        double halfMaximum = orbi_psf(1., 1.5);
        shouldEqual(fullMaximum, 2. * halfMaximum);

        // auto calibration
        MSTK_LOG(logINFO) << "Calibrating orbitrap peak shape function.";
	psf::Spectrum spectrum;
	loadSpectrumElements(spectrum,dirTestdata + "/PeakShapeFunctions/realistic_ms1.wsv");
	psf::MzExtractor get_mz;
	psf::IntensityExtractor get_int;
        
        orbi_psf.setA(0); // reset 
        orbi_psf.calibrateFor(get_mz, get_int, spectrum.begin(), spectrum.end());
        shouldEqualTolerance(orbi_psf.getA(), 1.19781e-06, 0.00001);       
    }

    void testGaussianPeakShapeFunction() {
        psf::GaussianPeakShapeFunction gaussian_psf;
        shouldEqual(gaussian_psf.getType().toEnum(), psf::gaussian);

        // construction
        psf::GaussianPeakShapeFunction psf(0.11442);
        shouldEqual(psf.getA(), 0.11442);

        // gaussian psf should be independent of the mass channel observed
        gaussian_psf.setA(3.);
        shouldEqualTolerance(gaussian_psf(400.,402.), 0.291632, 0.000001);
        shouldEqualTolerance(gaussian_psf(400., 397.64), 0.17982, 0.00001);
        shouldEqualTolerance(gaussian_psf(600.,602.), 0.291632, 0.000001);
        shouldEqualTolerance(gaussian_psf(600., 597.64), 0.17982, 0.00001);

        // test half maximum at full width
        gaussian_psf.setA(2.); // corresponds to const. fwhm of 2.0
        double fullMaximum = gaussian_psf(400., 400.);
        double halfMaximum = gaussian_psf(400., 401.);
        shouldEqual(fullMaximum, 2. * halfMaximum);

        // auto calibration
        MSTK_LOG(logINFO) << "Calibrating gaussian peak shape function.";
	psf::Spectrum spectrum;
	loadSpectrumElements(spectrum,dirTestdata + "/PeakShapeFunctions/realistic_ms1.wsv");
	psf::MzExtractor get_mz;
	psf::IntensityExtractor get_int;
        
        gaussian_psf.setA(0); // reset 
        gaussian_psf.calibrateFor(get_mz, get_int, spectrum.begin(), spectrum.end());
        shouldEqualTolerance(gaussian_psf.getA(), 0.031325, 0.000001);
    }

    void testOperator() {
        psf::PeakShapeFunctionTemplate<psf::GaussianPeakShape, psf::TofFwhm, psf::tof> gen;
        psf::TofFwhm fwhm;
        psf::GaussianPeakShape ps;

        gen.setA(0.43);
        gen.setB(0.76);
        fwhm.setA(0.43);
        fwhm.setB(0.76);

        ps.setFwhm(fwhm.at(400.));

        // ensure to be inside the threshold
        should(gen.getSupportThreshold(400.) > 5.0);

        shouldEqual(gen(400., 404.5), ps.at(404.5 - 400.));
        shouldEqual(gen(400., 397.2), ps.at(400. - 397.2));
        shouldEqual(gen(400., 400.), ps.at(0.));

        // now, test the threshold behaviour

        double threshold = gen.getSupportThreshold(400.);
        // Should be not zero at the threshold
        should(gen(400., 400. + threshold) > 0.);
        should(gen(400., 400. - threshold) > 0.);

        // zero after the threshold
        double delta = std::numeric_limits<double>().epsilon() + std::numeric_limits<double>().round_error();
        shouldEqual(gen(400., 400. + threshold + delta), 0.0);
        shouldEqual(gen(400., 400. - (threshold + delta)), 0.0);       
    }

    void testGetSupportThreshold() {
        psf::PeakShapeFunctionTemplate<psf::GaussianPeakShape, psf::TofFwhm, psf::tof> gen;
        psf::TofFwhm fwhm;
        psf::GaussianPeakShape ps;

        gen.setA(0.43);
        gen.setB(0.76);
        fwhm.setA(0.43);
        fwhm.setB(0.76);

        ps.setFwhm(fwhm.at(400.));
        
        double threshold = ps.getSupportThreshold();
        shouldEqual(gen.getSupportThreshold(400.), threshold);
    }

    void testSet_GetMinimalPeakHeightForCalibration() {
        psf::PeakShapeFunctionTemplate<psf::GaussianPeakShape, psf::OrbitrapFwhm, psf::orbi> psf;
       
        psf.setMinimalPeakHeightForCalibration(4.2);
        shouldEqual(psf.getMinimalPeakHeightForCalibration(), 4.2);

        psf.setMinimalPeakHeightForCalibration(0);
        shouldEqual(psf.getMinimalPeakHeightForCalibration(), 0);

        psf.setMinimalPeakHeightForCalibration(-0.87);
        shouldEqual(psf.getMinimalPeakHeightForCalibration(), -0.87);
    }

    void testOrbiFwhmLinearSqrtPeakShape() {
        // Note: This was the test for the first version of the orbi psf, but doesn't apply anymore.
        // Nevertheless, it is still a valid test for its combination of template parameters. Therefore,
        // we 'retired' the test here.

        psf::PeakShapeFunctionTemplate<psf::GaussianPeakShape, psf::OrbitrapFwhm, psf::orbi>  orbi_psf;
        shouldEqual(orbi_psf.getType().toEnum(), psf::orbi);

        // construction
        psf::PeakShapeFunctionTemplate<psf::GaussianPeakShape, psf::OrbitrapFwhm, psf::orbi> psf1(0.1214);
        shouldEqual(psf1.getA(), 0.1214);
        shouldEqual(psf1.getB(), orbi_psf.getB());

        psf::PeakShapeFunctionTemplate<psf::GaussianPeakShape, psf::OrbitrapFwhm, psf::orbi> psf2(0.1214, 1.20392);
        shouldEqual(psf2.getA(), 0.1214);
        shouldEqual(psf2.getB(), 1.20392);

        // operator()
        orbi_psf.setA(0.0123);
        orbi_psf.setB(0.0234);
        shouldEqualTolerance(orbi_psf(400., 402.), 0.998856, 0.000001); 
        shouldEqualTolerance(orbi_psf(400., 397.64), 0.998407, 0.000001);
        shouldEqualTolerance(orbi_psf(600., 602.), 0.999661, 0.000001);
        shouldEqualTolerance(orbi_psf(600., 597.64), 0.999528, 0.000001);

        // test half maximum at full width
        orbi_psf.setA(1.); // corresponds to const. fwhm of 2.0
        orbi_psf.setB(0.);
        double fullMaximum = orbi_psf(1., 1.);
        double halfMaximum = orbi_psf(1., 1.5);
        shouldEqual(fullMaximum, 2. * halfMaximum); 
    }

};


int main()
{
    PsfTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


