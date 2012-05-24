/*
 * SpectrumAlgorithm-test.cpp
 *
 * Copyright (C) 2011 Bernharc Kausler
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
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

#include <MSTK/common/Log.hpp>
#include <MSTK/common/Error.hpp>
#include <MSTK/psf/Error.hpp>
#include <MSTK/psf/SpectrumAlgorithm.hpp>
#include <MSTK/psf/types/Spectrum.hpp>

#include "testdata.hpp"

#include "unittest.hxx"

using namespace mstk::psf;
using namespace mstk;

namespace std {
ostream& operator<<(ostream& os, const pair<int*, int*>& p) {
    os << "<";
    if(p.first!=0)
        os << p.first;
    else
        os << "NULL";
    os << ", ";
    if (p.second!=0)
        os << p.second;
    else
        os << "NULL";
    os << ">";
    return os;
}

ostream& operator<<(ostream& os, const pair<double, double>& p) {
    os << "<" << p.first << ", " << p.second << ">";
    return os;
}
} /*namespace std*/

struct SpectrumAlgorithmTestSuite : vigra::test_suite {
    SpectrumAlgorithmTestSuite() : vigra::test_suite("SpectrumAlgorithm") {
        add( testCase(&SpectrumAlgorithmTestSuite::testFindBump));
        add( testCase(&SpectrumAlgorithmTestSuite::testMeasureFullWidths));
    }

    void testFindBump() {
    /* test data */

    // A sequence with two bumps. The maxima are at 69 and 21. A bump with a plateu is at 19 (shouldn't be found).
    // Note: The length of the array is 23.
    int twoBumps[] = {100, 81, 56, 56, 57, 69, 40, 13, 13, 9, 18, 21, 19, 15, 16, 19, 19, 18, 12, 11, 17, 22, 47};
    
    std::pair<int*, int*> firstBump(twoBumps + 3, twoBumps + 7); 
    shouldEqual(*(firstBump.first), 56);
    shouldEqual(*(firstBump.second), 13);
    
    std::pair<int*, int*> secondBump(twoBumps + 9, twoBumps + 13); 
    shouldEqual(*(secondBump.first), 9);
    shouldEqual(*(secondBump.second), 15);

    std::pair<int*, int*> twoBumps_noneFound(twoBumps + 23, twoBumps + 23);
    shouldEqual(*(twoBumps_noneFound.first-1), 47);
    shouldEqual(*(twoBumps_noneFound.second-1), 47);

    // no bumps; only a 'valley'
    // Note: The length of the array is 9.
    int noBumps[] = {99, 67, 98, 98, 98, 110, 112, 117, 121};

    std::pair<int*, int*> noBumps_noneFound(noBumps + 9, noBumps + 9);
    shouldEqual(*(noBumps_noneFound.first - 1), 121);
    shouldEqual(*(noBumps_noneFound.second - 1), 121);  
   

    std::pair<int*, int*> returned;

    // find both bumps
    returned = psf::findBump(twoBumps, twoBumps+23, std::less<int>());   
    shouldEqual(returned, firstBump);
    returned = psf::findBump(returned.second, twoBumps+23, std::less<int>());
    shouldEqual(returned, secondBump);

    // no bump is in this sequence
    returned = psf::findBump(noBumps, noBumps + 9, std::less<int>());
    shouldEqual(returned, noBumps_noneFound);
    }

    void testMeasureFullWidths() {
        using namespace psf;
        using namespace std;
	MzExtractor get_mz;
	IntensityExtractor get_int;

        // A spectrum with three pure peaks and some noise.
        // Peak# | lowness | full width 0.7 | full width 0.5 | full width 0.1
        // 1     | 30%     | 2              | nan            | nan
        // 2     | 50%     | 2              | 4              | nan
        // 3     | 90%     | 2              | 2              | 4
        Spectrum s1;
        s1.push_back(SpectrumElement(1, 9));
        s1.push_back(SpectrumElement(2, 8));
        s1.push_back(SpectrumElement(2.9, 6.8));//Peak 1 left
        s1.push_back(SpectrumElement(3, 7));  
        s1.push_back(SpectrumElement(4, 10)); // Peak 1 max
        s1.push_back(SpectrumElement(5, 7));
        s1.push_back(SpectrumElement(5.1, 6.8));  // Peak 1 right
        s1.push_back(SpectrumElement(6, 4.9));
        s1.push_back(SpectrumElement(6.9, 4.9));  // Peak 2 left  
        s1.push_back(SpectrumElement(7, 5));
        s1.push_back(SpectrumElement(8, 7));  // ...
        s1.push_back(SpectrumElement(9, 10)); // Peak 2 max
        s1.push_back(SpectrumElement(10, 7)); // ...
        s1.push_back(SpectrumElement(11, 5)); // ...
        s1.push_back(SpectrumElement(12, 1));
        s1.push_back(SpectrumElement(12.1, 0.9)); // Peak 2 right, Peak 3 left
        s1.push_back(SpectrumElement(12.2, 1));
        s1.push_back(SpectrumElement(13, 5)); // ...
        s1.push_back(SpectrumElement(12.5, 7));//...
        s1.push_back(SpectrumElement(14, 10)); // Peak 3
        s1.push_back(SpectrumElement(14.5, 7));//...
        s1.push_back(SpectrumElement(15, 5));  //...
        s1.push_back(SpectrumElement(16, 1));
        s1.push_back(SpectrumElement(16.1, 0.9)); // Peak 3 right
   
        typedef vector<pair< MzExtractor::result_type,
                             MzExtractor::result_type
                           > 
                      > MzWidthPair;
        MzWidthPair result;
    
        MSTK_LOG(logINFO) << "Testing fraction of 0.7 .";
        result = measureFullWidths(get_mz, get_int, s1.begin(), s1.end(), 0.7);
        shouldEqual(result.size(), (MzWidthPair::size_type)3);
        should(result.at(0).first == 4);
        shouldEqualTolerance(result.at(0).second, 2., 0.1);
        should(result.at(1).first == 9);
        shouldEqualTolerance(result.at(1).second, 2., 0.1);
        should(result.at(2).first == 14);
        shouldEqualTolerance(result.at(2).second, 2., 0.1);

        MSTK_LOG(logINFO) << "Testing fraction of 0.51 .";
        result = measureFullWidths(get_mz, get_int, s1.begin(), s1.end(), 0.51);
        shouldEqual(result.size(), (MzWidthPair::size_type)2);
        should(result.at(0).first == 9);
        shouldEqualTolerance(result.at(0).second, 4., 0.1);
        should(result.at(1).first == 14);
        shouldEqualTolerance(result.at(1).second, 2., 0.1);

        MSTK_LOG(logINFO) << "Testing fraction of 0.11 .";
        result = measureFullWidths(get_mz, get_int, s1.begin(), s1.end(), 0.11);
        shouldEqual(result.size(), (MzWidthPair::size_type)1);
        shouldEqual(result.at(0).first, 14);
        shouldEqualTolerance(result.at(0).second, 4., 0.1);

        MSTK_LOG(logINFO) << "Testing fraction of 0.001 .";
        result = measureFullWidths(get_mz, get_int, s1.begin(), s1.end(), 0.001);
        should(result.empty());

        // testing minimum intensity
        result = measureFullWidths(get_mz, get_int, s1.begin(), s1.end(), 0.7, 0);
        shouldEqual(result.size(), (MzWidthPair::size_type)3);
        result = measureFullWidths(get_mz, get_int, s1.begin(), s1.end(), 0.7, 11);
        shouldEqual(result.size(), (MzWidthPair::size_type)0);

        // test spectrum with no pure peaks
        MSTK_LOG(logINFO) << "Testing spectrum with no pure peaks.";
        Spectrum s_unpure;
        s_unpure.push_back(SpectrumElement(1, 9));
        s_unpure.push_back(SpectrumElement(2, 8));
        s_unpure.push_back(SpectrumElement(4, 8));
        s_unpure.push_back(SpectrumElement(5, 7));
        result = measureFullWidths(get_mz, get_int, s_unpure.begin(), s_unpure.end(), 0.5);
        should(result.empty());

        // test spectrum with duplicate mz values
        MSTK_LOG(logINFO) << "Testing spectrum with duplicate mz values.";
        Spectrum s_duplicate;
        s_duplicate.push_back(SpectrumElement(1, 9));
        s_duplicate.push_back(SpectrumElement(1, 9));
        s_duplicate.push_back(SpectrumElement(1, 9));
        result = measureFullWidths(get_mz, get_int, s_duplicate.begin(), s_duplicate.end(), 0.5);
        should(result.empty());

        s_duplicate.push_back(SpectrumElement(6.9, 4.9));  // Peak left  
        s_duplicate.push_back(SpectrumElement(7, 5));
        s_duplicate.push_back(SpectrumElement(8, 7));  // ...
        s_duplicate.push_back(SpectrumElement(9, 10)); // Peak max
        s_duplicate.push_back(SpectrumElement(10, 7)); // ...
        s_duplicate.push_back(SpectrumElement(11, 5)); // ...
        s_duplicate.push_back(SpectrumElement(12, 1));
        s_duplicate.push_back(SpectrumElement(12.1, 0.9)); // Peak right
        result = measureFullWidths(get_mz, get_int, s_duplicate.begin(), s_duplicate.end(), 0.51);
        shouldEqual(result.size(), (MzWidthPair::size_type)1);
        should(result.at(0).first == 9);
        shouldEqualTolerance(result.at(0).second, 4., 0.1);

        s_duplicate.push_back(SpectrumElement(12.1, 0.9));
        s_duplicate.push_back(SpectrumElement(12.1, 0.9));
        result = measureFullWidths(get_mz, get_int, s_duplicate.begin(), s_duplicate.end(), 0.51);
        shouldEqual(result.size(), (MzWidthPair::size_type)1);
        should(result.at(0).first == 9);
        shouldEqualTolerance(result.at(0).second, 4., 0.1);
        
        // test empty spectrum
        MSTK_LOG(logINFO) << "Testing with empty spectrum.";
        Spectrum s_empty;
        result = measureFullWidths(get_mz, get_int, s_empty.begin(), s_empty.end(), 0.5);
        should(result.empty());

        // test fraction range
        bool thrown = false;
        try {
            measureFullWidths(get_mz, get_int, s1.begin(), s1.end(), 0.);
        }
        catch(const PreconditionViolation& e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(!thrown);
        thrown = false;

        try {
            measureFullWidths(get_mz, get_int, s1.begin(), s1.end(), 1.);
        }
        catch(const PreconditionViolation& e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(!thrown);
        thrown = false;

        try {
            measureFullWidths(get_mz, get_int, s1.begin(), s1.end(), -0.3);
        }
        catch(const PreconditionViolation& e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        try {
            measureFullWidths(get_mz, get_int, s1.begin(), s1.end(), 1.3);
        }
        catch(const PreconditionViolation& e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
        thrown = false;

        MSTK_LOG(logINFO) << "Testing with realistic spectrum.";
	Spectrum spectrum;
	loadSpectrumElements(spectrum, dirTestdata + "/SpectrumAlgorithm/realistic_ms1.wsv");
        result = measureFullWidths(get_mz, get_int, spectrum.begin(), spectrum.end(), 0.5);
        /* result should be: (manually validated)
        879.98 0.0443713
        880.24 0.0195839
        880.45 0.0447405
        880.51 0.037717
        880.93 0.0320702
        881.12 0.0366487
        881.45 0.018084
        881.57 0.0425406
        881.64 0.0179091
        881.68 0.0195845
        */
        shouldEqual(result.at(0).first, 879.98);
        shouldEqualTolerance(result.at(0).second, 0.0443713, 0.000001); 
        shouldEqual(result.at(1).first, 880.24);
        shouldEqualTolerance(result.at(1).second, 0.0195839, 0.00001);
        shouldEqual(result.at(2).first, 880.45);
        shouldEqualTolerance(result.at(2).second, 0.0447405, 0.00001);
        shouldEqual(result.at(3).first, 880.51);
        shouldEqualTolerance(result.at(3).second, 0.037717, 0.00001);
        shouldEqual(result.at(4).first, 880.93);
        shouldEqualTolerance(result.at(4).second, 0.0320702, 0.00001);
        shouldEqual(result.at(5).first, 881.12);
        shouldEqualTolerance(result.at(5).second, 0.0366487, 0.00001);
        shouldEqual(result.at(6).first, 881.45);
        shouldEqualTolerance(result.at(6).second, 0.018084, 0.00001);
        shouldEqual(result.at(7).first, 881.57);
        shouldEqualTolerance(result.at(7).second, 0.0425406, 0.00001);
        shouldEqual(result.at(8).first, 881.64);
        shouldEqualTolerance(result.at(8).second, 0.0179091, 0.00001);
        shouldEqual(result.at(9).first, 881.68);
        shouldEqualTolerance(result.at(9).second, 0.0195845, 0.00001);   
    }
};

struct SpectralPeakTestSuite : vigra::test_suite {
    SpectralPeakTestSuite() : vigra::test_suite("SpectralPeak") {
        add( testCase(&SpectralPeakTestSuite::testHeight));
        add( testCase(&SpectralPeakTestSuite::testLowness));
        add( testCase(&SpectralPeakTestSuite::testFullWidthAtFractionOfMaximum));
    }

    void testHeight() {
        IntensityExtractor get_int;
        // A peak with height 3.1
        Spectrum s1;
        s1.push_back(SpectrumElement(1.1, 1.1));
        s1.push_back(SpectrumElement(1.2, 1.9));
        s1.push_back(SpectrumElement(1.4, 3.1));
        s1.push_back(SpectrumElement(1.5, 2.2));
        s1.push_back(SpectrumElement(1.69, 1.14));
        s1.push_back(SpectrumElement(1.76, 0.98));
	shouldEqual(SpectralPeak::height(get_int, s1.begin(), --(s1.end())), 3.1);
        
        // A empty spectrum
        Spectrum s2;
        bool thrown = false;
        try {
            SpectralPeak::height(get_int, s2.begin(), --(s2.end()));
        }
        catch(const PreconditionViolation& e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        should(thrown);
    }

    void testLowness() {
	IntensityExtractor get_int;
        // A quite normal spectral peak.
        //
        // The maximum abundance is 3.1
        // The lowest abundance on the 'left' is 1.1 and on the 'right' 0.98.
        // So, the lowness is (1 - 1.1/3.1) = 0.64516129.
        Spectrum s1;
        s1.push_back(SpectrumElement(1.1, 1.1));
        s1.push_back(SpectrumElement(1.2, 1.9));
        s1.push_back(SpectrumElement(1.4, 3.1));
        s1.push_back(SpectrumElement(1.5, 2.2));
        s1.push_back(SpectrumElement(1.69, 1.14));
        s1.push_back(SpectrumElement(1.76, 0.98));

        shouldEqual(SpectralPeak::lowness(get_int, s1.begin(), --(s1.end())), 1-(1.1/3.1));      

        // A peak with only one 'flank'.
        //
        // Lowness of a one flanked peak is 0.0.
        Spectrum s5;
        s5.push_back(SpectrumElement(1.1, 1.1));
        s5.push_back(SpectrumElement(1.2, 1.9));
        s5.push_back(SpectrumElement(1.4, 3.1));
        s5.push_back(SpectrumElement(1.5, 5.2));

        shouldEqual(SpectralPeak::lowness(get_int, s5.begin(), --(s5.end())), 0.0);     

        // An equiabundant sequence.
        //
        // The lowness of an equiabundant sequence is 0.0.
        Spectrum s2;
        s2.push_back(SpectrumElement(1.1, 1.1));
        s2.push_back(SpectrumElement(1.2, 1.1));
        s2.push_back(SpectrumElement(1.4, 1.1));
        s2.push_back(SpectrumElement(1.5, 1.1));

        shouldEqual(SpectralPeak::lowness(get_int, s2.begin(), --(s2.end())), 0.0);  

        // Quite unrealistic peak, with zero abundance elements
        //
        // Lowness is than a straight 1.0.
        Spectrum s6;
        s6.push_back(SpectrumElement(1.1, 0.1));
        s6.push_back(SpectrumElement(1.2, 0.0));
        s6.push_back(SpectrumElement(1.4, 1.1));
        s6.push_back(SpectrumElement(1.5, 1.2));
        s6.push_back(SpectrumElement(1.7, 0.0));
        s6.push_back(SpectrumElement(1.9, 1.1));
        s6.push_back(SpectrumElement(2.12, 0.9));

        shouldEqual(SpectralPeak::lowness(get_int, s6.begin(), --(s6.end())), 1.0);

        // Only one element.
        //
        // The lowness of one element is 0.0.
        Spectrum s3;
        s3.push_back(SpectrumElement(123.32, 89.1));

        shouldEqual(SpectralPeak::lowness(get_int, s3.begin(), --(s3.end())), 0.0);
    }

    void testFullWidthAtFractionOfMaximum() {
        MzExtractor get_mz;
        IntensityExtractor get_int;

        // A 'normal' peak.
        // Note the abundance twist in the last two elements.
        //
        // Fraction | Full width
        // 0.7      | 0.257459
        // 0.5      | 0.397029
        // 0.3      | lowness to small -> not defined     
        Spectrum s1;
        s1.push_back(SpectrumElement(0.4, 0.12));
        s1.push_back(SpectrumElement(1.1, 1.1));
        s1.push_back(SpectrumElement(1.2, 1.9));
        s1.push_back(SpectrumElement(1.4, 3.1));
        s1.push_back(SpectrumElement(1.5, 2.2));
        s1.push_back(SpectrumElement(1.6, 0.98));
        s1.push_back(SpectrumElement(1.69, 1.14));

        shouldEqualTolerance(SpectralPeak::fullWidthAtFractionOfMaximum(get_mz, get_int, s1.begin(), --(s1.end()), 0.7), 0.257459, 0.000001);
        shouldEqualTolerance(SpectralPeak::fullWidthAtFractionOfMaximum(get_mz, get_int, s1.begin(), --(s1.end()), 0.5), 0.397029, 0.000001);

        // not defined: full width at fraction of 0.3
        bool thrown = false;
        try {
            SpectralPeak::fullWidthAtFractionOfMaximum(get_mz, get_int, s1.begin(), --(s1.end()), 0.3);
        }
        catch(const Starvation& e) {
			MSTK_UNUSED(e);
            thrown = true;        
        }
        shouldMsg(thrown, "Starvation wasn't thrown despite of invalid input data.");
        thrown = false;

        // illegal fraction
        try {
            SpectralPeak::fullWidthAtFractionOfMaximum(get_mz, get_int, s1.begin(), --(s1.end()), 1.1);
        }
        catch(const PreconditionViolation& e) {
			MSTK_UNUSED(e);
            thrown = true;        
        }
        shouldMsg(thrown, "PreconditionViolation wasn't thrown despite of invalid fraction parameter.");
        thrown = false;

        try {
            SpectralPeak::fullWidthAtFractionOfMaximum(get_mz, get_int, s1.begin(), --(s1.end()), -0.3);
        }
        catch(const PreconditionViolation& e) {
			MSTK_UNUSED(e);
            thrown = true;        
        }
        shouldMsg(thrown, "PreconditionViolation wasn't thrown despite of invalid fraction parameter.");
        thrown = false;

        // legal border fraction values
        // no PreconditionViolation should be thrown here
        SpectralPeak::fullWidthAtFractionOfMaximum(get_mz, get_int, s1.begin(), --(s1.end()), 1.0);
        try {
            SpectralPeak::fullWidthAtFractionOfMaximum(get_mz, get_int, s1.begin(), --(s1.end()), 0.0);
        }
        catch(const Starvation& e) {
            MSTK_UNUSED(e);
        }

        // Special peak with elements exactly on target abundance
        Spectrum s_onTarget;
        s_onTarget.push_back(SpectrumElement(3, 7));
        s_onTarget.push_back(SpectrumElement(4, 10));
        s_onTarget.push_back(SpectrumElement(5, 7));

        shouldEqualTolerance(SpectralPeak::fullWidthAtFractionOfMaximum(get_mz, get_int, s_onTarget.begin(), --(s_onTarget.end()), 0.71), 2., 0.1);
    }
};

int main()
{
    SpectrumAlgorithmTestSuite test1;
    int failed = test1.run();
    std::cout << test1.report() << std::endl;

    SpectralPeakTestSuite test2;
    failed = test2.run();
    std::cout << test2.report() << std::endl;

    return failed;
}


