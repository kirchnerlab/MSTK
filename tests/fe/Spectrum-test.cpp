/*
 * Spectrum-test.cpp
 *
 * Copyright (C) 2012 Marc Kirchner
 * Copyright (c) 2008 Thorben Kroeger
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
#include <fstream>
#include <iostream>
#include <numeric>
// expose the class
#define private public
#define protected public
#include <MSTK/fe/types/Spectrum.hpp>
#undef private
#undef protected

using namespace mstk::fe;

class SpectrumTest {
public:

    void testMerging() {
        {
            Spectrum s1, s2;
            s1.push_back(Spectrum::Element(1, 1));
            s1.push_back(Spectrum::Element(2, 2));
            s1.push_back(Spectrum::Element(4, 4));
            s2.push_back(Spectrum::Element(3, 3));

            s1.merge(s2);
            shouldEqual(s1.size(), (unsigned int)4);
            for (int i = 0; i < 4; ++i) {
                shouldEqual(s1[i].mz, i + 1);
                shouldEqual(s1[i].abundance, i + 1);
            }
        }
        {
            Spectrum s1, s2;
            s1.push_back(Spectrum::Element(3, 3));
            s2.push_back(Spectrum::Element(1, 1));
            s2.push_back(Spectrum::Element(2, 2));
            s2.push_back(Spectrum::Element(4, 4));

            s1.merge(s2);
            shouldEqual(s1.size(), (unsigned int)4);
            for (int i = 0; i < 4; ++i) {
                shouldEqual(s1[i].mz, i + 1);
                shouldEqual(s1[i].abundance, i + 1);
            }
        }
        {
            Spectrum s1, s2;
            s1.push_back(Spectrum::Element(1, 3));
            s2.push_back(Spectrum::Element(1, 1));
            s1.merge(s2);
            shouldEqual(s1.size(), (unsigned int)1);
            shouldEqual(s1[0].mz, 1);
            shouldEqual(s1[0].abundance, 4);
        }
        {
            Spectrum s1, s2;
            s1.merge(s2);
            shouldEqual(s1.size(), (unsigned int)0);
        }
    }
    
    void testShiftMz() {
        //test all spectrum shifting functions

        Spectrum s, s2, s3, s4, s5;
        s.push_back(Spectrum::Element(1.0, 2.0));
        s.push_back(Spectrum::Element(1.5, 3.0));
        s2 = s;
        s3 = s;
        s4 = s;
        s5 = s;

        std::transform(s.begin(), s.end(), s.begin(), Spectrum::ShiftMz(0.5));
        s2.shiftBy(0.5);
        s3.shiftTo(1.5);
        s4.shiftMaxToMonoisotopicMass();
        s5.shiftBy(-0.5);

        shouldEqualTolerance(s[0].abundance, 2.0, 1E-6);
        shouldEqualTolerance(s[0].mz, 1.5, 1E-6);
        shouldEqualTolerance(s[1].abundance, 3.0, 1E-6);
        shouldEqualTolerance(s[1].mz, 2.0, 1E-6);
        
        shouldEqual(s, s2);
        shouldEqual(s, s3);
        shouldEqual(s2, s3);
        shouldEqual(s4, s5);
    }

    void testMaxAbundance() {
        //test Spectrum::getMaxAbundancePeak(), getTotalAbundance()

        int n=10;
        Spectrum s;
        int totalab=0;
        double mzab=0;
        for (int i=0; i<n; i++){
            s.push_back(Spectrum::Element(i, i+1));
            totalab+=i+1;
            mzab += i*(i+1);
        }
        int mzmax = 33;
        int abmax = 100;
        s.push_back(Spectrum::Element(mzmax, abmax));
        totalab+=abmax;
        mzab += mzmax*abmax;
        Spectrum::iterator iter;
        iter = s.getMaxAbundancePeak();
        shouldEqualTolerance((*iter).mz, mzmax, 1E-6);
        shouldEqualTolerance((*iter).abundance, abmax, 1E-6);

        //what if it's not at the end?
        mzmax = 50;
        abmax = abmax+1;
        s.insert(s.begin(), Spectrum::Element(mzmax, abmax));
        totalab+=abmax;
        mzab+=mzmax*abmax;
        iter = s.getMaxAbundancePeak();
        shouldEqualTolerance((*iter).mz, mzmax, 1E-6);
        shouldEqualTolerance((*iter).abundance, abmax, 1E-6);
        shouldEqualTolerance(s.getTotalAbundance(), totalab, 1E-6);

    }
    
    void test_LessThanMzScalar() {
        Spectrum s;
        s.push_back(Spectrum::Element(1,2)); //0
        s.push_back(Spectrum::Element(2,9)); //1
        s.push_back(Spectrum::Element(3,7)); //2
        s.push_back(Spectrum::Element(4,2)); //3
        
        Spectrum::iterator it;
        it = std::upper_bound(s.begin(), s.end(), 3.5, Spectrum::LessThanMz<double, Spectrum::Element>());
        shouldEqual(std::distance(s.begin(), it), 3);
    }
    void test_LessThanAbundanceScalar() {
        Spectrum s;
        s.push_back(Spectrum::Element(1,1)); //0
        s.push_back(Spectrum::Element(2,2)); //1
        s.push_back(Spectrum::Element(3,3)); //2
        s.push_back(Spectrum::Element(4,4)); //3
        
        Spectrum::iterator it;
        it = std::upper_bound(s.begin(), s.end(), 3.5, Spectrum::LessThanAbundance<double, Spectrum::Element>());
        shouldEqual(std::distance(s.begin(), it), 3);
    }
    void test_LessThanRt() {
        std::vector<Spectrum> ss;
        Spectrum s1; s1.setRetentionTime(1); ss.push_back(s1);
        Spectrum s2; s2.setRetentionTime(3); ss.push_back(s2);
        Spectrum s3; s3.setRetentionTime(2); ss.push_back(s3);
        
        std::sort(ss.begin(), ss.end(), Spectrum::LessThanRt<Spectrum, Spectrum>());
        shouldEqual(ss[0].getRetentionTime(), 1);
        shouldEqual(ss[1].getRetentionTime(), 2);
        shouldEqual(ss[2].getRetentionTime(), 3);
    }
    void test_LessThanRtScalar() {
        std::vector<Spectrum> ss;
        Spectrum s1; s1.setRetentionTime(1); ss.push_back(s1); //0
        Spectrum s2; s2.setRetentionTime(2); ss.push_back(s2); //1
        Spectrum s3; s3.setRetentionTime(3); ss.push_back(s3); //2
        
        std::vector<Spectrum>::iterator it;
        it = std::upper_bound(ss.begin(), ss.end(), 2.5, Spectrum::LessThanRt<double, Spectrum>());
        shouldEqual(std::distance(ss.begin(), it), 2);
    }
    void test_LessThanMsLevel() {
        std::vector<Spectrum> ss;
        Spectrum s1; s1.setMsLevel(1); ss.push_back(s1);
        Spectrum s2; s2.setMsLevel(3); ss.push_back(s2);
        Spectrum s3; s3.setMsLevel(2); ss.push_back(s3);
        
        std::sort(ss.begin(), ss.end(), Spectrum::LessThanMsLevel<Spectrum, Spectrum>());
        shouldEqual(ss[0].getMsLevel(), 1u);
        shouldEqual(ss[1].getMsLevel(), 2u);
        shouldEqual(ss[2].getMsLevel(), 3u);
    }
    void test_LessThanMsLevelScalar() {
        std::vector<Spectrum> ss;
        Spectrum s1; s1.setMsLevel(1); ss.push_back(s1); //0
        Spectrum s2; s2.setMsLevel(2); ss.push_back(s2); //1
        Spectrum s3; s3.setMsLevel(3); ss.push_back(s3); //2
        
        std::vector<Spectrum>::iterator it;
        it = std::upper_bound(ss.begin(), ss.end(), 2.5, Spectrum::LessThanMsLevel<double, Spectrum>());
        shouldEqual(std::distance(ss.begin(), it), 2);
    }
    void test_EqualMsLevel() {
        std::vector<Spectrum> ss;
        Spectrum s1; s1.setMsLevel(1); ss.push_back(s1); //0
        Spectrum s2; s2.setMsLevel(2); ss.push_back(s2); //1
        Spectrum s3; s3.setMsLevel(3); ss.push_back(s3); //2
        
        std::vector<Spectrum>::iterator it;
        it = std::find_if(ss.begin(), ss.end(), Spectrum::EqualMsLevel(2) );
        shouldEqual(std::distance(ss.begin(), it), 1);
    }
    void test_SumAbundance() {
        Spectrum s;
        s.push_back(Spectrum::Element(1,2)); //0
        s.push_back(Spectrum::Element(2,9)); //1
        
        double ab = std::accumulate(s.begin(), s.end(), 0.0, Spectrum::SumAbundance());
        shouldEqualTolerance(ab, 11, 1E-7);
    }
    void test_ConstructionFromIterators() {
        Spectrum s1;
        s1.push_back(Spectrum::Element(1,2));
        s1.push_back(Spectrum::Element(2,3));
        s1.push_back(Spectrum::Element(3,4));
        s1.push_back(Spectrum::Element(4,5));
        Spectrum::iterator it1 = s1.begin(); ++it1;
        Spectrum::iterator it2 = s1.begin(); it2 += 2;
        
        Spectrum s2(it1, it2);
        shouldEqual(s2.size(), 1u);
        shouldEqual(s2[0].abundance, 3);
        
        Spectrum::const_iterator cit1 = s1.begin(); ++cit1;
        Spectrum::const_iterator cit2 = s1.begin(); cit2 += 2;
        
        Spectrum s3(cit1, cit2);
        shouldEqual(s3.size(), 1u);
        shouldEqual(s3[0].abundance, 3);
    }

};

struct SpectrumTestSuite : public vigra::test_suite {
    SpectrumTestSuite() : vigra::test_suite("Spectrum Tests") {
        add( testCase(&SpectrumTest::testMerging ));
        add( testCase(&SpectrumTest::testShiftMz ));
        add( testCase(&SpectrumTest::testMaxAbundance ));
        add( testCase(&SpectrumTest::test_LessThanMzScalar ));
        add( testCase(&SpectrumTest::test_LessThanAbundanceScalar ));
        add( testCase(&SpectrumTest::test_LessThanRt ));
        add( testCase(&SpectrumTest::test_LessThanRtScalar ));
        add( testCase(&SpectrumTest::test_LessThanMsLevel ));
        add( testCase(&SpectrumTest::test_LessThanMsLevelScalar ));
        add( testCase(&SpectrumTest::test_EqualMsLevel ));
        add( testCase(&SpectrumTest::test_SumAbundance ));
        add( testCase(&SpectrumTest::test_ConstructionFromIterators ));
    }
};

int main()
{
    SpectrumTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}
