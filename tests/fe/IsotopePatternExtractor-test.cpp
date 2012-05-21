/*
 * IsotopePatternExtractor-test.cpp
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
#include "unittest.hxx"
#include "utilities.hpp"
#include <MSTK/fe/IsotopePatternExtractor.hpp>
#include <MSTK/fe/types/XicFbiTraits.hpp>
#include <MSTK/fe/UncenteredCorrelation.hpp>
#include <MSTK/fe/NopSplitter.hpp>
#include <MSTK/common/Types.hpp>

using namespace mstk::fe;
using namespace mstk;

struct IsotopePatternExtractorTestSuite : vigra::test_suite
{
    IsotopePatternExtractorTestSuite() :
            vigra::test_suite("IsotopePatternExtractor")
    {
        add(testCase(&IsotopePatternExtractorTestSuite::testSingleXic));
    }

    void testSingleXic()
    {
        std::vector < Xic > xics;
        // generate the input data
        double mz[] = { 500.001, 500.003, 500.002 };
        double rt[] = { 10.0, 11.0, 12.0 };
        unsigned int sn[] = { 2, 3, 4 };
        double ab[] = { 1.0, 2.0, 1.0 };
        Xic x = makeXic(3, mz, rt, sn, ab);
        xics.push_back(x);

        // charges
        std::set<int> charges;
        charges.insert(1);

        // mz shifts (the std MaxQuant distance)
        std::vector<double> mzShifts;
        mzShifts.push_back(1.00286864);

        // the box generators
        std::vector < XicBoxGenerator > boxGenerators;
        typedef std::vector<double>::const_iterator ShiftIT;
        for (ShiftIT s = mzShifts.begin(); s != mzShifts.end(); ++s) {
            typedef std::set<int>::const_iterator ChargeIT;
            for (ChargeIT c = charges.begin(); c != charges.end(); ++c) {
                boxGenerators.push_back(
                    XicBoxGenerator(*s, std::make_pair(2.0, 20.0),
                        std::make_pair(2.0, 10.0), *c));
            }
        }

        typedef IsotopePatternExtractor<UncenteredCorrelation,
                NopSplitter> MyIsotopePatternExtractor;

        MyIsotopePatternExtractor ipe;
        std::vector < IsotopePattern > ips;
        Size n = ipe(xics, boxGenerators, 0.6, 1, ips);

        // make sure we only get one isotope distribution
        shouldEqual(n, ips.size());
        shouldEqual(ips.size(), static_cast<size_t>(1));
        n = ipe(xics, boxGenerators, 0.6, 2, ips);

        // We do not expect to find anything.
        shouldEqual(ips.empty(), true);
    }
}
;

int main()
{
    IsotopePatternExtractorTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

