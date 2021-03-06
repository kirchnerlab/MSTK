/*
 * AminoAcid-test.cpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
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

#include <MSTK/aas/AminoAcid.hpp>

#include "unittest.hxx"

#include <iostream>

using namespace mstk;
using namespace aas;
using namespace aas::aminoAcids;
using namespace aas::stoichiometries;

/** Short description.
 * Long description.
 */
struct AminoAcidTestSuite : vigra::test_suite
{
    /** Constructor.
     * The AminoAcidTestSuite constructor adds all AminoAcid tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    AminoAcidTestSuite() :
            vigra::test_suite("AminoAcid")
    {
        add(testCase(&AminoAcidTestSuite::testAminoAcid));
        add(testCase(&AminoAcidTestSuite::testStoichiometryAminoAcid));

    }

    void testAminoAcid()
    {
        AminoAcid aa('Y');
        RawAminoAcidImpl raa('Y');
        RawAminoAcid raa_r('Y');

        shouldEqual(aa.getSymbol(), raa.getSymbol());
        shouldEqual(aa.getThreeLetterCode(), raa.getThreeLetterCode());
        shouldEqual(aa.getFullName(), raa.getFullName());
        shouldEqual(aa.getStoichiometry(), raa.getStoichiometry());
        shouldEqual(aa.getRawAminoAcidKey(), raa.getId());
        shouldEqual(aa.getRawAminoAcid(), raa);
        shouldEqual(&(aa.getRawAminoAcid().get()), &(raa_r.get()));
    }

    void testStoichiometryAminoAcid()
    {
        // setting up elements
        aas::elements::Element H(1);
        aas::elements::Element C(6);
        aas::elements::Element N(7);
        aas::elements::Element O(8);
        aas::elements::Element S(16);

        AminoAcid aa('C');

        shouldEqual(
            aa.getStoichiometryConfig(),
            StoichiometryConfig(StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG));
        Stoichiometry expected_s;
        expected_s.set(H, 5);
        expected_s.set(C, 3);
        expected_s.set(N, 1);
        expected_s.set(O, 1);
        expected_s.set(S, 1);
        shouldEqual(aa.getStoichiometry(), expected_s);

        // creating custom stoichiometry config
        std::vector<aas::elements::Isotope> is;
        is.push_back(aas::elements::Isotope(99.99, 1.11111111));
        Size freeID = aas::elements::ElementImpl::getNextId();
        aas::elements::addElement(freeID, "H", 1, is);

        StoichiometryConfigImpl::StoichiometryConfigImplKeyType sc_key =
                "Experiment 1";
        StoichiometryConfigImpl sc(sc_key);
        // inserting element "automatically" by fw<ElementImpl>
        sc.insertElement(aas::elements::Element(freeID));
        // free function to add stoichiometry config
        addStoichiometryConfig(sc);

        Stoichiometry expected_cs;
        expected_cs.set(aas::elements::Element(freeID), 5);
        expected_cs.set(C, 3);
        expected_cs.set(N, 1);
        expected_cs.set(O, 1);
        expected_cs.set(S, 1);

        aa.setStoichiometryConfig(sc_key);
        shouldEqual(aa.getStoichiometry(), expected_cs);
        aa.setStoichiometryConfig(
            StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);
        shouldEqual(aa.getStoichiometry(), expected_s);

        AminoAcid aa2('C', sc_key);
        shouldEqual(aa2.getStoichiometry(), expected_cs);

        shouldEqual(aa2 != aa, true);
    }

};

/** The main function that runs the tests for class AminoAcid.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    AminoAcidTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

