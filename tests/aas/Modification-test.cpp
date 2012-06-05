/*
 * Modification-test.cpp
 *
 * Copyright (c) 2011 Mathias Wilhelm
 * Copyright (c) 2011 Marc Kirchner
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

#include <MSTK/aas/Modification.hpp>
#include <MSTK/aas/Element.hpp>
#include <MSTK/aas/StoichiometryConfig.hpp>
#include <MSTK/aas/RawAminoAcid.hpp>

#include "unittest.hxx"

#include <iostream>

using namespace mstk;
using namespace aas;
using namespace aas::modifications;
using namespace aas::aminoAcids;
using namespace aas::stoichiometries;

/** Short description.
 * Long description.
 */
struct ModificationTestSuite : vigra::test_suite
{
    /** Constructor.
     * The ModificationTestSuite constructor adds all Modification tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    ModificationTestSuite() :
            vigra::test_suite("Modification")
    {
        add(testCase(&ModificationTestSuite::testModification));
        add(testCase(&ModificationTestSuite::testSpecificityModification));
        add(testCase(&ModificationTestSuite::testStoichiometryModification));
        add(testCase(&ModificationTestSuite::testModificationApplicable));
    }

    // testing getter/setter
    void testModification()
    {

        RawModificationImpl::RawModificationImplKeyType k = "Deamidated";
        Modification deamidated(k);

        shouldEqual(deamidated.getModificationId(), k);
        shouldEqual(deamidated.getName(), "Deamidated");
        shouldEqual(deamidated.getFullName(), "Deamidation");
        shouldEqual(deamidated.isVerified(), false);
        std::vector<String> altNames = deamidated.getAltNames();
        size_t altNameSize = 2;
        shouldEqual(altNames.size(), altNameSize);
        shouldEqual(altNames[0], "phenyllactyl from N-term Phe");
        shouldEqual(altNames[1], "Citrullination");

        aas::elements::Element H(1);
        aas::elements::Element N(7);
        aas::elements::Element O(8);
        Stoichiometry s = deamidated.getStoichiometry(), expected_s;
        expected_s.set(H, -1);
        expected_s.set(N, -1);
        expected_s.set(O, 1);
        shouldEqual(s, expected_s);
        shouldEqual(deamidated.getRawStoichiometry(), expected_s);

        std::vector<Specificity> specificities = deamidated.getSpecificities();
        size_t numberOfSpecificities = 4;
        shouldEqual(specificities.size(), numberOfSpecificities);
        Specificity spec0(RawAminoAcid('Q'), Specificity::ANYWHERE,
            Specificity::ARTEFACT);
        Specificity spec1(RawAminoAcid('R'), Specificity::ANYWHERE,
            Specificity::POST_TRANSLATIONAL);
        Specificity spec2(RawAminoAcid('N'), Specificity::ANYWHERE,
            Specificity::ARTEFACT);
        Specificity spec3(RawAminoAcid('F'), Specificity::PROTEIN_N_TERM,
            Specificity::POST_TRANSLATIONAL);
        shouldEqual(specificities[0], spec0);
        shouldEqual(specificities[1], spec1);
        shouldEqual(specificities[2], spec2);
        shouldEqual(specificities[3], spec3);

        Modification test;

        RawModification rm("Acetyl");

        test.setModification(rm);
        shouldEqual(test.getModification(), rm);
        test.setModification(rm.get_key());
        shouldEqual(test.getModification(), rm);

        shouldEqual(test.getName(), rm.get().getName());
        shouldEqual(test.getFullName(), rm.get(). getFullName());
        shouldEqual(test.isVerified(), rm.get().isVerified());
        shouldEqual(test.getStoichiometry(), rm.get().getStoichiometry());
        shouldEqual(test.getSpecificities(), rm.get().getSpecificities());
        shouldEqual(
            test.getStoichiometryConfig(),
            StoichiometryConfig(StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG));

        shouldEqual(test != deamidated, true);

        modifications::Modification tmp("Deamidated");
        tmp = test;
        shouldEqual(tmp, test);
    }

    // testing specificities of a modification (not rawmodification!)
    void testSpecificityModification()
    {

        Modification act("Acetyl");
        RawModificationImpl rawact("Acetyl");
        std::vector<Specificity> expectedCustom;

        shouldEqual(act.getSpecificities(), rawact.getSpecificities());
        shouldEqual(act.getRawSpecificities(), rawact.getSpecificities());
        shouldEqual(act.getCustomSpecificities(), expectedCustom);

        Specificity spec1(RawAminoAcid('A'), Specificity::ANYWHERE,
            Specificity::ARTEFACT);
        act.addCustomSpecificity(spec1);
        expectedCustom.push_back(spec1);

        shouldEqual(act.getSpecificities(), expectedCustom);
        shouldEqual(act.getCustomSpecificities(), expectedCustom);
        shouldEqual(act.getRawSpecificities(), rawact.getSpecificities());

        act.clearCustomSpecificities();
        expectedCustom.clear();

        shouldEqual(act.getSpecificities(), rawact.getSpecificities());
        shouldEqual(act.getRawSpecificities(), rawact.getSpecificities());
        shouldEqual(act.getCustomSpecificities(), expectedCustom);

        std::vector<Specificity> cSpecs;
        cSpecs.push_back(
            Specificity(RawAminoAcid('A'), Specificity::ANYWHERE,
                Specificity::POST_TRANSLATIONAL));
        act.setCustomSpecificities(cSpecs);
        shouldEqual(act.getSpecificities(), cSpecs);
        shouldEqual(act.getCustomSpecificities(), cSpecs);
    }

    // testing stoichiometry and stoichiometryconfig
    void testStoichiometryModification()
    {
        // setting up elements
        aas::elements::Element H(1);
        aas::elements::Element C(6);
        aas::elements::Element N(7);
        aas::elements::Element O(8);
        aas::elements::Element S(16);

        // creating default modification
        String k1 = "Acetyl";
        Modification m1(k1);
        shouldEqual(
            m1.getStoichiometryConfig(),
            StoichiometryConfig(StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG));
        Stoichiometry st1;
        st1.set(H, 2);
        st1.set(C, 2);
        st1.set(O, 1);
        shouldEqual(m1.getStoichiometry(), st1);

        // creating custom stoichiometry config
        std::vector<aas::elements::Isotope> is;
        is.push_back(aas::elements::Isotope(99.99, 1.11111111));
        Size freeID = aas::elements::ElementImpl::getNextId();
        aas::elements::addElement(freeID, "H", 1, is);

        StoichiometryConfigImpl sc("Experiment 1");
        // inserting element "automatically" by fw<ElementImpl>
        sc.insertElement(aas::elements::Element(freeID));
        // free function to add stoichiometry config
        addStoichiometryConfig(sc);

        // creating default modification
        Modification m2(k1);
        shouldEqual(m2.getStoichiometry(), st1);
        // changing stoichiometry config
        m2.setStoichiometryConfig("Experiment 1");
        shouldEqual(m2.getStoichiometryConfig(), sc);
        Stoichiometry st2;
        st2.set(aas::elements::Element(freeID), 2);
        st2.set(C, 2);
        st2.set(O, 1);
        shouldEqual(m2.getStoichiometry(), st2);

        // creating custom modification
        String k2 = "Amidated";
        Modification m3(k2, "Experiment 1");
        Stoichiometry st3;
        st3.set(aas::elements::Element(freeID), 1);
        st3.set(N, 1);
        st3.set(O, -1);
        shouldEqual(m3.getStoichiometry(), st3);
    }

    void testModificationApplicable()
    {
        Modification deam("Deamidated");
        Modification phos("Phosphopantetheine");
        RawAminoAcid prev('A');
        RawAminoAcid current('R');
        RawAminoAcid next('R');

        shouldEqual(deam.isApplicable(prev, current, next), true);
        shouldEqual(phos.isApplicable(prev, current, next), false);

        deam.addCustomSpecificity(
            Specificity(RawAminoAcid('A'), Specificity::ANYWHERE,
                Specificity::ARTEFACT));
        shouldEqual(deam.isApplicable(prev, current, next), false);
        deam.addCustomSpecificity(
            Specificity(RawAminoAcid('R'), Specificity::ANYWHERE,
                Specificity::PRE_TRANSLATIONAL));
        shouldEqual(deam.isApplicable(prev, current, next), true);
    }

};

/** The main function that runs the tests for class Modification.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    ModificationTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

