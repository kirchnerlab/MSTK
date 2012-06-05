/*
 * Residue-test.cpp
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

#include <MSTK/aas/Residue.hpp>
#include <MSTK/common/Error.hpp>

#include "unittest.hxx"

#include <iostream>

using namespace mstk;
using namespace aas;
using namespace aas::aminoAcids;
using namespace aas::modifications;
using namespace aas::stoichiometries;

/** Short description.
 * Long description.
 */
struct ResidueTestSuite : vigra::test_suite
{
    /** Constructor.
     * The ResidueTestSuite constructor adds all Residue tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    ResidueTestSuite() :
            vigra::test_suite("Residue")
    {
        add(testCase(&ResidueTestSuite::testResidue));
        add(testCase(&ResidueTestSuite::testResidueShared));
        add(testCase(&ResidueTestSuite::testResidueStoichiometry));
    }

    void testResidue()
    {

        RawAminoAcidImpl::RawAminoAcidImplKeyType aa_k = 'A';
        AminoAcid aa(aa_k);
        RawModificationImpl::RawModificationImplKeyType m_k = "Phospho";
        RawModificationImpl::RawModificationImplKeyType l_k = "ESP";
        Modification m(m_k);
        Modification l(l_k);

        Residue r(aa);

        shouldEqual(r.isCTerm(), false);
        shouldEqual(r.isNTerm(), false);
        shouldEqual(r.isModified(), false);
        shouldEqual(r.isLabeled(), false);

        r.changeType(aa);
        shouldEqual(r.getAminoAcid(), aa);

        r.changeType(aa_k);
        shouldEqual(r.getAminoAcid(), aa);

        Bool thrown = false;
        try {
            r.setModification(l);
        } catch (mstk::LogicError& e) {
            thrown = true;
        }shouldEqual(thrown, true);

        r.setModification(m);
        shouldEqual(r.getModification(), m);
        shouldEqual(r.isModified(), true);
        shouldEqual(r.hasModification(m), true);
        shouldEqual(r.hasModification(Modification("Oxidation")), false);

        thrown = false;
        try {
            r.setIsotopicLabel(m);
        } catch (mstk::LogicError& e) {
            thrown = true;
        }shouldEqual(thrown, true);

        r.setIsotopicLabel(l);
        shouldEqual(r.getIsotopicLabel(), l);
        shouldEqual(r.isLabeled(), true);
        shouldEqual(r.hasLabel(l), true);
        shouldEqual(r.hasLabel(l.getModificationId()), true);
        shouldEqual(r.hasLabel(Modification("Oxidation")), false);

        r.setIsotopicLabel(l.getModificationId());
        shouldEqual(r.getIsotopicLabel(), l);
        shouldEqual(r.isLabeled(), true);
        shouldEqual(r.hasLabel(l), true);
        shouldEqual(r.hasLabel(Modification("Oxidation")), false);

        r.removeIsotopicLabel();
        shouldEqual(r.getIsotopicLabel(), modifications::Modification(""));
        shouldEqual(r.isLabeled(), false);
        shouldEqual(r.hasLabel(l), false);

        r.changeType(RawAminoAcidImpl::PEPTIDE_C_TERM);
        shouldEqual(r.isCTerm(), true);
        r.changeType(RawAminoAcidImpl::PROTEIN_C_TERM);
        shouldEqual(r.isCTerm(), true);
        r.changeType(RawAminoAcidImpl::PEPTIDE_N_TERM);
        shouldEqual(r.isNTerm(), true);
        r.changeType(RawAminoAcidImpl::PROTEIN_N_TERM);
        shouldEqual(r.isNTerm(), true);

        Residue r1 = r;
        r1.changeType('C');
        r1.removeModification();
        r1.removeIsotopicLabel();

        shouldEqual(r != r1, true);
        shouldEqual(r == r1, false);
    }

    void testResidueShared()
    {
        Residue t('A', "Oxidation", "ICAT-G");
        shouldEqual(t.getAminoAcid(), AminoAcid('A'));
        shouldEqual(t.isModified(), true);
        shouldEqual(t.isLabeled(), true);

        Residue t_c(AminoAcid('A'), Modification("Phospho"),
            Modification("TMT"));
        t_c = t;

        StoichiometryConfigImpl::StoichiometryConfigImplKeyType k = "different";
        StoichiometryConfigImpl sci(k);
        addStoichiometryConfig(sci);

        t_c.applyIsotopicLabelStoichiometryConfig(k);
        t_c.applyModificationStoichiometryConfig(k);

        shouldEqual(t.getModification().getStoichiometryConfig().get_key(),
            StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);
        shouldEqual(t.getIsotopicLabel().getStoichiometryConfig().get_key(),
            StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);

        shouldEqual(t_c.getModification().getStoichiometryConfig().get_key(),
            k);
        shouldEqual(t_c.getIsotopicLabel().getStoichiometryConfig().get_key(),
            k);

        shouldEqual(t.getStoichiometry(), t_c.getStoichiometry());
    }

    void testResidueStoichiometry()
    {
        elements::Element H(1);
        elements::Element C(6);
        elements::Element N(7);
        elements::Element O(8);
        elements::Element S(16);

        Residue r('M');

        Stoichiometry expectedS;
        expectedS.set(H, 9);
        expectedS.set(C, 5);
        expectedS.set(N, 1);
        expectedS.set(O, 1);
        expectedS.set(S, 1);

        shouldEqual(r.getStoichiometry(), expectedS);
        r.setModification("Oxidation");
        expectedS.add(O, 1);
        shouldEqual(r.getStoichiometry(), expectedS);

        std::vector<elements::Isotope> cHi;
        elements::Element cH(
            elements::ElementImpl(elements::ElementImpl::getNextId(), "H", 1,
                cHi));
        StoichiometryConfigImpl sci("test1");
        sci.insertElement(cH);
        StoichiometryConfig sc(sci);

        r.applyAminoAcidStoichiometryConfig(
            StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);

        shouldEqual(r.getStoichiometry(), expectedS);

        r.applyAminoAcidStoichiometryConfig(sc);

        expectedS.set(H, 0);
        expectedS.set(cH, 9);

        shouldEqual(r.getStoichiometry(), expectedS);

        std::vector<elements::Isotope> cOi;
        elements::Element cO(
            elements::ElementImpl(elements::ElementImpl::getNextId(), "O", 8,
                cOi));
        StoichiometryConfigImpl scim("test2");
        scim.insertElement(cO);
        StoichiometryConfig scm(scim);
        r.applyModificationStoichiometryConfig("test2");

        expectedS.set(O, 1);
        expectedS.set(cO, 1);

        shouldEqual(r.getStoichiometry(), expectedS);

        r.setIsotopicLabel("ESP");
        expectedS.add(H, 26);
        expectedS.add(C, 16);
        expectedS.add(N, 4);
        expectedS.add(O, 2);
        expectedS.add(S, 1);

        shouldEqual(r.getStoichiometry(), expectedS);

        r.applyIsotopicLabelStoichiometryConfig(
            StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);
        shouldEqual(r.getStoichiometry(), expectedS);

        std::vector<elements::Isotope> cSi;
        elements::Element cS(
            elements::ElementImpl(elements::ElementImpl::getNextId(), "S", 16,
                cSi));
        StoichiometryConfigImpl scil("test3");
        scil.insertElement(cS);
        StoichiometryConfig scl(scil);
        r.applyIsotopicLabelStoichiometryConfig(scl);

        expectedS.add(S, -1);
        expectedS.add(cS, 1);
        shouldEqual(r.getStoichiometry(), expectedS);
    }

};

/** The main function that runs the tests for class Residue.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    ResidueTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

