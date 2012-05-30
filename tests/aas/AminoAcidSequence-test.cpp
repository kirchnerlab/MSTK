/*
 * AminoAcidSequence-test.cpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
 *
 */

#include <MSTK/aas/AminoAcidSequence.hpp>
#include <MSTK/common/Error.hpp>

#include "unittest.hxx"

#include <iostream>

using namespace mstk;
using namespace aas;
using namespace aas::stoichiometries;

/** Short description.
 * Long description.
 */
struct AminoAcidSequenceTestSuite : vigra::test_suite
{
    /** Constructor.
     * The AminoAcidSequenceTestSuite constructor adds all AminoAcidSequence tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    AminoAcidSequenceTestSuite() :
            vigra::test_suite("AminoAcidSequence")
    {
        add(testCase(&AminoAcidSequenceTestSuite::testAminoAcidSequence));
        add(
            testCase(&AminoAcidSequenceTestSuite::testAminoAcidSequenceCollection));
        add(
            testCase(&AminoAcidSequenceTestSuite::testAminoAcidSequenceSequenceAltering));
        add(
            testCase(&AminoAcidSequenceTestSuite::testAminoAcidSequenceApplyModifications));
        add(
            testCase(&AminoAcidSequenceTestSuite::testAminoAcidSequenceAminoAcidStoichiometry));
        add(testCase(&AminoAcidSequenceTestSuite::testCollection));
    }

    void testAminoAcidSequence()
    {
        String aass = "AACCCQ";
        String mods = "Phospho(C)@3; ICAT-G(C)@4; Oxidation(C)@5; ICAT-G(C)@5";
        String aassm = "AAC(Phospho)C(ICAT-G)C(Oxidation; ICAT-G)Q";

        AminoAcidSequence aas(aass);
        shouldEqual(aas.toString(), aass);
        shouldEqual(aas.toString(true), "0" + aass + "1");
        shouldEqual(aas.toUnmodifiedSequenceString(), aass);

        aas.applyModificationAtPosition("Phospho", 3);
        aas.applyModificationAtPosition("ICAT-G", 4);
        aas.applyModificationAtPosition("ICAT-G", 5);
        aas.applyModificationAtPosition("Oxidation", 5);

        shouldEqual(aas.toString(), aassm);
        shouldEqual(aas.toString(true), "0" + aassm + "1");
        shouldEqual(aas.getModificationString(), mods);

        shouldEqual(aas[3].isLabeled(), false);
        shouldEqual(aas[3].isModified(), true);
        shouldEqual(aas[3].hasModification("Phospho"), true);
        shouldEqual(aas[4].isModified(), false);
        shouldEqual(aas[4].isLabeled(), true);
        shouldEqual(aas[4].hasLabel("ICAT-G"), true);
        shouldEqual(aas[5].isLabeled(), true);
        shouldEqual(aas[5].hasLabel("ICAT-G"), true);
        shouldEqual(aas[5].isModified(), true);
        shouldEqual(aas[5].hasModification("Oxidation"), true);

        const AminoAcidSequence& aas_c = aas;
        const Residue& r1 = aas_c[1];
        shouldEqual(r1.isModified(), false);
        shouldEqual(r1.hasModification("Phospho"), false);
        shouldEqual(r1.isNTerm(), false);
        shouldEqual(r1.isCTerm(), false);

        AminoAcidSequence aasc(aas.begin(), aas.end());
        shouldEqual(aasc, aas);

        // remove oxidation
        shouldEqual(aasc[5].isModified(), true);
        aasc.remove("Oxidation");
        shouldEqual(aasc[5].isModified(), false);

        modifications::Modification cphospho("Phospho");
        StoichiometryConfig sc(StoichiometryConfigImpl("test"));
        cphospho.setStoichiometryConfig(sc);
        aasc.applyModificationAtPosition("Phospho", 4);
        aasc.applyModificationAtPosition(cphospho, 5);
        shouldEqual(aasc[4].isModified(), true);
        shouldEqual(aasc[5].isModified(), true);
        // remove default phospho but not the custom phospho since the stoichiometry config is different
        aasc.remove(modifications::Modification("Phospho"));
        shouldEqual(aasc[3].isModified(), false);
        shouldEqual(aasc[4].isModified(), false);
        shouldEqual(aasc[5].isModified(), true);

        // test apply stoichiometry config
        StoichiometryConfig dsc(
            StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);
        aas.applyAminoAcidStoichiometryConfig(sc);
        typedef AminoAcidSequence::const_iterator IT;
        for (IT it = aas.begin(); it != aas.end(); ++it) {
            shouldEqual(it->getAminoAcid().getStoichiometryConfig(), sc);
//            shouldEqual(it->getModification()->getStoichiometryConfig(), dsc);
        }
        aas.applyModificationStoichiometryConfig(sc);
        aas.applyAminoAcidStoichiometryConfig(dsc);
        typedef AminoAcidSequence::const_iterator IT;
        for (IT it = aas.begin(); it != aas.end(); ++it) {
            shouldEqual(it->getAminoAcid().getStoichiometryConfig(), dsc);
//            shouldEqual(it->getModification()->getStoichiometryConfig(), sc);
        }

        AminoAcidSequence tmp1("ASD");
        AminoAcidSequence tmp2(tmp1.begin() + 1, tmp1.begin() + 4);
        shouldEqual(tmp1, tmp2);
    }

    void testAminoAcidSequenceCollection()
    {
        // TODO we might want to implement those functions similar to push_back in aminoacidsequence, since these does not make sure that the peptide sequence consists of an n- and c-term
        String seq = "0ACGT1";
        AminoAcidSequence aas(seq);

        aas.assign(3, Residue('C'));
        shouldEqual(aas.toString(true), "CCC");

        aas.insert(aas.begin(), Residue('A'));
        shouldEqual(aas.toString(true), "ACCC");

        aas.insert(aas.begin() + 1, 2, Residue('D'));
        shouldEqual(aas.toString(true), "ADDCCC");

        aas.erase(aas.begin() + 1);
        aas.erase(aas.begin() + 1);
        shouldEqual(aas.toString(true), "ACCC");

        aas.erase(aas.begin() + 1, aas.end());
        shouldEqual(aas.toString(true), "A");

        AminoAcidSequence aas1(seq);
        Size i = seq.size() - 1;
        for (AminoAcidSequence::reverse_iterator it = aas1.rbegin();
                it != aas1.rend(); ++it) {
            shouldEqual(it->getAminoAcid().getSymbol(), seq[i--]);
        }

        shouldEqual(aas1.at(1u), Residue('A'));
        const AminoAcidSequence caas1 = aas1;
        shouldEqual(caas1.at(1u), Residue('A'));

        shouldEqual(aas1.max_size() >= aas1.size(), true);

        aas1.resize(100u);
        shouldEqual(aas1.size(), 100u);

        aas1.reserve(200u);
        shouldEqual(aas1.capacity(), 200u);

        AminoAcidSequence aasc = aas, aas1c = aas1;

        aas1.swap(aas);
        shouldEqual(aas1, aasc);
        shouldEqual(aas, aas1c);
    }

    void testAminoAcidSequenceSequenceAltering()
    {
        AminoAcidSequence ass("");
        size_t expectedSize = 2;
        shouldEqual(ass.size(), expectedSize);
        ass[0].changeType('A');
        ass[1].changeType('A');

        bool thrown = false;
        try {
            ass.makePeptideCTerm();
        } catch (RuntimeError& e) {
            thrown = true;
        }shouldEqual(thrown, true);

        thrown = false;
        try {
            ass.makePeptideNTerm();
        } catch (RuntimeError& e) {
            thrown = true;
        }shouldEqual(thrown, true);

        thrown = false;
        try {
            ass.makeProteinCTerm();
        } catch (RuntimeError& e) {
            thrown = true;
        }shouldEqual(thrown, true);

        thrown = false;
        try {
            ass.makeProteinNTerm();
        } catch (RuntimeError& e) {
            thrown = true;
        }shouldEqual(thrown, true);

        ass = AminoAcidSequence("");

        // testing different situation of push back and pop_back
        ass.push_back('A');
        expectedSize = 3;
        shouldEqual(ass[0], Residue('0'));
        shouldEqual(ass[1], Residue('A'));
        shouldEqual(ass[2], Residue('1'));
        shouldEqual(ass.size(), expectedSize);

        ass.makeProteinNTerm();
        shouldEqual(ass[0], Residue('2'));
        ass.makeProteinCTerm();
        shouldEqual(ass[2], Residue('3'));
        ass.makePeptideNTerm();
        shouldEqual(ass[0], Residue('0'));
        ass.makePeptideCTerm();
        shouldEqual(ass[2], Residue('1'));

        ass.push_back('C');
        expectedSize = 4;
        shouldEqual(ass[0], Residue('0'));
        shouldEqual(ass[1], Residue('A'));
        shouldEqual(ass[2], Residue('C'));
        shouldEqual(ass[3], Residue('1'));
        shouldEqual(ass.size(), expectedSize);

        ass.pop_back();
        expectedSize = 3;
        shouldEqual(ass[0], Residue('0'));
        shouldEqual(ass[1], Residue('A'));
        shouldEqual(ass[2], Residue('1'));
        shouldEqual(ass.size(), expectedSize);

        ass.push_back(Residue('3'));
        expectedSize = 3;
        shouldEqual(ass[0], Residue('0'));
        shouldEqual(ass[1], Residue('A'));
        shouldEqual(ass[2], Residue('3'));
        shouldEqual(ass.size(), expectedSize);

        // push back C term to empty aas
        AminoAcidSequence ass2("");
        ass2.push_back(Residue('3'));
        expectedSize = 2;
        shouldEqual(ass2[0], Residue('0'));
        shouldEqual(ass2[1], Residue('3'));
        shouldEqual(ass2.size(), expectedSize);

        ass2.pop_back();
        expectedSize = 2;
        shouldEqual(ass2[0], Residue('0'));
        shouldEqual(ass2[1], Residue('3'));
        shouldEqual(ass2.size(), expectedSize);

        // push back C term to empty ass
        AminoAcidSequence ass3("");
        ass3.push_back(Residue('3'));
        expectedSize = 2;
        shouldEqual(ass3[0], Residue('0'));
        shouldEqual(ass3[1], Residue('3'));
        shouldEqual(ass3.size(), expectedSize);

        // testing append in various situations
        AminoAcidSequence s1("ACA");
        AminoAcidSequence s2("GTG");
        s2.makeProteinCTerm();
        s2.makeProteinNTerm();
        AminoAcidSequence s3("");
        AminoAcidSequence s4("");

        s1.append(s2);
        shouldEqual(s1.toString(true), "0ACAGTG3");
        s2.append(s2);
        shouldEqual(s2.toString(true), "2GTGGTG3");

        s2.append(s3);
        shouldEqual(s2.toString(true), "2GTGGTG1");

        s3.append(s1);
        shouldEqual(s3.toString(true), "0ACAGTG3");

        s4.append(s4);
        shouldEqual(s4.toString(true), "01");

        s4.append(s2);
        shouldEqual(s4.toString(true), "0GTGGTG1");

        AminoAcidSequence tmp("");
        // TODO is clear supposed to do this?
        tmp.clear();
        tmp.push_back(Residue('A'));
        shouldEqual(tmp.toString(true), "0A1");
        tmp[2].changeType('C');
        tmp.push_back('D');
        shouldEqual(tmp.toString(true), "0ACD1");
        tmp[4].changeType('A');
        shouldEqual(tmp.toString(true), "0ACDA");
        tmp.pop_back();
        shouldEqual(tmp.toString(true), "0ACD1");

        tmp.clear();
        tmp.append(s2);
        shouldEqual(tmp.toString(true), "2GTGGTG1");
        tmp.clear();
        s2[0].changeType('A');
        tmp.append(s2);
        shouldEqual(tmp.toString(true), "0AGTGGTG1");
    }

    void testAminoAcidSequenceAminoAcidStoichiometry()
    {
        String aass = "ACQT";
        AminoAcidSequence ass(aass);

        elements::Element H(1);
        elements::Element C(6);
        elements::Element N(7);
        elements::Element O(8);
        elements::Element S(16);

        Stoichiometry expectedS;
        expectedS.set(H, 27);
        expectedS.set(C, 15);
        expectedS.set(N, 5);
        expectedS.set(O, 7);
        expectedS.set(S, 1);

        std::vector<elements::Isotope> cHi;
        elements::Element cH(
            elements::ElementImpl(elements::ElementImpl::getNextId(), "H", 1,
                cHi));

        shouldEqual(ass.getStoichiometry(), expectedS);

        StoichiometryConfigImpl sci("test1");
        sci.insertElement(cH);
        StoichiometryConfig sc(sci);

        expectedS.set(H, 0);
        expectedS.set(cH, 27);

        ass.applyAminoAcidStoichiometryConfig(sc);

        shouldEqual(ass.getStoichiometry(), expectedS);

        ass.applyAminoAcidStoichiometryConfig(
            StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);
        ass.applyAminoAcidStoichiometryConfig("test1");

        shouldEqual(ass.getStoichiometry(), expectedS);

        ass.applyModificationAtPosition("Oxidation", 2);
        expectedS.set(O, 8);

        shouldEqual(ass.getStoichiometry(), expectedS);

        std::vector<elements::Isotope> cOi;
        elements::Element cO(
            elements::ElementImpl(elements::ElementImpl::getNextId(), "O", 8,
                cOi));
        StoichiometryConfigImpl scim("test2");
        scim.insertElement(cO);
        StoichiometryConfig scm(scim);

        expectedS.set(O, 7);
        expectedS.set(cO, 1);

        ass.applyModificationStoichiometryConfig(scm);
        shouldEqual(ass.getStoichiometry(), expectedS);

        ass.applyModificationStoichiometryConfig(
            StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);
        ass.applyModificationStoichiometryConfig("test2");

        shouldEqual(ass.getStoichiometry(), expectedS);

        ass.applyModificationAtPosition(
            modifications::Modification("Acetyl:2H(3)"), 4);
        expectedS.add(H, -1);
        expectedS.add(
            elements::Element(
                elements::ElementImpl::getDefaultKeyForElementSymbol("2H")),
            3);
        expectedS.add(C, 2);
        expectedS.add(O, 1);

        shouldEqual(ass.getStoichiometry(), expectedS);

        std::vector<elements::Isotope> c2Hi;
        elements::Element c2H(
            elements::ElementImpl(elements::ElementImpl::getNextId(), "2H", 1,
                c2Hi));
        StoichiometryConfigImpl scil("test3");
        scil.insertElement(c2H);
        StoichiometryConfig scl(scil);

        expectedS.add(
            elements::Element(
                elements::ElementImpl::getDefaultKeyForElementSymbol("2H")),
            -3);
        expectedS.add(c2H, 3);

        ass.applyIsotopicLabelStoichiometryConfig(scl);

        shouldEqual(ass.getStoichiometry(), expectedS);

        ass.applyIsotopicLabelStoichiometryConfig(
            StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);
        ass.applyIsotopicLabelStoichiometryConfig("test3");

        shouldEqual(ass.getStoichiometry(), expectedS);

        Bool thrown = false;
        try {
            ass.applyModificationAtPosition("Acetyl:2H(3)", 4);
        } catch (RuntimeError& e) {
            thrown = true;
        }
        shouldEqual(thrown, true);
    }

    void testAminoAcidSequenceApplyModifications()
    {
        // test isotopic labels
        // test real n- and c-term mods
        String aass = "AACCGQQSSG";
        AminoAcidSequence aas(aass);

        modifications::RawModificationImpl::RawModificationImplKeyType k =
                "Oxidation";
        modifications::Modification mod(k);
        aas.applyModificationAtPosition(mod, 3);
        shouldEqual(aas[3].isModified(), true);
        shouldEqual(aas[3].hasModification(mod), true);

        bool thrown = false;
        try {
            aas.applyModificationAtPosition(mod, 1);
        } catch (Exception& e) {
            thrown = true;
        }shouldEqual(thrown, true);

        aas.remove(mod);
        aas.applyModificationAtPosition(k, 3);
        shouldEqual(aas[3].isModified(), true);
        shouldEqual(aas[3].hasModification(mod), true);

        aas.remove(mod);
        aas.applyModificationAtPosition(modifications::RawModification(k), 3);
        shouldEqual(aas[3].isModified(), true);
        shouldEqual(aas[3].hasModification(mod), true);

        AminoAcidSequence::ModificationList ml;
        modifications::RawModificationImpl::RawModificationImplKeyType k1 =
                "Phospho";
        modifications::RawModificationImpl::RawModificationImplKeyType k2 =
                "Trimethyl";

        std::vector<
                modifications::RawModificationImpl::RawModificationImplKeyType> mlk;
        mlk.push_back(k);
        mlk.push_back(k1);
        mlk.push_back(k2);

        std::vector<modifications::RawModification> mlrm;
        mlrm.push_back(modifications::RawModification(k));
        mlrm.push_back(modifications::RawModification(k1));
        mlrm.push_back(modifications::RawModification(k2));

        modifications::Modification mod1(k1);
        modifications::Modification mod2(k2);
        ml.push_back(mod);
        ml.push_back(mod1);
        ml.push_back(mod2);
        aas.makeProteinNTerm();

        AminoAcidSequence aas1 = aas, aas2 = aas;

        shouldEqual(aas1, aas);
        shouldEqual(aas2, aas);

        aas.applyFixedModifications(ml);
        aas1.applyFixedModifications(mlk);
        aas2.applyFixedModifications(mlrm);

        // N-term
        shouldEqual(aas[0].isModified(), false);
        // A
        shouldEqual(aas[1].hasModification(mod2), true);
        // A
        shouldEqual(aas[2].isModified(), false);
        // C
        shouldEqual(aas[3].hasModification(mod), true);
        // C
        shouldEqual(aas[4].hasModification(mod), true);
        // G
        shouldEqual(aas[5].isModified(), false);
        // Q
        shouldEqual(aas[6].isModified(), false);
        // Q
        shouldEqual(aas[7].isModified(), false);
        // S
        shouldEqual(aas[8].hasModification(mod1), true);
        // S
        shouldEqual(aas[9].hasModification(mod1), true);
        // G
        shouldEqual(aas[10].hasModification(mod), true);
        // C-term
        shouldEqual(aas[11].isModified(), false);

        shouldEqual(aas, aas1);
        shouldEqual(aas, aas2);
    }

    void testCollection()
    {
        std::vector<String> slv;
        Collection<String> sl;
        shouldEqual(slv.get_allocator() == sl.get_allocator(), true);
        shouldEqual(sl.size(), 0u);
        sl.push_back("ASD");
        shouldEqual(sl.size(), 1u);
        shouldEqual(sl[0], "ASD");
        sl.push_back("SDF");
        sl.pop_back();
        shouldEqual(sl.size(), 1u);
        shouldEqual(sl[0], "ASD");
        const Collection<String> sl_c = sl;
        shouldEqual(sl_c.size(), 1u);
        shouldEqual(sl_c[0], "ASD");
    }

};

/** The main function that runs the tests for class AminoAcidSequence.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    AminoAcidSequenceTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}
