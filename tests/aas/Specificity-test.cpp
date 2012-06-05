/*
 * Specificity-test.cpp
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

#include <MSTK/aas/Specificity.hpp>
#include <MSTK/common/Error.hpp>

#include "unittest.hxx"

#include <iostream>
#include <algorithm>

using namespace mstk;
using namespace aas::modifications;
using namespace aas::stoichiometries;

/** Short description.
 * Long description.
 */
struct SpecificityTestSuite : vigra::test_suite
{
    /** Constructor.
     * The SpecificityTestSuite constructor adds all Specificity tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    SpecificityTestSuite() :
            vigra::test_suite("Specificity")
    {
        add(testCase(&SpecificityTestSuite::testSpecificity));
        add(testCase(&SpecificityTestSuite::testStaticParser));
        add(testCase(&SpecificityTestSuite::testSpecificityApplicable));
    }

    // testing getter/setter
    void testSpecificity()
    {

        aas::aminoAcids::RawAminoAcid aa('A');
        Specificity::Position pos = Specificity::ANYWHERE;
        Specificity::Classification clas = Specificity::ARTEFACT;
        String comment = "Comment a";
        Specificity spec(aa, pos, clas);
        spec.setComment(comment);

        shouldEqual(spec.getSite(), aa);
        shouldEqual(spec.getClassification(), clas);
        shouldEqual(spec.getPosition(), pos);
        shouldEqual(spec.getComment(), comment);

        aas::aminoAcids::RawAminoAcid aan('C');
        Specificity::Position posn = Specificity::ANY_C_TERM;
        Specificity::Classification clasn = Specificity::POST_TRANSLATIONAL;
        spec.setSite(aan);
        spec.setPosition(posn);
        spec.setClassification(clasn);
        shouldEqual(spec.getSite(), aan);
        shouldEqual(spec.getPosition(), posn);
        shouldEqual(spec.getClassification(), clasn);

        aas::elements::Element H(1);
        aas::elements::Element C(6);
        aas::elements::Element N(7);
        aas::elements::Element O(8);
        aas::elements::Element S(16);

        Stoichiometry st1;
        st1.set(H, 1);
        st1.set(C, 2);
        st1.set(N, 3);

        Stoichiometry st2;
        st2.set(H, 1);
        st2.set(O, 2);
        st2.set(N, 3);

        Stoichiometry st3;
        st3.set(H, 1);
        st3.set(C, 2);
        st3.set(S, 3);

        Stoichiometry st4;
        st4.set(H, 1);
        st4.set(O, 2);
        st4.set(S, 3);

        std::vector<Stoichiometry> emptyLoss;
        std::vector<Stoichiometry> neutralLoss;
        neutralLoss.push_back(st1);
        neutralLoss.push_back(st2);
        std::vector<Stoichiometry> pepNeutralLoss;
        pepNeutralLoss.push_back(st3);
        pepNeutralLoss.push_back(st4);

        spec.addNeutralLoss(st1);
        spec.addNeutralLoss(st2);
        shouldEqual(spec.getNeutralLosses(), neutralLoss);

        spec.addPepNeutralLoss(st3);
        spec.addPepNeutralLoss(st4);
        shouldEqual(spec.getPepNeutralLosses(), pepNeutralLoss);

        spec.setNeutralLosses(neutralLoss);
        shouldEqual(spec.getNeutralLosses(), neutralLoss);

        spec.setPepNeutralLosses(pepNeutralLoss);
        shouldEqual(spec.getPepNeutralLosses(), pepNeutralLoss);

        spec.clearNeutralLosses();
        shouldEqual(spec.getNeutralLosses(), emptyLoss);

        spec.clearPepNeutralLosses();
        shouldEqual(spec.getPepNeutralLosses(), emptyLoss);

        Specificity spec1("C", "Any c-term", "post-translational");
        spec1.setComment(comment);
        shouldEqual(spec1, spec);

        Specificity spec2("G", "Any N-term", "artefact");
        Specificity& spec3 = spec1;
        spec3 = spec1;
        spec2 = spec1;
        shouldEqual(spec1, spec2);
        shouldEqual(spec1, spec3);
        shouldEqual(spec1 != spec2, false);
    }

    // testing static parser for position and classification
    void testStaticParser()
    {
        // both should fail in case of an invalid position/classification
        Bool thrown = false;
        try {
            Specificity::parsePositionString("ASD");
        } catch (mstk::LogicError& e) {
            thrown = true;
        }

        shouldEqual(thrown, true);

        thrown = false;
        try {
            Specificity::parseClassificationString("ASD");
        } catch (mstk::LogicError& e) {
            thrown = true;
        }

        shouldEqual(thrown, true);

        // testing all currently known positions + correct enum value
        String positions[] = { "Any N-term", "Any C-term",
                                    "Protein N-term", "Protein C-term",
                                    "Anywhere" };

        thrown = false;
        try {
            for (Size i = 0; i < 5; ++i) {
                shouldEqual(
                    (Size) Specificity::parsePositionString(positions[i]),
                    i);
                std::transform(positions[i].begin(), positions[i].end(),
                    positions[i].begin(), ::toupper);
                shouldEqual(
                    (Size) Specificity::parsePositionString(positions[i]),
                    i);
            }
        } catch (std::out_of_range& e) {
            thrown = true;
        }

        shouldEqual(thrown, false);

        // testing all currently known classifications + correct enum value
        String classifications[] = { "-", "Post-translational",
                                          "Co-translational",
                                          "Pre-translational",
                                          "Chemical derivative", "Artefact",
                                          "N-linked glycosylation",
                                          "O-linked glycosylation",
                                          "Other glycosylation",
                                          "Synth. pep. protect. gp.",
                                          "Isotopic label",
                                          "Non-standard residue", "Multiple",
                                          "Other" };

        thrown = false;
        try {
            for (Size i = 0; i < 13; ++i) {
                shouldEqual(
                    (Size) Specificity::parseClassificationString(classifications[i]),
                    i);
                std::transform(classifications[i].begin(),
                    classifications[i].end(), classifications[i].begin(),
                    ::toupper);
                shouldEqual(
                    (Size) Specificity::parseClassificationString(classifications[i]),
                    i);
            }
        } catch (std::out_of_range& e) {
            thrown = true;
        }

        shouldEqual(thrown, false);
    }

    void testSpecificityApplicable()
    {
        // testing all possible variations of specificities applied to all possible variations of amino acids
        aas::aminoAcids::RawAminoAcid c('C');
        aas::aminoAcids::RawAminoAcid a('A');
        aas::aminoAcids::RawAminoAcid e('\0');
        aas::aminoAcids::RawAminoAcid pep_n(
            aas::aminoAcids::RawAminoAcidImpl::PEPTIDE_N_TERM);
        aas::aminoAcids::RawAminoAcid prot_n(
            aas::aminoAcids::RawAminoAcidImpl::PROTEIN_N_TERM);
        aas::aminoAcids::RawAminoAcid pep_c(
            aas::aminoAcids::RawAminoAcidImpl::PEPTIDE_C_TERM);
        aas::aminoAcids::RawAminoAcid prot_c(
            aas::aminoAcids::RawAminoAcidImpl::PROTEIN_C_TERM);

        Specificity spec_c_any(c, Specificity::ANYWHERE,
            Specificity::ARTEFACT);

        Specificity spec_c_any_n(c, Specificity::ANY_N_TERM,
            Specificity::ARTEFACT);
        Specificity spec_c_prot_n(c, Specificity::PROTEIN_N_TERM,
            Specificity::ARTEFACT);
        Specificity spec_pepnterm_any_n(pep_n, Specificity::ANY_N_TERM,
            Specificity::ARTEFACT);
        Specificity spec_pepnterm_prot_n(pep_n, Specificity::PROTEIN_N_TERM,
            Specificity::ARTEFACT);
        Specificity spec_pepcterm_any_n(pep_c, Specificity::ANY_N_TERM,
            Specificity::ARTEFACT);
        Specificity spec_pepcterm_prot_n(pep_c, Specificity::PROTEIN_N_TERM,
            Specificity::ARTEFACT);

        Specificity spec_c_any_c(c, Specificity::ANY_C_TERM,
            Specificity::ARTEFACT);
        Specificity spec_c_prot_c(c, Specificity::PROTEIN_C_TERM,
            Specificity::ARTEFACT);
        Specificity spec_pepcterm_any_c(pep_c, Specificity::ANY_C_TERM,
            Specificity::ARTEFACT);
        Specificity spec_pepcterm_prot_c(pep_c, Specificity::PROTEIN_C_TERM,
            Specificity::ARTEFACT);
        Specificity spec_pepnterm_any_c(pep_n, Specificity::ANY_C_TERM,
            Specificity::ARTEFACT);
        Specificity spec_pepnterm_prot_c(pep_n, Specificity::PROTEIN_C_TERM,
            Specificity::ARTEFACT);

        shouldEqual(spec_c_any.isApplicable(c, a, c), false);
        shouldEqual(spec_c_any.isApplicable(c, c, c), true);
        shouldEqual(spec_c_any.isApplicable(pep_n, c, pep_c), true);
        shouldEqual(spec_c_any.isApplicable(pep_n, c, c), true);
        shouldEqual(spec_c_any.isApplicable(c, c, pep_c), true);
        shouldEqual(spec_c_any.isApplicable(prot_n, c, prot_c), true);
        shouldEqual(spec_c_any.isApplicable(prot_n, c, c), true);
        shouldEqual(spec_c_any.isApplicable(c, c, prot_c), true);
        shouldEqual(spec_c_any.isApplicable(prot_n, c, pep_c), true);
        shouldEqual(spec_c_any.isApplicable(pep_n, c, prot_c), true);
        shouldEqual(spec_c_any.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_c_any.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_c_any.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_c_any.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_c_any_n.isApplicable(c, a, c), false);
        shouldEqual(spec_c_any_n.isApplicable(c, c, c), false);
        shouldEqual(spec_c_any_n.isApplicable(pep_n, c, pep_c), true);
        shouldEqual(spec_c_any_n.isApplicable(pep_n, c, c), true);
        shouldEqual(spec_c_any_n.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_c_any_n.isApplicable(prot_n, c, prot_c), true);
        shouldEqual(spec_c_any_n.isApplicable(prot_n, c, c), true);
        shouldEqual(spec_c_any_n.isApplicable(c, c, prot_c), false);
        shouldEqual(spec_c_any_n.isApplicable(prot_n, c, pep_c), true);
        shouldEqual(spec_c_any_n.isApplicable(pep_n, c, prot_c), true);
        shouldEqual(spec_c_any_n.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_c_any_n.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_c_any_n.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_c_any_n.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_c_prot_n.isApplicable(c, a, c), false);
        shouldEqual(spec_c_prot_n.isApplicable(c, c, c), false);
        shouldEqual(spec_c_prot_n.isApplicable(pep_n, c, pep_c), false);
        shouldEqual(spec_c_prot_n.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_c_prot_n.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_c_prot_n.isApplicable(prot_n, c, prot_c), true);
        shouldEqual(spec_c_prot_n.isApplicable(prot_n, c, c), true);
        shouldEqual(spec_c_prot_n.isApplicable(c, c, prot_c), false);
        shouldEqual(spec_c_prot_n.isApplicable(prot_n, c, pep_c), true);
        shouldEqual(spec_c_prot_n.isApplicable(pep_n, c, prot_c), false);
        shouldEqual(spec_c_prot_n.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_c_prot_n.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_c_prot_n.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_c_prot_n.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_pepnterm_any_n.isApplicable(c, a, c), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(c, c, c), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(pep_n, c, pep_c), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(prot_n, c, prot_c),
            false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(prot_n, c, c), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(c, c, prot_c), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(prot_n, c, pep_c), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(pep_n, c, prot_c), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_pepnterm_any_n.isApplicable(e, prot_n, c), true);
        shouldEqual(spec_pepnterm_any_n.isApplicable(e, pep_n, c), true);

        shouldEqual(spec_pepnterm_prot_n.isApplicable(c, a, c), false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(c, c, c), false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(pep_n, c, pep_c), false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(prot_n, c, prot_c),
            false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(prot_n, c, c), false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(c, c, prot_c), false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(prot_n, c, pep_c),
            false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(pep_n, c, prot_c),
            false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_pepnterm_prot_n.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_pepcterm_any_n.isApplicable(c, a, c), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(c, c, c), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(pep_n, c, pep_c), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(prot_n, c, prot_c),
            false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(prot_n, c, c), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(c, c, prot_c), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(prot_n, c, pep_c), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(pep_n, c, prot_c), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_pepcterm_any_n.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_pepcterm_prot_n.isApplicable(c, a, c), false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(c, c, c), false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(pep_n, c, pep_c), false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(prot_n, c, prot_c),
            false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(prot_n, c, c), false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(c, c, prot_c), false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(prot_n, c, pep_c),
            false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(pep_n, c, prot_c),
            false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_pepcterm_prot_n.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_c_any_c.isApplicable(c, a, c), false);
        shouldEqual(spec_c_any_c.isApplicable(c, c, c), false);
        shouldEqual(spec_c_any_c.isApplicable(pep_n, c, pep_c), true);
        shouldEqual(spec_c_any_c.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_c_any_c.isApplicable(c, c, pep_c), true);
        shouldEqual(spec_c_any_c.isApplicable(prot_n, c, prot_c), true);
        shouldEqual(spec_c_any_c.isApplicable(prot_n, c, c), false);
        shouldEqual(spec_c_any_c.isApplicable(c, c, prot_c), true);
        shouldEqual(spec_c_any_c.isApplicable(prot_n, c, pep_c), true);
        shouldEqual(spec_c_any_c.isApplicable(pep_n, c, prot_c), true);
        shouldEqual(spec_c_any_c.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_c_any_c.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_c_any_c.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_c_any_c.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_c_prot_c.isApplicable(c, a, c), false);
        shouldEqual(spec_c_prot_c.isApplicable(c, c, c), false);
        shouldEqual(spec_c_prot_c.isApplicable(pep_n, c, pep_c), false);
        shouldEqual(spec_c_prot_c.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_c_prot_c.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_c_prot_c.isApplicable(prot_n, c, prot_c), true);
        shouldEqual(spec_c_prot_c.isApplicable(prot_n, c, c), false);
        shouldEqual(spec_c_prot_c.isApplicable(c, c, prot_c), true);
        shouldEqual(spec_c_prot_c.isApplicable(prot_n, c, pep_c), false);
        shouldEqual(spec_c_prot_c.isApplicable(pep_n, c, prot_c), true);
        shouldEqual(spec_c_prot_c.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_c_prot_c.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_c_prot_c.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_c_prot_c.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_pepcterm_any_c.isApplicable(c, a, c), false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(c, c, c), false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(pep_n, c, pep_c), false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(prot_n, c, prot_c),
            false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(prot_n, c, c), false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(c, c, prot_c), false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(prot_n, c, pep_c), false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(pep_n, c, prot_c), false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(c, prot_c, e), true);
        shouldEqual(spec_pepcterm_any_c.isApplicable(c, pep_c, e), true);
        shouldEqual(spec_pepcterm_any_c.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_pepcterm_any_c.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_pepcterm_prot_c.isApplicable(c, a, c), false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(c, c, c), false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(pep_n, c, pep_c), false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(prot_n, c, prot_c),
            false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(prot_n, c, c), false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(c, c, prot_c), false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(prot_n, c, pep_c),
            false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(pep_n, c, prot_c),
            false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_pepcterm_prot_c.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_pepnterm_any_c.isApplicable(c, a, c), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(c, c, c), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(pep_n, c, pep_c), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(prot_n, c, prot_c),
            false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(prot_n, c, c), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(c, c, prot_c), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(prot_n, c, pep_c), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(pep_n, c, prot_c), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_pepnterm_any_c.isApplicable(e, pep_n, c), false);

        shouldEqual(spec_pepnterm_prot_c.isApplicable(c, a, c), false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(c, c, c), false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(pep_n, c, pep_c), false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(pep_n, c, c), false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(c, c, pep_c), false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(prot_n, c, prot_c),
            false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(prot_n, c, c), false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(c, c, prot_c), false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(prot_n, c, pep_c),
            false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(pep_n, c, prot_c),
            false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(c, prot_c, e), false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(c, pep_c, e), false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(e, prot_n, c), false);
        shouldEqual(spec_pepnterm_prot_c.isApplicable(e, pep_n, c), false);
    }

};

/** The main function that runs the tests for class Specificity.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    SpecificityTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

