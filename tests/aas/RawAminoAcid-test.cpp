/*
 * RawRawAminoAcid-test.cpp
 *
 * Copyright (c) 2012 Mathias Wilhelm
 * Copyright (c) 2012 Marc Kirchner
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

#include <MSTK/aas/RawAminoAcid.hpp>
#include <MSTK/aas/Stoichiometry.hpp>
#include <MSTK/aas/Element.hpp>
#include <MSTK/common/Error.hpp>

#include "unittest.hxx"

#include <iostream>
#include <algorithm>

using namespace mstk;
using namespace aas::aminoAcids;
using namespace aas::stoichiometries;

/** Short description.
 * Long description.
 */
struct RawAminoAcidTestSuite : vigra::test_suite
{
    /** Constructor.
     * The RawRawAminoAcidTestSuite constructor adds all RawRawAminoAcid tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    RawAminoAcidTestSuite() :
            vigra::test_suite("RawRawAminoAcid")
    {
        add(testCase(&RawAminoAcidTestSuite::testRawAminoAcid));
        add(testCase(&RawAminoAcidTestSuite::testRawAminoAcidFW));
        add(testCase(&RawAminoAcidTestSuite::testStaticParser));
        add(testCase(&RawAminoAcidTestSuite::testRawAminoAcidRef));
        add(testCase(&RawAminoAcidTestSuite::testAddRawAminoAcid));
        add(testCase(&RawAminoAcidTestSuite::testAddRawAminoAcidRef));
        add(
            testCase(&RawAminoAcidTestSuite::testOverrideUninitializedRawAminoAcid));
        add(
            testCase(&RawAminoAcidTestSuite::testOverrideInitializedRawAminoAcid));
        add(testCase(&RawAminoAcidTestSuite::testCreateAminoAcid));
    }

    void testRawAminoAcid()
    {
        aas::elements::Element H(1);
        aas::elements::Element C(6);
        aas::elements::Element N(7);
        aas::elements::Element O(8);

        RawAminoAcidImpl::RawAminoAcidImplKeyType k = 'C';
        Char symbol = 'T';
        String three = "Cys";
        String full = "Cysteine";
        Stoichiometry ts;
        ts.set(H, 6);
        ts.set(C, 2);
        ts.set(O, 1);
        RawAminoAcidImpl aa_c(k, symbol, ts);
        aa_c.setThreeLetterCode(three);
        aa_c.setFullName(full);

        shouldEqual(aa_c.getId(), k);
        shouldEqual(aa_c.getSymbol(), symbol);
        shouldEqual(aa_c.getThreeLetterCode(), three);
        shouldEqual(aa_c.getFullName(), full);
        shouldEqual(aa_c.getStoichiometry(), ts);
        shouldEqual(aa_c.isCTerm(), false);
        shouldEqual(aa_c.isNTerm(), false);
        RawAminoAcidImpl aa_c_c = aa_c;
        shouldEqual(aa_c, aa_c_c);

        shouldEqual(
            RawAminoAcidImpl(RawAminoAcidImpl::PEPTIDE_C_TERM).isCTerm(),
            true);
        shouldEqual(
            RawAminoAcidImpl(RawAminoAcidImpl::PROTEIN_C_TERM).isCTerm(),
            true);
        shouldEqual(
            RawAminoAcidImpl(RawAminoAcidImpl::PEPTIDE_N_TERM).isNTerm(),
            true);
        shouldEqual(
            RawAminoAcidImpl(RawAminoAcidImpl::PROTEIN_N_TERM).isNTerm(),
            true);

        // test a standard amino acid
        RawAminoAcidImpl::RawAminoAcidImplKeyType k1 = 'A';
        RawAminoAcidImpl aa(k1);

        Stoichiometry s;
        s.set(H, 5);
        s.set(C, 3);
        s.set(N, 1);
        s.set(O, 1);

        shouldEqual(aa.getId(), k1);
        shouldEqual(aa.getSymbol(), 'A');
        shouldEqual(aa.getStoichiometry(), s);
        shouldEqual(aa.getThreeLetterCode(), "Ala");
        shouldEqual(aa.getFullName(), "Alanine");
        shouldEqual(aa_c == aa, false);
        shouldEqual(aa_c != aa, true);

        aa.setSymbol('a');
        shouldEqual(aa.getSymbol(), 'a');
        Stoichiometry ns;
        aa.setStoichiometry(ns);
        shouldEqual(aa.getStoichiometry(), ns);

        RawAminoAcidImpl tmp1 = aa;
        tmp1 = aa_c;
        shouldEqual(aa_c, tmp1);
        shouldEqual(aa_c != tmp1, false);

        shouldEqual(aa == tmp1, false);
        shouldEqual(aa != tmp1, true);
    }

    void testStaticParser()
    {
        String toparse[] = { "A", "a", "Ala", "ALA", "Alanine", "ALAninE",
                                  "N-term", "Peptide N-Term", "Protein n-term",
                                  "ProTEin C-TERM", "peptide C-TERM" };
        RawAminoAcidImpl::RawAminoAcidImplKeyType expectedRawAminoAcid[] = {
                'A', 'A', 'A', 'A', 'A', 'A', RawAminoAcidImpl::PEPTIDE_N_TERM,
                RawAminoAcidImpl::PEPTIDE_N_TERM,
                RawAminoAcidImpl::PROTEIN_N_TERM,
                RawAminoAcidImpl::PROTEIN_C_TERM,
                RawAminoAcidImpl::PEPTIDE_C_TERM };
        size_t n = 11;

        for (size_t i = 0; i < n; ++i) {
            shouldEqual(RawAminoAcidImpl::getKeyForAminoAcidString(toparse[i]),
                expectedRawAminoAcid[i]);
        }

        bool thrown = false;
        try {
            RawAminoAcidImpl::getKeyForAminoAcidString("unknown");
        } catch (mstk::LogicError& e) {
            thrown = true;
        }shouldEqual(thrown, true);

        thrown = false;
        try {
            RawAminoAcidImpl::getKeyForAminoAcidString("ttt");
        } catch (mstk::LogicError& e) {
            thrown = true;
        }
        shouldEqual(thrown, true);
    }

    void testRawAminoAcidRef()
    {
        RawAminoAcidImpl::RawAminoAcidImplKeyType k1 = 'L', k2 = 'G';

        // create two different amino acids without using the factory
        RawAminoAcidImpl e_2(k1);
        RawAminoAcidImpl e_3(k2);

        // create the same amino acids again, but using the flyweight implementation
        RawAminoAcid er_2(k1);
        RawAminoAcid er_3(k2);

        // testing equality
        shouldEqual(e_2, er_2.get());
        shouldEqual(e_3, er_3.get());

        // testing whether the id extractor works correctly
        shouldEqual(er_2.get_key(), er_2.get().getId());
        shouldEqual(er_3.get_key(), er_3.get().getId());

        // retrieve the amino acids again
        RawAminoAcid er_t2(k1);
        RawAminoAcid er_t3(k2);

        // test whether the old references and new references point to the same
        // object
        shouldEqual(&er_2.get(), &er_t2.get());
        shouldEqual(&er_3.get(), &er_t3.get());
    }

    // testing static function to add an amino acid
    void testAddRawAminoAcid()
    {
        // setting up test data
        aas::elements::Element H(1);
        aas::elements::Element C(6);
        aas::elements::Element O(8);

        RawAminoAcidImpl::RawAminoAcidImplKeyType k = 'Z';
        Char symbol = 'e';
        Char symbol2 = 't';
        String three = "Zet";
        String full = "Zetet";
        Stoichiometry ts;
        ts.set(H, 6);
        ts.set(C, 2);
        ts.set(O, 1);

        // adding not existing element
        shouldEqual(addRawAminoAcid(k, symbol, three, full, ts), true);
        RawAminoAcidImpl e(k, symbol, ts);
        e.setThreeLetterCode(three);
        e.setFullName(full);
        RawAminoAcid e_r(k);
        shouldEqual(e_r, e);

        // a second try to add an element with the same id should fail
        shouldEqual(addRawAminoAcid(k, symbol2, three, full, ts), false);
        // and the amino acid added in the frist place should stay the same
        shouldEqual(RawAminoAcid(k), e);
        RawAminoAcidImpl e2(k, symbol2, ts);
        shouldEqual(RawAminoAcid(k) == e2, false);
    }

    void testAddRawAminoAcidRef()
    {
        // first test should result in an out_of_range exceptions since the
        // default amino acid table does no contain an amino acid with char 'z'
        bool thrown = false;
        try {
            RawAminoAcid test('z');
        } catch (mstk::LogicError& e) {
            thrown = true;
        }shouldEqual(thrown, true);

        // create an arbitrary amino acid
        RawAminoAcidImpl::RawAminoAcidImplKeyType k1 = 'z';
        Char name = 'z';
        Stoichiometry ts;
        RawAminoAcidImpl t(k1, name, ts);
        // add amino acid to the flyweight table
        RawAminoAcid tr(t);

        // test equality of the original amino acid and the const ref
        shouldEqual(tr, t);
        shouldEqual(t.getId(), tr.get_key());

        // retrieve the amino acid directly from the flyweight table
        RawAminoAcid tr_t(k1);

        // test equality of two const refs of the same amino acids
        shouldEqual(tr, tr_t);
        // test reference pointer to ensure it is the same object
        shouldEqual(&tr.get(), &tr_t.get());

        RawAminoAcidImpl::RawAminoAcidImplKeyType kk = 'Q';
        RawAminoAcid wt;
        wt = kk;
        shouldEqual(RawAminoAcid('Q'), wt);
    }

    void testOverrideUninitializedRawAminoAcid()
    {
        // create different amino acid T
        RawAminoAcidImpl::RawAminoAcidImplKeyType k1 = 'S';
        Char symbol = 'Q';
        Stoichiometry ts;
        RawAminoAcidImpl t(k1, symbol, ts);
        // store it in flyweight table
        RawAminoAcid tr(t);
        // retrieve element id k1
        RawAminoAcid tr_t(k1);

        // getSymbol returns Q instead of T
        shouldEqual(tr_t.get().getSymbol(), symbol);
        shouldEqual(tr_t.get().getSymbol() != RawAminoAcidImpl(k1).getSymbol(),
            true);
    }

    void testOverrideInitializedRawAminoAcid()
    {
        RawAminoAcidImpl::RawAminoAcidImplKeyType k1 = 'A';
        RawAminoAcid tr_1(k1);

        RawAminoAcidImpl t(k1, 'T', Stoichiometry());
        // flyweight checks the id of the amino acid t and recognizes it, since it
        // was retrieved(initialized) earlier
        // as a consequence, the flyweight factory returns the known object
        // instead of overriding it with the given one!
        RawAminoAcid tr_2(t);

        shouldEqual(tr_2.get().getSymbol() != 'T', true);
    }

    void testCreateAminoAcid()
    {
        const Char stoi_chars[] = { 'A', 'C', 'D', 'E', 'F', 'G', 'H',
                                         'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                                         'R', 'S', 'T', 'V', 'W', 'Y',
                                         RawAminoAcidImpl::PEPTIDE_N_TERM,
                                         RawAminoAcidImpl::PEPTIDE_C_TERM,
                                         RawAminoAcidImpl::PROTEIN_N_TERM,
                                         RawAminoAcidImpl::PROTEIN_C_TERM };
        Size n = 24;
        for (Size i = 0; i < n; ++i) {
            RawAminoAcid a(stoi_chars[i]);
            shouldEqual(a.get_key(), stoi_chars[i]);
        }
    }

    void testRawAminoAcidFW()
    {
        RawAminoAcid a1('A');
        RawAminoAcid a2('G');
        RawAminoAcid a3('Q');

        shouldEqual(a1 < a1, false);
        shouldEqual(a1 < a2, true);
        shouldEqual(a1 < a3, true);

        shouldEqual(a1 <= a1, true);
        shouldEqual(a1 <= a2, true);
        shouldEqual(a1 <= a3, true);

        shouldEqual(a1 > a1, false);
        shouldEqual(a1 > a2, false);
        shouldEqual(a1 > a3, false);

        shouldEqual(a1 >= a1, true);
        shouldEqual(a1 >= a2, false);
        shouldEqual(a1 >= a3, false);
    }
};

/** The main function that runs the tests for class RawRawAminoAcid.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    RawAminoAcidTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

