/*
 * Stoichiometry-test.cpp
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

#include <MSTK/aas/Stoichiometry.hpp>
#include <MSTK/aas/Element.hpp>
#include <MSTK/common/Error.hpp>

#include "unittest.hxx"

#include <iostream>

using namespace mstk;
using namespace aas;
using namespace aas::stoichiometries;

/** Short description.
 * Long description.
 */
struct StoichiometryTestSuite : vigra::test_suite
{
    /** Constructor.
     * The StoichiometryTestSuite constructor adds all Stoichiometry tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    StoichiometryTestSuite() :
            vigra::test_suite("Stoichiometry")
    {
        add(testCase(&StoichiometryTestSuite::testStoichiometry));
        add(testCase(&StoichiometryTestSuite::testStoichiometryArithmetic));
        add(testCase(&StoichiometryTestSuite::testApplyStoichiometryConfig));
    }

    void testStoichiometry()
    {
        elements::Element H(1);
        elements::Element C(6);
        elements::Element N(7);
        elements::Element O(8);
        elements::Element S(16);

        Stoichiometry tmp;

        // testing size/empty/zero/nonNegative
        Size expectedSize = 0;
        shouldEqual(tmp.size(), expectedSize);
        shouldEqual(tmp.empty(), true);
        shouldEqual(tmp.nonNegative(), true);

        // testing set/get annotation id
        Int id = 123;
        tmp.setAnnotationId(id);
        shouldEqual(id, tmp.getAnnotationId());

        // testing set/add element count
        tmp.set(H, -3);
        tmp.add(H, 2);
        tmp.add(C, 2);
        tmp.set(N, 3);
        tmp.add(N, -3);
        tmp.add(O, 0);
        shouldEqual(tmp.get(H), -1);
        shouldEqual(tmp.get(C), 2);
        shouldEqual(tmp.get(N), 0);
        shouldEqual(tmp.get(O), 0);
        shouldEqual(tmp.get(S), 0);
        const Stoichiometry& tmpC = tmp;
        tmp = tmpC;
        shouldEqual(tmp, tmpC);

        // testing size/empty/zero/nonNegative
        expectedSize = 2;
        shouldEqual(tmp.size(), expectedSize);
        shouldEqual(tmp.empty(), false);
        shouldEqual(tmp.nonNegative(), false);

        // testing iterator
        // ?

        String expectedFormula = "H(-1)C(2)";
        String formula = tmp.toString();
        shouldEqual(expectedFormula, formula);

        // testing clear
        tmp.clear();
        shouldEqual(tmp.get(H), 0);
        shouldEqual(tmp.get(C), 0);
        shouldEqual(tmp.get(N), 0);
        shouldEqual(tmp.get(O), 0);
        shouldEqual(tmp.get(S), 0);

        Stoichiometry tmp1;
        elements::ElementImpl elem(aas::elements::ElementImpl::getNextId(),
            "Zz", 1002);
        aas::elements::addElement(elem);
        tmp1.set(elements::Element(elem.getId()), 5);
        Bool thrown = false;
        try {
            tmp1.applyStoichiometryConfiguration(
                StoichiometryConfig(
                    StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG));
        } catch (mstk::RuntimeError& e) {
            thrown = true;
        }
        shouldEqual(thrown, true);
    }

    void testStoichiometryArithmetic()
    {
        elements::Element H(1);
        elements::Element C(6);
        elements::Element N(7);
        elements::Element O(8);
        elements::Element S(16);

        Stoichiometry s1;
        s1.set(H, 2);
        s1.set(C, 4);
        s1.set(N, 8);
        Stoichiometry s2;
        s2.set(H, 16);
        s2.set(C, 8);
        s2.set(N, 4);
        s2.set(O, 2);

        Stoichiometry s1_minus_s2;
        s1_minus_s2.set(H, -14);
        s1_minus_s2.set(C, -4);
        s1_minus_s2.set(N, 4);
        s1_minus_s2.set(O, -2);

        Stoichiometry s2_minus_s1;
        s2_minus_s1.set(H, 14);
        s2_minus_s1.set(C, 4);
        s2_minus_s1.set(N, -4);
        s2_minus_s1.set(O, 2);

        Stoichiometry s1_plus_s2;
        s1_plus_s2.set(H, 18);
        s1_plus_s2.set(C, 12);
        s1_plus_s2.set(N, 12);
        s1_plus_s2.set(O, 2);

        Stoichiometry s2_plus_s1;
        s2_plus_s1.set(H, 18);
        s2_plus_s1.set(C, 12);
        s2_plus_s1.set(N, 12);
        s2_plus_s1.set(O, 2);

        Stoichiometry s1_si = s1;
        Stoichiometry s2_si = s2;
        {
            Stoichiometry tmp1 = s1 - s2;
            Stoichiometry tmp2 = s1;
            tmp2 -= s2;

            shouldEqual(tmp1, s1_minus_s2);
            shouldEqual(tmp1, tmp2);
            shouldEqual(tmp2 != s1, true);
            shouldEqual(s1, s1_si);
            shouldEqual(s2, s2_si);
        }

        {
            Stoichiometry tmp1 = s2 - s1;
            Stoichiometry tmp2 = s2;
            tmp2 -= s1;

            shouldEqual(tmp1, s2_minus_s1);
            shouldEqual(tmp1, tmp2);
            shouldEqual(tmp2 != s2, true);
            shouldEqual(s1, s1_si);
            shouldEqual(s2, s2_si);
        }

        {
            Stoichiometry tmp1 = s1 + s2;
            Stoichiometry tmp2 = s1;
            tmp2 += s2;

            shouldEqual(tmp1, s1_plus_s2);
            shouldEqual(tmp1, tmp2);
            shouldEqual(s1, s1_si);
            shouldEqual(s2, s2_si);
            shouldEqual(tmp2 != s1, true);
        }

        {
            Stoichiometry tmp1 = s2 + s1;
            Stoichiometry tmp2 = s2;
            tmp2 += s1;

            shouldEqual(tmp1, s2_plus_s1);
            shouldEqual(tmp1, tmp2);
            shouldEqual(s1, s1_si);
            shouldEqual(s2, s2_si);
            shouldEqual(tmp2 != s2, true);
        }

    }

    void testApplyStoichiometryConfig()
    {
        elements::Element H(1);
        elements::Element C(6);
        elements::Element N(7);
        elements::Element O(8);

        std::vector<elements::Isotope> cCi, cNi;
        elements::Element cC(
            elements::ElementImpl(elements::ElementImpl::getNextId(), "C", 13,
                cCi));
        elements::Element cN(
            elements::ElementImpl(elements::ElementImpl::getNextId(), "N", 14,
                cNi));

        Stoichiometry s;
        s.set(H, 10);
        s.set(C, 15);
        s.set(N, 20);
        s.set(O, 25);

        StoichiometryConfigImpl::StoichiometryConfigImplKeyType sck = "test";
        StoichiometryConfigImpl sci(sck);
        sci.insertElement(cC);
        sci.insertElement(cN);

        StoichiometryConfig sc(sci);

        Stoichiometry ns = s.recalculatesWithConfiguration(sc);
        s.applyStoichiometryConfiguration(sc);

        Stoichiometry ex_s;
        ex_s.set(H, 10);
        ex_s.set(cC, 15);
        ex_s.set(cN, 20);
        ex_s.set(O, 25);

        shouldEqual(s, ex_s);
        shouldEqual(ns, ex_s);
    }

};

/** The main function that runs the tests for class Stoichiometry.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    StoichiometryTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

