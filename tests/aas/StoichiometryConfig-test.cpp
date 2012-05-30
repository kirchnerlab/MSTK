/*
 * StoichiometryConfig-test.cpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
 *
 */

#include <MSTK/aas/StoichiometryConfig.hpp>

#include "unittest.hxx"

#include <iostream>

using namespace mstk;
using namespace aas;
using namespace aas::stoichiometries;

/** Short description.
 * Long description.
 */
struct StoichiometryConfigTestSuite : vigra::test_suite
{
    /** Constructor.
     * The StoichiometryConfigTestSuite constructor adds all StoichiometryConfig tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    StoichiometryConfigTestSuite() :
            vigra::test_suite("StoichiometryConfig")
    {
        add(testCase(&StoichiometryConfigTestSuite::testStoichiometryConfig));
        add(
            testCase(&StoichiometryConfigTestSuite::testAddStoichiometryConfig));
        add(
            testCase(&StoichiometryConfigTestSuite::testStoichiometryConfigFW));
    }

    void testStoichiometryConfig()
    {
        StoichiometryConfigImpl::StoichiometryConfigImplKeyType sck = "test";
        StoichiometryConfigImpl sc(sck);

        elements::Element H(1);
        elements::Element C(6);
        elements::Element N(7);
        elements::Element O(8);
        elements::Element S(16);

        sc.insertElement(H);
        sc.insertElement(C);
        sc.insertElement(N.get().getSymbol(), N.get_key());
        sc.insertElement(O.get().getSymbol(), O.get_key());
        sc.insertElement(S.get().getSymbol(), S.get_key());

        shouldEqual(sc.getId(), sck);
        shouldEqual(sc.getKeyForSymbol(H.get().getSymbol()), H.get_key());
        shouldEqual(sc.getKeyForSymbol(C.get().getSymbol()), C.get_key());
        shouldEqual(sc.getKeyForSymbol(N.get().getSymbol()), N.get_key());
        shouldEqual(sc.getKeyForSymbol(O.get().getSymbol()), O.get_key());
        shouldEqual(sc.getKeyForSymbol(S.get().getSymbol()), S.get_key());

        StoichiometryConfigImpl::DataType map = sc.getMapping();
        sc.setMapping(map);

        shouldEqual(sc.getKeyForSymbol(H.get().getSymbol()), H.get_key());
        shouldEqual(sc.getKeyForSymbol(C.get().getSymbol()), C.get_key());
        shouldEqual(sc.getKeyForSymbol(N.get().getSymbol()), N.get_key());
        shouldEqual(sc.getKeyForSymbol(O.get().getSymbol()), O.get_key());
        shouldEqual(sc.getKeyForSymbol(S.get().getSymbol()), S.get_key());

        typedef StoichiometryConfigImpl::const_iterator CIT;
        for (CIT it = sc.begin(); it != sc.end(); ++it) {

        }

        typedef StoichiometryConfigImpl::iterator IT;
        for (IT it = sc.begin(); it != sc.end(); ++it) {

        }

        StoichiometryConfigImpl sc1("asd");
        shouldEqual(sc1 != sc, true);
        shouldEqual(sc1 == sc, false);
        sc1 = sc;
        shouldEqual(sc, sc1);
        shouldEqual(sc1 != sc, false);

        StoichiometryConfigImpl sc2 = sc.clone("test1");
        shouldEqual(sc1.getMapping() == sc2.getMapping(), true);
        shouldEqual(sc2.getId(), "test1");

        StoichiometryConfigImpl& sc3 = sc;
        shouldEqual(sc3, sc);
    }

    void testAddStoichiometryConfig()
    {
        StoichiometryConfigImpl::StoichiometryConfigImplKeyType k = "other";
        elements::Element H(1);
        elements::Element C(6);
        elements::Element S(16);
        StoichiometryConfigImpl::DataType map;
        map.insert(
            StoichiometryConfigImpl::EntryType(H.get().getSymbol(), 10));
        map.insert(
            StoichiometryConfigImpl::EntryType(C.get().getSymbol(), 20));
        map.insert(StoichiometryConfigImpl::EntryType(S.get().getSymbol(), 5));

        StoichiometryConfigImpl sci(k);
        sci.setMapping(map);

        shouldEqual(addStoichiometryConfig(k, map), true);

        shouldEqual(StoichiometryConfig(k), sci);

        shouldEqual(addStoichiometryConfig(sci), true);
    }

    void testStoichiometryConfigFW()
    {
        StoichiometryConfigImpl::StoichiometryConfigImplKeyType sck1 = "test";
        StoichiometryConfigImpl::StoichiometryConfigImplKeyType sck2 = "stte";
        StoichiometryConfigImpl::StoichiometryConfigImplKeyType sck3 = "uvw";
        {
            StoichiometryConfigImpl sc1(sck1);
            StoichiometryConfigImpl sc2(sck1);
            StoichiometryConfigImpl sc3(sck1);

            elements::Element H(1);
            elements::Element C(6);
            elements::Element S(16);

            sc1.insertElement(H);
            sc2.insertElement(C);
            sc3.insertElement(S);

            addStoichiometryConfig(sc1);
            addStoichiometryConfig(sc2);
            addStoichiometryConfig(sc3);
        }

        StoichiometryConfig sc1(sck1);
        StoichiometryConfig sc2(sck2);
        StoichiometryConfig sc3(sck3);

        shouldEqual(sc1 < sc1, false);
        shouldEqual(sc1 < sc2, false);
        shouldEqual(sc1 < sc3, true);

        shouldEqual(sc1 <= sc1, true);
        shouldEqual(sc1 <= sc2, false);
        shouldEqual(sc1 <= sc3, true);

        shouldEqual(sc1 > sc1, false);
        shouldEqual(sc1 > sc2, true);
        shouldEqual(sc1 > sc3, false);

        shouldEqual(sc1 >= sc1, true);
        shouldEqual(sc1 >= sc2, true);
        shouldEqual(sc1 > sc3, false);

    }

};

/** The main function that runs the tests for class StoichiometryConfig.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    StoichiometryConfigTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

