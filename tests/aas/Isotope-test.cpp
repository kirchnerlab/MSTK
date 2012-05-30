/*
 * Isotope-test.cpp
 *
 * Copyright (c) 2011 Mathias Wilhelm
 * Copyright (c) 2011 Marc Kirchner
 *
 */

#include <MSTK/aas/Isotope.hpp>

#include "unittest.hxx"

#include <iostream>

using namespace mstk;
using namespace aas::elements;

/** Short description.
 * Long description.
 */
struct IsotopeTestSuite : vigra::test_suite
{
    /** Constructor.
     * The IsotopeTestSuite constructor adds all Isotope tests to
     * the test suite. If you write an additional test, add the test
     * case here.
     */
    IsotopeTestSuite() :
            vigra::test_suite("Isotope")
    {
        add(testCase(&IsotopeTestSuite::testIsotope));
    }

    void testIsotope()
    {
        double mass = 101.1;
        double frequency = 0.99;
        Isotope i1(mass, frequency);
        Isotope i2(12.4, 0.32);
        i2 = i1;
        Isotope i3(mass + 1, frequency);

        shouldEqual(i1.getMass(), mass);
        shouldEqual(i1.getFrequency(), frequency);
        shouldEqual(i1, i2);
        shouldEqual(i1 == i3, false);
        shouldEqual(i1 != i3, true);
    }

};

/** The main function that runs the tests for class Isotope.
 * Under normal circumstances you need not edit this.
 */
int main()
{
    IsotopeTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}

