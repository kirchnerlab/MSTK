/*
 * Digester-test.cpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
 *
 */

#include <MSTK/aas/tools/Digester.hpp>

#include "unittest.hxx"

#include <iostream>

using namespace mstk;
using namespace aas;
using namespace aas::tools;

/** Short description.
 * Long description.
 */
struct DigesterTestSuite: vigra::test_suite {
	/** Constructor.
	 * The DigesterTestSuite constructor adds all Digester tests to
	 * the test suite. If you write an additional test, add the test
	 * case here.
	 */
	DigesterTestSuite() :
			vigra::test_suite("Digester") {
		add(testCase(&DigesterTestSuite::testDigester));
		add(testCase(&DigesterTestSuite::testMissedCleavages));
	}

	void testDigester() {
		String regex = "(R|K)([^P])";
		tools::Digester d(regex);
		AminoAcidSequence aas("AAARCCCKDDDRPEEERKFFF");
		Digester::AminoAcidSequences frags;
		d.operator ()(aas, frags, 0);
		shouldEqual(frags.size(), 5);
		shouldEqual(frags[0].toUnmodifiedSequenceString(), "AAAR");
		shouldEqual(frags[1].toUnmodifiedSequenceString(), "CCCK");
		shouldEqual(frags[2].toUnmodifiedSequenceString(), "DDDRPEEER");
		shouldEqual(frags[3].toUnmodifiedSequenceString(), "FFF");
		shouldEqual(frags[4].toUnmodifiedSequenceString(), "K");
	}

	void testMissedCleavages() {
		String regex = "(R|K)([^P])";
		tools::Digester d(regex);
		AminoAcidSequence aas("AAARCCCKDDDRPEEERKFFF");
		Digester::AminoAcidSequences frags;
		d.operator ()(aas, frags, 1);
		shouldEqual(frags.size(), 9);
		d.operator ()(aas, frags, 2);
		shouldEqual(frags.size(), 12);
	}

};

/** The main function that runs the tests for class Digester.
 * Under normal circumstances you need not edit this.
 */
int main() {
	DigesterTestSuite test;
	int success = test.run();
	std::cout << test.report() << std::endl;
	return success;
}

