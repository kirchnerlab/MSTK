/*
 * Combination-test.hpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
 *
 */

#include "MSTK/aas/adapter/LibIPACA.hpp"

#include "MSTK/ipaca/Mercury7.hpp"
#include "MSTK/ipaca/Mercury7Impl.hpp"

#include "unittest.hxx"

#include <iostream>

using namespace mstk;

/** Test suite for the Mercury7 interface.
 */

typedef ipaca::detail::Stoichiometry MyStoichiometry;

//
// ipaca configuration starts here
//
struct SpectrumConverter {
	void operator()(const ipaca::detail::Spectrum& lhs,
			ipaca::detail::Spectrum& rhs) {
		rhs = lhs;
	}
};

struct StoichiometryConverter {
	void operator()(const MyStoichiometry& lhs,
			ipaca::detail::Stoichiometry& rhs) {
		rhs = lhs;
	}
};

namespace mstk {

namespace ipaca {

template<>
struct Traits<MyStoichiometry, ipaca::detail::Spectrum> {
	typedef SpectrumConverter spectrum_converter;
	typedef StoichiometryConverter stoichiometry_converter;
	static detail::Element getHydrogens(const Size n);
	static Bool isHydrogen(const detail::Element&);
	static Double getElectronMass();
};

detail::Element Traits<MyStoichiometry, ipaca::detail::Spectrum>::getHydrogens(
		const Size n) {
	return ipaca::detail::getHydrogens(n);
}

Bool Traits<MyStoichiometry, ipaca::detail::Spectrum>::isHydrogen(
		const detail::Element& e) {
	return ipaca::detail::isHydrogen(e);
}

Double Traits<MyStoichiometry, ipaca::detail::Spectrum>::getElectronMass() {
	return ipaca::detail::getElectronMass();
}

} // namespace ipaca

} // namespace mstk

using namespace aas;
using namespace ipaca;

/** <+Short description of the test suite+>
 * <+Longer description of the test suite+> 
 */
struct CombinationTestSuite: vigra::test_suite {
	/** Constructor.
	 * The CombinationTestSuite constructor adds all Combination tests to
	 * the test suite. If you write an additional test, add the test
	 * case here.
	 */
	CombinationTestSuite() :
			vigra::test_suite("Combination") {
		add(testCase(&CombinationTestSuite::test));
	}

	MyStoichiometry createIntegerH2O() {
		detail::Stoichiometry h2o;
		detail::Isotope i;
		detail::Element h;
		double massesH[] = { 1.007825, 2.01410178 };
		double freqsH[] = { 0.99985, 0.00015 };
		for (size_t k = 0; k < 2; ++k) {
			i.mz = massesH[k];
			i.ab = freqsH[k];
			h.isotopes.push_back(i);
		}
		h.count = 2.0;
		detail::Element o;
		double massesO[] = { 15.9949, 16.9991, 17.9992 };
		double freqsO[] = { 0.99759, 0.000374, 0.002036 };
		for (size_t k = 0; k < 3; ++k) {
			i.mz = massesO[k];
			i.ab = freqsO[k];
			o.isotopes.push_back(i);
		}
		o.count = 1.0;
		h2o.push_back(h);
		h2o.push_back(o);
		return h2o;
	}

	void test() {
		typedef ipaca::Mercury7<aas::stoichiometries::Stoichiometry,
				ipaca::detail::Spectrum> MyMercury7;
		MyMercury7 m;

		aas::stoichiometries::Stoichiometry s;
		s.add(aas::elements::Element(1), 2);
		s.add(aas::elements::Element(8), 1);

		ipaca::detail::Spectrum spectrum = m(s, 0, MyMercury7::PROTON);

		typedef Mercury7<MyStoichiometry, ipaca::detail::Spectrum> My2Mercury7;
		My2Mercury7 m2;
		MyStoichiometry s2 = createIntegerH2O();
		ipaca::detail::Spectrum spectrum2 = m2(s2, 0, My2Mercury7::PROTON);

		shouldEqual(spectrum.size(), spectrum2.size());
		ipaca::detail::Spectrum::const_iterator i1 = spectrum.begin(), i2 =
				spectrum2.begin(), e1 = spectrum.end(), e2 = spectrum2.end();
		for (; i1 != e1 && i2 != e2; ++i1, ++i2) {
			shouldEqualTolerance(i1->mz, i2->mz, 0.0001);
			shouldEqualTolerance(i1->ab, i2->ab, 0.0001);
		}
	}

};

/** The main function that runs the tests for class Combination.
 * Under normal circumstances you need not edit this.
 */
int main() {
	CombinationTestSuite test;
	int success = test.run();
	std::cout << test.report() << std::endl;
	return success;
}

