/*
 * Combination-test.hpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
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

#include <MSTK/aas/adapter/LibIPACA.hpp>

#include <MSTK/ipaca/Mercury7.hpp>
#include <MSTK/ipaca/Mercury7Impl.hpp>

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

/** The main function that runs the tests for class Combination.
 * Under normal circumstances you need not edit this.
 */
int main() {
	typedef ipaca::Mercury7<aas::stoichiometries::Stoichiometry,
			ipaca::detail::Spectrum> MyMercury7;
	MyMercury7 m;

	aas::stoichiometries::Stoichiometry s;
	s.add(aas::elements::Element(1), 2);
	s.add(aas::elements::Element(8), 1);

	ipaca::detail::Spectrum spectrum = m(s, 0, MyMercury7::PROTON);

	std::cout << spectrum << std::endl;

	return 0;
}

