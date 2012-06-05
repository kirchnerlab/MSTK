/*
 * AminoAcid.cpp
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

#include <MSTK/aas/AminoAcid.hpp>

int main() {

	// ========================================================================
	std::cout << "Create an amino acid" << std::endl;

	mstk::aas::aminoAcids::AminoAcid aa('C');

	// ========================================================================
	std::cout << "Set a stoichiometry configuration" << std::endl;

	// start creating a custom stoichiometry configuration
	mstk::aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType custom_key =
			"My_Stoichiometry_Configuration 1";
	mstk::aas::stoichiometries::StoichiometryConfigImpl customConfig(custom_key);
	// create custom element ...
	mstk::aas::elements::ElementImpl customElement(
			mstk::aas::elements::ElementImpl::getNextId(), "C", 6);
	customElement.addIsotope(13.0034, 1.0);
	mstk::aas::elements::Element customElementRef(customElement);
	// ... and add by instance of an element
	customConfig.insertElement(
			mstk::aas::elements::Element(customElement.getId()));
	addStoichiometryConfig(customConfig);

	aa.setStoichiometryConfig(custom_key);

	// ========================================================================
	std::cout << "Retrieve the stoichiometry" << std::endl;

	aa.getStoichiometry();

	return 0;
}

