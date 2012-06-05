/*
 * RawAminoAcid.cpp
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

#include <MSTK/aas/RawAminoAcid.hpp>

int main()
{

    // ========================================================================
    std::cout << "Retrieve a standard raw amino acid" << std::endl;
    mstk::aas::aminoAcids::RawAminoAcid cystein('C');

    // ========================================================================
    std::cout << "Create a custom raw amino acid" << std::endl;
    mstk::aas::aminoAcids::RawAminoAcidImpl::RawAminoAcidImplKeyType customAminoAcidKey =
            'U';
    mstk::aas::stoichiometries::Stoichiometry customAAStoichiometry;
    customAAStoichiometry.set(
        mstk::aas::elements::Element(
            mstk::aas::elements::ElementImpl::getDefaultKeyForElementSymbol("C")),
        3);
    customAAStoichiometry.set(
        mstk::aas::elements::Element(
            mstk::aas::elements::ElementImpl::getDefaultKeyForElementSymbol("H")),
        7);
    customAAStoichiometry.set(
        mstk::aas::elements::Element(
            mstk::aas::elements::ElementImpl::getDefaultKeyForElementSymbol("N")),
        1);
    customAAStoichiometry.set(
        mstk::aas::elements::Element(
            mstk::aas::elements::ElementImpl::getDefaultKeyForElementSymbol("O")),
        2);
    customAAStoichiometry.set(
        mstk::aas::elements::Element(
            mstk::aas::elements::ElementImpl::getDefaultKeyForElementSymbol(
                "Se")), 1);
    mstk::aas::aminoAcids::RawAminoAcidImpl customAminoAcid(customAminoAcidKey,
        'U', customAAStoichiometry);
    customAminoAcid.setThreeLetterCode("Sec");
    customAminoAcid.setFullName("Selenocysteine");

    // ========================================================================
    std::cout << "Add a custom raw amino acid" << std::endl;
    mstk::aas::aminoAcids::addRawAminoAcid(customAminoAcid);

    // ========================================================================
    std::cout << "Override a standard raw amino acid" << std::endl;
    mstk::aas::aminoAcids::RawAminoAcidImpl serine('S');
    serine.setFullName("Serine is a proteinogenic amino acid");
    mstk::aas::aminoAcids::addRawAminoAcid(serine);

    // ========================================================================
    std::cout << "Retrieve the stoichiometry" << std::endl;
    serine.getStoichiometry();

    return 0;
}

