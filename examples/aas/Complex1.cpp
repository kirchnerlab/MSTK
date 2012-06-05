/*
 * Complex1.cpp
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

#include <MSTK/aas/AminoAcidSequence.hpp>

int main()
{

    // ========================================================================
    std::cout
            << "Create an amino acid with a custom modification using a different element configuration"
            << std::endl;
    // ========================================================================
    std::cout << "  create/add custom element" << std::endl;

    mstk::aas::elements::ElementImpl custom13C(
        mstk::aas::elements::ElementImpl::getNextId(), "13C", 6);
    custom13C.addIsotope(mstk::aas::elements::Isotope(12.0, 0.01));
    custom13C.addIsotope(mstk::aas::elements::Isotope(13.0033554, 0.99));
    mstk::aas::elements::addElement(custom13C);

    // ========================================================================
    std::cout << "  create/add stoichiometry configuration" << std::endl;
    mstk::aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType configKey =
            "Label Config";
    mstk::aas::stoichiometries::StoichiometryConfigImpl customConfig(configKey);
    customConfig.insertElement(custom13C.getSymbol(), custom13C.getId());
    mstk::aas::stoichiometries::addStoichiometryConfig(customConfig);

    // ========================================================================
    std::cout << "  create/add custom raw modification" << std::endl;
    mstk::aas::modifications::RawModificationImpl::RawModificationImplKeyType modKey =
            "Label:13C(7)";
    mstk::aas::modifications::RawModificationImpl customRawMod(modKey,
        "Label:13C(7)", "13C(7) custom label", false);
    mstk::aas::stoichiometries::Stoichiometry customModStoichiometry;
    customModStoichiometry.set(
        mstk::aas::elements::Element(
            mstk::aas::elements::ElementImpl::getDefaultKeyForElementSymbol("C")),
        7);
    customModStoichiometry.set(
        mstk::aas::elements::Element(
            mstk::aas::elements::ElementImpl::getDefaultKeyForElementSymbol(
                "13C")), 7);
    customRawMod.setStoichiometry(customModStoichiometry);
    customRawMod.addSpecificity(
        mstk::aas::modifications::Specificity(
            mstk::aas::aminoAcids::RawAminoAcid('Y'),
            mstk::aas::modifications::Specificity::ANYWHERE,
            mstk::aas::modifications::Specificity::ISOTOPIC_LABEL));
    customRawMod.addSpecificity(
        mstk::aas::modifications::Specificity(
            mstk::aas::aminoAcids::RawAminoAcid('F'),
            mstk::aas::modifications::Specificity::ANYWHERE,
            mstk::aas::modifications::Specificity::ISOTOPIC_LABEL));
    mstk::aas::modifications::addRawModification(customRawMod);

    // ========================================================================
    std::cout
            << "  create a modification using the custom raw modification and custom stoichiometry configuration"
            << std::endl;
    mstk::aas::modifications::Modification customMod(modKey, configKey);

    // ========================================================================
    std::cout << "  create/add raw amino acid" << std::endl;
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
    mstk::aas::aminoAcids::addRawAminoAcid(customAminoAcid);

    // ========================================================================
    std::cout << "  create custom modification which binds to U" << std::endl;
    mstk::aas::modifications::RawModificationImpl::RawModificationImplKeyType acetylKey =
            "Acetyl";
    mstk::aas::modifications::RawModificationImpl rawAcetyl(acetylKey);

    mstk::aas::modifications::Modification customAcetyl(acetylKey);
    customAcetyl.setCustomSpecificities(customAcetyl.getSpecificities());
    customAcetyl.addCustomSpecificity(
        mstk::aas::modifications::Specificity(
            mstk::aas::aminoAcids::RawAminoAcid(customAminoAcidKey),
            mstk::aas::modifications::Specificity::ANYWHERE,
            mstk::aas::modifications::Specificity::POST_TRANSLATIONAL));

    // ========================================================================
    std::cout << "  create residue and set amino acid and custom modification"
            << std::endl;
    mstk::aas::AminoAcidSequence aas("ACQUY");
    aas.applyModificationAtPosition(customMod, 5);
    aas.applyModificationAtPosition(customAcetyl, 4);

    // ========================================================================
    std::cout << "Human readable stoichiometry: "
            << aas.getStoichiometry().toString() << std::endl;
    std::cout << "Stoichiometry including isotopic information: "
            << aas.getStoichiometry() << std::endl;
    std::cout << "Amino acid sequence including its modificaitons: "
            << aas.toString(false) << std::endl;

    return 0;
}
