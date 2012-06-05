/*
 * AminoAcidSequence.cpp
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
    std::cout << "	Create an amino acid sequence" << std::endl;

    mstk::aas::AminoAcidSequence aas("DCGTQS");

    // ========================================================================
    std::cout << "	Add a standard modification" << std::endl;

    mstk::aas::modifications::Modification mod1("Acetyl");
    aas.applyModificationAtPosition(mod1, 2);

    mstk::aas::modifications::Modification mod2("Oxidation");
    mstk::aas::AminoAcidSequence::ModificationList ml;
    ml.push_back(mod2);
    aas.applyFixedModifications(ml);

    // ========================================================================
    std::cout << "	Add a custom modification" << std::endl;

    // start creating custom raw modification
    mstk::aas::modifications::RawModificationImpl::RawModificationImplKeyType id =
            "My Mod";
    mstk::String name = "My mod";
    mstk::String fullName = "My mod is special";
    mstk::Bool verified = true;
    mstk::aas::modifications::RawModificationImpl customMod(id, name, fullName,
        verified);

    mstk::String altName = "My mod is super special";
    customMod.addAltName(altName);

    mstk::aas::stoichiometries::Stoichiometry myStoichiometry;
    myStoichiometry.set(mstk::aas::elements::Element(8), 1);
    myStoichiometry.set(mstk::aas::elements::Element(1), 1);
    customMod.setStoichiometry(myStoichiometry);

    customMod.addSpecificity(
        mstk::aas::modifications::Specificity(
            mstk::aas::aminoAcids::RawAminoAcid('D'),
            mstk::aas::modifications::Specificity::ANYWHERE,
            mstk::aas::modifications::Specificity::ARTEFACT));

    // add custom raw modification
    addRawModification(customMod);

    aas.applyModificationAtPosition(id, 1);

    // ========================================================================
    std::cout << "	Remove a modifications" << std::endl;

    // ========================================================================
    std::cout << "    a remove mod by key" << std::endl;

    aas.remove("Oxidation");

    // ========================================================================
    std::cout << "    b remove mod by ref" << std::endl;

    aas.remove(mod1);

    // ========================================================================
    std::cout << "  Apply stoichiometry config" << std::endl;
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

    aas.applyAminoAcidStoichiometryConfig(custom_key);

    aas.applyModificationStoichiometryConfig(custom_key);

    // ========================================================================
    std::cout << "	Retrieve the stoichiometry" << std::endl;

    aas.getStoichiometry();

    return 0;
}
