/*
 * Modifications.cpp
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

#include <MSTK/aas/Modification.hpp>

#include <iostream>

using namespace mstk::aas;
using namespace mstk::aas::modifications;
using namespace mstk::aas::aminoAcids;
using namespace mstk::aas::elements;

int main()
{

    // ========================================================================
    std::cout << "Use a standard modification" << std::endl;
    Modification acetyl("Acetyl");

    std::cout << "  Standard modification Acetyl: " << acetyl << std::endl;

    // ========================================================================
    std::cout << "Create a custom modification" << std::endl;

    // start creating custom raw modification
    RawModificationImpl::RawModificationImplKeyType id = "My Mod";
    mstk::String name = "My mod";
    mstk::String fullName = "My mod is special";
    mstk::Bool verified = true;
    RawModificationImpl customMod(id, name, fullName, verified);

    mstk::String altName = "My mod is super special";
    customMod.addAltName(altName);

    stoichiometries::Stoichiometry myStoichiometry;
    myStoichiometry.set(Element(8), 1);
    myStoichiometry.set(Element(1), 1);
    customMod.setStoichiometry(myStoichiometry);

    // add custom raw modification
    addRawModification(customMod);

    Modification custom(id);

    std::cout << "  Custom modification " << id << ": " << custom << std::endl;

    // ========================================================================
    std::cout << "Add a custom specificity" << std::endl;

    Specificity spec1(RawAminoAcid('K'), Specificity::ANYWHERE,
        Specificity::CHEMICAL_DERIVATIVE);
    spec1.setComment("Some comment");
    acetyl.addCustomSpecificity(spec1);

    std::cout << " Custom specificity modification acetyl: " << acetyl
            << std::endl;

    // ========================================================================
    std::cout << "Use a custom modification" << std::endl;

    // ========================================================================
    std::cout << "Retrieve the stoichiometry" << std::endl;

    stoichiometries::Stoichiometry stoichiometry = acetyl.getStoichiometry();
    std::cout << "  Retrieved stoichiometry of acetyl: " << stoichiometry
            << std::endl;

    // ========================================================================
    std::cout << "Set a stoichiometry configuration" << std::endl;

    // start creating a custom stoichiometry configuration
    stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType custom_key =
            "My_Stoichiometry_Configuration 1";
    stoichiometries::StoichiometryConfigImpl customConfig(custom_key);
    // create custom element ...
    ElementImpl customElement(ElementImpl::getNextId(), "C", 6);
    customElement.addIsotope(13.0034, 1.0);
    Element customElementRef(customElement);
    // ... and add by instance of an element
    customConfig.insertElement(Element(customElement.getId()));

    stoichiometries::addStoichiometryConfig(customConfig);

    acetyl.setStoichiometryConfig(custom_key);
    std::cout << "  Retrieved stoichiometry with custom config of acetyl: "
            << acetyl.getStoichiometry() << std::endl;

    return 0;
}
