/*
 * RawModifications.cpp
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

#include <MSTK/aas/RawModification.hpp>
#include <MSTK/aas/Stoichiometry.hpp>

#include <iostream>

using namespace mstk::aas;
using namespace mstk::aas::elements;
using namespace mstk::aas::modifications;
using namespace mstk::aas::aminoAcids;

int main()
{
    // ========================================================================
    std::cout << "Retrieve a standard modification" << std::endl;

    RawModification mod("Acetyl");
    std::cout << "  Standard modification Acetyl: " << mod << std::endl;

    // ========================================================================
    std::cout << "Create a custom modification" << std::endl;

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

    std::cout << "  Custom modification " << id << ": " << customMod
            << std::endl;

    // ========================================================================
    std::cout << "Add a custom specificity" << std::endl;

    Specificity spec1(RawAminoAcid('K'), Specificity::ANYWHERE,
        Specificity::CHEMICAL_DERIVATIVE);
    spec1.setComment("Some comment");
    customMod.addSpecificity(spec1);

    std::cout << "  Custom modification " << id << ": " << customMod
            << std::endl;

    // ========================================================================
    std::cout << "Add a custom modification" << std::endl;

    // ------------------------------------------------------------------------
    std::cout << " a) manually" << std::endl;

    RawModification customModRef(customMod);

    std::cout << "  Added custom modification " << id << ": " << customModRef
            << std::endl;

    // ------------------------------------------------------------------------
    std::cout << " b) by convenience functions" << std::endl;

    if (!addRawModification(customMod)) {
        std::cout << "  Custom raw modification not added correctly"
                << std::endl;
    }

    std::vector<Specificity> specVector;
    specVector.push_back(spec1);
    std::vector<mstk::String> altNames;
    altNames.push_back(altName);
    if (!addRawModification(id, name, fullName, altNames, myStoichiometry,
        specVector, verified)) {
        std::cout << "  Custom raw modification not added correctly"
                << std::endl;
    }

    // ========================================================================
    std::cout << "Retrieve a custom modification" << std::endl;

    RawModification retrievedCustomMod(id);
    std::cout << "  Retrieved custom raw modification " << id << ": "
            << retrievedCustomMod << std::endl;

    // ========================================================================
    std::cout << "Retrieve the stoichiometry" << std::endl;

    stoichiometries::Stoichiometry stoichiometry = retrievedCustomMod.get().getStoichiometry();
    std::cout << "  Stoichiometry: " << stoichiometry << std::endl;

    // ========================================================================
    std::cout << "Override a standard modification" << std::endl;

    // ------------------------------------------------------------------------
    std::cout << " a) from scratch" << std::endl;
    RawModificationImpl mod1("Amidated", "Amidated", "Full Name", true);
    if (addRawModification(mod1)) {
        std::cout << "  Overridden standard modification Amidated: "
                << RawModification("Amidated") << std::endl;
        std::cout << "  Standard modification Amidated: "
                << RawModificationImpl("Amidated") << std::endl;
    } else {
        std::cout << "  Cannot override standard modification Amidated"
                << std::endl;
    }

    // ------------------------------------------------------------------------
    std::cout << " b) using standard modification" << std::endl;
    RawModificationImpl mod2("Biotin");
    mod2.setFullName("Custom Biotin");
    std::vector<Specificity> specs;
    mod2.setSpecificities(specs);

    if (addRawModification(mod2)) {
        std::cout << "  Overridden standard modification Biotin: "
                << RawModification("Biotin") << std::endl;
        std::cout << "  Standard modification Biotin: "
                << RawModificationImpl("Biotin") << std::endl;
    } else {
        std::cout << "  Cannot override standard modification Biotin"
                << std::endl;
    }

    return 0;
}
