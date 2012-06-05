/*
 * Residue.cpp
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

#include <MSTK/aas/Residue.hpp>

int main()
{

    // ========================================================================
    std::cout << "	Create a Residue" << std::endl;

    mstk::aas::Residue res('A');

    // ========================================================================
    std::cout << "	Set/Change the amino acid" << std::endl;

    res.changeType('C');

    // ========================================================================
    std::cout << "	Set/Change the modification" << std::endl;

    res.setModification("Acetyl");

    // ========================================================================
    std::cout << "  Apply stoichiometry configurations" << std::endl;
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

    res.applyAminoAcidStoichiometryConfig(custom_key);
    res.applyModificationStoichiometryConfig(
        mstk::aas::stoichiometries::StoichiometryConfig(custom_key));

    // ========================================================================
    std::cout << "	Retrieve the stoichiometry" << std::endl;

    res.getStoichiometry();

    return 0;
}
