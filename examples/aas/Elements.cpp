/*
 * Elements.cpp
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

#include <MSTK/aas/Element.hpp>

#include <iostream>

using namespace mstk::aas::elements;

int main()
{

    // ========================================================================
    std::cout << "Retrieve a standard element" << std::endl;

    Element C(6);
    std::cout << "  Standard element C: " << C << std::endl;

    // ========================================================================
    std::cout << "Create a custom element" << std::endl;

    ElementImpl::ElementImplKeyType custom_key = ElementImpl::getNextId();
    mstk::String symbol = "C";
    mstk::Size atomicNumber = 6;
    ElementImpl customElement(custom_key, symbol, atomicNumber);
    customElement.addIsotope(13.0034, 1.0);
    std::cout << "  Custom element " << symbol << ": " << customElement
            << std::endl;

    // ========================================================================
    std::cout << "Add a custom element" << std::endl;

    // ------------------------------------------------------------------------
    std::cout << " a) manually" << std::endl;

    Element customElementRef(customElement);
    std::cout << "  Added custom element " << symbol << ": "
            << customElementRef << std::endl;

    // ------------------------------------------------------------------------
    std::cout << " b) by convenience functions" << std::endl;

    if (!addElement(customElement)) {
        std::cout << "  Custom element not added correctly" << std::endl;
    }

    if (addElement(customElement.getId(), customElement.getSymbol(),
        customElement.getAtomicNumber(), customElement.getIsotopes())) {
        std::cout << "  Custom element added correctly" << std::endl;
    } else {
        std::cout << "  Custom element not added correctly" << std::endl;
    }

    // ========================================================================
    std::cout << "Retrieve a custom element" << std::endl;

    Element retrievedCustomElementRef(custom_key);
    std::cout << "  Retrieved custom element " << symbol << ": "
            << retrievedCustomElementRef << std::endl;

    // ========================================================================
    std::cout << "Override a standard element" << std::endl;

    // ------------------------------------------------------------------------
    std::cout << " a) from scratch" << std::endl;

    ElementImpl::ElementImplKeyType N_key = 7;
    mstk::String new_N_symbol = "C";
    mstk::Size new_N_atomicNumber = 6;
    ElementImpl new_N(N_key, new_N_symbol, new_N_atomicNumber);
    new_N.addIsotope(14.0031, 0.5);
    new_N.addIsotope(15.0001, 0.5);
    if (addElement(new_N)) {
        std::cout << "  Overridden standard element N: " << Element(N_key)
                << std::endl;
        std::cout << "  Non-overridden standard element N: "
                << ElementImpl(N_key) << std::endl;
    }

    // ------------------------------------------------------------------------
    std::cout << " b) using the standard element" << std::endl;

    ElementImpl O(8);
    O.clearIsotopes();
    O.addIsotope(15.9994, 1.0);
    if (addElement(O)) {
        std::cout << "  Overridden standard element O: " << Element(8)
                << std::endl;
        std::cout << "  Non-overridden standard element O: " << ElementImpl(8)
                << std::endl;
    }
    // but the default O is no longer available as flyweight

    return 0;
}
