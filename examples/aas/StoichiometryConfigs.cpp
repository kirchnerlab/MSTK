/*
 * StoichiometryConfigs.cpp
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

#include <MSTK/aas/StoichiometryConfig.hpp>

#include <iostream>

using namespace mstk::aas;
using namespace mstk::aas::elements;

int main()
{

    // ========================================================================
    std::cout << "Retrieve default stoichiometry configuration" << std::endl;

    // in order to show how to override the default configuration, we are not
    // allowed to retrieve the default config within this application before we
    // override it.

    /*
     StoichiometryConfig defaultStoichiometry(
     StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);
     std::cout << "  Default stoichiometry configuration: "
     << defaultStoichiometry << std::endl;
     */

    // ========================================================================
    std::cout << "Create a custom stoichiometry configuration" << std::endl;

    stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType custom_key =
            "My_Stoichiometry_Configuration 1";
    stoichiometries::StoichiometryConfigImpl customConfig(custom_key);
    // add element mapping manually
    ElementImpl::ElementImplSymbolType es1 = "O";
    ElementImpl::ElementImplKeyType ek1 = 8;
    customConfig.insertElement(es1, ek1);
    // add element mapping by instance of an element
    Element C(6);
    customConfig.insertElement(C);
    // create custom element ...
    ElementImpl customElement(ElementImpl::getNextId(), "my O 18", 8);
    customElement.addIsotope(17.9992, 1.0);
    Element customElementRef(customElement);
    // ... and add by instance of an element
    customConfig.insertElement(Element(customElement.getId()));

    std::cout << "  Custom stoichiometry configuration " << customConfig
            << std::endl;

    // ========================================================================
    std::cout << "Add a custom stoichiometry configuration" << std::endl;

    // ------------------------------------------------------------------------
    std::cout << " a) manually" << std::endl;

    stoichiometries::StoichiometryConfig customConfigRef(customConfig);
    std::cout << "  Added custom stoichiometry configuration: "
            << customConfigRef << std::endl;

    // ------------------------------------------------------------------------
    std::cout << " b) by convenience functions" << std::endl;
    if (!addStoichiometryConfig(customConfig)) {
        std::cout
                << "  Custom stoichiometry configuration not added correctly."
                << std::endl;
    }

    if (!stoichiometries::addStoichiometryConfig(custom_key, customConfig.getMapping())) {
        std::cout
                << "  Custom stoichiometry configuration not added correctly."
                << std::endl;
    }

    // ========================================================================
    std::cout << "Retrieve a custom stoichiometry configuration" << std::endl;

    stoichiometries::StoichiometryConfig retrievedCustomConfig(custom_key);
    std::cout << "  Retrieved custom stoichiometry configuration: "
            << retrievedCustomConfig << std::endl;

    // ========================================================================
    std::cout << "Override the default stoichiometry configuration"
            << std::endl;

    // ------------------------------------------------------------------------
    std::cout << " b) using standard mapping" << std::endl;
    stoichiometries::StoichiometryConfigImpl defaultConfig(
        stoichiometries::StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);
    defaultConfig.insertElement(Element(customElement.getId()));
    stoichiometries::StoichiometryConfig defaultConfigRef(defaultConfig);
    if (addStoichiometryConfig(defaultConfig)) {
        std::cout << "  New default stoichiometry configuration: "
                << defaultConfig << std::endl;
        std::cout << "  Current default stoichiometry configuration: "
                << stoichiometries::StoichiometryConfig(
                    stoichiometries::StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG)
                << std::endl;
    }

    return 0;
}
