/*
 * ElementImpl.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_ELEMENTIMPL_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_ELEMENTIMPL_HPP__

#include "MSTK/aas/Isotope.hpp"
#include "MSTK/common/Types.hpp"

#include <string>
#include <vector>
#include <map>
#include <iostream>

namespace mstk {
namespace aas {
namespace elements {

/** @addtogroup mstk_aas
 * @{
 */

/**Representation of an element.
 *
 * This class holds all information about an element.
 *
 * TODO do we need a function which checks whether all isotopes add up to 1.0?
 */
class ElementImpl
{

public:

    /**Convenience typedef of the symbol type of an element.
     */
    typedef mstk::String ElementImplSymbolType;
    /**Convenience typedef of the key type of an element.
     */
    typedef mstk::Size ElementImplKeyType;

    /**Default constructor to create a standard element.
     *
     * This constructor has access to static variables containing all
     * information necessary to create and fill the internal variables.
     *
     * @param[in] id Key/Id of the element
     * @throws Throws an exceptions of the given id is not recognized in the
     * list of standard elements.
     */
    ElementImpl(const ElementImplKeyType& id);

    /**Constructor to create a custom element.
     * @param[in] id Key/Id of the element
     * @param[in] symbol Symbol of the element
     * @param[in] atomicNumber Atomic number of the element
     */
    ElementImpl(const ElementImplKeyType& id, const mstk::String& symbol,
        const mstk::Size& atomicNumber);

    /**Constructor to create a custom element.
     * @param[in] id Key/Id of the element
     * @param[in] symbol Symbol of the element
     * @param[in] atomicNumber Atomic number of the element
     * @param[in] isotopes List of isotopes
     */
    ElementImpl(const ElementImplKeyType& id, const mstk::String& symbol,
        const mstk::Size& atomicNumber, const std::vector<Isotope>& isotopes);

    /**Returns the key/id of the element.
     * @returns Key/Id of the element.
     */
    const ElementImplKeyType& getId() const;

    /**Returns the element symbol.
     * @returns Symbol of the element.
     */
    const mstk::String& getSymbol() const;

    /**Returns the atomic number of the element.
     * @return Atomic number of the element
     */
    const mstk::Size& getAtomicNumber() const;

    /**Returns the list of isotopes of the element.
     * @returns List of isotopes of the element.
     */
    const std::vector<Isotope>& getIsotopes() const;

    /**Adds an isotope to the element.
     * @param[in] isotope Isotope
     */
    void addIsotope(const Isotope& isotope);

    /**Adds an isotope to the element.
     * @param[in] mass Mass
     * @param[in] frequency Frequency of the
     */
    void addIsotope(const mstk::Double& mass, const mstk::Double& frequency);

    /**Clears the list of isotopes.
     */
    void clearIsotopes();

    /**Sets the list of isotopes of the element.
     * @param[in] isotopes List of isotopes
     */
    void setIsotopes(const std::vector<Isotope>& isotopes);

    /**Sets a copy of the argument as the new content for the element object.
     * The previous content is dropped.
     * @param[in] rhs Element to copy
     * @returns *this
     */
    ElementImpl& operator=(const ElementImpl& rhs);

    /**Compares the element against another.
     * @param[in] e Element object to compare *this with
     * @returns true if both elements are the same, false otherwise
     */
    bool operator==(const ElementImpl& e) const;

    /**Compares the element against another, with opposite result of
     * ElementImpl::operator==.
     * @param[in] e Element object to compare *this with
     * @returns true if the elements are different, false otherwise.
     */
    bool operator!=(const ElementImpl& e) const;

    /**Returns the number of standard elements.
     * @returns The number of standard elements.
     */
    static mstk::Size getNumberOfStandardElements();

    /**Returns a free key/id which is available to use for a custom element.
     * @returns Free key/id
     */
    static ElementImplKeyType getNextId();

    /**Returns the default mapping of an element symbol to an element key.
     * This mapping is used to create the default stoichiometry configuration.
     * @returns Default mapping of an element symbol to the element key.
     */
    static const std::map<ElementImplSymbolType, ElementImplKeyType>&
    getDefaultMapping();

    /**Returns the key for a given element symbol. This functions only keeps
     * track of standard elements and will not be able to return a key for a
     * custom element.
     * @param[in] symbol Symbol of an element
     * @returns The element key
     * @throws Throws an exception in case the given sybol is not present in
     * the default element mapping.
     */
    static ElementImplKeyType getDefaultKeyForElementSymbol(
        const ElementImplSymbolType& symbol);

private:

    /** Key/Id of the element.
     */
    ElementImplKeyType id_;
    /** Symbol of the element.
     */
    mstk::String symbol_;
    /** Atomic number of the element.
     */
    mstk::Size atomicNumber_;
    /** Isotopes possible for this element.
     */
    std::vector<Isotope> isotopes_;

    /**Indicates a possible free Id.
     */
    static ElementImplKeyType freeId;

    /**Default element mapping.
     */
    static std::map<ElementImplSymbolType, ElementImplKeyType> elementMapping;

};
// class ElementImpl

std::ostream& operator<<(std::ostream&, const ElementImpl&);

inline const ElementImpl::ElementImplKeyType& ElementImpl::getId() const
{
    return id_;
}

/** @\ */

} // namespace elements
} // namespace aas
} // nameapace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_ELEMENTIMPL_HPP__ */
