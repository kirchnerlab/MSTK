/*
 * StoichiometryConfigImpl.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_STOICHIOMETRYCONFIGIMPL_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_STOICHIOMETRYCONFIGIMPL_HPP__

#include "MSTK/aas/Element.hpp"
#include "MSTK/common/Types.hpp"

#include <map>
#include <iostream>

namespace mstk {
namespace aas {
namespace stoichiometries {

/** @addtogroup mstk_aas
 * @{
 */

/** Representation of a stoichiometry configuration.
 *
 */
class StoichiometryConfigImpl
{

public:

    /**Convenience typedef of the stoichiometry configuration's data type.
     */
    typedef std::map<elements::ElementImpl::ElementImplSymbolType,
            elements::ElementImpl::ElementImplKeyType> DataType;
    /**Convenience typedef of the stoichiometry configuration's entry type.
     */
    typedef std::pair<elements::ElementImpl::ElementImplSymbolType,
            elements::ElementImpl::ElementImplKeyType> EntryType;

    /**Convenience typedef of the stoichiometry configuration's const_iterator.
     */
    typedef DataType::const_iterator const_iterator;
    /**Convenience typedef of the stoichiometry configuration's iterator.
     */
    typedef DataType::iterator iterator;
    /**Convenience typedef of the stoichiometry configuration's value_type.
     */
    typedef DataType::value_type value_type;

    /**Convenience typedef of the stoichiometry configuration's  key type.
     */
    typedef mstk::String StoichiometryConfigImplKeyType;

    /** Default constructor.
     * Creates an empty configuration in case the id does not match
     * DEFAULT_ELEMENT_CONFIG.
     * @param[in] id Key/Id of the stoichiometry configuration
     */
    StoichiometryConfigImpl(const StoichiometryConfigImplKeyType& id);

    /**Returns the key/id of the configuration.
     * @returns id of the stoichiometry configuration
     */
    const StoichiometryConfigImplKeyType& getId() const;

    /**Return const iterator to beginning.
     * @returns Const Iterator to beginning
     */
    const_iterator begin() const;

    /**Return iterator to beginning.
     * @returns Iterator to beginning
     */
    iterator begin();

    /**Return iterator to const end.
     * @returns Const iterator to end
     */
    const_iterator end() const;

    /**Return iterator to end.
     * @returns Iterator to end
     */
    iterator end();

    /**Inserts a mapping of element.getSymbol() to element.getId().
     * @param[in] element
     */
    void insertElement(const elements::Element& element);

    /**Inserts a mapping of and element symbol to the element.
     * @param[in] symbol Element symbol
     * @param[in] key Key/Id of the element
     */
    void insertElement(
        const elements::ElementImpl::ElementImplSymbolType& symbol,
        const elements::ElementImpl::ElementImplKeyType& key);

    /**Returns the id of the element which matches the given symbol.
     * @param[in] symbol Element symbol
     * @returns Id of the element
     * @throws Throws an exception if the given symbol is not present in the
     * configuration.
     */
    const elements::ElementImpl::ElementImplKeyType& getKeyForSymbol(
        const elements::ElementImpl::ElementImplSymbolType& symbol) const;

    /**Return mapping of symbol to element ids.
     * @returns Mapping of symbol to element id.
     */
    const DataType& getMapping() const;

    /**Sets element mapping.
     * @param[in] mapping Mapping of an element symbol to an element id
     */
    void setMapping(const DataType& mapping);

    /**Clones the properties of the object and sets the id to the given id.
     * @param[in] id Key/Id of the new stoichiometry configuration
     * @returns A new stoichiometry configuration
     */
    StoichiometryConfigImpl clone(
        const StoichiometryConfigImplKeyType& id) const;

    /**Sets a copy of the argument as the new content for the stoichiometry
     * configuration object.
     * The previous content is dropped.
     * @param[in] rhs Stoichiometry configuration to copy
     * @returns *this
     */
    StoichiometryConfigImpl& operator=(const StoichiometryConfigImpl& rhs);

    /**Compares the stoichiometry configuration against another.
     * @param[in] s Stoichiometry configuration object to compare *this with
     * @returns true if both stoichiometry configurations are the same, false
     * otherwise
     */
    bool operator==(const StoichiometryConfigImpl& s) const;

    /**Compares the stoichiometry configuration against another, with opposite
     * result of StoichiometryConfigurationImpl::operator==.
     * @param[in] s Stoichiometry configuration object to compare *this with
     * @returns true if the stoichiometry configuration are different, false
     * otherwise
     */
    bool operator!=(const StoichiometryConfigImpl& s) const;

    /**The Key/Id of the default stoichiometry configuration. Using this key
     * will result in a stoichiometry configuration holding the default mapping
     * of all standard elements.
     */
    static StoichiometryConfigImplKeyType DEFAULT_ELEMENT_CONFIG;

private:

    /**Key/Id of the stoichiometry configuration
     */
    StoichiometryConfigImplKeyType id_;
    /**Mapping of ElementImpl::ElementSymbolType to ElementImpl::ElementKeyType
     */
    DataType map_;

};
// class StoichiometryConfigImpl

std::ostream& operator<<(std::ostream&, const StoichiometryConfigImpl&);

/** @\ */

} // namespace stoichiometries
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_STOICHIOMETRYCONFIGIMPL_HPP__ */
