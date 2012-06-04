/*
 * Element.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_ELEMENT_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_ELEMENT_HPP__

#include "MSTK/aas/ElementImpl.hpp"
#include "MSTK/common/Types.hpp"

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>

namespace mstk {
namespace aas {
namespace elements {

/** @addtogroup mstk_aas
 * @{
 */

/**The class ElementIdExtractor is used allow the instantiation of
 * flyweight<ElementImpl>(Key). This simplifies the instantiation and provides
 * a speed up, since the ElementImpl is instantiated only once.
 */
struct ElementIdExtractor
{
    /**Returns the key of the element.
     * @param[in] e Instance of an element implementation
     * @returns The key of the element
     */
    const ElementImpl::ElementImplKeyType& operator()(
        const ElementImpl& e) const
    {
        return e.getId();
    }
};

/**Typedef to simplify the data type flyweight<ElementImpl>
 */
typedef boost::flyweight<
        boost::flyweights::key_value<ElementImpl::ElementImplKeyType,
                ElementImpl, ElementIdExtractor>,
        boost::flyweights::no_tracking> Element;

/**Convenience function to add a custom element to this list of known elements.
 * This methods calls addElement(ElementImpl)
 * @param[in] id Key/Id of the element
 * @param[in] symbol Symbol of the element
 * @param[in] atomicNumber Atomic number of the element
 * @param[in] isotopes List of isotopes of the element
 * @returns True if the given element is added correctly, false otherwise.
 */
mstk::Bool addElement(const ElementImpl::ElementImplKeyType& id,
    const mstk::String& symbol, const mstk::Size& atomicNumber,
    const std::vector<Isotope>& isotopes);

/**Convenience function to add a custom element to the list of known elements.
 * Note: Once an element is added it is not possible to alter its properties.
 * In order to change an internal property you have to create a new element
 * with the correct values and add it to the list. In order to add it the
 * element must have a unique key.
 * IT IS NOT POSSIBLE TO OVERRIDE A REFERENCE (flyweight).
 * @param[in] element Instance of an ElementImpl
 * @returns True if the given element is added correctly, false otherwise.
 */
mstk::Bool addElement(const ElementImpl& element);

bool operator<(const Element& lhs, const Element& rhs);
bool operator<=(const Element& lhs, const Element& rhs);
bool operator>(const Element& lhs, const Element& rhs);
bool operator>=(const Element& lhs, const Element& rhs);

/** @\ */

} // namespace elements
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_ELEMENT_HPP__ */
