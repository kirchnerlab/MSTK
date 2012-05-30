/*
 * Element.cpp
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

#include "MSTK/aas/Element.hpp"

namespace mstk {
namespace aas {
namespace elements {

bool operator<(const Element& lhs, const Element& rhs)
{
    return lhs.get_key() < rhs.get_key();
}

bool operator<=(const Element& lhs, const Element& rhs)
{
    return lhs.get_key() <= lhs.get_key();
}

bool operator>(const Element& lhs, const Element& rhs)
{
    return lhs.get_key() > lhs.get_key();
}

bool operator>=(const Element& lhs, const Element& rhs)
{
    return lhs.get_key() >= rhs.get_key();
}

Bool addElement(const ElementImpl::ElementImplKeyType& id,
    const String& symbol, const Size& atomicNumber,
    const std::vector<Isotope>& isotopes)
{
    return addElement(ElementImpl(id, symbol, atomicNumber, isotopes));
}

Bool addElement(const ElementImpl& element)
{
    Element element_ref(element);
    return element_ref == element;
}

} // namespace elements
} // namespace aas
} // namespace mstk
