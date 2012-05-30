/*
 * StoichiometryConfigImpl.cpp
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

#include "MSTK/aas/StoichiometryConfigImpl.hpp"
#include "MSTK/common/Error.hpp"

namespace mstk {
namespace aas {
namespace stoichiometries {

StoichiometryConfigImpl::StoichiometryConfigImplKeyType StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG =
        "DEFAULT_STOICHIOMETRY_CONFIG";

StoichiometryConfigImpl::StoichiometryConfigImpl(
    const StoichiometryConfigImplKeyType& id) :
        id_(id)
{
    if (id == DEFAULT_ELEMENT_CONFIG) {
        std::map < elements::ElementImpl::ElementImplSymbolType, elements::ElementImpl::ElementImplKeyType
                > map = elements::ElementImpl::getDefaultMapping();
        typedef std::map<elements::ElementImpl::ElementImplSymbolType,
                elements::ElementImpl::ElementImplKeyType>::const_iterator IT;
        for (IT it = map.begin(); it != map.end(); ++it) {
            map_.insert(EntryType(it->first, it->second));
        }
    }
}

const StoichiometryConfigImpl::StoichiometryConfigImplKeyType& StoichiometryConfigImpl::getId() const
{
    return id_;
}

StoichiometryConfigImpl::const_iterator StoichiometryConfigImpl::begin() const
{
    return map_.begin();
}

StoichiometryConfigImpl::iterator StoichiometryConfigImpl::begin()
{
    return map_.begin();
}

StoichiometryConfigImpl::const_iterator StoichiometryConfigImpl::end() const
{
    return map_.end();
}

StoichiometryConfigImpl::iterator StoichiometryConfigImpl::end()
{
    return map_.end();
}

void StoichiometryConfigImpl::insertElement(const elements::Element& element)
{
    map_.insert(EntryType(element.get().getSymbol(), element.get().getId()));
}

void StoichiometryConfigImpl::insertElement(
    const elements::ElementImpl::ElementImplSymbolType& symbol,
    const elements::ElementImpl::ElementImplKeyType& key)
{
    map_.insert(EntryType(symbol, key));
}

const elements::ElementImpl::ElementImplKeyType& StoichiometryConfigImpl::getKeyForSymbol(
    const elements::ElementImpl::ElementImplSymbolType& symbol) const
{
    typedef const_iterator IT;
    IT tmp = map_.find(symbol);
    if (tmp != end()) {
        return tmp->second;
    } else {
        throw mstk::LogicError(
            "StoichiometryConfigImpl::getKeyForSymbol(): Cannot find symbol in stoichiometry config.");
    }
}

const StoichiometryConfigImpl::DataType& StoichiometryConfigImpl::getMapping() const
{
    return map_;
}

void StoichiometryConfigImpl::setMapping(const DataType& mapping)
{
    map_ = mapping;
}

StoichiometryConfigImpl StoichiometryConfigImpl::clone(
    const StoichiometryConfigImplKeyType& id) const
{
    StoichiometryConfigImpl tmp(id);
    typedef StoichiometryConfigImpl::const_iterator IT;
    for (IT it = begin(); it != end(); ++it) {
        tmp.map_.insert(EntryType(it->first, it->second));
    }
    return tmp;
}

bool StoichiometryConfigImpl::operator==(
    const StoichiometryConfigImpl& s) const
{
    return id_ == s.id_ && map_ == s.map_;
}

bool StoichiometryConfigImpl::operator!=(
    const StoichiometryConfigImpl& s) const
{
    return !(operator ==(s));
}

StoichiometryConfigImpl& StoichiometryConfigImpl::operator=(
    const StoichiometryConfigImpl& rhs)
{
    if (this != &rhs) {
        id_ = rhs.id_;
        map_ = rhs.map_;
    }
    return *this;
}

std::ostream& operator<<(std::ostream& os, const StoichiometryConfigImpl& o)
{
    os << o.getId() << ":";
    typedef StoichiometryConfigImpl::const_iterator IT;
    for (IT it = o.begin(); it != o.end(); ++it) {
        os << it->first << " " << it->second << ", ";
    }
    return os;
}

} // namespace stoichiometries
} // namespace aas
} // namespace mstk

