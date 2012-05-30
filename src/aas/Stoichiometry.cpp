/*
 * Stoichiometry.cpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
 * Copyright (c) 2010 Nathan Hueksen
 * Copyright (c) 2009,2010 Marc Kirchner
 * Copyright (c) 2008 Thorben Kroeger
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

#include "MSTK/aas/Stoichiometry.hpp"
#include "MSTK/common/Error.hpp"

#include <sstream>

namespace mstk {
namespace aas {
namespace stoichiometries {

Stoichiometry::Stoichiometry() :
        annotationId_(0), counts_()
{
}

void Stoichiometry::clear()
{
    counts_.clear();
    annotationId_ = 0;
}

Stoichiometry::const_iterator Stoichiometry::begin() const
{
    return counts_.begin();
}

Stoichiometry::iterator Stoichiometry::begin()
{
    return counts_.begin();
}

Stoichiometry::const_iterator Stoichiometry::end() const
{
    return counts_.end();
}

Stoichiometry::iterator Stoichiometry::end()
{
    return counts_.end();
}

Bool Stoichiometry::empty() const
{
    return counts_.empty();
}

Size Stoichiometry::size() const
{
    return counts_.size();
}

Double Stoichiometry::get(const aas::elements::Element& element) const
{
    const_iterator elem = counts_.find(element);
    if (elem == counts_.end()) {
        return 0.0; // Not found, so zero atoms present
    } else {
        return elem->second;
    }
}

void Stoichiometry::set(const aas::elements::Element& element,
    const Double& count)
{
    if (count == 0.0) {
        counts_.erase(element); // Not in here has the same meaning as count =0
    } else {
        counts_[element] = count;
    }
}

void Stoichiometry::add(const aas::elements::Element& element,
    const Double& count)
{
    Double& c = counts_[element];
    c += count;
    if (c == 0.0) {
        counts_.erase(element); // Not present has the same meaning as count = 0
    }
}

Bool Stoichiometry::nonNegative() const
{
    const_iterator iter;
    for (iter = counts_.begin(); iter != counts_.end(); ++iter) {
        if (iter->second < 0.0) {
            return false;
        }
    }
    return true;
}

String Stoichiometry::toString() const
{
    std::ostringstream oss;
    typedef DataType::const_iterator IT;
    for (IT it = counts_.begin(); it != counts_.end(); ++it) {
        oss << it->first.get().getSymbol() << "(" << it->second << ")";
    }
    return oss.str();
}

// TODO try to change StoichiometryConfig to const StoichiometryConfigImpl& to avoid all this .get() stuff where possible.
// this applies also for all methods which receive a fw but acutally do not need it
// the conversion between both types should be done implicitly

void Stoichiometry::applyStoichiometryConfiguration(
    const StoichiometryConfig& config)
{
    typedef Stoichiometry::iterator IT;

    Stoichiometry ret;

    StoichiometryConfig defaultConfig = StoichiometryConfig(
        StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);

    // iterate over all elements in rawStoichiometry
    for (IT it = begin(); it != end(); ++it) {
        elements::ElementImpl::ElementImplKeyType elementId = 0;
        const String& symbol = it->first.get().getSymbol();
        try {
            elementId = config.get().getKeyForSymbol(symbol);
        } catch (mstk::LogicError& e) {
            // cannot find symbol in custom map
            // searching symbol in default map
            try {
                elementId = defaultConfig.get().getKeyForSymbol(symbol);
            } catch (mstk::LogicError& e) {
                // cannot find symbol in default map
                mstk_fail(
                    "Stoichiometry::applyStoichiometryConfiguration(): Cannot find element symbol.");
            }
        }
        // MAYBE optimize by changing values directly
        // remove ++it in for
        //		if (elementId != it->first.get_key()) {
        //			Double count = it->second;
        //			elements::Element e = it->first;
        //			++it;
        //			counts_.erase(e);
        //			set(aas::elements::Element(elementId), count);
        //		} else {
        //			++it;
        //		}
        if (elementId != it->first.get_key()) {
            ret.set(aas::elements::Element(elementId), it->second);
            ret.set(it->first, -it->second);
        }
    }
    *this += ret;
}

Stoichiometry Stoichiometry::recalculatesWithConfiguration(
    const StoichiometryConfig& config) const
{
    Stoichiometry retStoichiometry = *this;
    retStoichiometry.applyStoichiometryConfiguration(config);
    return retStoichiometry;
}

Stoichiometry& Stoichiometry::operator=(const Stoichiometry& s)
{
    if (this != &s) {
        annotationId_ = s.annotationId_;
        counts_ = s.counts_;
    }
    return *this;
}

bool Stoichiometry::operator==(const Stoichiometry& s) const
{
    // Since it is assured, that no zero elements exists,
    // the stoichiometries are equal if the maps are equal
    return counts_ == s.counts_;
}

bool Stoichiometry::operator!=(const Stoichiometry& s) const
{
    return !operator==(s);
}

Stoichiometry& Stoichiometry::operator+=(const Stoichiometry& s)
{
    const_iterator other_iter = s.counts_.begin();
    for (other_iter = s.begin(); other_iter != s.end(); ++other_iter) {
        add(other_iter->first, other_iter->second);
    }
    return *this;
}

Stoichiometry Stoichiometry::operator+(const Stoichiometry& s)
{
    Stoichiometry stoi = *this; // Make a copy
    stoi += s;
    return stoi;
}

Stoichiometry& Stoichiometry::operator-=(const Stoichiometry& s)
{
    const_iterator other_iter = s.begin();
    for (other_iter = s.begin(); other_iter != s.end(); ++other_iter) {
        add(other_iter->first, -other_iter->second);
    }
    return *this;
}

Stoichiometry Stoichiometry::operator-(const Stoichiometry& s)
{
    Stoichiometry stoi = *this; //Make a copy
    stoi -= s;
    return stoi;
}

std::ostream& operator<<(std::ostream& o, const Stoichiometry& s)
{
    for (Stoichiometry::const_iterator it = s.begin(); it != s.end(); ++it) {
        o << "(" << it->first << ")" << it->second << " ";
    }
    o << '[' << s.getAnnotationId() << ']';
    return o;
}

std::ostream& operator<<(std::ostream& o,
    const std::vector<Stoichiometry>& s)
{
    for (std::vector<Stoichiometry>::const_iterator it = s.begin();
            it != s.end(); ++it) {
        o << *it << std::endl;
    }
    return o;
}

} // namespace stoichiometries
} // namespace aas
} // namespace mstk
