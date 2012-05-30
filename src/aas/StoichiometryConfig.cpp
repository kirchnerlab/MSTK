/*
 * StoichiometryConfig.cpp
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

#include "MSTK/aas/StoichiometryConfig.hpp"

namespace mstk {
namespace aas {
namespace stoichiometries {

bool operator<(const StoichiometryConfig& lhs, const StoichiometryConfig& rhs)
{
    return lhs.get_key() < rhs.get_key();
}

bool operator<=(const StoichiometryConfig& lhs, const StoichiometryConfig& rhs)
{
    return lhs.get_key() <= rhs.get_key();
}

bool operator>(const StoichiometryConfig& lhs, const StoichiometryConfig& rhs)
{
    return lhs.get_key() > rhs.get_key();
}

bool operator>=(const StoichiometryConfig& lhs, const StoichiometryConfig& rhs)
{
    return lhs.get_key() >= rhs.get_key();
}

Bool addStoichiometryConfig(
    const StoichiometryConfigImpl::StoichiometryConfigImplKeyType& id,
    const StoichiometryConfigImpl::DataType& map)
{
//    ElementConfig dec(ElementConfigImpl::DEFAULT_ELEMENT_CONFIG);
    StoichiometryConfigImpl ec(id);
    typedef StoichiometryConfigImpl::DataType::const_iterator IT;
    for (IT it = map.begin(); it != map.end(); ++it) {
        ec.insertElement(it->first, it->second);
    }
    return addStoichiometryConfig(ec);
}

Bool addStoichiometryConfig(
    const StoichiometryConfigImpl& stoichiometryConfig)
{
    return StoichiometryConfig(stoichiometryConfig) == stoichiometryConfig;
}

} // namespace stoichiometries
} // namespace aas
} // namespace mstk
