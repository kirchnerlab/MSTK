/*
 * RawModification.cpp
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

#include "MSTK/aas/RawModification.hpp"

namespace mstk {
namespace aas {
namespace modifications {

Bool addRawModification(
    const RawModificationImpl::RawModificationImplKeyType& id,
    const String& name, const String& fullName,
    const std::vector<String>& altNames, const aas::stoichiometries::Stoichiometry& stoichiometry,
    const std::vector<Specificity>& specificities, const Bool& verified)
{
    RawModificationImpl rm(id, name, fullName, verified);
    rm.setAltNames(altNames);
    rm.setStoichiometry(stoichiometry);
    rm.setSpecificities(specificities);
    return addRawModification(rm);
}

Bool addRawModification(const RawModificationImpl& rawModification)
{
    return RawModification(rawModification) == rawModification;
}

bool operator<(const RawModification& lhs, const RawModification& rhs)
{
    return lhs.get_key() < rhs.get_key();
}

bool operator<=(const RawModification& lhs, const RawModification& rhs)
{
    return lhs.get_key() <= rhs.get_key();
}

bool operator>(const RawModification& lhs, const RawModification& rhs)
{
    return lhs.get_key() > rhs.get_key();
}

bool operator>=(const RawModification& lhs, const RawModification& rhs)
{
    return lhs.get_key() >= rhs.get_key();
}

} // namespace modifications
} // namespace aas
} // namespace mstk
