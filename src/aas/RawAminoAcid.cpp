/*
 * RawAminoAcid.cpp
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

#include "MSTK/aas/RawAminoAcid.hpp"

namespace mstk {
namespace aas {
namespace aminoAcids {

bool operator<(const RawAminoAcid& lhs, const RawAminoAcid& rhs)
{
    return lhs.get_key() < rhs.get_key();
}

bool operator<=(const RawAminoAcid& lhs, const RawAminoAcid& rhs)
{
    return lhs.get_key() <= rhs.get_key();
}

bool operator>(const RawAminoAcid& lhs, const RawAminoAcid& rhs)
{
    return lhs.get_key() > rhs.get_key();
}

bool operator>=(const RawAminoAcid& lhs, const RawAminoAcid& rhs)
{
    return lhs.get_key() >= rhs.get_key();
}

Bool addRawAminoAcid(const RawAminoAcidImpl::RawAminoAcidImplKeyType& id,
    const Char symbol, const String& threeLetterCode, const String& fullName,
    const aas::stoichiometries::Stoichiometry& stoichiometry)
{
    RawAminoAcidImpl aa(id, symbol, stoichiometry);
    aa.setThreeLetterCode(threeLetterCode);
    aa.setFullName(fullName);
    return addRawAminoAcid(aa);
}

Bool addRawAminoAcid(const RawAminoAcidImpl& aminoAcid)
{
    RawAminoAcid aminoAcid_r(aminoAcid);
    // in case the key of the amino acid aminoAcid was already added or is a
    // standard amino acid which was retrieved earlier, the reference aminoAcid_r
    // will not contain the information as given in aminoAcid
    return aminoAcid_r == aminoAcid;
}

} // namespace aminoAcids
} // namespace aas
} // namespace mstk
