/*
 * AminoAcid.hpp
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

#include "MSTK/aas/AminoAcid.hpp"

namespace mstk {
namespace aas {
namespace aminoAcids {

AminoAcid::AminoAcid(
    const RawAminoAcidImpl::RawAminoAcidImplKeyType& aminoAcidKey,
    const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configid) :
        rawAminoAcid_(aminoAcidKey), stoichiometryConfig_(configid)
{
}

AminoAcid::AminoAcid(const RawAminoAcid& aminoAcid,
    const aas::stoichiometries::StoichiometryConfig& config) :
        rawAminoAcid_(aminoAcid), stoichiometryConfig_(config)
{
}

Char AminoAcid::getSymbol() const
{
    return rawAminoAcid_.get().getSymbol();
}

const RawAminoAcidImpl::RawAminoAcidImplKeyType& AminoAcid::getRawAminoAcidKey() const
{
    return rawAminoAcid_.get_key();
}

const RawAminoAcid& AminoAcid::getRawAminoAcid() const
{
    return rawAminoAcid_;
}

const String& AminoAcid::getThreeLetterCode() const
{
    return rawAminoAcid_.get().getThreeLetterCode();
}

const String& AminoAcid::getFullName() const
{
    return rawAminoAcid_.get().getFullName();
}

Bool AminoAcid::isNTerm() const
{
    return rawAminoAcid_.get().isNTerm();
}

Bool AminoAcid::isCTerm() const
{
    return rawAminoAcid_.get().isCTerm();
}

aas::stoichiometries::Stoichiometry AminoAcid::getStoichiometry() const
{
    if (stoichiometryConfig_.get_key()
            == aas::stoichiometries::StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG) {
        return rawAminoAcid_.get().getStoichiometry();
    }
    return rawAminoAcid_.get().getStoichiometry().recalculatesWithConfiguration(
        stoichiometryConfig_);
}

void AminoAcid::setStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfig& config)
{
    if (&config != &stoichiometryConfig_) {
        stoichiometryConfig_ = config;
    }
}

void AminoAcid::setStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configid)
{
    aas::stoichiometries::StoichiometryConfig config(configid);
    if (&config != &stoichiometryConfig_) {
        stoichiometryConfig_ = config;
    }
}

const aas::stoichiometries::StoichiometryConfig& AminoAcid::getStoichiometryConfig() const
{
    return stoichiometryConfig_;
}

AminoAcid& AminoAcid::operator=(const AminoAcid& a)
{
    if (this != &a) {
        rawAminoAcid_ = a.rawAminoAcid_;
        stoichiometryConfig_ = a.stoichiometryConfig_;
    }
    return *this;
}

bool AminoAcid::operator==(const AminoAcid& a) const
{
    return rawAminoAcid_ == a.rawAminoAcid_
            && stoichiometryConfig_ == a.stoichiometryConfig_;
}

bool AminoAcid::operator!=(const AminoAcid& a) const
{
    return !(operator ==(a));
}

std::ostream& operator<<(std::ostream& os, const AminoAcid& a)
{
    os << a.getRawAminoAcid();
    os << a.getStoichiometry();
    return os;
}

} // namespace aminoAcids
} // namespace aas
} // namespace mstk
