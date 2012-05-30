/*
 * Modification.cpp
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

#include "MSTK/aas/Modification.hpp"

namespace mstk {
namespace aas {
namespace modifications {

Modification::Modification(const RawModification& modification,
    const aas::stoichiometries::StoichiometryConfig& config) :
        rawModification_(modification), stoichiometryConfig_(config), customSpecificities_()
{
}

Modification::Modification(
    const RawModificationImpl::RawModificationImplKeyType& modid,
    const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configid) :
        rawModification_(RawModification(modid)), stoichiometryConfig_(
            aas::stoichiometries::StoichiometryConfig(configid)), customSpecificities_()
{
}

void Modification::setModification(const RawModification& modification)
{
    rawModification_ = modification;
}

void Modification::setModification(
    const RawModificationImpl::RawModificationImplKeyType& modid)
{
    rawModification_ = RawModification(modid);
}

const RawModification& Modification::getModification() const
{
    return rawModification_;
}

void Modification::setStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfig& config)
{
    if (&config != &stoichiometryConfig_) {
        stoichiometryConfig_ = config;
    }
}

void Modification::setStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configid)
{
    aas::stoichiometries::StoichiometryConfig config(configid);
    if (&config != &stoichiometryConfig_) {
        stoichiometryConfig_ = config;
    }
}

const aas::stoichiometries::StoichiometryConfig& Modification::getStoichiometryConfig() const
{
    return stoichiometryConfig_;
}

aas::stoichiometries::Stoichiometry Modification::getStoichiometry() const
{
    if (stoichiometryConfig_.get_key()
            == aas::stoichiometries::StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG) {
        return rawModification_.get().getStoichiometry();
    }
    return rawModification_.get().getStoichiometry().recalculatesWithConfiguration(
        stoichiometryConfig_);
}

void Modification::addCustomSpecificity(const Specificity& specificity)
{
    customSpecificities_.push_back(specificity);
}

void Modification::setCustomSpecificities(
    const std::vector<Specificity>& specificities)
{
    customSpecificities_ = specificities;
}

const std::vector<Specificity>& Modification::getSpecificities() const
{
    if (customSpecificities_.empty()) {
        return rawModification_.get().getSpecificities();
    }
    return customSpecificities_;
}

const std::vector<Specificity>& Modification::getCustomSpecificities() const
{
    return customSpecificities_;
}

void Modification::clearCustomSpecificities()
{
    customSpecificities_.clear();
}

const RawModificationImpl::RawModificationImplKeyType& Modification::getModificationId() const
{
    return rawModification_.get().getId();
}

const String& Modification::getName() const
{
    return rawModification_.get().getName();
}

const String& Modification::getFullName() const
{
    return rawModification_.get().getFullName();
}

const std::vector<String>& Modification::getAltNames() const
{
    return rawModification_.get().getAltNames();
}

const aas::stoichiometries::Stoichiometry& Modification::getRawStoichiometry() const
{
    return rawModification_.get().getStoichiometry();
}

const std::vector<Specificity>& Modification::getRawSpecificities() const
{
    return rawModification_.get().getSpecificities();
}

Bool Modification::isIsotopicLabel() const
{
    // TODO this is true for most of the isotopic labels, but not "TMT", "Succinyl" and "Propionyl"
    typedef std::vector<Specificity>::const_iterator IT;
    for (IT it = getSpecificities().begin(); it != getSpecificities().end();
            ++it) {
        if (it->getClassification() == Specificity::ISOTOPIC_LABEL) {
            return true;
        }
    }
    return false;
}

Bool Modification::isVerified() const
{
    return rawModification_.get().isVerified();
}

Bool Modification::isApplicable(const aminoAcids::AminoAcid& prev,
    const aminoAcids::AminoAcid& current,
    const aminoAcids::AminoAcid& next) const
{
    // TODO maybe it is easier to get the specs from the modification and test directly
    // currently it is forwarded to the instance which is able to decide this
    if (customSpecificities_.empty()) {
        // custom specificities are empty, so check against raw modification
        return rawModification_.get().isApplicable(prev, current, next);
    } else {
        // custom specificites are not empty, so check against those
        typedef std::vector<Specificity>::const_iterator IT;
        IT end = customSpecificities_.end();
        for (IT it = customSpecificities_.begin(); it != end; ++it) {
            if (it->isApplicable(prev.getRawAminoAcid(),
                current.getRawAminoAcid(), next.getRawAminoAcid())) {
                return true;
            }
        }
        return false;
    }
}

bool Modification::operator==(const Modification& m) const
{
    return rawModification_ == m.rawModification_
            && stoichiometryConfig_ == m.stoichiometryConfig_
            && customSpecificities_ == m.customSpecificities_;
}

bool Modification::operator!=(const Modification& m) const
{
    return !(operator ==(m));
}

Modification& Modification::operator=(const Modification& rhs)
{
    if (this != &rhs) {
        rawModification_ = rhs.getModification();
        stoichiometryConfig_ = rhs.getStoichiometryConfig();
        customSpecificities_ = rhs.getCustomSpecificities();
    }
    return *this;
}

std::ostream& operator<<(std::ostream& os, const Modification& o)
{
    os << o.getModification() << "\t" << o.getStoichiometryConfig() << "\t"
            << o.getStoichiometry() << "\t" << o.getCustomSpecificities();
    return os;
}

} // namespace modifications
} // namespace aas
} // namespace mstk
