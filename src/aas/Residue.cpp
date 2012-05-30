/*
 * Residue.cpp
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

#include "MSTK/aas/Residue.hpp"
#include "MSTK/common/Error.hpp"

#include <boost/make_shared.hpp>

#include <stdio.h>
#include <sstream>

namespace mstk {
namespace aas {

Residue::ModificationPtr Residue::EMPTY_MOD;

Residue::Residue(
    const aas::aminoAcids::RawAminoAcidImpl::RawAminoAcidImplKeyType& aminoAcidKey,
    const aas::modifications::RawModificationImpl::RawModificationImplKeyType& modificationKey,
    const aas::modifications::RawModificationImpl::RawModificationImplKeyType& labelKey) :
        aminoacid_(aminoAcidKey), modification_(
            boost::make_shared<modifications::Modification>()), isotopicLabel_(
            boost::make_shared<modifications::Modification>())
{
    if (!EMPTY_MOD.get()) {
        EMPTY_MOD = ModificationPtr(new modifications::Modification(""));
    }
    if (modificationKey != "") {
        setModification(modificationKey);
    } else {
        modification_ = EMPTY_MOD;
    }
    if (labelKey != "") {
        setIsotopicLabel(labelKey);
    } else {
        isotopicLabel_ = EMPTY_MOD;
    }
}

Residue::Residue(const aas::aminoAcids::AminoAcid& aminoAcid,
    const aas::modifications::Modification& mod,
    const aas::modifications::Modification& label) :
        aminoacid_(aminoAcid), modification_(
            boost::make_shared<modifications::Modification>()), isotopicLabel_(
            boost::make_shared<modifications::Modification>())
{
    if (!EMPTY_MOD.get()) {
        EMPTY_MOD = ModificationPtr(new modifications::Modification(""));
    }
    if (mod.getModificationId() != "") {
        setModification(mod);
    } else {
        modification_ = EMPTY_MOD;
    }
    if (label.getModificationId() != "") {
        setIsotopicLabel(label);
    } else {
        isotopicLabel_ = EMPTY_MOD;
    }
}

void Residue::changeType(
    const aminoAcids::RawAminoAcidImpl::RawAminoAcidImplKeyType& aminoAcidKey)
{
    aminoacid_ = aas::aminoAcids::AminoAcid(aminoAcidKey);
}

void Residue::changeType(const aminoAcids::AminoAcid& aminoAcid)
{
    aminoacid_ = aminoAcid;
}

Bool Residue::hasModification(
    const modifications::RawModificationImpl::RawModificationImplKeyType& modificationKey) const
{
    return modification_->getModificationId() == modificationKey;
}

Bool Residue::hasModification(
    const modifications::Modification& modification) const
{
    return *modification_ == modification;
}

Bool Residue::isModified() const
{
    return modification_->getModificationId() != "";
}

Bool Residue::hasLabel(
    const modifications::RawModificationImpl::RawModificationImplKeyType& labelKey) const
{
    return isotopicLabel_->getModificationId() == labelKey;
}

Bool Residue::hasLabel(const modifications::Modification& label) const
{
    return *isotopicLabel_ == label;
}

Bool Residue::isLabeled() const
{
    return isotopicLabel_->getModificationId() != "";
}

void Residue::removeModification()
{
    modification_ = EMPTY_MOD;
}

void Residue::removeIsotopicLabel()
{
    isotopicLabel_ = EMPTY_MOD;
}

void Residue::applyAminoAcidStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configKey)
{
    aminoacid_.setStoichiometryConfig(configKey);
}

void Residue::applyAminoAcidStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfig& config)
{
    aminoacid_.setStoichiometryConfig(config);
}

void Residue::applyModificationStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configKey)
{
    applyModificationStoichiometryConfig(aas::stoichiometries::StoichiometryConfig(configKey));
}

void Residue::applyModificationStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfig& config)
{
    if (isModified()) {
        if (modification_.use_count() > 1) {
            modification_ = ModificationPtr(
                new modifications::Modification(*modification_));
        }
        modification_->setStoichiometryConfig(config);
    }
}

void Residue::applyIsotopicLabelStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configKey)
{
    applyIsotopicLabelStoichiometryConfig(aas::stoichiometries::StoichiometryConfig(configKey));
}

void Residue::applyIsotopicLabelStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfig& config)
{
    if (isLabeled()) {
        if (isotopicLabel_.use_count() > 1) {
            isotopicLabel_ = ModificationPtr(
                new modifications::Modification(*isotopicLabel_));
        }
        isotopicLabel_->setStoichiometryConfig(config);
    }
}

aas::stoichiometries::Stoichiometry Residue::getStoichiometry() const
{
    aas::stoichiometries::Stoichiometry s = aminoacid_.getStoichiometry();
    s += modification_->getStoichiometry();
    s += isotopicLabel_->getStoichiometry();
    return s;
}

String Residue::toString() const
{
    std::ostringstream oss;
    oss << aminoacid_.getSymbol();
    Bool labeled = isLabeled(), modified = isModified();
    if (labeled || modified) {
        oss << "(";
    }
    if (isModified()) {
        oss << modification_->getModificationId();
    }
    if (isLabeled()) {
        if (modified) {
            oss << "; ";
        }
        oss << isotopicLabel_->getModificationId();
    }
    if (labeled || modified) {
        oss << ")";
    }

    return oss.str();
}

Bool Residue::isNTerm() const
{
    return aminoacid_.isNTerm();
}

Bool Residue::isCTerm() const
{
    return aminoacid_.isCTerm();
}

const aas::aminoAcids::AminoAcid& Residue::getAminoAcid() const
{
    return aminoacid_;
}

aas::aminoAcids::AminoAcid& Residue::getAminoAcid()
{
    return aminoacid_;
}

void Residue::setModification(
    const aas::modifications::RawModificationImpl::RawModificationImplKeyType& modificationKey)
{
    setModification(modifications::Modification(modificationKey));
}

void Residue::setModification(
    const aas::modifications::Modification& modification)
{
    if (!modification.isIsotopicLabel()) {
        modification_ = ModificationPtr(
            new modifications::Modification(modification));
    } else {
        std::ostringstream oss;
        oss << "Residue::setModification(): Given modification '";
        oss << modification.getModificationId();
        oss
                << "' is an isotopic label since modification.isIsotopicLabel() returned true. Use setIsotopicLabel() instead.";
        throw mstk::LogicError(oss.str());
    }
}

const modifications::Modification& Residue::getModification() const
{
    return *modification_;
}

void Residue::setIsotopicLabel(
    const aas::modifications::RawModificationImpl::RawModificationImplKeyType& isotopicLabelKey)
{
    setIsotopicLabel(modifications::Modification(isotopicLabelKey));
}

void Residue::setIsotopicLabel(
    const aas::modifications::Modification& isotopicLabel)
{
    if (isotopicLabel.isIsotopicLabel()) {
        isotopicLabel_ = ModificationPtr(
            new modifications::Modification(isotopicLabel));
    } else {
        std::ostringstream oss;
        oss << "Residue::setIsotopicLabel(): Given isotopic label '";
        oss << isotopicLabel.getModificationId();
        oss
                << "' is a standard modification since isotopicLabel.isIsotopicLabel() returned false. Use setModification() instead.";
        throw mstk::LogicError(oss.str());
    }
}

const modifications::Modification& Residue::getIsotopicLabel() const
{
    return *isotopicLabel_;
}

Residue& Residue::operator=(const Residue& rhs)
{
    if (this != &rhs) {
        aminoacid_ = rhs.aminoacid_;
        modification_ = rhs.modification_;
        isotopicLabel_ = rhs.isotopicLabel_;
    }
    return *this;
}

bool Residue::operator==(const Residue& r) const
{
    return aminoacid_ == r.aminoacid_ && *modification_ == *r.modification_
            && *isotopicLabel_ == *r.isotopicLabel_;
}

bool Residue::operator!=(const Residue& r) const
{
    return !(operator ==(r));
}

std::ostream& operator<<(std::ostream& os, const Residue& o)
{
    os << o.getAminoAcid() << "\t" << o.getModification() << "\t"
            << o.getIsotopicLabel();
    return os;
}

} // namespace aas
} // namespace mstk
