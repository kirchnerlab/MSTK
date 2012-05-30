/*
 * Specificity.cpp
 *
 * Copyright (c) 2011 Mathias Wilhelm
 * Copyright (c) 2011 Marc Kirchner
 * Copyright (c) 2010 Nathan Huesken
 * Copyright (c) 2010 Marc Kirchner
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
#include "MSTK/aas/Specificity.hpp"
#include "MSTK/common/Error.hpp"

#include <algorithm>

namespace mstk {
namespace aas {
namespace modifications {

Specificity::Specificity(const aas::aminoAcids::RawAminoAcid& site,
    const Position position, const Classification classification) :
        site_(site), position_(position), classification_(classification), neutralLosses_(), pepNeutralLosses_(), comment_(
            "")
{
}

Specificity::Specificity(const String& site, const String& position,
    const String& classification) :
        neutralLosses_(), pepNeutralLosses_(), comment_("")
{
    site_ = aas::aminoAcids::RawAminoAcid(
        aas::aminoAcids::RawAminoAcidImpl::getKeyForAminoAcidString(site));
    position_ = Specificity::parsePositionString(position);
    classification_ = Specificity::parseClassificationString(classification);
}

void Specificity::setSite(const aas::aminoAcids::RawAminoAcid& aminoAcid)
{
    site_ = aminoAcid;
}

const aas::aminoAcids::RawAminoAcid& Specificity::getSite() const
{
    return site_;
}

void Specificity::setClassification(const Classification& classification)
{
    classification_ = classification;
}

const Specificity::Classification& Specificity::getClassification() const
{
    return classification_;
}

void Specificity::setPosition(const Position& position)
{
    position_ = position;
}
const Specificity::Position& Specificity::getPosition() const
{
    return position_;
}

void Specificity::addNeutralLoss(
    const aas::stoichiometries::Stoichiometry& stoichiometry)
{
    neutralLosses_.push_back(stoichiometry);
}

void Specificity::setNeutralLosses(
    const std::vector<aas::stoichiometries::Stoichiometry>& stoichiometries)
{
    neutralLosses_ = stoichiometries;
}

const std::vector<aas::stoichiometries::Stoichiometry>& Specificity::getNeutralLosses() const
{
    return neutralLosses_;
}

void Specificity::clearNeutralLosses()
{
    neutralLosses_.clear();
}

void Specificity::addPepNeutralLoss(
    const aas::stoichiometries::Stoichiometry& stoichiometry)
{
    pepNeutralLosses_.push_back(stoichiometry);
}

void Specificity::setPepNeutralLosses(
    const std::vector<aas::stoichiometries::Stoichiometry>& stoichiometries)
{
    pepNeutralLosses_ = stoichiometries;
}

const std::vector<aas::stoichiometries::Stoichiometry>& Specificity::getPepNeutralLosses() const
{
    return pepNeutralLosses_;
}

void Specificity::clearPepNeutralLosses()
{
    pepNeutralLosses_.clear();
}

void Specificity::setComment(const String& comment)
{
    comment_ = comment;
}

const String& Specificity::getComment() const
{
    return comment_;
}

Bool Specificity::isApplicable(const aas::aminoAcids::RawAminoAcid& prev,
    const aas::aminoAcids::RawAminoAcid& current,
    const aas::aminoAcids::RawAminoAcid& next) const
{
    Bool matchingSites = site_.get_key() == current.get_key();
    // MAYBE optimize by evaluating the boolean statements?
    switch (position_) {
        case ANY_N_TERM:
            // either site match current and prev is n-term or site is n-term and current is n-term
            if (!((matchingSites && prev.get().isNTerm())
                    || (site_.get().isNTerm() && current.get().isNTerm()))) {
                return false;
            }
            break;
        case ANY_C_TERM:
            // either site match current and next is c-term or site is c-term and current is c-term
            if (!((matchingSites && next.get().isCTerm())
                    || (site_.get().isCTerm() && current.get().isCTerm()))) {
                return false;
            }
            break;
        case PROTEIN_N_TERM:
            // either site match current and prev is prot n-term or site and current is prot n-term
            if (!((matchingSites
                    && prev.get_key()
                            == aminoAcids::RawAminoAcidImpl::PROTEIN_N_TERM)
                    || (site_.get_key()
                            == aminoAcids::RawAminoAcidImpl::PROTEIN_N_TERM
                            && current.get_key()
                                    == aminoAcids::RawAminoAcidImpl::PROTEIN_N_TERM))) {
                return false;
            }
            break;
        case PROTEIN_C_TERM:
            // either site match current and next is prot c-term or site and current is prot c-term
            if (!((matchingSites
                    && next.get_key()
                            == aminoAcids::RawAminoAcidImpl::PROTEIN_C_TERM)
                    || (site_.get_key()
                            == aminoAcids::RawAminoAcidImpl::PROTEIN_C_TERM
                            && current.get_key()
                                    == aminoAcids::RawAminoAcidImpl::PROTEIN_C_TERM))) {
                return false;
            }
            break;
        case ANYWHERE:
            if (!matchingSites) {
                return false;
            }
            break;
    }

    return true;
}

bool Specificity::operator==(const Specificity& s) const
{
    return site_ == s.site_ && classification_ == s.classification_
            && position_ == s.position_ && neutralLosses_ == s.neutralLosses_
            && pepNeutralLosses_ == s.pepNeutralLosses_
            && comment_ == s.comment_;
}

bool Specificity::operator!=(const Specificity& s) const
{
    return !(operator ==(s));
}

Specificity& Specificity::operator=(const Specificity& rhs)
{
    if (this != &rhs) {
        site_ = rhs.site_;
        classification_ = rhs.classification_;
        position_ = rhs.position_;
        neutralLosses_ = rhs.neutralLosses_;
        pepNeutralLosses_ = rhs.pepNeutralLosses_;
        comment_ = rhs.comment_;
    }
    return *this;
}

Specificity::Classification Specificity::parseClassificationString(
    const String& classification)
{
    String classification_tmp = classification;
    std::transform(classification_tmp.begin(), classification_tmp.end(),
        classification_tmp.begin(), ::tolower);

    if (classification_tmp == "-") {
        return Specificity::NONE;
    } else
        if (classification_tmp == "post-translational") {
            return Specificity::POST_TRANSLATIONAL;
        } else
            if (classification_tmp == "co-translational") {
                return Specificity::CO_TRANSLATIONAL;
            } else
                if (classification_tmp == "pre-translational") {
                    return Specificity::PRE_TRANSLATIONAL;
                } else
                    if (classification_tmp == "chemical derivative") {
                        return Specificity::CHEMICAL_DERIVATIVE;
                    } else
                        if (classification_tmp == "artefact") {
                            return Specificity::ARTEFACT;
                        } else
                            if (classification_tmp
                                    == "n-linked glycosylation") {
                                return Specificity::N_LINKED_GLYCOSYLATION;
                            } else
                                if (classification_tmp
                                        == "o-linked glycosylation") {
                                    return Specificity::O_LINKED_GLYCOSYLATION;
                                } else
                                    if (classification_tmp
                                            == "other glycosylation") {
                                        return Specificity::OTHER_GLYCOSYLATION;
                                    } else
                                        if (classification_tmp
                                                == "synth. pep. protect. gp.") {
                                            return Specificity::SYNTH_PEP_PROTECT_GP;
                                        } else
                                            if (classification_tmp
                                                    == "isotopic label") {
                                                return Specificity::ISOTOPIC_LABEL;
                                            } else
                                                if (classification_tmp
                                                        == "non-standard residue") {
                                                    return Specificity::NON_STANDARD_RESIDUE;
                                                } else
                                                    if (classification_tmp
                                                            == "multiple") {
                                                        return Specificity::MULTIPLE;
                                                    } else
                                                        if (classification_tmp
                                                                == "other") {
                                                            return Specificity::OTHER;
                                                        } else
                                                            if (classification_tmp
                                                                    == "aa substitution") {
                                                                return Specificity::AA_SUBSTITUTION;
                                                            }
    throw mstk::LogicError(
        "Specificity::parseClassificationString(): Given classification string does not represent a known classification.");
}

Specificity::Position Specificity::parsePositionString(
    const String& position)
{
    String position_tmp = position;
    std::transform(position_tmp.begin(), position_tmp.end(),
        position_tmp.begin(), ::tolower);
    if (position_tmp == "any n-term") {
        return Specificity::ANY_N_TERM;
    } else
        if (position_tmp == "any c-term") {
            return Specificity::ANY_C_TERM;
        } else
            if (position_tmp == "protein n-term") {
                return Specificity::PROTEIN_N_TERM;
            } else
                if (position_tmp == "protein c-term") {
                    return Specificity::PROTEIN_C_TERM;
                } else
                    if (position_tmp == "anywhere") {
                        return Specificity::ANYWHERE;
                    }
    throw mstk::LogicError(
        "Specificity::parsePositionString(): Given position string does not represent a known position.");
}

std::ostream& operator<<(std::ostream& os, const Specificity& s)
{
    os << s.getSite() << "\t" << s.getClassification() << "\t"
            << s.getPosition() << "\t" << s.getPosition() << "\t"
            << s.getNeutralLosses() << "\t" << s.getPepNeutralLosses() << "\t"
            << s.getComment();
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<Specificity>& s)
{
    typedef std::vector<Specificity>::const_iterator IT;
    for (IT it = s.begin(); it != s.end(); ++it) {
        os << *it << "|";
    }
    return os;
}

} // namespace modifications
} // namespace aas
} // namespace mstk
