/*
 * AminoAcidSequence.cpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2011,2012 Marc Kirchner
 * Copyright (c) 2010 Nathan Huesken
 * Copyright (c) 2009 Marc Kirchner
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

#include "MSTK/aas/AminoAcidSequence.hpp"
#include "MSTK/common/Error.hpp"

#include <sstream>

namespace mstk {
namespace aas {

AminoAcidSequence::AminoAcidSequence(const String& aminoAcidSequence,
    const aas::stoichiometries::StoichiometryConfig& aminoAcidConfig)
{
    // TODO do we really want an empty aas in case the string is empty?
    if (!aminoAcidSequence.empty()) {
        // prepend peptide n-term if sequence starts without a n-term
        if (!aminoAcids::AminoAcid(aminoAcidSequence[0]).isNTerm()) {
            c_.push_back(aminoAcids::RawAminoAcidImpl::PEPTIDE_N_TERM);
        }

        typedef String::const_iterator IT;
        IT end = aminoAcidSequence.end();
        Residue r;
        for (IT it = aminoAcidSequence.begin(); it != end; ++it) {
            c_.push_back(*it);
            c_.back().applyAminoAcidStoichiometryConfig(aminoAcidConfig);
        }

        // append a peptide c-term if the sequence starts without a c-term
        if (!aminoAcids::AminoAcid(
            aminoAcidSequence[aminoAcidSequence.size() - 1]).isCTerm()) {
            c_.push_back(aminoAcids::RawAminoAcidImpl::PEPTIDE_C_TERM);
        }
    } else {
        c_.push_back(aminoAcids::RawAminoAcidImpl::PEPTIDE_N_TERM);
        c_.push_back(aminoAcids::RawAminoAcidImpl::PEPTIDE_C_TERM);
    }
}

AminoAcidSequence::AminoAcidSequence(const_iterator first, const_iterator last)
{
    // When the original sequence has no N and/or C term, we add peptide C and N terms as default
    if (!(*first).isNTerm()) {
        c_.push_back(aminoAcids::RawAminoAcidImpl::PEPTIDE_N_TERM);
    }
    std::copy(first, last, std::back_inserter(c_));
    if (!c_[c_.size() - 1].isCTerm()) {
        c_.push_back(aminoAcids::RawAminoAcidImpl::PEPTIDE_C_TERM);
    }
}

void AminoAcidSequence::push_back(const Residue& value)
{
    // MAYBE optimize by storing n and c-term as member variables (will result in checked [], might be to slow)
    //If a C-terminal is passed, the previous C-terminal is replaced by it
    Residue last(
        aminoAcids::AminoAcid(aminoAcids::RawAminoAcidImpl::PEPTIDE_C_TERM));
    if (size() == 0) {
        if (!value.isNTerm()) {
            c_.push_back(aminoAcids::RawAminoAcidImpl::PEPTIDE_N_TERM);
        }
        c_.push_back(value);
    } else {
        if (c_[size() - 1].isCTerm()) {
            //copy the C-terminal, in case there are mods
            last = c_[size() - 1];
            c_.pop_back();
            c_.push_back(value);
        } else {
            c_.push_back(value);
        }
    }
    if (!value.isCTerm())
        c_.push_back(last);
}

void AminoAcidSequence::pop_back()
{
    // TODO shall we prepend and append an N and C term in case they are missing?
    // any more boundary conditions?
    Size csize = c_.size();
    if (csize > 0) {
        if (c_[csize - 1].isCTerm()) {
            if (csize > 1 && !c_[csize - 2].isNTerm()) {
                Residue last = c_[csize - 1];
                c_.pop_back();
                c_.pop_back();
                c_.push_back(last);
            }
        } else
            if (!c_[csize - 1].isNTerm()) {
                c_.pop_back();
                c_.push_back(aminoAcids::RawAminoAcidImpl::PEPTIDE_C_TERM);
            }
    }
}

void AminoAcidSequence::makePeptideCTerm()
{
    if (size() == 0 || !c_.back().isCTerm()) {
        // TODO shall we append a peptide c-term in this case?
        mstk_fail("Unable to change amino acid sequence C-term"
        "to peptide C-term, because there is no C-term.");
    }
    if (c_.back().getAminoAcid().getRawAminoAcidKey()
            == aminoAcids::RawAminoAcidImpl::PEPTIDE_C_TERM) {
        return;
    }
    c_.back().changeType(aminoAcids::RawAminoAcidImpl::PEPTIDE_C_TERM);
    // TODO keep stoich config of c term?
}

void AminoAcidSequence::makePeptideNTerm()
{
    if (size() == 0 || !c_.front().isNTerm()) {
        // TODO shall we prepend a peptide N-term in this case?
        mstk_fail("Unable to change amino acid sequence N-term"
        "to protein N-term, because there is no N-term.");
    }
    if (c_.front().getAminoAcid().getRawAminoAcidKey()
            == aminoAcids::RawAminoAcidImpl::PEPTIDE_N_TERM) {
        return;
    }
    c_.front().changeType(aminoAcids::RawAminoAcidImpl::PEPTIDE_N_TERM);
    // TODO keep stoich config of n term?
}

void AminoAcidSequence::makeProteinCTerm()
{
    if (size() == 0 || !c_.back().isCTerm()) {
        // TODO shall we append a protein c-term in this case?
        mstk_fail("Unable to change amino acid sequence C-term"
        "to protein C-term, because there is no C-term.");
    }
    if (c_.back().getAminoAcid().getRawAminoAcidKey()
            == aminoAcids::RawAminoAcidImpl::PROTEIN_C_TERM) {
        return;
    }
    c_.back().changeType(aminoAcids::RawAminoAcidImpl::PROTEIN_C_TERM);
    // TODO keep stoich config of c term?
}

void AminoAcidSequence::makeProteinNTerm()
{
    if (size() == 0 || !c_.front().isNTerm()) {
        // TODO shall we prepend a protein n-term in this case?
        mstk_fail("Unable to change amino acid sequence N-term"
        "to protein N-term, because there is no N-term.");
    }
    if (c_.front().getAminoAcid().getRawAminoAcidKey()
            == aminoAcids::RawAminoAcidImpl::PROTEIN_N_TERM) {
        return;
    }
    c_.front().changeType(aminoAcids::RawAminoAcidImpl::PROTEIN_N_TERM);
    // TODO keep stoich config of n term?
}

void AminoAcidSequence::remove(
    const modifications::RawModificationImpl::RawModificationImplKeyType& modKey)
{
    for (AminoAcidSequence::iterator it = begin(); it != end(); ++it) {
        if (it->hasModification(modKey)) {
            it->removeModification();
        }
    }
}

void AminoAcidSequence::remove(const modifications::Modification& mod)
{
    for (AminoAcidSequence::iterator it = begin(); it != end(); ++it) {
        if (it->hasModification(mod)) {
            it->removeModification();
        }
    }
}

void AminoAcidSequence::append(const AminoAcidSequence& sequence)
{
    if (!sequence.empty()) {
        if (this != &sequence) {
            //The C-terminal is taken from the second sequence, so the
            //possible modifications of the C-terminal of the first sequence
            //get lost here
            if (c_.size() == 0) {
                // if we have a n-term mod at sequence, should we just use the sequence[0] instead of a new peptide n-term here?
                if (sequence.c_.front().isNTerm()) {
                    c_.push_back(sequence.c_.front());
                } else {
                    c_.push_back(aminoAcids::RawAminoAcidImpl::PEPTIDE_N_TERM);
                    c_.push_back(sequence[0]);
                }
            }
            if (c_[c_.size() - 1].isCTerm()) {
                c_.pop_back();
            }
            std::copy(sequence.begin() + 1, sequence.end(),
                std::back_inserter(c_));
        } else {
            Residue last(
                aminoAcids::AminoAcid(
                    aminoAcids::RawAminoAcidImpl::PEPTIDE_C_TERM));
            if (c_[c_.size() - 1].isCTerm()) {
                last = c_.back();
                c_.pop_back();
            }
            std::copy(sequence.begin() + 1, sequence.end(),
                std::back_inserter(c_));
            c_.push_back(last);
        }
    }
}

void AminoAcidSequence::applyFixedModifications(
    const std::vector<
            modifications::RawModificationImpl::RawModificationImplKeyType>& mods)
{
    ModificationList mmods;
    typedef std::vector<
            modifications::RawModificationImpl::RawModificationImplKeyType>::const_iterator IT;
    IT end = mods.end();
    for (IT iter = mods.begin(); iter != end; ++iter) {
        mmods.push_back(modifications::Modification(*iter));
    }
    applyFixedModifications(mmods);
}

void AminoAcidSequence::applyFixedModifications(
    const std::vector<modifications::RawModification>& mods)
{
    ModificationList mmods;
    typedef std::vector<modifications::RawModification>::const_iterator IT;
    IT end = mods.end();
    for (IT iter = mods.begin(); iter != end; ++iter) {
        mmods.push_back(modifications::Modification(*iter));
    }
    applyFixedModifications(mmods);
}

void AminoAcidSequence::applyFixedModifications(const ModificationList& mods)
{
    typedef ModificationList::const_iterator IT;
    IT end = mods.end();
    for (IT it = mods.begin(); it != end; ++it) {
        // iterate through all modifications and apply them to the complete peptide
        Size end = size() - 1;
        const modifications::Modification& mod = *it;
        for (Size pos = 1; pos < end; ++pos) {
            try {
                applyModificationAtPosition(mod, pos);
            } catch (mstk::Exception& e) {
                // nothing to do here
            	// TODO or is it?
            }
        }
    }
}

void AminoAcidSequence::applyModificationAtPosition(
    const modifications::Modification& mod, const Size& pos)
{
    // TODO we implicitly set a label or modification if the modification "thinks" it is a label. shall we use the return code of the getSpecificities instead?
    mstk_precondition(
        pos < c_.size() && pos > 0,
        "AminoAcidSequence::applyModificationAtPosition(): Trying to apply modification at position out of bound.");

    if (mod.isIsotopicLabel() && c_[pos].isLabeled()) {
        mstk_fail(
            "AminoAcidSequence::applyModificationAtPosition(): Trying to apply label modification at position which is already labeled.");
    }
    if (!mod.isIsotopicLabel() && c_[pos].isModified()) {
        mstk_fail(
            "AminoAcidSequence::applyModificationAtPosition(): Trying to apply modification at position which is already modified.");
    }

    // get prev amino acid
    aminoAcids::AminoAcid prev('\0');
    if (pos != 0) {
        prev = operator[](pos - 1).getAminoAcid();
    }

    // get next amino acid
    aminoAcids::AminoAcid next('\0');
    if (pos != size() - 1) {
        next = operator[](pos + 1).getAminoAcid();
    }

    // get pos amino acid
    Residue current = operator[](pos);
    // check if mod is applicable to this position
    if (!mod.isApplicable(prev, current.getAminoAcid(), next)) {
        mstk_fail(
            "AminoAcidSequence::applyModificationAtPosition(): Cannot apply mod to this position.");
    }

    // in case the modification is an isotopic label, we will set the isotopic label
    if (mod.isIsotopicLabel()) {
        operator[](pos).setIsotopicLabel(mod);
    } else {
        operator[](pos).setModification(mod);
    }
}

void AminoAcidSequence::applyAminoAcidStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& aminoAcidConfigKey)
{
    applyAminoAcidStoichiometryConfig(aas::stoichiometries::StoichiometryConfig(aminoAcidConfigKey));
}

void AminoAcidSequence::applyAminoAcidStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfig& aminoAcidConfig)
{
    for (iterator it = begin(); it != end(); ++it) {
        it->applyAminoAcidStoichiometryConfig(aminoAcidConfig);
    }
}

void AminoAcidSequence::applyModificationStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& modificationConfigKey)
{
    applyModificationStoichiometryConfig(
        aas::stoichiometries::StoichiometryConfig(modificationConfigKey));
}

void AminoAcidSequence::applyModificationStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfig& modificationConfig)
{
    for (iterator it = begin(); it != end(); ++it) {
        it->applyModificationStoichiometryConfig(modificationConfig);
    }
}

void AminoAcidSequence::applyIsotopicLabelStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& labelConfigKey)
{
    applyIsotopicLabelStoichiometryConfig(aas::stoichiometries::StoichiometryConfig(labelConfigKey));
}

void AminoAcidSequence::applyIsotopicLabelStoichiometryConfig(
    const aas::stoichiometries::StoichiometryConfig& labelConfig)
{
    for (iterator it = begin(); it != end(); ++it) {
        it->applyIsotopicLabelStoichiometryConfig(labelConfig);
    }
}

aas::stoichiometries::Stoichiometry AminoAcidSequence::getStoichiometry() const
{
    // MAYBE optimize by storing the stoichiometry
    aas::stoichiometries::Stoichiometry ret;
    for (const_iterator iter = begin(); iter != end(); ++iter) {
        ret += iter->getStoichiometry();
    }
    return ret;
}

String AminoAcidSequence::toUnmodifiedSequenceString() const
{
    std::string s;
    for (const_iterator it = begin(); it != end(); ++it) {
        if (!it->isNTerm() && !it->isCTerm())
            s += it->getAminoAcid().getSymbol();
    }
    return s;
}

String AminoAcidSequence::toString(const Bool showTerminals) const
{
    std::string seq;
    for (const_iterator it = begin(); it != end(); ++it) {
        if (showTerminals || (!(*it).isNTerm() && !(*it).isCTerm()))
            seq += (*it).toString();
    }
    return seq;
}
String AminoAcidSequence::getModificationString() const
{
    std::ostringstream modoss;
    int pos = 0;
    for (const_iterator it = begin(); it != end(); ++it, ++pos) {
        if (it->isModified()) {
            if (!modoss.str().empty()) {
                modoss << "; ";
            }
            modoss << it->getModification().getModificationId() << "("
                    << it->getAminoAcid().getSymbol() << ")@" << pos;
        }
        if (it->isLabeled()) {
            if (!modoss.str().empty()) {
                modoss << "; ";
            }
            modoss << it->getIsotopicLabel().getModificationId() << "("
                    << it->getAminoAcid().getSymbol() << ")@" << pos;
        }
    }
    return modoss.str();
}

std::ostream& operator<<(std::ostream& os, const AminoAcidSequence& o)
{
    typedef AminoAcidSequence::const_iterator IT;
    IT end = o.end();
    for (IT it = o.begin(); it != end; ++it) {
        os << *it << "\t";
    }
    return os;
}

} // namespace aas
} // namespace mstk
