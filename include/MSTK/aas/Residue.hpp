/*
 * Residue.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_RESIDUE_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_RESIDUE_HPP__

#include "MSTK/aas/AminoAcid.hpp"
#include "MSTK/aas/Modification.hpp"
#include "MSTK/common/Types.hpp"

#include <boost/shared_ptr.hpp>

#include <iostream>

namespace mstk {
namespace aas {

/**Representation of a residue.
 *
 * A residue holds an amino acid and its modification.
 */
class Residue
{

public:

    /** Creates a residue with a modification
     *
     * Note: the constructor does not check whether the modification is
     * applicable to this amino acid.
     *
     * @param[in] aminoAcidKey The key/id of an amino acid
     * @param[in] modificationKey The key/id of the modification
     * @param[in] labelKey The key/id of the label
     * @throws Throws an aas::errors::LogicError exception if given modification or label key are no modifications
     * which can be added as modification or label.
     */
    Residue(
        const aas::aminoAcids::RawAminoAcidImpl::RawAminoAcidImplKeyType& aminoAcidKey =
                '\0',
        const aas::modifications::RawModificationImpl::RawModificationImplKeyType& modificationKey =
                "",
        const aas::modifications::RawModificationImpl::RawModificationImplKeyType& labelKey =
                "");

    /**Creates a residue with a modification.
     *
     * Note: the constructor does not check whether the modification is
     * applicable to this amino acid.
     *
     * @param[in] aminoAcid The amino acid
     * @param[in] mod The modification
     * @param[in] label The isotopic label
     * @throws Throws an aas::errors::LogicError exception if given modification or label key are no modifications
     * which can be added as modification or label.
     */
    Residue(const aas::aminoAcids::AminoAcid& aminoAcid,
        const aas::modifications::Modification& mod =
                modifications::Modification(),
        const aas::modifications::Modification& label =
                modifications::Modification());

    /**Change type of the amino acid.
     * @param[in] aminoAcidKey The key of the amino acid
     *
     * TODO the modifications stays the same. shall we remove it?
     */
    void changeType(
        const aminoAcids::RawAminoAcidImpl::RawAminoAcidImplKeyType& aminoAcidKey);

    /**Change type of the amino acid.
     * @param[in] aminoAcid The amino acid
     *
     * TODO the modifications stays the same. shall we remove it?
     */
    void changeType(const aminoAcids::AminoAcid& aminoAcid);

    /**Returns the amino acid
     * @returns The amino acid
     */
    const aas::aminoAcids::AminoAcid& getAminoAcid() const;

    /**Returns a modifiable reference to the amino acid of this residue.
     * @return A reference to the amino acid.
     */
    aas::aminoAcids::AminoAcid& getAminoAcid();

    /**Checks whether the amino acid is N-terminal.
     * @returns True if the amino acid is AminoAcidImpl::PEPTIDE_N_TERM or
     * AminoAcidImpl::PROTEIN_N_TERM.
     */
    mstk::Bool isNTerm() const;

    /**Checks whether the amino acid is C-terminal.
     * @returns True if the amino acid is AminoAcidImpl::PEPTIDE_C_TERM or
     * AminoAcidImpl::PEPTIDE_C_TERM.
     */
    mstk::Bool isCTerm() const;

    /**Sets the modification.
     *
     * Calls aas::Residue::setModification(Modification)
     *
     * @param[in] modificationKey The key of a modification
     */
    void setModification(
        const aas::modifications::RawModificationImpl::RawModificationImplKeyType& modificationKey);

    /**Sets the modification.
     *
     * Note: this method does not check whether the modification is applicable
     * to this position.
     *
     * @param[in] modification The modifiaction
     * @throws Throws an aas::error::LogicError if the given modification is an isotopic label.
     */
    void setModification(const aas::modifications::Modification& modification);

    /**Returns the modification of this residue
     * @returns The modification
     */
    const modifications::Modification& getModification() const;

    /**Sets the isotopic label.
     *
     * calls aas::Residue::setIsotopicLabel(Modification)
     *
     * @param[in] isotopicLabelKey The key/id of the isotopic label.
     */
    void setIsotopicLabel(
        const aas::modifications::RawModificationImpl::RawModificationImplKeyType& isotopicLabelKey);

    /**Sets the isotopic label.
     *
     * Note: this method does not check whether the isotopic label is applicable
     * to this position.
     *
     * @param[in] isotopicLabel The isotopic label
     * @throws Throws an aas::errors::LogicError exception if the given isotopicLabel
     * is a standard modification.
     */
    void setIsotopicLabel(
        const aas::modifications::Modification& isotopicLabel);

    /**Returns the isotopic label of this residue
     * @returns The isotopic label
     */
    const modifications::Modification& getIsotopicLabel() const;

    /**Checks whether the modification id is equal to the given key.
     *
     * This check does not include internal properties of the modification,
     * such as custom specificities or the stoichiometry configuration
     *
     * @param[in] modificationKey
     * @returns true if modification_.get_key() == modificationKey, false otherwise
     */
    mstk::Bool hasModification(
        const modifications::RawModificationImpl::RawModificationImplKeyType& modificationKey) const;

    /**Checks whether the modification is equal to the given one.
     * @param[in] modification Modification
     * @returns true if modification_ == modification, false otherwise
     */
    mstk::Bool hasModification(
        const modifications::Modification& modification) const;

    /**Checks whether the residue is modified.
     * @returns true if the residue is modified, false otherwise
     */
    mstk::Bool isModified() const;

    /**Checks whether the isotopic label id is equal to the given key.
     *
     * This check does not include internal properties of the isotopic label,
     * such as custom specificities or the stoichiometry configuration
     *
     * @param[in] labelKey
     * @returns true if isotopicLabel_.get_key() == labelKey, false otherwise
     */
    mstk::Bool hasLabel(
        const modifications::RawModificationImpl::RawModificationImplKeyType& labelKey) const;

    /**Checks whether the isotopic label is equal to the given one.
     * @param[in] label Isotopic label
     * @returns true if isotopicLabel_ == label, false otherwise
     */
    mstk::Bool hasLabel(const modifications::Modification& label) const;

    /**Checks whether the residue is modified.
     * @returns true if the residue is modified, false otherwise
     */
    mstk::Bool isLabeled() const;

    /**Removes the modification.
     */
    void removeModification();

    /**Removes the isotopic label.
     */
    void removeIsotopicLabel();

    /**Sets the stoichiometry configuration key of the amino acid.
     * @param[in] configKey Stoichiometry configuration key
     */
    void applyAminoAcidStoichiometryConfig(
        const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configKey);

    /**Sets the stoichiometry configuration of the amino acid.
     * @param[in] config Stoichiometry configuration
     */
    void applyAminoAcidStoichiometryConfig(const aas::stoichiometries::StoichiometryConfig& config);

    /**Sets the stoichiometry configuration key of the modification.
     * @param[in] configKey Stoichiometry configuration key
     */
    void applyModificationStoichiometryConfig(
        const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configKey);

    /**Sets the stoichiometry configuration of the modification.
     * @param[in] config Stoichiometry configuration
     */
    void applyModificationStoichiometryConfig(
        const aas::stoichiometries::StoichiometryConfig& config);

    /**Sets the stoichiometry configuration key of the isotopic label.
     * @param[in] configKey Stoichiometry configuration key
     */
    void applyIsotopicLabelStoichiometryConfig(
        const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configKey);

    /**Sets the stoichiometry configuration of the isotopic label.
     * @param[in] config Stoichiometry configuration
     */
    void applyIsotopicLabelStoichiometryConfig(
        const aas::stoichiometries::StoichiometryConfig& config);

    /**Returns the stoichiometry of the residue.
     * Note: This method actually calculates the stoichiometry by summing up the
     * stoichiometry of the amino acid, the modification and the isotopic label.
     *
     * @returns Stoichiometry of the amino acid plus its modification
     */
    aas::stoichiometries::Stoichiometry getStoichiometry() const;

    /**Returns a human readable string containing the symbol of the amino acid
     * the, if present, its modification in parenthesis, i.e. C(Oxidation)
     * @returns human readable string containing the amino acid symbol and its
     * modification
     */
    mstk::String toString() const;

    /**Sets a copy of the argument as the new content for the residue object.
     * The previous content is dropped.
     * @param[in] rhs Residue to copy
     * @returns *this
     */
    Residue& operator=(const Residue& rhs);

    /**Compares the residue against another.
     * @param[in] r Residue object to compare *this with
     * @returns true if both residues are the same, false otherwise
     */
    bool operator==(const Residue& r) const;

    /**Compares the residue against another, with opposite result of
     * Residue::operator==.
     * @param[in] r Residue object to compare *this with
     * @returns true if the residues are different, false otherwise.
     */
    bool operator!=(const Residue& r) const;

private:

    /**Convenience typedef.
     */
    typedef boost::shared_ptr<modifications::Modification> ModificationPtr;

    /**Private reference to an empty modification, used to keep memory consumption
     * minimal in case of an unmodified and unlabeled residue.
     */
    static ModificationPtr EMPTY_MOD;

    /** The amino acid.
     */
    aas::aminoAcids::AminoAcid aminoacid_;
    /** The modification of the amino acid.
     */
    ModificationPtr modification_;
    /** The isotopic label of the amino acid.
     */
    ModificationPtr isotopicLabel_;

};
// class Residue

std::ostream& operator<<(std::ostream&, const Residue&);

} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_RESIDUE_HPP__ */
