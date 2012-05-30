/*
 * Specificity.hpp
 *
 * Copyright (c) 2011 Mathias Wilhelm
 * Copyright (c) 2010,2011 Marc Kirchner
 * Copyright (c) 2010 Nathan Huesken
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_SPECIFICITY_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_SPECIFICITY_HPP__

#include "MSTK/aas/RawAminoAcid.hpp"
#include "MSTK/common/Types.hpp"

#include <vector>
#include <iostream>

namespace mstk {
namespace aas {
namespace modifications {

/** Representation of the specificity.
 */
class Specificity
{

public:

    /**Available positions.
     */
    enum Position
    {
        ANY_N_TERM = 0, ANY_C_TERM, PROTEIN_N_TERM, PROTEIN_C_TERM, ANYWHERE
    };

    /** Available classifications.
     */
    enum Classification
    {
        NONE = 0,
        POST_TRANSLATIONAL,
        CO_TRANSLATIONAL,
        PRE_TRANSLATIONAL,
        CHEMICAL_DERIVATIVE,
        ARTEFACT,
        N_LINKED_GLYCOSYLATION,
        O_LINKED_GLYCOSYLATION,
        OTHER_GLYCOSYLATION,
        SYNTH_PEP_PROTECT_GP,
        ISOTOPIC_LABEL,
        NON_STANDARD_RESIDUE,
        MULTIPLE,
        OTHER,
        AA_SUBSTITUTION
    };

    /**Default constructor to create a specificity.
     * @param[in] site Site at which an event can happen
     * @param[in] position Position of the specificity
     * @param[in] classification Classification of the specificity
     */
    Specificity(const aas::aminoAcids::RawAminoAcid& site,
        const Position position, const Classification classification);

    /**Convenience constructor.
     *
     * @param[in] site Site at which an event can happen
     * @param[in] position Position of the specificity
     * @param[in] classification Classification of the specificity
     * @throws Throws an exception if one of the given strings cannot be parsed
     * correctly.
     */
    Specificity(const mstk::String& site, const mstk::String& position,
        const mstk::String& classification);

    /**Sets the site.
     * @param[in] aminoAcid
     */
    void setSite(const aas::aminoAcids::RawAminoAcid& aminoAcid);

    /**Returns the site.
     * @returns
     */
    const aas::aminoAcids::RawAminoAcid& getSite() const;

    /**Sets the classification of the specificity.
     * @param[in] classification Classification
     */
    void setClassification(const Classification& classification);

    /**Returns the classification of the specificity
     * @returns The classification
     */
    const Classification& getClassification() const;

    /**Sets the position of the specificity.
     * @param[in] position Position
     */
    void setPosition(const Position& position);

    /**Returns the position of the specificity
     * @returns The position
     */
    const Position& getPosition() const;

    /**Adds a neutral loss.
     * @param[in] stoichiometry Stoichiometry of the neutral loss
     */
    void addNeutralLoss(const aas::stoichiometries::Stoichiometry& stoichiometry);

    /**Sets the list of neutral losses.
     * @param[in] stoichiometry List of stoichiometries
     */
    void setNeutralLosses(const std::vector<aas::stoichiometries::Stoichiometry>& stoichiometry);

    /**Returns the list of neutral losses.
     * @returns List of neutral losses
     */
    const std::vector<aas::stoichiometries::Stoichiometry>& getNeutralLosses() const;

    /**Clears the list of neutral losses.
     */
    void clearNeutralLosses();

    /**Adds a peptide neutral loss
     * @param[in] stoichiometry Stoichiometry of the peptide neutral loss
     */
    void addPepNeutralLoss(const aas::stoichiometries::Stoichiometry& stoichiometry);

    /**Sets the list of peptide neutral losses
     * @param[in] stoichiometries List of stoichiometries
     */
    void setPepNeutralLosses(
        const std::vector<aas::stoichiometries::Stoichiometry>& stoichiometries);

    /**Returns all peptide neutral losses.
     * @returns List of stoichiometries
     */
    const std::vector<aas::stoichiometries::Stoichiometry>& getPepNeutralLosses() const;

    /**Clears the list of peptide neutral losses.
     */
    void clearPepNeutralLosses();

    /**Sets the comment.
     * @param[in] comment Comment
     */
    void setComment(const mstk::String& comment);

    /**Returns the comment.
     * @returns The comment
     */
    const mstk::String& getComment() const;

    /**Checks whether this specificity matches the surrounding conditions of the
     * amino acid.
     * @param[in] prev Previous amino acid
     * @param[in] current Amino acid which should be modified with this
     * modification
     * @param[in] next Next amino acid
     * @returns true if the modification is applicable to the current amino acid,
     * false otherwise
     */
    mstk::Bool isApplicable(const aas::aminoAcids::RawAminoAcid& prev,
        const aas::aminoAcids::RawAminoAcid& current,
        const aas::aminoAcids::RawAminoAcid& next) const;

    /**Sets a copy of the argument as the new content for the specificity object.
     * The previous content is dropped.
     * @param[in] rhs Specificity to copy
     * @returns *this
     */
    Specificity& operator=(const Specificity& rhs);

    /**Compares the specificity against another.
     * @param[in] s Specificity object to compare *this with
     * @returns true if both specificites are the same, false otherwise
     */
    bool operator==(const Specificity& s) const;

    /**Compares the specificity against another, with opposite result of
     * Specificity::operator==.
     * @param[in] s Specificity object to compare *this with
     * @returns true if the specificities are different, false otherwise.
     */
    bool operator!=(const Specificity& s) const;

    /**Converts the given string to lower case and tries to find the matching
     * enumeration.
     * @param[in] classification Classification string
     * @returns Classification
     * @throws Throws an exception if the given classification string does not
     * match any of the classification enumeration.
     */
    static Classification parseClassificationString(
        const mstk::String& classification);

    /**Converts the given string to lower case and tries to find the matching
     * position.
     * @param[in] position Position string
     * @returns Positions
     * @throws Throws an exception if the given position string does not match
     * any of the position enumerations.
     */
    static Position parsePositionString(const mstk::String& position);

private:

    /** Site at which an event can happen.
     */
    aas::aminoAcids::RawAminoAcid site_;
    /** Position of the site.
     */
    Position position_;
    /** Classification of the specificity.
     */
    Classification classification_;
    /** List of possible neutral losses.
     */
    std::vector<aas::stoichiometries::Stoichiometry> neutralLosses_;
    /** List of possible peptide neutral losses.
     */
    std::vector<aas::stoichiometries::Stoichiometry> pepNeutralLosses_;
    /** A comment.
     */
    mstk::String comment_;

};
// class Specificity

std::ostream& operator<<(std::ostream&, const Specificity&);
std::ostream& operator<<(std::ostream&, const std::vector<Specificity>&);

} // namespace modifications
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_SPECIFICITY_HPP__ */
