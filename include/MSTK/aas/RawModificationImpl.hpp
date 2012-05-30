/*
 * RawModificationImpl.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_RAWMODIFICATIONIMPL_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_RAWMODIFICATIONIMPL_HPP__

#include "MSTK/aas/Stoichiometry.hpp"
#include "MSTK/aas/Specificity.hpp"
#include "MSTK/common/Types.hpp"
#include "MSTK/aas/AminoAcid.hpp"

#include <vector>
#include <iostream>

namespace mstk {
namespace aas {
namespace modifications {

/** Representation of a modification.
 *
 * This class is used as an internal data holder to allow the usage of flyweights.
 */
class RawModificationImpl
{

public:

    /**Convenience typedef.
     */
    typedef mstk::String RawModificationImplKeyType;

    /**Default constructor to instantiate a standard modification.
     * @param[in] id Key/Id of the modification
     * @throws Throws an aas::errors::LogicError exception if the given id is not in
     * the list of standard modifications.
     */
    RawModificationImpl(const RawModificationImplKeyType& id);

    /**Constructor to create custom modifications.
     * @param[in] id Key/Id of the modification
     * @param[in] name Name of the modification
     * @param[in] fullName Full name of the modification (description)
     * @param[in] verified Verified
     */
    RawModificationImpl(
        const RawModificationImpl::RawModificationImplKeyType& id,
        const mstk::String& name, const mstk::String& fullName,
        const mstk::Bool& verified);

    /**Retruns the id of the modification.
     * @returns Id of the modification
     */
    const RawModificationImplKeyType& getId() const;

    /**Sets the name of the modification.
     * @param[in] name Name
     */
    void setName(const mstk::String& name);

    /**Returns the name of the modification.
     * @returns Name
     */
    const mstk::String& getName() const;

    /**Sets the full name (description) of the modification.
     * @param[in] fullName Full name
     */
    void setFullName(const mstk::String& fullName);

    /**Returns the full name (description) of the modification.
     * @returns Full name
     */
    const mstk::String& getFullName() const;

    /**Adds an alternative name.
     * @param[in] altName Alternative name
     */
    void addAltName(const mstk::String& altName);

    /**Sets alternative names.
     * @param[in] altNames Alternative names
     */
    void setAltNames(const std::vector<mstk::String>& altNames);

    /**Returns all alternative names.
     * @returns List of alternative names
     */
    const std::vector<mstk::String>& getAltNames() const;

    /**Sets the stoichiometry of the modification
     * @param[in] stoichiometry Stoichiometry
     */
    void setStoichiometry(const aas::stoichiometries::Stoichiometry& stoichiometry);

    /**Returns the stoichiometry of the modification.
     *
     * In case this is a standard modification, the stoichiometry is calculated
     * by using the default stoichiometry configuration.
     *
     * @returns The stoichiometry of the modification
     */
    const aas::stoichiometries::Stoichiometry& getStoichiometry() const;

    /**Adds a specificity to the modification.
     * @param[in] specificity Specificity
     */
    void addSpecificity(const Specificity& specificity);

    /**Sets specificities.
     * @param[in] specificities List of specificities
     */
    void setSpecificities(const std::vector<Specificity>& specificities);

    /**Returns all specificities of the modification.
     * @returns List of specificities
     */
    const std::vector<Specificity>& getSpecificities() const;

    /**Set the status of the modification.
     * @param[in] verified Verfified
     */
    void setVerified(mstk::Bool verified);

    /**Returns the status of the modification.
     * @returns True if the modification is verified, otherwise false.
     */
    mstk::Bool isVerified() const;

    /**Checks whether this raw modification is applicable to the amino acid
     * current.
     * This function iterates over all specificities and calls
     * Specificity::isApplicable.
     * @param[in] prev Previous amino acid
     * @param[in] current Amino acid which should be modified with this
     * modification
     * @param[in] next Next amino acid
     * @returns true if the modification is applicable to the current amino acid,
     * false otherwise
     */
    mstk::Bool isApplicable(const aminoAcids::AminoAcid& prev,
        const aminoAcids::AminoAcid& current,
        const aminoAcids::AminoAcid& next) const;

    /**Sets a copy of the argument as the new content for the raw modification
     * object.
     * The previous content is dropped.
     * @param[in] rhs Raw modification to copy
     * @returns *this
     */
    RawModificationImpl& operator=(const RawModificationImpl& rhs);

    /**Compares the raw modification against another.
     * @param[in] m Raw moficifation object to compare *this with
     * @returns true if both raw modifications are the same, false otherwise
     */
    bool operator==(const RawModificationImpl& m) const;

    /**Compares the raw modification against another, with opposite result of
     * RawModificationImpl::operator==.
     * @param[in] m Raw modification object to compare *this with
     * @returns true if the raw modifications are different, false otherwise.
     */
    bool operator!=(const RawModificationImpl& m) const;

private:

    /** Key/Id of the raw modification.
     */
    RawModificationImplKeyType id_;
    /** Name of the raw modification.
     */
    mstk::String name_;
    /** Full name or description of the raw modification.
     */
    mstk::String fullName_;
    /** List of alternative names of the raw modification.
     */
    std::vector<mstk::String> altNames_;
    /** Stoichiometry of the raw modification
     */
    aas::stoichiometries::Stoichiometry stoichiometry_;
    /** List of specificities for the raw modification.
     */
    std::vector<Specificity> specificities_;
    /** Status of the raw modification.
     */
    mstk::Bool verified_;

    // unimod also contains references and a note
};
// class RawModificationImpl

std::ostream& operator<<(std::ostream&, const RawModificationImpl&);

inline const RawModificationImpl::RawModificationImplKeyType& RawModificationImpl::getId() const
{
    return id_;
}

} // namespace modifications
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_RAWMODIFICATIONIMPL_HPP__ */
