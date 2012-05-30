/*
 * RawAminoAcidImpl.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_RAWAMINOACIDIMPL_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_RAWAMINOACIDIMPL_HPP__

#include "MSTK/aas/Stoichiometry.hpp"
#include "MSTK/common/Types.hpp"

#include <iostream>

namespace mstk {
namespace aas {
namespace aminoAcids {

/**Representation of an amino acid.
 *
 * This class only represents the amino acid without its modifications.
 * Furthermore, this class is usually not used directly, since the appropriate
 * flyweight AminoAcid provides the necessary constructors and functions to
 * use amino acids in a memory friendly manner.
 */
class RawAminoAcidImpl
{

public:

    /**Convenience typedef for the key type.
     */
    typedef mstk::Char RawAminoAcidImplKeyType;

    /**Default constructor to create a standard amino acid.
     * In order to create a custom amino acid, use
     * RawAminoAcidImpl(id, symbol, stoichiometry).
     *
     * This constructor has access to static variables containing all
     * information necessary to create and fill the internal variables.
     *
     * @param[in] id Key/Id of the amino acid.
     * @throws Throws an aas::errors::LogicError exception in case the id is not in
     * the list of standard amino acids.
     */
    RawAminoAcidImpl(const RawAminoAcidImplKeyType& id = '\0');

    /**Cosntructor to create a custom amino acid.
     * @param[in] id Key/Id of the amino acid
     * @param[in] symbol Symbol of the amino acid
     * @param[in] stoichiometry Stoichiometry of the amino acid
     */
    RawAminoAcidImpl(const RawAminoAcidImplKeyType& id, const char symbol,
        const aas::stoichiometries::Stoichiometry& stoichiometry);

    /**Returns the key of the amino acid.
     * @returns Key of the amino acid.
     */
    const RawAminoAcidImplKeyType& getId() const;

    /**Sets the symbol of the amino acid.
     * @param[in] symbol Symbol of the amino acid
     */
    void setSymbol(const mstk::Char& symbol);

    /**Returns the symbol of the amino acid
     * @returns Symbol of the amino acid.
     */
    mstk::Char getSymbol() const;

    /**Sets the stoichiometry of the amino acid.
     * @param[in] stoichiometry Stoichiometry of the amino acid
     */
    void setStoichiometry(const aas::stoichiometries::Stoichiometry& stoichiometry);

    /**Returns the stoichiometry of the amino acid.
     * @returns stoichiometry of the amino acid.
     */
    const aas::stoichiometries::Stoichiometry& getStoichiometry() const;

    /**Sets the three letter code for the amino acid.
     * @param[in] threeLetterCode The three letter code of the amino acid
     */
    void setThreeLetterCode(const mstk::String& threeLetterCode);

    /**Returns the three letter code of the amino acid.
     * @return The three letter code of the amino acid
     */
    const mstk::String& getThreeLetterCode() const;

    /**Sets the full name of the amino acid.
     * @param[in] fullName Full name of the amino acid
     */
    void setFullName(const mstk::String& fullName);

    /**Returns the full name of the amino acid.
     * @returns The full name of the amino acid
     */
    const mstk::String& getFullName() const;

    /**Checks whether the amino acid is N-terminal.
     * @returns True if the amino acid is RawAminoAcidImpl::PEPTIDE_N_TERM or
     * RawAminoAcidImpl::PROTEIN_N_TERM.
     */
    mstk::Bool isNTerm() const;

    /**Checks whether the amino acid is C-terminal.
     * @returns True if the amino acid is RawAminoAcidImpl::PEPTIDE_C_TERM or
     * RawAminoAcidImpl::PEPTIDE_C_TERM.
     */
    mstk::Bool isCTerm() const;

    /**Sets a copy of the argument as the new content for the raw amino acid object.
     * The previous content is dropped.
     * @param[in] a Raw amino acid to copy
     * @returns *this
     */
    RawAminoAcidImpl& operator=(const RawAminoAcidImpl& a);

    /**Compares the raw amino acid against another.
     * @param[in] a Raw amino acid object to compare *this with
     * @returns true if both raw amino acids are the same, false otherwise
     */
    bool operator==(const RawAminoAcidImpl& a) const;

    /**Compares the raw amino acid against another, with opposite result of
     * RawAminoAcidImpl::operator==.
     * @param[in] a Raw amino acid object to compare *this with
     * @returns true if the raw amino acids are different, false otherwise.
     */
    bool operator!=(const RawAminoAcidImpl& a) const;

    /** Key for the artificial peptide N-term amino acid.
     */
    static const mstk::Char PEPTIDE_N_TERM;
    /** Key for the artificial protein N-term amino acid.
     */
    static const mstk::Char PROTEIN_N_TERM;
    /** Key for the artificial peptide C-term amino acid.
     */
    static const mstk::Char PEPTIDE_C_TERM;
    /** Key for the artificial protein C-term amino acid.
     */
    static const mstk::Char PROTEIN_C_TERM;

    /**Returns the corresponding key of the amino acid string aminoAcid. The
     * string can contain the one-letter code, three-letter code and full
     * name of the amino acid.
     *
     * Note: n-term will translate to PEPTIDE_N_TERM
     *       c-term will translate to PEPTIDE_C_TERM
     *
     * @param[in] aminoAcid String representing the amino acid.
     * @returns Corresponding key of the amino acid with the string aminoAcid
     * @throws Throws an exception if the given amino acid string is not in the
     * list of standard amino acids.
     */
    static RawAminoAcidImplKeyType getKeyForAminoAcidString(
        const mstk::String& aminoAcid);

private:

    /** Key/Id of the amino acid.
     */
    RawAminoAcidImplKeyType id_;
    /** One letter symbol of the amino acid.
     */
    mstk::Char symbol_;
    /** Three letter code of the amino acid.
     */
    mstk::String threeLetterCode_;
    /** Full name of the amino acid.
     */
    mstk::String fullName_;
    /** Stoichiometry of the amino acid.
     */
    aas::stoichiometries::Stoichiometry stoichiometry_;

};
// class RawAminoAcidImpl

std::ostream& operator<<(std::ostream&, const RawAminoAcidImpl&);

inline const RawAminoAcidImpl::RawAminoAcidImplKeyType& RawAminoAcidImpl::getId() const
{
    return id_;
}

} // namespace aminoAcids
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_RAWAMINOACIDIMPL_HPP__ */
