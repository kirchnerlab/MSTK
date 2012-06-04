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

#ifndef __MSTK_INCLUDE_MSTK_AAS_AMINOACID_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_AMINOACID_HPP__

#include "MSTK/aas/RawAminoAcid.hpp"
#include "MSTK/aas/StoichiometryConfig.hpp"
#include "MSTK/aas/Stoichiometry.hpp"
#include "MSTK/common/Types.hpp"

#include <iostream>

namespace mstk {
namespace aas {
namespace aminoAcids {

/** @addtogroup mstk_aas
 * @{
 */

/**Representation of an amino acid using aas::aminoAcids::RawAminoAcid.
 */
class AminoAcid
{
public:

    /**Constructor.
     * @param[in] aminoAcidKey The key/id of the raw amino acid
     * @param[in] configid The stoichiometry config key/id
     * @throws Throwsn an aas::errors::LogicError exception in case the given amino
     * acid key cannot be resolved.
     */
    AminoAcid(
        const RawAminoAcidImpl::RawAminoAcidImplKeyType& aminoAcidKey = '\0',
        const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configid =
                aas::stoichiometries::StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG);

    /**Constructor.
     * @param[in] aminoAcid The raw amino acid
     * @param[in] config The stoichiometry config which is used to calculate the
     * stoichiometry
     */
    AminoAcid(
        const RawAminoAcid& aminoAcid,
        const aas::stoichiometries::StoichiometryConfig& config =
                aas::stoichiometries::StoichiometryConfig(
                    aas::stoichiometries::StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG));

    /**Returns the symbol of the amino acid.
     * @returns The symbol of the amino acid.
     */
    mstk::Char getSymbol() const;

    /**Returns the key/id of the raw amino acid.
     * @returns The key/id of the raw amino acid.
     */
    const RawAminoAcidImpl::RawAminoAcidImplKeyType& getRawAminoAcidKey() const;

    /**Returns the raw amino acid.
     * @returns The raw amino acid.
     */
    const RawAminoAcid& getRawAminoAcid() const;

    /**Returns the three letter code of the amino acid.
     * @return The three letter code of the amino acid
     */
    const mstk::String& getThreeLetterCode() const;

    /**Returns the full name of the amino acid.
     * @returns The full name of the amino acid
     */
    const mstk::String& getFullName() const;

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

    /**Calculates and returns a copy of the stoichiometry of this amino acid.
     * The calculation is skipped in case the present stoichiometry configuration is
     * StoichiometryConfigImpl::DEFAULT_ELEMENT_CONFIG.
     * This method calls recalculatesWithConfiguration() on the stoichiometry retrieved
     * by the raw amino acid.
     * @returns The stoichiometry calculated with the present stoichiometry configuration
     * @throws Throws an aas::errors::RuntimeError in case one or more elements cannot
     * be resolved by the stochiometry config.
     */
    aas::stoichiometries::Stoichiometry getStoichiometry() const;

    /**Sets the stoichiometry configuration.
     * Note: This will trigger the recalculation of the stoichiometry using the
     * given configuration.
     * @param[in] config Stoichiometry configuration
     * @throws Throws an exception in case one ore more elements cannot be
     * resolved by the stoichiometry configuration.
     */
    void setStoichiometryConfig(
        const aas::stoichiometries::StoichiometryConfig& config);

    /**Sets the stoichiometry configuration.
     * Note: This will trigger the recalculation of the stoichiometry using the
     * given configuration.
     * @param[in] configid Id of the stoichiometry configuration
     * @throws Throws an exception in case one ore more elements cannot be
     * resolved by the stoichiometry configuration.
     */
    void
    setStoichiometryConfig(
        const aas::stoichiometries::StoichiometryConfigImpl::StoichiometryConfigImplKeyType& configid);

    /**Returns the stoichiometry configuration.
     * @returns Stoichiometry configuration
     */
    const aas::stoichiometries::StoichiometryConfig& getStoichiometryConfig() const;

    /**Sets a copy of the argument as the new content for the amino acid object.
     * The previous content is dropped.
     * @param[in] a Amino acid to copy
     * @returns *this
     */
    AminoAcid& operator=(const AminoAcid& a);

    /**Compares the amino acid against another.
     * @param[in] a Amino acid object to compare *this with
     * @returns true if both amino acids are the same, false otherwise
     */
    bool operator==(const AminoAcid& a) const;

    /**Compares the amino acid against another, with opposite result of
     * RawAminoAcidImpl::operator==.
     * @param[in] a Amino acid object to compare *this with
     * @returns true if the amino acids are different, false otherwise.
     */
    bool operator!=(const AminoAcid& a) const;

private:

    /**The raw amino acid.
     */
    RawAminoAcid rawAminoAcid_;
    /**The stoichiometry configuration used to calculate the stoichiometry.
     */
    aas::stoichiometries::StoichiometryConfig stoichiometryConfig_;

};

std::ostream& operator<<(std::ostream& os, const AminoAcid& a);

/** @} */

} // namespace aminoAcids
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_AMINOACID_HPP__ */
