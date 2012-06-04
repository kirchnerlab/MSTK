/*
 * RawAminoAcid.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_RAWAMINOACID_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_RAWAMINOACID_HPP__

#include "MSTK/aas/RawAminoAcidImpl.hpp"
#include "MSTK/common/Types.hpp"

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>

namespace mstk {
namespace aas {
namespace aminoAcids {

/** @addtogroup mstk_aas
 * @{
 */

/**@brief ID extractor for raw amino acids.
 *
 * The class RawAminoAcidIdExtractor is used allow the instantiation of
 * flyweight<RawAminoAcidImpl>(Key) in order to simplify the access.
 */
struct RawAminoAcidIdExtractor
{
    /**Returns the key of the raw amino acid.
     * @param[in] a instance of an raw amino acid implementation
     * @returns The key of the raw amino acid
     */
    const RawAminoAcidImpl::RawAminoAcidImplKeyType& operator()(
        const RawAminoAcidImpl& a) const
    {
        return a.getId();
    }
};

/**Typedef to simplify the data type flyweight<RawAminoAcidImpl>
 */
typedef boost::flyweight<
        boost::flyweights::key_value<RawAminoAcidImpl::RawAminoAcidImplKeyType,
                RawAminoAcidImpl, RawAminoAcidIdExtractor>
        , boost::flyweights::no_tracking> RawAminoAcid;

/**Convenience function to add a custom raw amino acid to this list of known amino
 * acids.
 * This methods calls addAminoAcid(RawAminoAcidImpl)
 * @param[in] id Key/Id of the amino acid
 * @param[in] symbol Symbol of the amino acid
 * @param[in] threeLetterCode Three letter code of the amino acid
 * @param[in] fullName Full name of the amino acid
 * @param[in] stoichiometry Stoichiometry of the amino acid
 * @returns True if the given amino acid is added correctly, false otherwise.
 */
mstk::Bool addRawAminoAcid(const RawAminoAcidImpl::RawAminoAcidImplKeyType& id,
    const mstk::Char symbol, const mstk::String& threeLetterCode, const mstk::String& fullName,
    const aas::stoichiometries::Stoichiometry& stoichiometry);

/**Convenience function to add a custom raw amino acid to the list of known amino
 * acids.
 * Note: Once an amino acid is added it is not possible to alter its properties.
 * In order to change an internal property, you have to create a new amino acid
 * with the correct values.
 * IT IS NOT POSSIBLE TO OVERRIDE A REFERENCE (flyweight).
 * @param[in] aa Instance of an RawAminoAcidImpl
 * @returns True if the given amino acid is added correctly, false otherwise.
 */
mstk::Bool addRawAminoAcid(const RawAminoAcidImpl& aa);

bool operator<(const RawAminoAcid& lhs, const RawAminoAcid& rhs);
bool operator<=(const RawAminoAcid& lhs, const RawAminoAcid& rhs);
bool operator>(const RawAminoAcid& lhs, const RawAminoAcid& rhs);
bool operator>=(const RawAminoAcid& lhs, const RawAminoAcid& rhs);

/** @\ */

} // namespace aminoAcids
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_RAWAMINOACID_HPP__ */
