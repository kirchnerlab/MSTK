/*
 * RawModification.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_RAWMODIFICATION_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_RAWMODIFICATION_HPP__

#include "MSTK/aas/RawModificationImpl.hpp"
#include "MSTK/common/Types.hpp"

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>

namespace mstk {
namespace aas {
namespace modifications {

/**@brief ID extrator for raw modifications.
 *
 * The class RawModificationIdExtractor is used allow the instantiation of
 * flyweight<RawModificationImpl>(Key) in order to simplify the access.
 */
struct RawModificationIdExtractor
{
    /**Returns the key of the raw modification.
     * @param[in] m instance of a raw modification implementation
     * @returns The key of the raw modification
     */
    const RawModificationImpl::RawModificationImplKeyType& operator()(
        const RawModificationImpl& m) const
    {
        return m.getId();
    }
};

/**Typedef to simplify the data type flyweight<RawModificationImpl>
 */
typedef boost::flyweight<
        boost::flyweights::key_value<
                RawModificationImpl::RawModificationImplKeyType,
                RawModificationImpl, RawModificationIdExtractor>
        , boost::flyweights::no_tracking> RawModification;

/**Convenience function to add a custom raw modification to the list of known
 * raw modifications.
 * @param[in] id Key/Id of the raw modification
 * @param[in] name Name of the raw modification
 * @param[in] fullname Full name (description) of the raw modification
 * @param[in] altNames Alternative names of the raw modification
 * @param[in] stoichiometry Stoichiometry of the modification
 * @param[in] specificities specificities of the modification
 * @param[in] verified Raw modification is verified
 * @returns True if the given raw modification is added correctly, false otherwise.
 */
mstk::Bool addRawModification(
    const RawModificationImpl::RawModificationImplKeyType& id,
    const mstk::String& name, const mstk::String& fullname,
    const std::vector<mstk::String>& altNames, const aas::stoichiometries::Stoichiometry& stoichiometry,
    const std::vector<Specificity>& specificities, const mstk::Bool& verified);

/**Convenience function to add a custom raw modification to the list of known
 * raw modifications.
 * Note: Once an amino acid is added it is not possible to alter its properties.
 * In order to change an internal property, you have to create a new amino acid
 * with the correct values.
 * IT IS NOT POSSIBLE TO OVERRIDE A REFERENCE (flyweight).
 * @param[in] rawModification Instance of a RawModificationImpl
 * @returns True if the given raw modification is added correctly, false otherwise.
 */
mstk::Bool addRawModification(const RawModificationImpl& rawModification);

bool operator<(const RawModification& lhs, const RawModification& rhs);
bool operator<=(const RawModification& lhs, const RawModification& rhs);
bool operator>(const RawModification& lhs, const RawModification& rhs);
bool operator>=(const RawModification& lhs, const RawModification& rhs);

} // namespace modifications
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_RAWMODIFICATION_HPP__ */
