/*
 * StoichiometryConfig.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_STOICHIOMETRYCONFIG_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_STOICHIOMETRYCONFIG_HPP__

#include "MSTK/aas/StoichiometryConfigImpl.hpp"
#include "MSTK/common/Types.hpp"

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>

namespace mstk {
namespace aas {
namespace stoichiometries {

/** @addtogroup mstk_aas
 * @{
 */

/**@brief ID extractor for stoichiometry configurations.
 *
 * The class StoichiometryConfigIdExtractor is used allow the instantiation of
 * flyweight<StoichiometryConfigImpl>(Key) in order to simplify the access.
 */
struct StoichiometryConfigIdExtractor
{
    /**Returns the key of the stoichiometry configuration.
     * @param[in] e instance of a stoichiometry configuration implementation
     * @returns The key of the stoichiometry configuration
     */
    const StoichiometryConfigImpl::StoichiometryConfigImplKeyType& operator()(
        const StoichiometryConfigImpl& e) const
    {
        return e.getId();
    }
};

/**Typedef to simplify the data type flyweight<StoichiometryConfigImpl>
 */
typedef boost::flyweight<
        boost::flyweights::key_value<
                StoichiometryConfigImpl::StoichiometryConfigImplKeyType,
                StoichiometryConfigImpl, StoichiometryConfigIdExtractor>
        , boost::flyweights::no_tracking> StoichiometryConfig;

/**Convenience function to add a custom stoichiometry configuration to the list of
 * known stoichiometry configurations.
 * This method creates a StoichiometryConfigImpl and calls
 * addStoichiometryConfig(stoichiometryConfig) to add it.
 * @param[in] id Key/Id of the stoichiometry configuration
 * @param[in] map Mapping of an element symbol to an element id
 * @returns True if the given amino acid is added correctly, false otherwise.
 */
mstk::Bool
addStoichiometryConfig(
    const StoichiometryConfigImpl::StoichiometryConfigImplKeyType& id,
    const StoichiometryConfigImpl::DataType& map);

/**Convenience function to add a custom stoichiometry configuration to the list of
 * known stoichiometry configurations.
 * Note: Once a configuration is added it is not possible to alter its properties.
 * In order to change an internal property, you have to create a new amino acid
 * with the correct values.
 * @param[in] stoichiometryConfig Stoichiometry configuration
 * @returns True if the given amino acid is added correctly, false otherwise.
 */
mstk::Bool addStoichiometryConfig(
    const StoichiometryConfigImpl& stoichiometryConfig);

bool operator<(const StoichiometryConfig& lhs, const StoichiometryConfig& rhs);
bool operator<=(const StoichiometryConfig& lhs,
    const StoichiometryConfig& rhs);
bool operator>(const StoichiometryConfig& lhs, const StoichiometryConfig& rhs);
bool operator>=(const StoichiometryConfig& lhs,
    const StoichiometryConfig& rhs);

/** @\ */

} // namespace stoichiometries
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_STOICHIOMETRYCONFIG_HPP__ */
