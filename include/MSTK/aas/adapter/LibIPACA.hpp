/*
 * LibIPACA.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_ADAPTER_LIBIPACA_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_ADAPTER_LIBIPACA_HPP__

#include "MSTK/aas/Stoichiometry.hpp"
#include "MSTK/aas/Element.hpp"

#include "MSTK/common/Types.hpp"

#include "MSTK/ipaca/Spectrum.hpp"
#include "MSTK/ipaca/Stoichiometry.hpp"
#include "MSTK/ipaca/Traits.hpp"

#include <vector>

namespace mstk {
namespace aas {
namespace adapter {

/** @addtogroup mstk_aas
 * @{
 */

typedef ipaca::detail::Spectrum LibaasSpectrum;
typedef aas::stoichiometries::Stoichiometry LibaasStoichiometry;

struct SpectrumConverter
{
    void operator()(const ipaca::detail::Spectrum& lhs, LibaasSpectrum& rhs)
    {
        rhs = lhs;
    }
};

struct ElementConverter
{
    void operator()(const aas::elements::Element& lhs,
        ipaca::detail::Element& rhs)
    {
        rhs.count = 0.0;
        rhs.isotopes.clear();
        typedef std::vector<aas::elements::Isotope> AASIsotopeList;
        const AASIsotopeList& eis = lhs.get().getIsotopes();
        // TODO make sure isotopes are added sorted by mass!
        for (AASIsotopeList::const_iterator it = eis.begin(); it != eis.end();
                ++it) {
            ipaca::detail::Isotope i;
            i.ab = it->getFrequency();
            i.mz = it->getMass();
            rhs.isotopes.push_back(i);
        }
    }
};

struct StoichiometryConverter
{
    void operator()(const LibaasStoichiometry& lhs,
        ipaca::detail::Stoichiometry& rhs)
    {
        typedef LibaasStoichiometry::const_iterator SIT;
        typedef std::vector<aas::elements::Isotope>::const_iterator IIT;
        ElementConverter ec;
        for (SIT it = lhs.begin(); it != lhs.end(); ++it) {
            ipaca::detail::Element e;
            ec(it->first, e);
            e.count = it->second;
            rhs.push_back(e);
        }
    }
};

/** @\ */

} // namespace adapter
} // namespace aas

namespace ipaca {

/** @addtogroup mstk_aas
 * @{
 */

template<>
struct Traits<mstk::aas::adapter::LibaasStoichiometry, mstk::aas::adapter::LibaasSpectrum>
{
    typedef aas::adapter::SpectrumConverter spectrum_converter;
    typedef aas::adapter::StoichiometryConverter stoichiometry_converter;
    static ipaca::detail::Element getHydrogens(const Size n);
    static Bool isHydrogen(const ipaca::detail::Element&);
    static Double getElectronMass();
};

ipaca::detail::Element Traits<aas::adapter::LibaasStoichiometry,
        aas::adapter::LibaasSpectrum>::getHydrogens(const Size n)
{
    ipaca::detail::Element h;
    aas::adapter::ElementConverter ec;
    ec(aas::elements::Element(1), h);
    h.count = static_cast<Double>(n);
    return h;
}

Bool Traits<aas::adapter::LibaasStoichiometry,
        aas::adapter::LibaasSpectrum>::isHydrogen(
    const ipaca::detail::Element& e)
{
    ipaca::detail::Element h;
    aas::adapter::ElementConverter ec;
    ec(aas::elements::Element(1), h);
    typedef std::vector<ipaca::detail::Isotope>::const_iterator IT;
    for (IT et = e.isotopes.begin(), ht = h.isotopes.begin();
            et != e.isotopes.end() && ht != h.isotopes.end(); ++et, ++ht) {
        if (et->ab != ht->ab || et->mz != ht->mz) {
            return false;
        }
    }
    return true;
}

Double Traits<aas::adapter::LibaasStoichiometry,
        aas::adapter::LibaasSpectrum>::getElectronMass()
{
    return ipaca::detail::getElectronMass();
}

/** @\ */

} // namespace ipaca
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_ADAPTER_LIBIPACA_HPP__ */
