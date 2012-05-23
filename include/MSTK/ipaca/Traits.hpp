/*
 * Traits.hpp
 *
 *  Copyright (C) 2012 Marc Kirchner
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
#ifndef __MSTK_INCLUDE_MSTK_IPACA_TRAITS_HPP__
#define __MSTK_INCLUDE_MSTK_IPACA_TRAITS_HPP__

#include <MSTK/config.hpp>
#include <MSTK/ipaca/Stoichiometry.hpp>
#include <MSTK/common/Types.hpp>
#include <MSTK/common/Error.hpp>

#include <algorithm>

namespace mstk { 

namespace ipaca {

namespace detail {

// default convenience implementations
detail::Element getHydrogens(const Size n);
Bool isHydrogen(const detail::Element&);
Double getElectronMass();

}

template<typename StoichiometryType, typename SpectrumType>
struct Traits
{
    typedef typename StoichiometryType::converter stoichiometry_converter;
    typedef typename SpectrumType::converter spectrum_converter;
    static detail::Element getHydrogens(const Size n);
    static Bool isHydrogen(const detail::Element&);
    static Double getElectronMass();
};

/** Adjust the stoichiometry according to charge and charge type.
 * This adjusts the stoichiometry to reflect the changes induced by adding
 * charges and allows \c libipaca to keep all charge calculations out of the
 * distribution calculation logic. The function is a no-op for electron-based
 * charges, but adjusts the number of hydrogens in protonation-based charge
 * situations. Hence the isotopic distribution of the charge protons
 * contributes to the isotope pattern (as it should be).
 * @param s The uncharged stoichiometry.
 * @param charge The charge, signed.
 * @param p The type of charge mechanics used: either PROTON or ELECTRON.
 */
template<typename StoichiometryType, typename SpectrumType>
void adjustStoichiometryForProtonation(detail::Stoichiometry& s, const Int charge);

}

//
// template implementations
//
namespace ipaca {

namespace detail {

template<typename StoichiometryType, typename SpectrumType>
void adjustStoichiometryForProtonation(detail::Stoichiometry& s,
    const Int charge)
{
    // find the hydrogen entry
    typedef Stoichiometry::iterator IT;
    IT h = std::find_if(s.begin(), s.end(),
        Traits<StoichiometryType, SpectrumType>::isHydrogen);
    Double c = static_cast<Double>(charge);
    if (h != s.end()) {
        h->count = charge > 0 ? h->count + c : h->count - c;
        if (h->count < 0) {
            throw RuntimeError("Requested deprotonation but number of "
                "hydrogens is insufficient.");
        }
    } else {
        if (charge > 0) {
        s.push_back(Traits<StoichiometryType, SpectrumType>::getHydrogens(charge));
        } else {
            throw RuntimeError("Requested deprotonation but no hydrogens present.");
        }
    }
}

} // namespace detail

} // namespace ipaca

} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_IPACA_TRAITS_HPP__ */
