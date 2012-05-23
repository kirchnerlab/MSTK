/*
 * Mercury7.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_IPACA_MERCURY7_HPP__
#define __MSTK_INCLUDE_MSTK_IPACA_MERCURY7_HPP__
#include <MSTK/config.hpp>
#include <MSTK/ipaca/Mercury7Impl.hpp>
#include <MSTK/ipaca/Spectrum.hpp>
#include <MSTK/ipaca/Stoichiometry.hpp>
#include <MSTK/common/Types.hpp>
#include <MSTK/ipaca/Traits.hpp>
#include <boost/shared_ptr.hpp>

namespace mstk {
    
namespace ipaca {

/** Calculates a theoretical isotope distribution from an
 *  elemental composition (stoichiometry).
 *
 *  The client is required to provide type information for the
 *  stoichiometry and the resulting spectra as well as two conversion
 *  function in \c ipaca::Traits<U,V> that specifiy how the user-
 *  defined stoichiometry and spectrum types are converted into
 *  the internale representation used by ipaca (i.e. \c ipaca::Spectrum
 *  and \c ipaca::Stoichiometry).
 *  See \c test/Mercury7-test.cpp for a simple example.
 */
template<typename StoichiometryType, typename SpectrumType>
class Mercury7
{
public:
    /** The type of particle that carries the charge.
     */
    enum Particle
    {
        ELECTRON, PROTON
    };

    /** Constructor.
     */
    Mercury7();

    /** Functor method to calculate the theoretical isotope
     *         distribution of a compound.
     * @param stoichiometry The stoichiometry for which the isotope
     *                      distribution should be calculated.
     * @param limit The abundance limit below which peaks are pruned
     *              during the processing
     *
     * The procedure is based on Perttu Haimi's and Alan Rockwood's
     * sparse/binary convolution algorithm.
     */
    SpectrumType
    operator()(const StoichiometryType& stoichiometry, const int charge,
        const Particle particle, const Double limit = 1e-26) const;

    /** calculate the monoisotopic mass of a given stoichiometry
     *  @param stoichiometry The stoichiometry to calculate the mass for.
     *  @param charge The charge at which the monoisotopic mass is desired
     *                (zero for theoretical but unobservable mass).
     */
    Double
    getMonoisotopicMass(const StoichiometryType& stoichiometry) const;

    /** calculate the average mass of a given stoichiometry
     *  @param stoichiometry The stoichiometry to calculate the mass for.
     *  @param charge The charge at which the average mass is desired
     *                (zero for theoretical but unobservable mass).
     */
    Double getAverageMass(const StoichiometryType& stoichiometry) const;
private:
    boost::shared_ptr<detail::Mercury7Impl> pImpl_;
};

//
// template implementation
//

template<typename StoichiometryType, typename SpectrumType>
Mercury7<StoichiometryType, SpectrumType>::Mercury7() :
    pImpl_(new detail::Mercury7Impl)
{
}

template<typename StoichiometryType, typename SpectrumType>
SpectrumType Mercury7<StoichiometryType, SpectrumType>::operator()(
    const StoichiometryType& stoichiometry, const int charge,
    const Particle particle, const Double limit) const
{
    // convert the user type to our internal type
    detail::Stoichiometry s;
    typename Traits<StoichiometryType, SpectrumType>::stoichiometry_converter
            stoi_conv;
    stoi_conv(stoichiometry, s);
    // adjust the stoichiometry for charge and particle type


    // Adjust the number of hydrogens.
    if (charge != 0 && particle == PROTON) {
        detail::adjustStoichiometryForProtonation<StoichiometryType, SpectrumType>(s, charge);
    }
    detail::Spectrum result = pImpl_->operator()(s, limit);
    // Do the charge adjustment. This is the same for all types of charges
    // because we adjusted the number of hydrogens earlier.
    if (charge != 0) {
        Int absCharge = (abs)(charge);
        Double e = Traits<StoichiometryType, SpectrumType>::getElectronMass();
        typedef detail::Spectrum::iterator IT;
        for (IT i = result.begin(); i != result.end(); ++i) {
            i->mz = (i->mz - (charge * e)) / absCharge;
        }
    }
    SpectrumType spectrum;
    typename Traits<StoichiometryType, SpectrumType>::spectrum_converter
            spec_conv;
    spec_conv(result, spectrum);
    return spectrum;
}

template<typename StoichiometryType, typename SpectrumType>
Double Mercury7<StoichiometryType, SpectrumType>::getMonoisotopicMass(
    const StoichiometryType& stoichiometry) const
{
    detail::Stoichiometry s;
    typename Traits<StoichiometryType, SpectrumType>::stoichiometry_converter
            stoi_conv;
    stoi_conv(stoichiometry, s);
    return pImpl_->getMonoisotopicMass(s);
}

template<typename StoichiometryType, typename SpectrumType>
Double Mercury7<StoichiometryType, SpectrumType>::getAverageMass(
    const StoichiometryType& stoichiometry) const
{
    detail::Stoichiometry s;
    typename Traits<StoichiometryType, SpectrumType>::stoichiometry_converter
            stoi_conv;
    stoi_conv(stoichiometry, s);
    return pImpl_->getAverageMass(s);
}

} // namespace ipaca

} // namespace mstk

#endif
