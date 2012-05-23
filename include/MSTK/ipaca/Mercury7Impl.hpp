/*
 * Mercury7Impl.hpp
 *
 * Copyright (c) 2007-2011 Marc Kirchner
 * Copyright (c) 2008 Thorben Kroeger
 * Copyright (c) 2007 Xinghua Lou
 * Copyright (c) 2007 Bjoern Voss
 * Based on the emass implementation of P. Haimi (see Copyright notice below)
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
/*
 * Copyright (c) 2005 Perttu Haimi and Alan L. Rockwood
 *
 * All rights reserved.
 *
 * Based on an algorithm developed by Alan L. Rockwood.
 *
 * Published in
 * Rockwood, A.L. and Haimi, P.: "Efficent calculation of
 * Accurate Masses of Isotopic Peaks",
 * Journal of The American Society for Mass Spectrometry
 * JASMS 03-2263, 2006
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted provided
 * that the following conditions are met:
 *
 *    * Redistributions of source code must retain the
 *      above copyright notice, this list of conditions
 *      and the following disclaimer.
 *    * Redistributions in binary form must reproduce
 *      the above copyright notice, this list of conditions
 *      and the following disclaimer in the documentation
 *      and/or other materials provided with the distribution.
 *    * Neither the author nor the names of any contributors
 *      may be used to endorse or promote products derived
 *      from this software without specific prior written
 *      permission.
 */

/*
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __MSTK_INCLUDE_MSTK_IPACA_MERCURY7IMPL_HPP__
#define __MSTK_INCLUDE_MSTK_IPACA_MERCURY7IMPL_HPP__
#include <MSTK/config.hpp>
#include <MSTK/ipaca/Spectrum.hpp>
#include <MSTK/ipaca/Stoichiometry.hpp>
#include <MSTK/common/Types.hpp>
#include <vector>
#include <exception>
#include <stdexcept>

namespace mstk {

namespace ipaca {

namespace detail {
/** Calculates a theoretical isotope distribution from an
 *  elemental composition (stoichiometry).
 *  @ingroup asap
 */
class Mercury7Impl
{
public:
    /** Functor method to calculate the theoretical isotope
     *         distribution of a compound.
     * @param stoichiometry The stoichiometry for which the isotope
     *                      distribution should be calculated.
     * @param charge The charge of the isotope distribution
     * @param limit The abundance limit below which peaks are pruned
     *              during the processing
     * @param includeMonoisotopicMass For larger molecules, the monoisotopic
     *              mass may have a normalized abundance (probability) that
     *              falls under the pruning threshold. Thes flag controls if the
     *              monoisotopic mass will always be present in the result of
     *              not. If set to true (the default), Mercury will re-add the
     *              monoisotopic mass if it is lost on the way.
     * @throws mstk::OutOfRange Tried to access an incomplete isotope table.
     *
     * The procedure is based on Perttu Haimi's and Alan Rockwood's
     * sparse/binary convolution algorithm.
     */
    detail::Spectrum
    operator()(const detail::Stoichiometry& stoichiometry,
        const Double limit = 1e-26) const;

    /** calculate the monoisotopic mass of a given stoichiometry
     *  @param stoichiometry The stoichiometry to calculate the mass for.
     *  @param charge The charge at which the monoisotopic mass is desired
     *                (zero for theoretical but unobservable mass).
     */
    Double
    getMonoisotopicMass(const detail::Stoichiometry& stoichiometry) const;

    /** calculate the average mass of a given stoichiometry
     *  @param stoichiometry The stoichiometry to calculate the mass for.
     *  @param charge The charge at which the average mass is desired
     *                (zero for theoretical but unobservable mass).
     */
    Double getAverageMass(const detail::Stoichiometry& stoichiometry) const;

private:
    /** Calculate the theoretical isotope distribution of a compound
     * of integer stoichiometries.
     */
    void integerMercury(const detail::Stoichiometry& stoichiometry,
        const Double limit, detail::Spectrum& spectrum) const;

    /** Calculate the theoretical isotope distribution of a compound
     * of fractional stoichiometries.
     */
    void fractionalMercury(const detail::Stoichiometry& stoichiometry,
        const Double limit, detail::Spectrum& spectrum) const;

    /** Convolves two isotope distributions.
     * @param s1 Spectrum on the left hand side of the convolution.
     * @param s2 Spectrum on the right hand side of the convolution.
     * @param result The result of the convolution.
     */
    void convolve(const detail::Spectrum& s1, const detail::Spectrum& s2,
        detail::Spectrum& result) const;

    /** Prunes sparse isotope distributions based on the observed intensities.
     * @param spectrum A \c detail::Spectrum object.
     * @param limit The (relative) abundance limit below which isotope peaks
     *              are discarded.
     *
     * Discards all entries in the mass and abundance vectors whose abundance is
     * below the abundance limit.
     */
    void prune(detail::Spectrum& spectrum, const Double limit) const;
};

} // namespace detail

} // namespace ipaca

} // namespace mstk

#endif
