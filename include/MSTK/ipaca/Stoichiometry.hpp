/*
 * Stoichiometry.hpp
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
#ifndef __LIBIPACAP_INCLUDE_MSTK_IPACA_STOICHIOMETRY_HPP__
#define __LIBIPACAP_INCLUDE_MSTK_IPACA_STOICHIOMETRY_HPP__

#include <MSTK/config.hpp>
#include <vector>
#include <iosfwd>
#include <MSTK/common/Types.hpp>

namespace mstk {

namespace ipaca {

namespace detail {

/*
 * Keep things as simple as possible.
 */

/** An isotope.
 */
struct Isotope
{
    Double mz, ab;
};

/** An isotope distribution.
 */
typedef std::vector<Isotope> Isotopes;

/** A mass spectrum is the same as an isotope distribution.
 */
typedef Isotope SpectrumElement;
typedef std::vector<SpectrumElement> Spectrum;

/** The relevant element information in a stoichiometry. We only need
 * the isotopic distribution and the number of occurences.
 */
struct Element
{
    Isotopes isotopes;
    Double count;
};

/** A stoichiometry is simply a list of elements with their isotopic
 * distributions and number of occurences.
 */
typedef std::vector<Element> Stoichiometry;

//
// a few free functions
//

/** Check if all entries in a stoichiometry are non-negative and if there is
 * at least one entry that is positive. Here, 'plausible' means that
 * a meaningful calculation is possible, not necessarily that the stoichiometry
 * is chemically feasible (e.g O_32 would be considered plausible, but O_{-3}
 * and H_{0}C_{0} would not).
 * @param s The \c detail::Stoichiometry object to test.
 * @return A boolean, true if the stoichiometry is plausible.
 */
Bool isPlausibleStoichiometry(const Stoichiometry& s);

/** Split a stoichiometry into integer and fractional contributions.
 *
 * @param s The \c Stoichiometry to split.
 * @param intStoi All integer contributions to the stoichiometry.
 * @param fracStoi All fractional contributions.
 */
void splitStoichiometry(const Stoichiometry& s, Stoichiometry& intStoi,
    Stoichiometry& fracStoi);

/** Stream operator.
 *
 * @param os A reference to the outstream.
 * @param s The Stoichiometry to stream.
 * @return A reference to the same outstream.
 */
std::ostream& operator<<(std::ostream& os, const Stoichiometry& s);

} // namespace detail

} // namespace ipaca

} // namespace mstk

#endif /* __LIBIPACAP_INCLUDE_MSTK_IPACA_STOICHIOMETRY_HPP__ */
