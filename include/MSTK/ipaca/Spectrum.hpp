/*
 * Spectrum.hpp
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
#ifndef __LIBIPACAP_INCLUDE_MSTK_IPACA_SPECTRUM_HPP__
#define __LIBIPACAP_INCLUDE_MSTK_IPACA_SPECTRUM_HPP__

#include <MSTK/config.hpp>
#include <vector>
#include <MSTK/ipaca/Stoichiometry.hpp>

namespace mstk {

namespace ipaca {

namespace detail {

/** A mass spectrum is the same as an isotope distribution.
 */
typedef Isotope SpectrumElement;
typedef std::vector<SpectrumElement> Spectrum;

/** A stream operator for the Spectrum class.
 *
 * @param os The stream.
 * @param s The spectrum.
 * @return The stream after adding the spectrum info.
 */
std::ostream& operator<<(std::ostream& os, const Spectrum& s);

} // namespace detail

} // namespace ipaca

} // namespace mstk

#endif /* __LIBIPACAP_INCLUDE_MSTK_IPACA_SPECTRUM_HPP__ */
