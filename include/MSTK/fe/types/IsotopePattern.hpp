/*
 * IsotopePattern.hpp
 *
 * Copyright (C) 2011 Marc Kirchner
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
#ifndef __MSTK_INCLUDE_MSTK_FE_TYPES_ISOTOPEPATTERN_HPP__
#define __MSTK_INCLUDE_MSTK_FE_TYPES_ISOTOPEPATTERN_HPP__

#include <MSTK/config.hpp>
#include <MSTK/common/Collection.hpp>
#include <MSTK/fe/types/Xic.hpp>
#include <MSTK/fe/types/Spectrum.hpp>

#include <functional>
#include <vector>
#include <iostream>
#include <set>

namespace mstk {

namespace fe {

/* A representation for an isotope pattern.
 */
class IsotopePattern : public Collection<Xic>
{
public:
    /** Constructor
     */
    IsotopePattern();

    /** Constructor.
     * Constructs an isotope pattern from a set of Xics.
     */
    IsotopePattern(const_iterator first, const_iterator last);

    /** Set the isotope pattern charge state information.
     *
     * @param[in] z The charge states for the isotope pattern.
     */
    void setCharges(const std::set<int>& z);

    /** Get the isotope pattern charge state information.
     *
     * @return The set of charge states associated with the isotope pattern.
     */
    const std::set<int>& getCharges() const;

    /** Convert the isotope pattern in a Spectrum representation.
     *
     * @param[out] ss The isotope pattern in the form of a Spectrum.
     */
    void asSpectrum(Spectrum& ss) const;

    /** Split the current isotope pattern into multiple isotope patterns,
     *  i.e. unmix by charge state.
     * @param isotopePatterns A vector of single-charge (not singly-charged)
     *                        isotope patterns.
     */
    void split(std::vector<IsotopePattern>& isotopePatterns);

    /** Calculate the overall abundance of the isotope pattern
     * @return The abundance (the sum of all XIC abundances).
     */
    double getAbundance() const;

private:
    /** The set of charges associated with the isotope pattern.
     */
    std::set<int> charges_;
};

/** Output operator for IsotopePatterns.
 *
 * @param os Reference to the output stream.
 * @param p The isotope pattern to be streamed.
 * @return Reference to the output stream into which p has been streamed.
 */
std::ostream& operator<<(std::ostream& os, const IsotopePattern& p);

} // namespace fe

} // namespace mstk

#endif

