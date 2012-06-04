/*
 * Isotope.hpp
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_ISOTOPE_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_ISOTOPE_HPP__

#include "MSTK/common/Types.hpp"

#include <vector>
#include <iostream>

namespace mstk {
namespace aas {
namespace elements {

/** @addtogroup mstk_aas
 * @{
 */

/**Representation of an isotope.
 */
class Isotope
{

public:
    /** Default constructor.
     * @param[in] mass Mass
     * @param[in] frequency Frequency
     */
    Isotope(const mstk::Double& mass, const mstk::Double& frequency);

    /**Returns the mass of the isotope.
     * @returns Mass of the isotope.
     */
    const mstk::Double& getMass() const;

    /**Returns the frequency of the isotope.
     * @returns Frequency of the isotope.
     */
    const mstk::Double& getFrequency() const;

    /**Sets a copy of the argument as the new content for the isotope object.
     * The previous content is dropped.
     * @param[in] rhs Isotope to copy
     * @returns *this
     */
    Isotope& operator=(const Isotope& rhs);

    /**Compares the isotope against another.
     * @param[in] i Isotope object to compare *this with
     * @returns true if both isotopes are the same, false otherwise
     */
    bool operator==(const Isotope& i) const;

    /**Compares the isotope against another, with opposite result of
     * Isotope::operator==.
     * @param[in] i Isotope object to compare *this with
     * @returns true if the isotopes are different, false otherwise.
     */
    bool operator!=(const Isotope& i) const;

private:

    /** Mass of the isotope.
     */
    mstk::Double mass_;
    /** Frequency of the isotope.
     */
    mstk::Double frequency_;

};
// class Isotope

inline const mstk::Double& Isotope::getMass() const
{
    return mass_;
}

inline const mstk::Double& Isotope::getFrequency() const
{
    return frequency_;
}

std::ostream& operator<<(std::ostream&, const Isotope&);
std::ostream& operator<<(std::ostream&, const std::vector<Isotope>&);

/** @\ */

} // namespace elements
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_ISOTOPE_HPP__ */
