/*
 * Traits.cpp
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
#include <MSTK/ipaca/Traits.hpp>
#include <MSTK/ipaca/Stoichiometry.hpp>

namespace mstk {

namespace ipaca {

namespace detail {

Element getHydrogens(const Size n)
{
    // function static object is initialized on first pass
    Element e;
    if (e.isotopes.empty()) {
        detail::Isotope i1, i2;
        i1.mz = 1.007825;
        i1.ab = 0.99985;
        i2.mz = 2.01410178;
        i2.ab = 0.00015;
        detail::Isotopes i;
        i.push_back(i1);
        i.push_back(i2);
        e.isotopes = i;
    }
    e.count = static_cast<Double> (n);
    return e;
}

Double getElectronMass()
{
    return 0.00054857990946;
}

Bool isHydrogen(const detail::Element& e)
{
    return e.isotopes[0].mz == 1.007825 && e.isotopes[1].mz == 2.01410178;
}

} // namespace detail

} // namespace ipaca

} // namespace mstk
