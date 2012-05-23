/*
 * Stoichiometry.cpp
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
#include <MSTK/ipaca/Stoichiometry.hpp>
#include <cmath>
#include <iostream>

namespace mstk {

namespace ipaca {

Bool detail::isPlausibleStoichiometry(const detail::Stoichiometry& s)
{
    typedef detail::Stoichiometry::const_iterator CI;
    bool allZero = true;
    for (CI i = s.begin(); i < s.end(); ++i) {
        Double c = i->count;
        if (c < 0) {
            return false;
        }
        if (allZero && c > 0) {
            allZero = false;
        }
    }
    if (!allZero) {
        return true;
    } else {
        return false;
    }
}

void detail::splitStoichiometry(const detail::Stoichiometry& s,
    detail::Stoichiometry& intStoi, detail::Stoichiometry& fracStoi)
{
    intStoi.clear();
    fracStoi.clear();
    typedef detail::Stoichiometry::const_iterator CI;
    for (CI i = s.begin(); i != s.end(); ++i) {
        Double integer = trunc(i->count);
        Double fractional = i->count - integer;
        if (integer > 0.0) {
            intStoi.push_back(*i);
            intStoi.back().count = integer;
        }
        if (fractional > 0.0) {
            fracStoi.push_back(*i);
            fracStoi.back().count = fractional;
        }
    }
}

std::ostream& detail::operator<<(std::ostream& os, const detail::Stoichiometry& s)
{
    typedef detail::Stoichiometry::const_iterator CI;
    os << "(";
    for (CI i = s.begin(); i != s.end(); ++i) {
        os << "(" << i->isotopes[0].mz << ", " << i->count << ")";
    }
    os << ")";
    return os;
}

} // namespace ipaca

} // namespace mstk

