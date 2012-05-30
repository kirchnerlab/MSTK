/*
 * Isotope.cpp
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

#include "MSTK/aas/Isotope.hpp"

namespace mstk {
namespace aas {
namespace elements {

Isotope::Isotope(const mstk::Double& mass, const mstk::Double& frequency) :
        mass_(mass), frequency_(frequency)
{
}

bool Isotope::operator==(const Isotope& i) const
{
    return i.mass_ == mass_ && i.frequency_ == frequency_;
}

bool Isotope::operator!=(const Isotope& i) const
{
    return !(operator ==(i));
}

Isotope& Isotope::operator=(const Isotope& rhs)
{
    if (this != &rhs) {
        mass_ = rhs.mass_;
        frequency_ = rhs.frequency_;
    }
    return *this;
}

std::ostream& operator<<(std::ostream& os, const Isotope& o)
{
    os << o.getMass() << "; " << o.getFrequency();
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<Isotope>& o)
{
    typedef std::vector<Isotope>::const_iterator IT;
    for (IT it = o.begin(); it != o.end(); ++it) {
        os << "(" << *it << ")";
    }
    return os;
}

} // namespace elements
} // namespace aas
} // namespace mstk
