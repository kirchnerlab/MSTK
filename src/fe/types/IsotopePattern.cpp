/*
 * IsotopePattern.cpp
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
#include <MSTK/fe/types/IsotopePattern.hpp>
#include <MSTK/common/Error.hpp>

#include <functional>

// anonymous namespace for local stuff
namespace {

struct Xic2SpectrumElement : public std::unary_function<
        mstk::fe::Spectrum::Element, mstk::fe::Xic>
{
    mstk::fe::Spectrum::Element operator()(const mstk::fe::Xic& xic)
    {
        return mstk::fe::Spectrum::Element(xic.getMz(), xic.getAbundance());
    }
};

} // anonymous namespace

namespace mstk {

namespace fe {

IsotopePattern::IsotopePattern() :
    Collection<Xic> ()
{
}

IsotopePattern::IsotopePattern(const_iterator first, const_iterator last) :
    Collection<Xic>(first, last)
{
}

void IsotopePattern::setCharges(const std::set<int>& z)
{
    charges_ = z;
}

const std::set<int>& IsotopePattern::getCharges() const
{
    return charges_;
}

void IsotopePattern::asSpectrum(Spectrum& ss) const
{
    // clear the spectrum
    ss.clear();
    // make sure the XIC vector is sorted
    std::vector<Xic>::const_iterator first = c_.begin();
    std::vector<Xic>::const_iterator last = c_.end();
    // and transform
    std::transform(c_.begin(), c_.end(), std::back_inserter(ss),
        Xic2SpectrumElement());
    std::sort(
        ss.begin(),
        ss.end(),
        Spectrum::LessThanMz<Spectrum::Element,
                Spectrum::Element>());
    // Spectrum only supports a single precursor charge; hence, only
    // specify a charge if there is no uncertainty.
    if (charges_.size() == 1) {
        ss.setPrecursorCharge(*(charges_.begin()));
    } else {
        // TODO: extend Spectrum to support multiple charges
        ss.setPrecursorCharge(0);
    }
}

void IsotopePattern::split(std::vector<IsotopePattern>& pips)
{
    // TODO: implement!
    mstk_fail("not implemented");
    return;
}

double IsotopePattern::getAbundance() const
{
    double ab = 0.0;
    for (const_iterator i = this->begin(); i != this->end(); ++i) {
        ab += i->getAbundance();
    }
    return ab;
}

std::ostream& operator<<(std::ostream& os, const IsotopePattern& p)
{
    typedef std::set<int>::const_iterator SIT;
    os << p.size() << '\t';
    for (SIT c = p.getCharges().begin(); c != p.getCharges().end(); ++c) {
        os << *c << ',';
    }
    os << '\t' << p.getAbundance();
    /*
    typedef Collection<Xic>::const_iterator IT;
    for (IT i = p.begin(); i != p.end(); ++i) {
        os << *i << '\n';
    }
    */
    return os;
}

} // namespace fe

} // namespace mstk

