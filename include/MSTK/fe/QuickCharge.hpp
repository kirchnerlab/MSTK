/* 
 * QuickCharge.hpp
 *
 * Copyright (C) 2008 Marc Kirchner
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
#ifndef __MSTK_INCLUDE_MSTK_COMMON_QUICKCHARGE_HPP__
#define __MSTK_INCLUDE_MSTK_COMMON_QUICKCHARGE_HPP__
#include <MSTK/config.hpp>

namespace mstk {

namespace fe {

/**
 * Implementation of the QuickCharge algorithm, as proposed by Michael Hoopman,
 * MacCoss lab, Seattle.
 */
class QuickCharge
{
public:
    /**
     * operator()
     * @param first Iterator pointing at the starting element of a sequence.
     * @param last Iterator pointing at the last element of a sequence.
     * @param out Output iterator for the detected charges.
     */
    template<typename InputIterator, typename OutputIterator>
    void operator()(InputIterator first, InputIterator last,
        OutputIterator out);
};

} // namespace fe

} // namespace mstk

#include <vector>
#include <set>
#include <MSTK/common/Types.hpp>
#include <cmath>
#include <iterator>
#include <MSTK/common/Log.hpp>
#include <MSTK/fe/XicTraits.hpp>

namespace mstk {

namespace fe {

template<typename InputIterator, typename OutputIterator>
void QuickCharge::operator()(InputIterator first, InputIterator last,
    OutputIterator out)
{
    // calculate vector of mass distances
    Size n = std::distance(first, last);
    if (n < 1) {
        return;
    }
    std::vector<Double> diff(n - 1);
    InputIterator l = first;
    InputIterator r = l; std::advance(r, 1);
    Size i = 0;
    typedef typename std::iterator_traits<InputIterator>::value_type ValueType;
    typedef typename XicTraits<ValueType>::MzAccessor MzAccessor;
    MzAccessor accMz;
    while (r != last) {
        diff[i] = accMz(*r) - accMz(*l);
        ++i;
        ++r;
        ++l;
    }
    // walk through all observed masses and calculate
    // potential charges. The rationale is to look at
    // all following peaks in a 1.1 Da window and to determine
    // the inverse distances
    std::set<Int> charges;
    charges.clear();
    for (Size i = 0; i < (n - 1); ++i) {
        Size j = i;
        Double delta(0.0);
        Int oldCharge(0);
        while (j < (n - 1)) {
            delta = delta + diff[j];
            if (delta > 1.1) {
                break;
            }
            Int charge = static_cast<Int>(std::floor(1.0 / delta + 0.5));
            if (charge != oldCharge) {
                oldCharge = charge;
                if (!charges.count(charge)) {
                    *out = charge;
                    ++out;
                    charges.insert(charge);
                }
            }
            ++j;
        }
    }
}

} // namespace fe

} // namespace mstk

#endif
