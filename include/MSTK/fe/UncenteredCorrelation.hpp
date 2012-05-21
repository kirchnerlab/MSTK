/*
 * UncenteredCorrelation.hpp
 *
 * Copyright (C) 2012 Marc Kirchner
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
#ifndef __MSTK_INCLUDE_MSTK_FE_UNCENTEREDCORRELATION_HPP__
#define __MSTK_INCLUDE_MSTK_FE_UNCENTEREDCORRELATION_HPP__

#include <MSTK/config.hpp>

namespace mstk {

namespace fe {

class UncenteredCorrelation
{
public:
    typedef double ThresholdType;

protected:
    virtual ~UncenteredCorrelation() = 0;

    template<typename InputIterator>
    ThresholdType correlate(InputIterator lhsFirst, InputIterator lhsLast,
        InputIterator rhsFirst, InputIterator rhsLast);
};

}

}

//
// template implementation
//
#include <iterator>
#include <limits>
#include <cmath>
#include <MSTK/fe/CentroidTraits.hpp>
#include <MSTK/common/Log.hpp>

namespace mstk {

namespace fe {

template<typename InputIterator>
UncenteredCorrelation::ThresholdType UncenteredCorrelation::correlate(
    InputIterator lhsFirst, InputIterator lhsLast, InputIterator rhsFirst,
    InputIterator rhsLast)
{
    // get the accessors
    typedef typename std::iterator_traits<InputIterator>::value_type ValueType;
    typedef typename CentroidTraits<ValueType>::AbundanceAccessor AbundanceAccessor;
    typedef typename CentroidTraits<ValueType>::RtAccessor RtAccessor;
    AbundanceAccessor accAb;
    RtAccessor accRt;
    typedef typename AbundanceAccessor::value_type AbType;
    typedef typename RtAccessor::value_type RtType;

    InputIterator lhsIt = lhsFirst;
    InputIterator rhsIt = rhsFirst;
    AbType l, r;
    AbType lr(0.0), lsq(0.0), rsq(0.0);
    RtType lhsRt, rhsRt;
    while ((lhsIt != lhsLast) || (rhsIt != rhsLast)) {
        if (lhsIt != lhsLast) {
            lhsRt = accRt(*lhsIt);
            l = accAb(*lhsIt);
        } else {
            lhsRt = std::numeric_limits<RtType>::max();
            l = 0.0;
        }
        if (rhsIt != rhsLast) {
            rhsRt = accRt(*rhsIt);
            r = accAb(*rhsIt);
        } else {
            rhsRt = std::numeric_limits<RtType>::max();
            r = 0.0;
        }
        if (lhsRt == rhsRt) {
            // same scan
            lr += l * r;
            lsq += l * l;
            rsq += r * r;
            ++lhsIt;
            ++rhsIt;
        } else
            if (lhsRt < rhsRt) {
                // no rhs measurement at lhs rt pos, hence add lhs only
                lsq += l * l;
                ++lhsIt;
            } else {
                // no lhs measurement at rhs rt pos, hence add rhs only
                rsq += r * r;
                ++rhsIt;
            }
    }
    return static_cast<ThresholdType>((lr != 0) ? (lr / std::sqrt(lsq * rsq)) : 0.0);
}

} // namespace fe

} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_FE_UNCENTEREDCORRELATION_HPP__ */
