/*
 * OnePassCentroidMerger.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_CENTROIDWEIGHTEDMEANDISAMBIGUATOR_HPP__
#define __MSTK_INCLUDE_MSTK_FE_CENTROIDWEIGHTEDMEANDISAMBIGUATOR_HPP__

#include <MSTK/config.hpp>
#include <MSTK/fe/CentroidTraits.hpp>
#include <MSTK/common/Error.hpp>
#include <MSTK/common/Log.hpp>
#include <iterator>

namespace mstk {

namespace fe {

class CentroidWeightedMeanDisambiguator
{
protected:
    virtual ~CentroidWeightedMeanDisambiguator() = 0;

    template<typename InputIterator>
    InputIterator disambiguate(InputIterator first, InputIterator last);
};

} // namespace fe

} // namespace mstk

//
// template implementation
//
namespace mstk {

namespace fe {

template<typename InputIterator>
InputIterator CentroidWeightedMeanDisambiguator::disambiguate(
    InputIterator first, InputIterator last)
{
    // Single-pass implementation that merges (consecutive) duplicate scan
    // number entries (ie. requires an scanNumber-sorted XIC, which should be
    // the same as an rt-sorting).
    // The code uses the weighted mean of the accurate masses of the
    // merged centroids for the new centroid. It does not drop zero abundance
    // centroids unless they are merged with non-zero abundance centroid. In
    // the latter case, the weighted mean will automatically ignore the m/z
    // contribution of the zero abundance centroid.
    // In cases where there are zero abundance centroids to be merged, the
    // will determine their arithmetic mean for the merged centroid.

    MSTK_LOG(logDEBUG) << "disambiguate: starting with "
            << std::distance(first, last) << " centroids.";

    // get the accessors
    typedef typename std::iterator_traits<InputIterator>::value_type ValueType;
    typedef typename CentroidTraits<ValueType>::MzAccessor MzAccessor;
    MzAccessor accMz;
    typedef typename CentroidTraits<ValueType>::AbundanceAccessor AbundanceAccessor;
    AbundanceAccessor accAb;
    typedef typename CentroidTraits<ValueType>::RtAccessor RtAccessor;
    RtAccessor accRt;

    InputIterator cur = first;
    InputIterator p = cur;
    while (p != last) {
        if (cur != p) {
            *cur = *p;
            ++p;
        }
        // note: this assumes that MzAccessor::value_type and
        // AbundanceAccessor::value_type can multiplied and converted back
        // to MzAccessor::value_type.
        typename MzAccessor::value_type wMz = accMz(*cur) * accAb(*cur);
        typename MzAccessor::value_type sMz = accMz(*cur);
        typename AbundanceAccessor::value_type sAb = accAb(*cur);
        // while we have not reached the end and while the rt measurements match
        while (p != last && accRt(*cur) == accRt(*p)) {
            if (cur != p) {
                sAb += accAb(*p);
                wMz += accMz(*p) * accAb(*p);
                sMz += accMz(*p);
            }
            ++p;
        }
        mstk_assert(sAb >= 0.0, "Negative total abundance not allowed.");
        if (sAb > 0.0) {
            // Use weighted mean if there was at least one abundance
            // measurement above zero.
            accMz(*cur) = wMz / sAb;
        } else {
            // Use normal mean for all-zero measurements.
            // We must not delete zeros because they may influence splitting.
            accMz(*cur) = sMz / std::distance(cur, p);
        }
        accAb(*cur) = sAb;
        ++cur;
    }
    MSTK_LOG(logDEBUG) << "disambiguate: finishing with "
            << std::distance(first, cur) << " centroids.";
    return cur;
}

} // namespace fe

} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_FE_CENTROIDWEIGHTEDMEANDISAMBIGUATOR_HPP__ */
