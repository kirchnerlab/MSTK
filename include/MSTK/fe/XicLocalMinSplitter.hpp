/*
 * XicLocalMinSplitter.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_XICLOCALMINSPLITTER_HPP__
#define __MSTK_INCLUDE_MSTK_FE_XICLOCALMINSPLITTER_HPP__

#include <MSTK/config.hpp>
#include <MSTK/common/Types.hpp>
#include <MSTK/fe/CentroidTraits.hpp>

namespace mstk {

namespace fe {

template<typename XicType>
class XicLocalMinSplitter
{
public:
    typedef typename XicType::const_iterator IteratorType;
    typedef std::pair<IteratorType, IteratorType> value_type;
    typedef std::vector<value_type> IteratorContainer;
    typedef typename IteratorContainer::const_iterator const_iterator;
    typedef typename CentroidTraits<typename XicType::value_type>::AbundanceAccessor::value_type ThresholdType;

    template <typename ForwardIterator1, typename ForwardIterator2>
    Size split(ForwardIterator1 rawFirst, ForwardIterator1 rawLast,
        ForwardIterator2 smoothFirst, ForwardIterator2 smoothLast,
        ThresholdType minDepth=0.76);
    const_iterator begin() const;
    const_iterator end() const;
    Size size() const;

private:
    IteratorContainer iterators_;
};

} // namespace fe

} // namespace mstk

//
// template implementation
//

#include <MSTK/fe/utils.hpp>
#include <MSTK/common/Log.hpp>

namespace mstk {

namespace fe {

template<typename XicType>
template <typename ForwardIterator1, typename ForwardIterator2>
Size XicLocalMinSplitter<XicType>::split(ForwardIterator1 rawFirst, ForwardIterator1 rawLast,
    ForwardIterator2 smoothFirst, ForwardIterator2 smoothLast, ThresholdType minDepth)
{
    // get the accessors
    typedef typename std::iterator_traits<ForwardIterator2>::value_type ValueType;
    typedef typename CentroidTraits<ValueType>::AbundanceAccessor AbundanceAccessor;
    typedef typename CentroidTraits<ValueType>::RtAccessor RtAccessor;
    AbundanceAccessor accAb;
    RtAccessor accRt;
    typedef typename AbundanceAccessor::value_type AbType;
    typedef typename RtAccessor::value_type RtType;

    MSTK_LOG(logDEBUG3) << "assigning XIC of length: "
            << std::distance(smoothFirst, smoothLast);
    // clean up
    iterators_.clear();
    // if the assignment is empty, we're done
    if (smoothFirst == smoothLast)
        return 0;
    // with less than 4 measurements, there is no point in splitting because we
    // will always generate a single measurement XIC
    if (std::distance(smoothFirst, smoothLast) < 4) {
        iterators_.push_back(std::make_pair(smoothFirst, smoothLast));
        MSTK_LOG(logDEBUG3) << __FUNCTION__ << ": size too small (" << size()
                << "<4)";
        return 1;
    }

    // iterate through the smoothed XIC, find minima and maxima
    ForwardIterator2 l = smoothFirst;
    ForwardIterator2 i = l;
    std::advance(i, 1);
    ForwardIterator2 r = l;
    std::advance(r, 2);
    ForwardIterator2 currentMin = smoothLast;
    ForwardIterator2 previousMax = smoothLast;
    ForwardIterator2 nextMax = smoothLast;
    // deal with the left boundary: the left boundary is the first
    // lastMin; it can also be the first previousMax, depending if we first hit
    // a local max or a local min later.
    ForwardIterator2 lastMin = smoothFirst;
    if (accAb(*i) <= accAb(*l)) {
        // the first smoothed XIC value is a local max
        previousMax = l;
    }

    ForwardIterator2 smoothLastValid = smoothFirst;
    std::advance(smoothLastValid, std::distance(smoothFirst, smoothLast) - 1);
    while (r != smoothLastValid) {
        double lAb = accAb(*l);
        double iAb = accAb(*i);
        double rAb = accAb(*r);
        // check if we hit a minimum
        if (lAb >= iAb && iAb < rAb) {
            // when we hit a minimum, we just need to record it.
            currentMin = i;
            MSTK_LOG(logDEBUG3) << __FUNCTION__
                    << ": setting currentMin to rt=" << accRt(*i);
        }
        // check if we hit a maximum
        if (lAb < iAb && iAb >= rAb) {
            MSTK_LOG(logDEBUG3) << __FUNCTION__ << ": local max at rt="
                    << accRt(*i);
            // smoothFirst thing we need to check is if this is the smoothFirst time we find
            // a maximum:
            if (previousMax != smoothLast) {
                // this is the second max in the max-min-max combination
                nextMax = i;
                MSTK_LOG(logDEBUG3) << __FUNCTION__
                        << ": setting nextMax to rt=" << accRt(*nextMax);

                // Check if the minimum is deep enough. We use the same criterion as
                // MaxQuant (see Supplementary Info for Cox & Mann, 2008)
                if (accAb(*currentMin)
                        < minDepth
                                * (std::min)(accAb(*previousMax),
                                    accAb(*nextMax))) {
                    // store the range in the original (raw) XIC: all iterators
                    // address the smoothed XIC - hence, we need to find the position
                    // of the minimum in the true XIC. This can be done by simple
                    // distance calculations (the XICs have the same length).
                    ForwardIterator1 f = rawFirst;
                    ForwardIterator1 l = f;
                    std::advance(f, std::distance(smoothFirst, lastMin));
                    // the minimum itself goes to the right
                    std::advance(l, std::distance(smoothFirst, currentMin));
                    // copy, but make sure that we do not copy singles
                    if (std::distance(f, l) > 1) {
                        // push the range
                        iterators_.push_back(std::make_pair(f, l));
                        // advance the iterators
                        lastMin = currentMin;
                    }
                    previousMax = nextMax;
                    MSTK_LOG(logDEBUG3) << __FUNCTION__
                            << ": shifting previousMax to rt="
                            << accRt(*previousMax);

                } else {
                    MSTK_LOG(logDEBUG3) << __FUNCTION__ << "failed criterion: "
                            << accAb(*previousMax) << " > " << iAb << ", < "
                            << accAb(*nextMax) << ", minDepth: " << minDepth
                            << ", critval: "
                            << minDepth
                                    * (std::min)(accAb(*previousMax),
                                        accAb(*nextMax));
                    if (accAb(*nextMax) > accAb(*previousMax)) {
                        previousMax = nextMax;
                        MSTK_LOG(logDEBUG3) << __FUNCTION__
                                << "Not shifting max";
                    }
                }
                // keep the max (no matter what)
                //previousMax = nextMax;
                MSTK_LOG(logDEBUG3) << __FUNCTION__
                        << ": shifting previousMax to rt="
                        << accRt(*previousMax);
            } else {
                // the smoothFirst max we found. Record it and continue.
                previousMax = i;
                MSTK_LOG(logDEBUG3) << __FUNCTION__
                        << ": previousMax set to rt=" << accRt(*i);
            }
        }
        ++l;
        ++i;
        ++r;
    }
    // now deal with what is left
    ForwardIterator1 f = rawFirst;
    std::advance(f, std::distance(smoothFirst, lastMin));
    // the criterion in the loop guarantees that we never copy singles
    iterators_.push_back(std::make_pair(f, rawLast));
    if (std::distance(f, rawLast) <= 1) {
        MSTK_LOG(logWARNING) << "Splitting generated a size "
                << std::distance(f, rawLast) << " XIC.";
    }
    return iterators_.size();
}

template<typename XicType>
typename XicLocalMinSplitter<XicType>::const_iterator XicLocalMinSplitter<XicType>::begin() const
{
    return iterators_.begin();
}

template<typename XicType>
typename XicLocalMinSplitter<XicType>::const_iterator XicLocalMinSplitter<XicType>::end() const
{
    return iterators_.end();
}

template<typename XicType>
Size XicLocalMinSplitter<XicType>::size() const
{
    return iterators_.size();
}

} // namespace fe

} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_FE_XICLOCALMINSPLITTER_HPP__ */
