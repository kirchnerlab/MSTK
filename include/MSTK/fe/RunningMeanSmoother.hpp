/*
 * RunningMeanSmoother.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_RUNNINGMEANSMOOTHER_HPP__
#define __MSTK_INCLUDE_MSTK_FE_RUNNINGMEANSMOOTHER_HPP__

#include <MSTK/config.hpp>

namespace mstk {

namespace fe {

class RunningMeanSmoother
{
protected:
    virtual ~RunningMeanSmoother() = 0;

    template<typename InputIterator>
    void smooth(InputIterator first, InputIterator last);
};

}

}

//
// template implementation
//
#include <iterator>
#include <MSTK/fe/CentroidTraits.hpp>
#include <MSTK/common/Log.hpp>

namespace mstk {

namespace fe {

template<typename InputIterator>
void RunningMeanSmoother::smooth(InputIterator first,
    InputIterator last)
{
    // check that the XIC has at least the size of the structuring element
    if (std::distance(first, last) < 3) {
        MSTK_LOG(logWARNING)
                << "RunningMeanSmoother: input shorter than structuring element.";
        return;
    }
    // get the accessors
    typedef typename std::iterator_traits<InputIterator>::value_type ValueType;
    typedef typename CentroidTraits<ValueType>::AbundanceAccessor AbundanceAccessor;
    typedef typename CentroidTraits<ValueType>::RtAccessor RtAccessor;
    AbundanceAccessor accAb;
    RtAccessor accRt;
    typedef typename AbundanceAccessor::value_type AbType;
    typedef typename RtAccessor::value_type RtType;

    // position the structuring element
    InputIterator l = first;
    InputIterator i = first;
    std::advance(i, 1);
    InputIterator r = first;
    std::advance(r, 2);

    // Set up the variables; this must take into account the assignments at
    // the top of the loop below, hence it looks strange.
    // Setting variables up that way enables us to do in-place smoothing.
    RtType lRt;
    RtType iRt = accRt(*l);
    RtType rRt = accRt(*i);
    AbType lAb;
    AbType iAb = accAb(*l);
    AbType rAb = accAb(*i);

    // Apply a running mean on the observed vector: we use a sparse running
    // mean in which the contributions are weighted be the relative distances
    // of the left and right neighbor. The first and last entry in the vector
    // are always left untouched.
    while (r != last) {
        lRt = iRt;
        iRt = rRt;
        rRt = accRt(*r);
        // get the span between l and r and left/right distances from i
        RtType span = rRt - lRt;
        RtType dl = iRt - lRt;
        RtType dr = rRt - iRt;

        lAb = iAb;
        iAb = rAb;
        rAb = accAb(*r);
        AbType abundance = (2
                * ((1.0 - dl / span) * lAb + (1.0 - dr / span) * rAb) + iAb)
                / 3.0;
        // in-place smoothing
        accAb(*i) = abundance;
        // advance the iterators
        ++i;
        ++r;
    }
}

} // namespace fe

} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_FE_RUNNINGMEANSMOOTHER_HPP__ */
