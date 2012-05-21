/*
 * Splitter.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_SPLITTER_HPP__
#define __MSTK_INCLUDE_MSTK_FE_SPLITTER_HPP__

#include <MSTK/config.hpp>
#include <MSTK/common/Types.hpp>

namespace mstk {

namespace fe {

template<typename T, typename PeakShapeFunction>
class Splitter
{
public:
    typedef typename T::const_iterator IteratorType;
    typedef std::pair<IteratorType, IteratorType> value_type;
    typedef std::vector<value_type> IteratorContainer;
    typedef typename IteratorContainer::const_iterator const_iterator;
    Splitter(const PeakShapeFunction& psf);
    void assign(IteratorType first, IteratorType last);
    const_iterator begin() const;
    const_iterator end() const;
    Size size() const;

private:
    const PeakShapeFunction& psf_;
    IteratorContainer iterators_;
};

}

}

//
// template implementation
//
#include <MSTK/fe/SpectrumTraits.hpp>
#include <MSTK/fe/utils.hpp>
#include <MSTK/common/Log.hpp>

namespace mstk {

namespace fe {

template<typename T, typename PeakShapeFunction>
Splitter<T, PeakShapeFunction>::Splitter(const PeakShapeFunction& psf) :
    psf_(psf)
{
}

template<typename T, typename PeakShapeFunction>
void Splitter<T, PeakShapeFunction>::assign(IteratorType first, IteratorType last)
{
    MSTK_LOG(logDEBUG2)
            << "assigning spectrum of length: " << std::distance(first, last);

    // clean up
    iterators_.clear();
    // if the assignment is empty, we're done
    if (first == last)
        return;
    IteratorType i = first;
    IteratorType j = first;
    ++j;
    typename SpectrumValueTraits<typename SpectrumTraits<T>::Value>::MzAccessor
            mzAcc;
    DiffByMz<typename SpectrumTraits<T>::Value> mzDelta;
    while (j != last) {
        // ask the PSF for the max allowed distance between peaks
        double mzThreshold = psf_.getSupportThreshold(mzAcc(*j));

        if (mzDelta(*j, *i) > mzThreshold) {
            MSTK_LOG(logDEBUG2)
                    << "Splitting at d: " << mzDelta(*j, *i) << " > thresh="
                            << mzThreshold;
            // split if the m/z position are further apart
            iterators_.push_back(std::make_pair(first, j));
            first = i = j;
        } else {
            MSTK_LOG(logDEBUG2)
                    << "Found d: " << mzDelta(*j, *i) << " < thresh="
                            << mzThreshold;
            ++i;
        }
        ++j;
    }
    iterators_.push_back(std::make_pair(first, j));
}

template<typename T, typename PeakShapeFunction>
typename Splitter<T, PeakShapeFunction>::const_iterator Splitter<T, PeakShapeFunction>::begin() const
{
    return iterators_.begin();
}

template<typename T, typename PeakShapeFunction>
typename Splitter<T, PeakShapeFunction>::const_iterator Splitter<T, PeakShapeFunction>::end() const
{
    return iterators_.end();
}

template<typename T, typename PeakShapeFunction>
Size Splitter<T, PeakShapeFunction>::size() const
{
    return iterators_.size();
}

} // namespace fe

} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_FE_SPLITTER_HPP__ */
