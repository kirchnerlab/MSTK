/*
 * GaussianMeanAccumulator.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_GAUSSIANMEANACCUMULATOR_HPP__
#define __MSTK_INCLUDE_MSTK_FE_GAUSSIANMEANACCUMULATOR_HPP__

#include <MSTK/config.hpp>
#include <MSTK/common/Error.hpp>
#include <MSTK/common/Log.hpp>
#include <MSTK/fe/utils.hpp>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <utility>

namespace mstk {

namespace fe {

class GaussianMeanAccumulator
{
protected:
    virtual ~GaussianMeanAccumulator() = 0;

    template<typename InputIterator>
    double mean(InputIterator first, InputIterator last);

    template<typename InputIterator>
    void trimAndMax(InputIterator first, InputIterator last,
        InputIterator& start, InputIterator& stop, InputIterator& maxElement);
};

}

}

//
// template implementation
//
namespace mstk {

namespace fe {

template<typename InputIterator>
void GaussianMeanAccumulator::trimAndMax(InputIterator first,
    InputIterator last, InputIterator& start, InputIterator& stop,
    InputIterator& maxElement)
{
    // the position of the maximum abundance element
    maxElement = last;
    // the beginning of the non-zero data
    start = last;
    // the end of the non-zero data
    stop = last;

    // empty input?
    if (first == last) {
        return;
    }

    // and a simple running iterator
    InputIterator i = first;

    // get an accessor to access the abundances
    typedef typename std::iterator_traits<InputIterator>::value_type ValueType;
    typename SpectrumValueTraits<ValueType>::AbundanceAccessor accAb;

    // skip leading zeros
    while (accAb(*i) == 0.0) {
        ++i;
    }
    if (i != last) {
        // find the maximum (until we hit a zero)
        start = i;
        maxElement = start;
        double maxValue = accAb(*maxElement);
        while (i != last) {
            // update the max abundance position and value, if necessary
            if (maxValue < accAb(*i)) {
                maxElement = i;
                maxValue = accAb(*maxElement);
            }
            // check if we his a zero; once this happens, there will only
            // be trailing zeros (because our input is a bump).
            if (accAb(*i) == 0.0) {
                stop = i;
                break;
            }
            ++i;
        }
    }
    // NOTE: if there are only zeros, then first, last and maxElement will
    //       all point to last (which is why start and stop are initialized
    //       to last).
}

template<typename InputIterator>
double GaussianMeanAccumulator::mean(InputIterator first, InputIterator last)
{
    // There is no way of returning a reasonable value if the input is
    // empty. Hence, we make non-empty input a precondition.
    mstk_precondition(first != last,
            "GaussianMeanAccumulator::mean: cannot calculate mean of empty input.");

    typedef typename std::iterator_traits<InputIterator>::value_type ValueType;
    typename SpectrumValueTraits<ValueType>::MzAccessor accMz;
    typename SpectrumValueTraits<ValueType>::AbundanceAccessor accAb;

    InputIterator m;
    // Get rid of leading and trailing zeros, and get the position of the
    // maximum abundance in a single pass.
    trimAndMax(first, last, first, last, m);
    if (first == last) {
        // There were only zero abundances in [first, last).
        // Return the mean of the m/z values as a best guess.
        double mass = 0.0;
        for (InputIterator i = first; i != last; ++i) {
            mass += accAb(*i);
        }
        return mass / static_cast<double> (std::distance(first, last));
    }
    // calculate the accurate mass
    double mass = 0.0;
    // depending on the location of the maximum relative to the beginning
    // and end of the bump, we need to choose between simple averaging
    // and a three-point Gaussian estimate.
    typedef typename std::iterator_traits<InputIterator>::difference_type
            DiffType;
    DiffType nl = std::distance(first, m);
    DiffType nr = std::distance(m, last);
    MSTK_LOG(logDEBUG2)
            << "Got nl=" << nl << ", nr=" << nr;
    if ((nl > 0) && (nr > 1)) {
        // We have (at least) one non-zero measurement to the left and one
        // measurement to the right of the maximum; hence, we can do a
        // three-point Gaussian fit.
        InputIterator l = first;
        std::advance(l, nl - 1);
        InputIterator r = m;
        std::advance(r, 1);
        double lAb = accAb(*l);
        double lMz = accMz(*l);
        double mAb = accAb(*m);
        double mMz = accMz(*m);
        double rAb = accAb(*r);
        double rMz = accMz(*r);
        double numerator = (std::log(mAb) - std::log(rAb)) * (lMz * lMz)
                + (std::log(rAb) - std::log(lAb)) * (mMz * mMz) + (std::log(
            lAb) - std::log(mAb)) * (rMz * rMz);
        double denominator = 2 * ((std::log(mAb) - std::log(rAb)) * lMz
                + (std::log(rAb) - std::log(lAb)) * mMz + (std::log(lAb)
                - std::log(mAb)) * rMz);
        mass = numerator / denominator;
    } else {
        if ((nl > 0) && (nr <= 1)) {
            // average to the left
            InputIterator l = first;
            std::advance(l, nl - 1);
            double lAb = accAb(*l);
            double lMz = accMz(*l);
            double mAb = accAb(*m);
            double mMz = accMz(*m);
            double denominator = (lAb + mAb);
            if (denominator > 0.0) {
                // return the weighted mean position
                mass = (lMz * lAb + mMz * mAb) / denominator;
            } else {
                // return the unweighted mean position
                mass = (lMz + mMz) / 2.0;
            }
        } else {
            if ((nl == 0) && (nr > 1)) {
                // average to the right
                InputIterator r = m;
                std::advance(r, 1);
                double mAb = accAb(*m);
                double mMz = accMz(*m);
                double rAb = accAb(*r);
                double rMz = accMz(*r);
                double denominator = (mAb + rAb);
                if (denominator > 0.0) {
                    // return the weighted mean position
                    mass = (mMz * mAb + rMz * rAb) / denominator;
                } else {
                    // return the unweighted mean position
                    mass = (mMz + rMz) / 2.0;
                }
            } else {
                if ((nl == 0) && (nr == 1)) {
                    // just take it
                    MSTK_LOG(mstk::logDEBUG)
                            << "Single-value centroid.";
                    mass = accMz(*m);
                } else {
                    // Something is very wrong. Die.
                    assert(false && "Spectrum length is zero.");
                    return 0;
                }
            }
        }
    }
    assert(!(mass < 0.0) && "Bogus mass value.");
    MSTK_LOG(logDEBUG2)
            << "mass: " << mass;
    return mass;
}

}

}

#endif /* __MSTK_INCLUDE_MSTK_FE_GAUSSIANMEANACCUMULATOR_HPP__ */
