/*
 * Centroider.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_CENTROIDER_HPP__
#define __MSTK_INCLUDE_MSTK_FE_CENTROIDER_HPP__

#include <MSTK/common/Types.hpp>
#include <MSTK/common/Log.hpp>
#include <iterator>
#include <utility>

namespace mstk {

namespace fe {

template<class CentroidType, class BumpFinder, class MeanAccumulator,
        class AbundanceAccumulator>
class Centroider : public BumpFinder,
                   public MeanAccumulator,
                   public AbundanceAccumulator
{
public:
    template<typename InputIterator, typename OutputIterator>
    void operator()(InputIterator first, InputIterator last,
        Double retentionTime, UnsignedInt scanNumber, OutputIterator out);
};

//
// template implementation
//

template<class CentroidType, class BumpFinder, class MeanAccumulator,
        class AbundanceAccumulator>
template<class InputIterator, typename OutputIterator>
void Centroider<CentroidType, BumpFinder, MeanAccumulator,
        AbundanceAccumulator>::operator()(InputIterator first,
    InputIterator last, Double retentionTime, UnsignedInt scanNumber,
    OutputIterator out)
{
    // Implement the Centroider workhorse function using the Template Method
    // Design Pattern: the base class simply calls the functions provided
    // by the different algorithm objects. Hence, the following lines implicitly
    // define the interface that must be supported by the different policies.
    MSTK_LOG(logDEBUG2)
            << "Got Spectrum of size: " << std::distance(first, last);
    std::pair<InputIterator, InputIterator> bumpBounds = std::make_pair(first,
        first);
    // find bumps until the right bump bound equals the right range bound
    while (bumpBounds.second != last) {
        // always start the next bump from the right bound of the previous
        bumpBounds = this->findBump(bumpBounds.second, last);
        MSTK_LOG(logDEBUG2)
                << "+- Got bump: [" << std::distance(first, bumpBounds.first)
                        << ", " << std::distance(first, bumpBounds.second)
                        << ").";
        double mz = this->mean(bumpBounds.first, bumpBounds.second);
        double ab = this->abundance(bumpBounds.first, bumpBounds.second);
        // centroids must have a valid abundance measurement
        if (ab > 0.0) {
            CentroidType c(retentionTime, mz, scanNumber, ab,
                bumpBounds.first, bumpBounds.second);
            *out = c;
            ++out;
        }
    }
}

} // namespace fe

} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_FE_CENTROIDER_HPP__ */
