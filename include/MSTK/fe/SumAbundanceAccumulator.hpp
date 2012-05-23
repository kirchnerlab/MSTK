/*
 * SumAbundanceAccumulator.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_SUMABUNDANCEACCUMULATOR_HPP__
#define __MSTK_INCLUDE_MSTK_FE_SUMABUNDANCEACCUMULATOR_HPP__

#include <MSTK/config.hpp>
#include <MSTK/fe/SpectrumTraits.hpp>
#include <MSTK/fe/utils.hpp>
#include <MSTK/common/Error.hpp>
#include <algorithm>
#include <iterator>

namespace mstk {

namespace fe {

class SumAbundanceAccumulator
{
protected:
    virtual ~SumAbundanceAccumulator() = 0;
    template<typename InputIterator>
    double abundance(InputIterator first, InputIterator last);
};

}

}

//
// template implementation
//
namespace mstk {

namespace fe {

template<typename InputIterator>
double SumAbundanceAccumulator::abundance(InputIterator first, InputIterator last)
{
    // summation is invariant towards the order of the elements
    typedef typename std::iterator_traits<InputIterator>::value_type ValueType;
    typename SpectrumValueTraits<ValueType>::AbundanceAccessor abAcc;
    double abundance = 0.0;
    for (InputIterator i = first; i != last; ++i) {
        abundance += abAcc(*i);
    }
    mstk_assert(!(abundance < 0.0), "Bogus abundance value: must be non-negative.");
    return abundance;
}

}

}

#endif /* __MSTK_INCLUDE_MSTK_FE_SUMABUNDANCEACCUMULATOR_HPP__ */
