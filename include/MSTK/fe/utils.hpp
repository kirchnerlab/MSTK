/*
 * utils.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_UTILS_HPP__
#define __MSTK_INCLUDE_MSTK_FE_UTILS_HPP__

#include <functional>
#include <MSTK/fe/SpectrumTraits.hpp>

namespace mstk {

namespace fe {

template<typename T>
class LessByAbundance : public std::binary_function<bool, T, T>
{
public:
    LessByAbundance() :
        abAcc_()
    {
    }

    bool operator()(const T& lhs, const T& rhs) const
    {
        return abAcc_(lhs) < abAcc_(rhs);
    }
private:
    typename SpectrumValueTraits<T>::AbundanceAccessor abAcc_;
};

template<typename T>
struct DiffByMz : public std::binary_function<double, T, T>
{
public:
    DiffByMz() :
        mzAcc_()
    {
    }

    double operator()(const T& lhs, const T& rhs)
    {
        return mzAcc_(lhs) - mzAcc_(rhs);
    }
private:
    typename SpectrumValueTraits<T>::MzAccessor mzAcc_;
};

} // namespace fe

} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_FE_UTILS_HPP__ */
