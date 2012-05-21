/*
 * constrainToRange.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_COMMON_CONSTRAINTORANGE_HPP__
#define __MSTK_INCLUDE_MSTK_COMMON_CONSTRAINTORANGE_HPP__

namespace mstk {

/** Returns a value, constrained to the specified range.
 * Returns the value \c value if value lies inside (\c minVal, \c maxVal ) or
 * the minimum/maximum boundary otherwise. The implementation only requires a
 * less than-operator. The return values of \c constrainToRange is guaranteed to
 * lie in (\c minVal, \c maxVal).
 * @param[in] value The value that is compared against the boundaries.
 * @param[in] minVal Lower boundary of the valid range.
 * @param[in] maxVal Upper boundary of the valid range.
 * @return The value \c value, or one of the boundary values.
 */
template <typename T>
inline T constrainToRange(const T& value, const T& minVal, const T& maxVal)
{
    if (value < minVal) {
        return minVal;
    }
    if (maxVal < value) {
        return maxVal;
    }
    return value;
}

/** Returns a value, constrained to the specified range.
 * Returns the value \c value if \c comp(\c minVal, \c value) \c && \c comp(\c
 * value, \c maxValue). If the first test fails, the function returns the
 * minimum, if the second test fails, it returns the maximum.
 * @param[in] value The value that is compared against the boundaries.
 * @param[in] minVal Lower boundary of the valid range.
 * @param[in] maxVal Upper boundary of the valid range.
 * @return The value \c value if it lies inside the boundaries, the lower of
 *         upper boundary value otherwise. 
 */
template <typename T, typename Comp>
inline T constrainToRange(const T& value, const T& minVal, const T& maxVal, Comp comp)
{
    if (comp(value, minVal)) {
        return minVal;
    }
    if (comp(maxVal, value)) {
        return maxVal;
    }
    return value;
}

} // namespace mstk

#endif


