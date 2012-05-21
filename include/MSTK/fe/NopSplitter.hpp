/*
 * NopSplitter.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_NOPSPLITTER_HPP__
#define __MSTK_INCLUDE_MSTK_FE_NOPSPLITTER_HPP__

namespace mstk {

namespace fe {

/** No-operation (NOP) isotope pattern splitter implementation.
 * Use this as a policy class for mstk::fe::IsotopePatternExtractor if you do 
 * not want to carry out isotope pattern splitting.
 */
class NopSplitter
{
public:
    template <typename XicSets>
    void split(XicSets& xicsets);
};

} // namespace fe

} // namespace mstk

//
// template implementation
//

namespace mstk {

namespace fe {

template<typename XicSets>
void NopSplitter::split(XicSets& xicSets)
{
    return;
}

} // namespace fe

} // namespace mstk

#endif

