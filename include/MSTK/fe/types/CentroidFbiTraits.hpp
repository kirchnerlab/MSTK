/*
 * CentroidFbiTraits.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_TYPES_CENTROIDFBITRAITS_HPP__
#define __MSTK_INCLUDE_MSTK_FE_TYPES_CENTROIDFBITRAITS_HPP__

#include <MSTK/config.hpp>
#include <MSTK/fe/types/Centroid.hpp>
#include <fbi/fbi.h>
#include <utility>

namespace fbi {

template<>
struct Traits<mstk::fe::Centroid> : mpl::TraitsGenerator<double, double>
{
};

}

namespace mstk {

namespace fe {

struct CentroidBoxGenerator
{
    template<size_t N>
    typename std::tuple_element<N,
            typename fbi::Traits<Centroid>::key_type>::type
    get(const Centroid &) const;

    int snTolerance_;
    double mzTolerance_; // in ppm

    CentroidBoxGenerator(int snTolerance, double mzTolerance) :
        snTolerance_(snTolerance), mzTolerance_(mzTolerance)
    {
    }
};

}

}

//
// template implementation
//

namespace mstk {

namespace fe {

template<>
std::pair<double, double> CentroidBoxGenerator::get<0>(const Centroid & centroid) const
{
    double sn = centroid.getScanNumber();
    return std::make_pair(sn - snTolerance_ - 0.1, sn + snTolerance_ + 0.1);
}

template<>
std::pair<double, double> CentroidBoxGenerator::get<1>(const Centroid & centroid) const
{
    double shift = mzTolerance_ * 1e-6;
    double mz = centroid.getMz();
    return std::make_pair(mz * (1 - shift), mz * (1 + shift));
}

}

}


#endif
