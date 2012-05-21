/*
 * XicFbiTraits.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_TYPES_XICFBITRAITS_HPP__
#define __MSTK_INCLUDE_MSTK_FE_TYPES_XICFBITRAITS_HPP__

#include <MSTK/config.hpp>
#include <MSTK/fe/types/Xic.hpp>
#include <fbi/fbi.h>
#include <utility>

namespace fbi {

template<>
struct Traits<mstk::fe::Xic> : mpl::TraitsGenerator<double, double>
{
};

} // namespace fbi

namespace mstk {

namespace fe {

struct XicBoxGenerator
{
    template<size_t N>
    typename std::tuple_element<N, typename fbi::Traits<Xic>::key_type>::type
    get(const Xic &) const;

    XicBoxGenerator(double mzShift, std::pair<double, double> rtToleranceRange,
        std::pair<double, double> mzToleranceRange, int charge);

    double mzShift_; // in Da
    double minRtTolerance_; // in s
    double maxRtTolerance_; // in s
    double minMzTolerance_; // in ppm
    double maxMzTolerance_; // in ppm
    int charge_;
    static const double deltaS_; // in Da
};

} // namespace fe

} // namespace mstk

//
// template implementation
//
#include <MSTK/common/Log.hpp>
#include <MSTK/common/constrainToRange.hpp>

namespace mstk {

namespace fe {

const double XicBoxGenerator::deltaS_ = 0.0109135;

XicBoxGenerator::XicBoxGenerator(double mzShift,
    std::pair<double, double> rtToleranceRange,
    std::pair<double, double> mzToleranceRange, int charge) :
        mzShift_(mzShift), minRtTolerance_(rtToleranceRange.first), maxRtTolerance_(
            rtToleranceRange.second), minMzTolerance_(mzToleranceRange.first), maxMzTolerance_(
            mzToleranceRange.second), charge_(std::abs(charge))
{
    if (charge_ == 0) {
        mstk_fail("XicBoxGenerator: no zero charges allowed.");
    }
}

template<>
std::pair<double, double> XicBoxGenerator::get<0>(const Xic & xic) const
{
    double rt = xic.getRetentionTime();
    double rtTol = constrainToRange(xic.getRetentionTimeTolerance(),
        minRtTolerance_, maxRtTolerance_);
    return std::make_pair(rt - rtTol, rt + rtTol);
}

template<>
std::pair<double, double> XicBoxGenerator::get<1>(const Xic & xic) const
{
    double mz = xic.getMz() + mzShift_ / charge_; // add shift immediately
    double minMzTol = mz * minMzTolerance_ * 1e-6;
    double maxMzTol = mz * maxMzTolerance_ * 1e-6;
    double mzTolEstimate = xic.getMzTolerance();
    double deltaS = deltaS_ / charge_;
    double mzTol = constrainToRange(
        std::sqrt(mzTolEstimate * mzTolEstimate + deltaS * deltaS), minMzTol,
        maxMzTol);
    MSTK_LOG(logDEBUG3)
            << "m/z tolerance: " << mzTol << " at " << mz << ", " << 1e6
                    * mzTol / mz << "ppm.";
    return std::make_pair(mz - mzTol, mz + mzTol);
}

} // namespace fe

} // namespace mstk

#endif
