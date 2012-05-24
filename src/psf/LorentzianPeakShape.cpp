/*
 * LorentzianPeakShape.hpp
 *
 * Copyright (C) 2011 Bernhard Kausler
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
#include <cmath>

#include <MSTK/common/Error.hpp>
#include <MSTK/psf/PeakShape.hpp>

namespace mstk {

namespace psf {

double LorentzianPeakShape::at(const double xCoordinate) const {
    return fwhm_ / ((xCoordinate * xCoordinate) + (fwhm_*fwhm_));
}

double LorentzianPeakShape::getSupportThreshold() const {
    return this->getFwhm() * this->getFwhmFactorForSupportThreshold();
}

// construction
LorentzianPeakShape::LorentzianPeakShape(const double fwhm, const double fwhmFactorForSupportThreshold) 
    : fwhm_(fwhm), fwhmFactorForSupportThreshold_(fwhmFactorForSupportThreshold) {
    mstk_precondition(fwhm > 0, "LorentzianPeakShape::LorentzianPeakShape(): Parameter fwhm has to be positive.");
    mstk_precondition(fwhmFactorForSupportThreshold > 0, "LorentzianPeakShape::LorentzianPeakShape(): fwhmFactorForSupportThreshold has to be positive.");
}


// setter/getter
void LorentzianPeakShape::setFwhm(const double fwhm) {
    mstk_precondition(fwhm > 0, "LorentzianPeakShape::LorentzianPeakShape(): Parameter fwhm has to be positive.");
    fwhm_ = fwhm;
}
double LorentzianPeakShape::getFwhm() const {
    return fwhm_;
}

void LorentzianPeakShape::setFwhmFactorForSupportThreshold(const double factor) {
    mstk_precondition(factor > 0, "LorentzianPeakShape::LorentzianPeakShape(): Parameter fwhmFactorForSupportThreshold has to be positive.");
    fwhmFactorForSupportThreshold_ = factor;
}
double LorentzianPeakShape::getFwhmFactorForSupportThreshold() const {
    return fwhmFactorForSupportThreshold_;
}

} // namespace psf

} // namespace mstk

