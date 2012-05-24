/*
 * GaussianPeakShape.cpp
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

using namespace mstk::psf;
using namespace mstk;

double GaussianPeakShape::at(const double xCoordinate) const {
    return std::exp(-(xCoordinate * xCoordinate) / (2 * sigma_ * sigma_));
}

double GaussianPeakShape::getSupportThreshold() const {
    return this->getSigma() * this->getSigmaFactorForSupportThreshold();
}

// construction
GaussianPeakShape::GaussianPeakShape(const double sigma, const double sigmaFactorForSupportThreshold) 
    : sigma_(sigma), sigmaFactorForSupportThreshold_(sigmaFactorForSupportThreshold) {
    mstk_precondition(sigma > 0, "GaussianPeakShape::GaussianPeakShape(): sigma has to be positive.");
    mstk_precondition(sigmaFactorForSupportThreshold > 0, "GaussianPeakShape::GaussianPeakShape(): sigmaFactorForSupportThreshold has to be positive.");
}


// setter/getter
void GaussianPeakShape::setSigma(const double sigma) {
    mstk_precondition(sigma > 0, "GaussianPeakShape::GaussianPeakShape(): Parameter sigma has to be positive."); 
    sigma_ = sigma; 
}


void GaussianPeakShape::setFwhm(const double fwhm) {
    mstk_precondition(fwhm > 0, "GaussianPeakShape::GaussianPeakShape(): Parameter fwhm has to be positive.");
    sigma_ = fwhm / sigmaToFwhmConversionFactor();
}
double GaussianPeakShape::getFwhm() const {
    return sigma_ * sigmaToFwhmConversionFactor();
}

void GaussianPeakShape::setSigmaFactorForSupportThreshold(const double factor) {
    mstk_precondition(factor > 0, "GaussianPeakShape::GaussianPeakShape(): sigmaFactorForSupportThreshold has to be positive.");
    sigmaFactorForSupportThreshold_ = factor;
}
double GaussianPeakShape::getSigmaFactorForSupportThreshold() const {
    return sigmaFactorForSupportThreshold_;
}

double GaussianPeakShape::sigmaToFwhmConversionFactor() const {
    return 2 * sqrt(2 * std::log(2.));
}
