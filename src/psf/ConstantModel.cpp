/*
 * ConstantModel.cpp
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
#include <MSTK/common/Error.hpp>
#include <MSTK/psf/PeakParameter.hpp>

namespace mstk {

namespace psf {

unsigned int ConstantModel::numberOfParameters() {
    return numberOfParameters_;
}

void ConstantModel::setParameter(unsigned index, double value) {
    mstk_precondition(index < numberOfParameters(), "ConstantModel::setParameter(): Parameter index out-of-range.");
    a_ = value;
}
double ConstantModel::getParameter(unsigned index) {
    mstk_precondition(index < numberOfParameters(), "ConstantModel::getParameter(): Parameter index out-of-range.");
    return a_;
}

double ConstantModel::at(const double x) const {
    return a_;
}

GeneralizedSlope ConstantModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {1., 0.};
    return GeneralizedSlope(slope, slope + 2);
}

// setter / getter
void ConstantModel::setA(const double a) {
    a_ = a;
}
double ConstantModel::getA() const {
    return a_;
}

} /* namespace psf */

} /* namespace mstk */

