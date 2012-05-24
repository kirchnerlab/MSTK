/*
 * SqrtModel.cpp
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
#include <MSTK/psf/PeakParameter.hpp>

namespace mstk {

namespace psf {

unsigned int SqrtModel::numberOfParameters() {
    return numberOfParameters_;
}

void SqrtModel::setParameter(unsigned index, double value) {
    mstk_precondition(index < numberOfParameters(), "SqrtModel::setParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        a_ = value;
    }
    else {
        b_ = value;    
    }
}
double SqrtModel::getParameter(unsigned index) {
    mstk_precondition(index < numberOfParameters(), "SqrtModel::getParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        return a_;
    }
    else {
        return b_;    
    }
}

double SqrtModel::at(const double x) const {
    mstk_precondition(x >= 0, "SqrtModel::at(): Parameter x hast to be >= 0.");
    return a_ * std::sqrt(x) + b_;
}

GeneralizedSlope SqrtModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {std::sqrt(x), 1., 0.};
    return GeneralizedSlope(slope, slope + 3);
}

// setter / getter
void SqrtModel::setA(const double a) {
    a_ = a;
}
double SqrtModel::getA() const {
    return a_;
}

void SqrtModel::setB(const double b) {
    b_ = b;
}
double SqrtModel::getB() const {
    return b_;
}

} // namespace psf

} // namespace mstk

