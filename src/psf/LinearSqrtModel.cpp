/*
 * LinearSqrtModel.hpp
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

unsigned int LinearSqrtModel::numberOfParameters() {
    return numberOfParameters_;
}

void LinearSqrtModel::setParameter(unsigned index, double value) {
    mstk_precondition(index < numberOfParameters(), "LinearSqrtModel::setParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        a_ = value;
    }
    else {
        b_ = value;    
    }
}
double LinearSqrtModel::getParameter(unsigned index) {
    mstk_precondition(index < numberOfParameters(), "LinearSqrtModel::getParameter(): Parameter index out-of-range.");
    if(index == 0) {    
        return a_;
    }
    else {
        return b_;    
    }
}

double LinearSqrtModel::at(const double x) const {
    mstk_precondition(x >= 0, "LinearSqrtModel::at(): Parameter x has to be >= 0.");
    return a_ * x * std::sqrt(x) + b_;
}

GeneralizedSlope LinearSqrtModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {x * std::sqrt(x), 1., 0.};
    return GeneralizedSlope(slope, slope + 3);
}

// setter / getter
void LinearSqrtModel::setA(const double a) {
    a_ = a;
}
double LinearSqrtModel::getA() const {
    return a_;
}

void LinearSqrtModel::setB(const double b) {
    b_ = b;
}
double LinearSqrtModel::getB() const {
    return b_;
}




unsigned int LinearSqrtOriginModel::numberOfParameters() {
    return numberOfParameters_;
}

void LinearSqrtOriginModel::setParameter(unsigned index, double value) {
    mstk_precondition(index < numberOfParameters(), "LinearSqrtModel::setParameter(): Parameter index out-of-range.");    
    a_ = value;
}
double LinearSqrtOriginModel::getParameter(unsigned index) {
    mstk_precondition(index < numberOfParameters(), "LinearSqrtModel::getParameter(): Parameter index out-of-range.");
    return a_;
}

double LinearSqrtOriginModel::at(const double x) const {
    mstk_precondition(x >= 0, "LinearSqrtOriginModel::at(): Parameter x has to be >= 0.");
    return a_ * x * std::sqrt(x);
}

GeneralizedSlope LinearSqrtOriginModel::slopeInParameterSpaceFor(double x) const {
    double slope[] = {x * std::sqrt(x), 0.};
    return GeneralizedSlope(slope, slope + 2);
}

// setter / getter
void LinearSqrtOriginModel::setA(const double a) {
    a_ = a;
}
double LinearSqrtOriginModel::getA() const {
    return a_;
}

} /* namespace psf */

} /* namespace mstk */

