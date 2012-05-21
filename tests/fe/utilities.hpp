/*
 * utilities.hpp
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
#ifndef __MSTK_TESTS_FP_UTILITIES_HPP__
#define __MSTK_TESTS_FP_UTILITIES_HPP__

#include <MSTK/config.hpp>
#include <MSTK/fe/types/IsotopePattern.hpp>
#include <MSTK/fe/types/Centroid.hpp>
#include <MSTK/fe/types/Xic.hpp>
#include <vector>

mstk::fe::IsotopePattern makeIsotopePattern(const std::vector<double>& mz,
    double rt, int charge);

std::vector<mstk::fe::Centroid> makeCentroids(size_t n, const double* mz,
    const double* rt, const unsigned int* sn, const double *ab);

mstk::fe::Xic makeXic(size_t n, const double* mz,
    const double* rt, const unsigned int* sn, const double *ab);

#endif /* __MSTK_TESTS_FP_UTILITIES_HPP__ */
