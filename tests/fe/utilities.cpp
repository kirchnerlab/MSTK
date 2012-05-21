/*
 * utilities.cpp
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
#include "utilities.hpp"
#include <set>
#include <MSTK/fe/types/IsotopePattern.hpp>

using namespace mstk::fe;

IsotopePattern makeIsotopePattern(const std::vector<double>& mz,
    double rt, int charge)
{
    // generates three centroids for every requested XIC position
    // add a bit of variance to the m/z measurements
    double deltaMz[] = { -0.002, 0.0, 0.002 };
    double deltaRt[] = { -2.0, 0.0, 2.0 };
    double abundances[] = { 500.0, 1000.0, 500.0 };
    IsotopePattern ip;
    typedef std::vector<double>::const_iterator IT;
    for (IT m = mz.begin(); m != mz.end(); ++m) {
        Xic x;
        for (size_t i = 0; i < 3; ++i) {
            Centroid c;
            c.setMz(*m + deltaMz[i]);
            c.setRetentionTime(rt + deltaRt[i]);
            c.setAbundance(abundances[i]);
            x.push_back(c);
        }
        x.recalculate();
        ip.push_back(x);
    }
    std::set<int> charges;
    charges.insert(charge);
    ip.setCharges(charges);
    return ip;
}

std::vector<Centroid> makeCentroids(size_t n, const double* mz,
    const double* rt, const unsigned int* sn, const double *ab)
{
    std::vector<Centroid> cs;
    for (size_t i = 0; i < n; ++i) {
        Centroid c;
        c.setMz(mz[i]);
        c.setRetentionTime(rt[i]);
        c.setScanNumber(sn[i]);
        c.setAbundance(ab[i]);
        cs.push_back(c);
    }
    return cs;
}

Xic makeXic(size_t n, const double* mz,
    const double* rt, const unsigned int* sn, const double *ab)
{
    std::vector<Centroid> cs = makeCentroids(n, mz, rt, sn, ab);
    Xic xic;
    xic.insert(xic.end(), cs.begin(), cs.end());
    xic.recalculate();
    return xic;
}
