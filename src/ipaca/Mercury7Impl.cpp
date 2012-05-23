/*
 * Mercury.cpp
 *
 * Copyright (c) 2007-2012 Marc Kirchner
 * Copyright (c) 2007 Xinghua Lou
 * Copyright (c) 2007 Bjoern Voss
 * Copyright (c) 2008 Thorben Kroeger
 *
 * Based on the emass implementation of Perttu Haimi.
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
#include <MSTK/ipaca/Mercury7Impl.hpp>
#include <MSTK/common/Error.hpp>
#include <cmath>

using namespace mstk::ipaca;
using namespace mstk;

void detail::Mercury7Impl::convolve(const detail::Spectrum& s1,
    const detail::Spectrum& s2, detail::Spectrum& result) const
{
    // Check if the input is non-empty. We use size() instead of
    // empty() because we need the values later.
    Size n1 = s1.size();
    Size n2 = s2.size();
    if ((n1 + n2) == 0) {
        result.clear();
        return;
    }
    // If one of the inputs is empty, keep the m/z values but set the
    // abundances to zero.
    if (n1 == 0 || n2 == 0) {
        result = (n1 == 0 ? s2 : s1);
        typedef detail::Spectrum::iterator IT;
        for (IT i = result.begin(); i != result.end(); ++i) {
            i->ab = 0.0;
        }
        return;
    }
    // No need to clear out the return values, we will overwrite them
    // anyways. Hence, simply make sure the elements exist and that we
    // will not need to reallocate inside the loop.
    result.resize(n1 + n2 - 1);
    for (size_t k = 0; k < n1 + n2 - 1; k++) {
        double totalAbundance = 0.0;
        double massExpectation = 0.0;
        size_t start = k < (n2 - 1) ? 0 : k - n2 + 1; // max(0, k-n2+1)
        size_t end = k < (n1 - 1) ? k : n1 - 1; // min(n1-1, k)
        // calculate the convolution of the k-th peak with everything else
        for (size_t i = start; i <= end; i++) {
            double ithAbundance = s1[i].ab * s2[k - i].ab;
            if (ithAbundance > 0.0) {
                // calculate the expected mass position
                totalAbundance += ithAbundance;
                double ithMass = s1[i].mz + s2[k - i].mz;
                massExpectation += ithAbundance * ithMass;
            }
        }
        // We cannot simply throw away isotopes with zero probability, as
        // this would mess up the isotope count k.
        result[k].mz = totalAbundance > 0 ? (massExpectation / totalAbundance)
                : 0;
        result[k].ab = totalAbundance;
    }
}

void detail::Mercury7Impl::prune(detail::Spectrum& s, const double limit) const
{
    // This is a private function, hence any call to prune with
    // a non-positve limit is a programming error. Parameter validity
    // must be checked in operator() (which is where it comes in).
    mstk_assert(limit > 0.0, "pruning limit must be strictly positive");
    // prune from the left
    typedef detail::Spectrum::const_iterator CI;
    CI l;
    for (l = s.begin(); l != s.end(); ++l) {
        if (l->ab > limit) {
            break;
        }
    }
    // prune from the right
    typedef detail::Spectrum::const_reverse_iterator RCI;
    RCI r;
    for (r = s.rbegin(); r != s.rend(); ++r) {
        if (r->ab > limit || r.base() == l) {
            break;
        }
    }
    // trim down using the swap trick; should be faster than two copies...
    detail::Spectrum(l, r.base()).swap(s);
}

void detail::Mercury7Impl::integerMercury(
    const detail::Stoichiometry& stoichiometry, const double limit,
    detail::Spectrum& msa) const
{
    mstk_assert(limit > 0.0, "pruning limit must be strictly positive");
    msa.clear();
    detail::Spectrum tmp, esa;
    Bool msa_initialized = false;

    // walk through the elements
    typedef detail::Stoichiometry::const_iterator SCI;
    for (SCI iter = stoichiometry.begin(); iter != stoichiometry.end(); ++iter) {
        // number of atoms at iterator position
        mstk_assert(iter->count >= 0.0, "expect at least one atom");
        Size n = static_cast<Size>(iter->count);
        // if the element is present in the composition,
        // then calculate ESA and update MSA
        if (n) {
            // initialize ESA
            esa.assign(iter->isotopes.begin(), iter->isotopes.end());
            mstk_assert(!esa.empty(), "expect non-empty ESA after assignment");
            while (1) {
                // check if we need to do the MSA update
                if (n & 1) {
                    // MSA update
                    if (msa_initialized) {
                        // normal update
                        convolve(msa, esa, tmp);
                        msa = tmp;
                    } else {
                        // initialize MSA=ESA
                        msa = esa;
                        msa_initialized = true;
                    }
                    prune(msa, limit);
                }
                // the ESA update is always carried out (with the exception of
                // the last time, i.e. when n==1)
                if (n == 1) {
                    break;
                }
                convolve(esa, esa, tmp);
                esa = tmp;
                prune(esa, limit);
                n = n >> 1;
            }
        }
    }
}

void detail::Mercury7Impl::fractionalMercury(const detail::Stoichiometry& s,
    double limit, detail::Spectrum& frac) const
{
    mstk_assert(limit > 0.0, "pruning limit must be strictly positive");
    frac.clear();
    typedef detail::Stoichiometry::const_iterator SCI;
    for (SCI i = s.begin(); i != s.end(); ++i) {
        // the last spectrum
        detail::Spectrum temp(frac);
        // initialize ESA
        detail::Spectrum esa;
        for (Size u = 0; u < i->isotopes.size(); ++u) {
            if (i->isotopes[u].ab <= 0.0) {
                continue;
            }
            detail::SpectrumElement se;
            if (u > 0) {
                se.mz = i->isotopes[u].mz - i->isotopes[0].mz + esa.front().mz;
                se.ab = i->isotopes[u].ab * i->count;
                esa.push_back(se);
            } else {
                se.mz = i->isotopes[0].mz * i->count;
                se.ab = (1 - i->count) + i->isotopes[0].ab * i->count;
            }
            esa.push_back(se);
        }
        if (i == s.begin()) {
            frac = esa;
        } else {
            convolve(esa, temp, frac);
        }
    }
}

detail::Spectrum detail::Mercury7Impl::operator()(
    const detail::Stoichiometry& stoichiometry, const double limit) const
{
    // check the parameters
    mstk_precondition(limit > 0.0, "require positive pruning limit.");
    // split the stoichiometry into integer and fractional parts
    detail::Stoichiometry intStoi;
    detail::Stoichiometry fracStoi;
    detail::splitStoichiometry(stoichiometry, intStoi, fracStoi);
    // check if there is any integer contribution, and calculate the mz and
    // abundance vectors if yes
    detail::Spectrum intSpec;
    bool hasValidIntegerStoichiometry = detail::isPlausibleStoichiometry(
        intStoi);
    if (hasValidIntegerStoichiometry) {
        integerMercury(intStoi, limit, intSpec);
    }
    // check if there is any fractional contribution and calculate the mz and
    // abundance vectors if yes
    detail::Spectrum fracSpec;
    bool hasValidFractionalStoichiometry = detail::isPlausibleStoichiometry(
        fracStoi);
    if (hasValidFractionalStoichiometry) {
        fractionalMercury(fracStoi, limit, fracSpec);
    }
    // if we have integer and fractional contributions, we need to convolve the
    // two; otherwise assign the resepctive non-zero contribution.
    detail::Spectrum result;
    if (hasValidIntegerStoichiometry && hasValidFractionalStoichiometry) {
        Mercury7Impl::convolve(intSpec, fracSpec, result);
        Mercury7Impl::prune(result, limit);
    } else {
        if (hasValidIntegerStoichiometry) {
            result = intSpec;
        } else {
            result = fracSpec;
        }
    }
    return result;
}

double detail::Mercury7Impl::getMonoisotopicMass(
    const detail::Stoichiometry& stoichiometry) const
{
    double mass = 0.0;
    // multiply the number of atoms for each element with the element's
    // monoisotopic mass
    typedef detail::Stoichiometry::const_iterator CI;
    for (CI i = stoichiometry.begin(); i != stoichiometry.end(); ++i) {
        mass += i->count * i->isotopes[0].mz;
    }
    return mass;
}

double detail::Mercury7Impl::getAverageMass(
    const detail::Stoichiometry& stoichiometry) const
{
    double avg = 0.0;
    typedef detail::Stoichiometry::const_iterator CI;

    for (CI i = stoichiometry.begin(); i != stoichiometry.end(); ++i) {
        // get the average mass of the stoichiometry entry
        double entryAvg = 0.0;
        typedef detail::Isotopes::const_iterator ICI;
        for (ICI j = i->isotopes.begin(); j != i->isotopes.end(); ++j) {
            entryAvg += j->mz * j->ab;
        }
        // add to the overall average
        avg += entryAvg;
    }
    return avg;
}

