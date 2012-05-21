/*
 * Xic.cpp
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
#include <MSTK/fe/types/Xic.hpp>

#include <cmath>
#include <algorithm>
#include <limits>
#include <cassert>

#include <MSTK/common/Log.hpp>
#include <MSTK/common/Error.hpp>
#include <MSTK/fe/triangularIntegration.hpp>

namespace mstk {

namespace fe {

bool Xic::LessThanAbundance::operator()(const Xic& lhs, const Xic& rhs) const
{
    return lhs.abundance_ < rhs.abundance_;
}

bool Xic::LessThanRt::operator()(const Xic& lhs, const Xic& rhs) const
{
    return lhs.rt_ < rhs.rt_;
}

bool Xic::LessThanMz::operator()(const Xic& lhs, const Xic& rhs) const
{
    return lhs.mz_ < rhs.mz_;
}

double& Xic::MzAccessor::operator()(Xic& c)
{
    return c.mz_;
}

double Xic::MzAccessor::operator()(const Xic& c) const
{
    return c.mz_;
}

double& Xic::RtAccessor::operator()(Xic& c)
{
    return c.rt_;
}

double Xic::RtAccessor::operator()(const Xic& c) const
{
    return c.rt_;
}

double& Xic::AbundanceAccessor::operator()(Xic& c)
{
    return c.abundance_;
}

double Xic::AbundanceAccessor::operator()(const Xic& c) const
{
    return c.abundance_;
}

Xic::Xic() :
    rt_(0.0), rtSigma_(0.0), mz_(0.0), mzSigma_(0.0), abundance_(0.0)
{
}

Xic::Xic(Xic::const_iterator first, Xic::const_iterator last) :
    rt_(0.0), rtSigma_(0.0), mz_(0.0), mzSigma_(0.0), abundance_(0.0)
{
    this->assign(first, last);
    this->recalculate();
}

bool Xic::operator==(const Xic& rhs)
{
    return rt_ == rhs.rt_ && rtSigma_ == rhs.rtSigma_ && mz_ == rhs.mz_
            && mzSigma_ == rhs.mzSigma_ && abundance_ == rhs.abundance_;
}

double Xic::getAbundance() const
{
    return abundance_;
}

double Xic::getMz() const
{
    return mz_;
}

double Xic::getMzTolerance() const
{
    return mzSigma_;
}

double Xic::getRetentionTime() const
{
    return rt_;
}

double Xic::getRetentionTimeTolerance() const
{
    return rtSigma_;
}

struct Xic::RtAbAccessor
{
    double x(const Centroid& l)
    {
        return l.getRetentionTime();
    }
    double y(const Centroid& l)
    {
        return l.getAbundance();
    }
};

void Xic::mergeDuplicates()
{
    // Single-pass implementation that merges (consecutive) duplicate scan
    // number entries (ie. requires an scanNumber-sorted XIC, which should be
    // the same as an rt-sorting).
    // The code uses the weighted mean of the accurate masses of the
    // merged centroids for the new centroid. It does not drop zero abundance
    // centroids unless they are merged with non-zero abundance centroid. In
    // the latter case, the weighted mean will automatically ignore the m/z
    // contribution of the zero abundance centroid.
    // In cases where there are zero abundance centroids to be merged, the
    // will determine their arithmetic mean for the merged centroid.
    iterator cur = begin();
    iterator p = cur;
    while (p != end()) {
        if (cur != p) {
            *cur = *p;
            ++p;
        }
        double wMz = cur->getMz() * cur->getAbundance();
        double sMz = cur->getMz();
        double sAb = cur->getAbundance();
        while (p != end() && cur->getScanNumber() == p->getScanNumber()) {
            if (cur != p) {
                sAb += p->getAbundance();
                wMz += p->getMz() * p->getAbundance();
                sMz += p->getMz();
            }
            ++p;
        }
        assert(sAb >= 0.0 && "Xic::mergeDiplicates: negative total abundance not allowed.");
        if (sAb > 0.0) {
            // Use weighted mean if there was at least one abundance
            // measurement above zero.
            cur->setMz(wMz/sAb);
        } else {
            // Use normal mean for all-zero measurements.
            // We must not delete zeros because they may influence splitting.
            cur->setMz(sMz/std::distance(cur, p));
        }
        cur->setAbundance(sAb);
        ++cur;
    }
    erase(cur, end());
}

void Xic::recalculate()
{
    // make sure the XIC is sorted by rt
    std::sort(begin(), end(), Centroid::LessThanRt());

    // make sure there are no rt/sn duplicates
    mergeDuplicates();

    // calculate the overall abundance of the XIC; use triangulation over rt
    abundance_ = triangularIntegration<double, const_iterator,
            Xic::RtAbAccessor> (begin(), end(), RtAbAccessor());
    // calculate abundance-weighted mean and variance of all centroids in XIC
    double swm(0.0); // sum of weighted masses
    double swsm(0.0); // sum of weighted squared masses
    double swr(0.0); // sum of weighted retention times
    double swsr(0.0); // sum of weighted retention times
    double sw(0.0); // sum of weights
    double ssw(0.0); // sum of squared weights
    for (const_iterator i = begin(); i != end(); ++i) {
        double mz = i->getMz();
        double ab = i->getAbundance();
        double rt = i->getRetentionTime();
        double swmTmp = ab * mz;
        swm += swmTmp;
        swsm += swmTmp * mz;
        double swrTmp = ab * rt;
        swr += swrTmp;
        swsr += swrTmp * rt;
        sw += ab;
        ssw += ab * ab;
    }
    // maximum position via empirical weighted mean
    mz_ = swm / sw;
    rt_ = swr / sw;
    // empirical weighted variance running sums
    // https://stat.ethz.ch/pipermail/r-help/2008-July/168762.html
    if (size() > 1) {
        // The procedure is not exactly numerically robust and can yield
        // small negative values, which will cause silent NaNs in the call
        // to std::sqrt. From a practical point of view, this means that
        // the calculated variance is very close to zero. We simply use
        // std::abs to avoid the NaNs.
        double m = ((swsm * sw) - (swm * swm)) / ((sw * sw) - ssw);
        double r = ((swsr * sw) - (swr * swr)) / ((sw * sw) - ssw);
        if (m < 0) {
            MSTK_LOG(logDEBUG2)
                    << "XIC::recalculate: stabilizing m/z "
                        "variance estimate: m=" << m;
            m = std::abs(m);
        }
        if (r < 0) {
            MSTK_LOG(logDEBUG2)
                    << "XIC::recalculate: stabilizing rt "
                        "variance estimate: r=" << r;
            r = std::abs(r);
        }
        mzSigma_ = std::sqrt(m);
        rtSigma_ = std::sqrt(r);
    } else {
        MSTK_LOG(logWARNING)
                << "Xic with " << size() << " entries!";
        mzSigma_ = 0.0;
        rtSigma_ = 0.0;
    }
    // MaxQuant argues that the measurements are not independent and
    // consequently uses a bootstrap estimate here.
    // TODO: their arguments are compelling. We should do the same.
}

void Xic::getSmoothedXic(Xic& xic)
{
    // make sure the XIC is sorted by rt
    std::sort(begin(), end(), Centroid::LessThanRt());
    // check that the XIC has at least the size of the structuring element
    if (size() < 3) {
        MSTK_LOG(logWARNING)
                << "Xic::getSmoothedXic: XIC too short to smooth: n=" << size();
        xic = *this;
        return;
    }
    // push the first element
    xic.push_back(*(begin()));
    // position the structuring element
    const_iterator l = begin();
    const_iterator i = begin();
    std::advance(i, 1);
    const_iterator r = begin();
    std::advance(r, 2);
    // apply a running mean on the observed XIC: we use a sparse running mean in
    // which the contributions are weighted be the relative distances of the
    // left and right neighbor 
    while (r != end()) {
        double lRt = l->getRetentionTime();
        double iRt = i->getRetentionTime();
        double rRt = r->getRetentionTime();
        // get the span between l and r and left/right distances from i
        double span = rRt - lRt;
        double dl = iRt - lRt;
        double dr = rRt - iRt;
        Centroid tmp(*i);
        double abundance = (2 * ((1.0 - dl / span) * l->getAbundance() + (1.0
                - dr / span) * r->getAbundance()) + i->getAbundance()) / 3.0;
        tmp.setAbundance(abundance);
        xic.push_back(tmp);
        ++l;
        ++i;
        ++r;
    }
    // push the last element
    xic.push_back(*(rbegin()));
}

void Xic::split(std::vector<Xic>& splitXics, double mindepth)
{
    // with less than 4 measurements, there is no point in splitting because we
    // will always generate a single measurement XIC
    if (size() < 4) {
        splitXics.push_back(*this);
        MSTK_LOG(logDEBUG3)
                << __FUNCTION__ << ": size too small (" << size() << "<4)";
        return;
    }
    // get a smoothed version of our own XIC
    Xic smoothXic;
    getSmoothedXic(smoothXic);
    // iterate through the smoothed XIC, find minima and maxima
    const_iterator l = smoothXic.begin();
    const_iterator i = l + 1;
    const_iterator r = l + 2;
    const_iterator currentMin = smoothXic.end();
    const_iterator previousMax = smoothXic.end();
    const_iterator nextMax = smoothXic.end();
    // deal with the left boundary: the left boundary is the first
    // lastMin; it can also be the first previousMax, depending if we first hit
    // a local max or a local min later.
    const_iterator lastMin = smoothXic.begin();
    if (i->getAbundance() <= l->getAbundance()) {
        // the first XIC value is a local max
        previousMax = l;
    }

    while (r != smoothXic.end() - 1) {
        // minimize number of function calls
        double lAb = l->getAbundance();
        double iAb = i->getAbundance();
        double rAb = r->getAbundance();
        // check if we hit a minimum
        if (lAb >= iAb && iAb < rAb) {
            // when we hit a minimum, we just need to record it.
            currentMin = i;
            MSTK_LOG(logDEBUG3)
                    << __FUNCTION__ << ": setting currentMin to rt="
                            << i->getRetentionTime();
        }
        // check if we hit a maximum
        if (lAb < iAb && iAb >= rAb) {
            MSTK_LOG(logDEBUG3)
                    << __FUNCTION__ << ": local max at rt="
                            << i->getRetentionTime();
            // first thing we need to check is if this is the first time we find
            // a maximum: 
            if (previousMax != smoothXic.end()) {
                // this is the second max in the max-min-max combination
                nextMax = i;
                MSTK_LOG(logDEBUG3)
                        << __FUNCTION__ << ": setting nextMax to rt="
                                << nextMax->getRetentionTime();

                // Check if the minimum is deep enough. We use the same criterion as
                // MaxQuant (see Supplementary Info for Cox & Mann, 2008)
                if (currentMin->getAbundance() < mindepth * (std::min)(
                    previousMax->getAbundance(), nextMax->getAbundance())) {
                    // copy the original XIC to a subXIC: all XXXXminIt iterators
                    // address the smoothedXic - hence, we need to find the position
                    // of the minimum in the true XIC. This can be done by simple
                    // distance calculations (the XICs have the same length).
                    iterator first = begin();
                    iterator last = first;
                    std::advance(
                        first,
                        std::distance(
                            static_cast<const_iterator> (smoothXic.begin()),
                            lastMin));
                    // the minimum itself goes to the right
                    std::advance(
                        last,
                        std::distance(
                            static_cast<const_iterator> (smoothXic.begin()),
                            currentMin));
                    // copy, but make sure that we do not copy singles
                    if (std::distance(first, last) > 1) {
                        Xic subXic;
                        subXic.assign(first, last);
                        subXic.recalculate();
                        splitXics.push_back(subXic);
                        // advance the iterators
                        lastMin = currentMin;
                    }
                    previousMax = nextMax;
                    MSTK_LOG(logDEBUG3)
                            << __FUNCTION__ << ": shifting previousMax to rt="
                                    << previousMax->getRetentionTime();

                } else {
                    MSTK_LOG(logDEBUG3)
                            << __FUNCTION__ << "failed criterion: "
                                    << previousMax->getAbundance() << " > "
                                    << iAb << ", < "
                                    << nextMax->getAbundance()
                                    << ", mindepth: " << mindepth
                                    << ", critval: " << mindepth * (std::min)(
                                previousMax->getAbundance(),
                                nextMax->getAbundance());
                    if (nextMax->getAbundance() > previousMax->getAbundance()) {
                        previousMax = nextMax;
                        MSTK_LOG(logDEBUG3)
                                << __FUNCTION__ << "Not shifting max";
                    }
                }
                // keep the max (no matter what)
                //previousMax = nextMax;
                MSTK_LOG(logDEBUG3)
                        << __FUNCTION__ << ": shifting previousMax to rt="
                                << previousMax->getRetentionTime();
            } else {
                // the first max we found. Record it and continue.
                previousMax = i;
                MSTK_LOG(logDEBUG3)
                        << __FUNCTION__ << ": previousMax set to rt="
                                << i->getRetentionTime();
            }
        }
        ++l;
        ++i;
        ++r;
    }
    // now deal with what is left 
    iterator first = begin();
    std::advance(
        first,
        std::distance(static_cast<const_iterator> (smoothXic.begin()), lastMin));
    // the criterion in the loop guarantees that we never copy singles
    Xic subXic;
    subXic.assign(first, end());
    if (subXic.size() <= 1) {
        MSTK_LOG(logWARNING)
                << "Splitting generated a size " << subXic.size() << " Xic.";
    }
    subXic.recalculate();
    splitXics.push_back(subXic);
}

double Xic::correlate(Xic& rhs)
{
    // make sure the Xics are sorted in rt
    // FIXME: Xic should not be derived from Collection<T> and
    //        should instead guarantee the sorting of its elements
    std::sort(begin(), end(), Centroid::LessThanScanNumber());
    std::sort(rhs.begin(), rhs.end(), Centroid::LessThanScanNumber());
    // get pointers to the centroids of the XICs
    iterator lhsCentroidIt = begin();
    iterator rhsCentroidIt = rhs.begin();
    double lr(0.0), lsq(0.0), rsq(0.0);
    unsigned int lhsSn, rhsSn;
    double l, r;
    while ((lhsCentroidIt != end()) || (rhsCentroidIt != rhs.end())) {
        if (lhsCentroidIt != end()) {
            lhsSn = lhsCentroidIt->getScanNumber();
            l = lhsCentroidIt->getAbundance();
        } else {
            lhsSn = std::numeric_limits<unsigned int>::max();
            l = 0.0;
        }
        if (rhsCentroidIt != rhs.end()) {
            rhsSn = rhsCentroidIt->getScanNumber();
            r = rhsCentroidIt->getAbundance();
        } else {
            rhsSn = std::numeric_limits<unsigned int>::max();
            r = 0.0;
        }
        if (lhsSn == rhsSn) {
            // same scan number
            lr += l * r;
            lsq += l * l;
            rsq += r * r;
            ++lhsCentroidIt;
            ++rhsCentroidIt;
        } else
            if (lhsSn < rhsSn) {
                // no rhs measurement at lhs scan number, hence add lhs only
                lsq += l * l;
                ++lhsCentroidIt;
            } else {
                // no lhs measurement at rhs scan number, hence add rhs only
                rsq += r * r;
                ++rhsCentroidIt;
            }
    }
    return (lr != 0) ? (lr / std::sqrt(lsq * rsq)) : 0.0;
}

std::ostream& operator<<(std::ostream& os, const Xic& x)
{
    os << x.getRetentionTime() << "\t" << x.getRetentionTimeTolerance()
            << "\t" << x.getMz() << "\t" << x.getMzTolerance() << "\t"
            << x.getAbundance();
    return os;
}

} // namespace fe

} // namespace mstk

