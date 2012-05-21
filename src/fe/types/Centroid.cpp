/*
 * Centroid.cpp
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
#include <MSTK/fe/types/Centroid.hpp>
#include <MSTK/common/Error.hpp>

namespace mstk {

namespace fe {

double& Centroid::MzAccessor::operator()(Centroid& c)
{
    return c.mz_;
}

double Centroid::MzAccessor::operator()(const Centroid& c) const
{
    return c.mz_;
}

double& Centroid::RtAccessor::operator()(Centroid& c)
{
    return c.rt_;
}

double Centroid::RtAccessor::operator()(const Centroid& c) const
{
    return c.rt_;
}

double& Centroid::AbundanceAccessor::operator()(Centroid& c)
{
    return c.ab_;
}

double Centroid::AbundanceAccessor::operator()(const Centroid& c) const
{
    return c.ab_;
}

Centroid::Centroid() :
    rt_(0.0), mz_(0.0), sn_(0), ab_(0.0), raw_()
{
}

Centroid::Centroid(const Centroid& rhs) :
        rt_(rhs.rt_), mz_(rhs.mz_), sn_(rhs.sn_), ab_(rhs.ab_), raw_(rhs.raw_)
{
}

Centroid::Centroid(Double retentionTime, Double mz, UnsignedInt sn, Double ab,
    Spectrum::const_iterator first, Spectrum::const_iterator last) :
    rt_(retentionTime), mz_(mz), sn_(sn), ab_(ab), raw_(first, last)
{
    mstk_precondition(retentionTime >= 0.0,
            "mstk::Centroid retention times cannot be negative.");
    mstk_precondition(mz >= 0.0,
            "mstk::Centroid m/z ratios cannot be negative.");
    mstk_precondition(ab >= 0.0,
            "mstk::Centroid abundance cannot be negative.");
}

bool Centroid::operator==(const Centroid& rhs) const
{
    return (rt_ == rhs.rt_ && mz_ == rhs.mz_ && sn_ == rhs.sn_ && ab_
            == rhs.ab_ && raw_ == rhs.raw_);
}

Centroid::~Centroid()
{
}

Double Centroid::getRetentionTime() const
{
    return rt_;
}

void Centroid::setRetentionTime(const Double rt)
{
    mstk_precondition(rt >= 0.0,
            "mstk::Centroid retention times cannot be negative.");
    rt_ = rt;
}

Double Centroid::getMz() const
{
    return mz_;
}

void Centroid::setMz(const Double mz)
{
    mstk_precondition(mz >= 0.0,
            "mstk::Centroid m/z ratios cannot be negative.");
    mz_ = mz;
}

UnsignedInt Centroid::getScanNumber() const
{
    return sn_;
}

void Centroid::setScanNumber(const UnsignedInt sn)
{
    sn_ = sn;
}

Double Centroid::getAbundance() const
{
    return ab_;
}

void Centroid::setAbundance(const Double ab)
{
    mstk_precondition(ab >= 0.0,
            "mstk::Centroid abundance cannot be negative.");
    ab_ = ab;
}

const Spectrum& Centroid::getRawData() const
{
    return raw_;
}

void Centroid::setRawData(const Spectrum& s)
{
    raw_ = s;
}

bool Centroid::LessThanRt::operator()(const Centroid& lhs, const Centroid& rhs)
{
    return lhs.getRetentionTime() < rhs.getRetentionTime();
}

bool Centroid::LessThanMz::operator()(const Centroid& lhs, const Centroid& rhs)
{
    return lhs.getMz() < rhs.getMz();
}

bool Centroid::LessThanScanNumber::operator()(const Centroid& lhs,
    const Centroid& rhs)
{
    return lhs.getScanNumber() < rhs.getScanNumber();
}

bool Centroid::LessThanAbundance::operator()(const Centroid& lhs,
    const Centroid& rhs)
{
    return lhs.getAbundance() < rhs.getAbundance();
}

std::ostream& operator<<(std::ostream& os, const Centroid& c)
{
    os << c.getRetentionTime() << '\t' << c.getScanNumber() << '\t'
            << c.getMz() << '\t' << c.getAbundance();
    return os;
}

} // namespace fe

} // namespace mstk

