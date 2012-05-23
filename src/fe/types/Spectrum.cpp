/*
 * Spectrum.cpp
 *
 * Copyright (c) 2007-2011 Marc Kirchner
 * Copyright (c) 2007 Bjoern Voss
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

#include <MSTK/fe/types/Spectrum.hpp>
#include <MSTK/common/Error.hpp>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

namespace mstk {

namespace fe {

Spectrum::Spectrum()
  : Collection<SpectrumElement>(), rt_(0.0), msLevel_(0), 
  scanNumber_(0), totalIonCurrent_(0.0), precursorScanNumber_(0),
  precursorMz_(0.0), precursorCharge_(0), precursorAbundance_(0.0)
{}

/** Range constructor.
 */
Spectrum::Spectrum(iterator first, iterator last) 
  : Collection<SpectrumElement>(first, last),
  rt_(0.0), msLevel_(0), 
  scanNumber_(0), totalIonCurrent_(0.0), precursorScanNumber_(0),
  precursorMz_(0.0), precursorCharge_(0), precursorAbundance_(0.0)
{} 

/** Range constructor.
 */
Spectrum::Spectrum(const_iterator first, const_iterator last) 
  : Collection<SpectrumElement>(first, last),
  rt_(0.0), msLevel_(0), 
  scanNumber_(0), totalIonCurrent_(0.0), precursorScanNumber_(0),
  precursorMz_(0.0), precursorCharge_(0), precursorAbundance_(0.0)
{}

Spectrum::Spectrum(const std::vector<double>& mz, 
  const std::vector<double>& abundances)
  : Collection<SpectrumElement>()
{
    // call default constructor
    clear();
    // read from vectors
    if (mz.size() == abundances.size()) {
        for (unsigned int i = 0; i < mz.size(); i++) {
            c_.push_back(Element(mz[i], abundances[i]));
        }
    }
}

// destructor
Spectrum::~Spectrum()
{
}

Spectrum& Spectrum::operator=(const Spectrum& rhs)
{
    if (this != &rhs) {
        this->clear();
        rt_ = rhs.rt_;
        msLevel_ = rhs.msLevel_;
        scanNumber_ = rhs.scanNumber_;
        totalIonCurrent_ = rhs.totalIonCurrent_;
        precursorScanNumber_ = rhs.precursorScanNumber_;
        precursorMz_ = rhs.precursorMz_;
        precursorCharge_ = rhs.precursorCharge_;
        precursorAbundance_ = rhs.precursorAbundance_;
        c_ = rhs.c_;    
    }
    return *this;
}
   
bool Spectrum::operator==(const Spectrum& s) const {
    return c_ == s.c_ && msLevel_ == s.msLevel_ && rt_ == s.rt_
      && s.scanNumber_ == scanNumber_ && s.totalIonCurrent_ == totalIonCurrent_
      && s.precursorScanNumber_ == precursorScanNumber_
      && s.precursorMz_ == precursorMz_
      && s.precursorCharge_ == precursorCharge_
      && s.precursorAbundance_ == precursorAbundance_;
}

void Spectrum::clear() {
    c_.clear();
    rt_ = 0.0;
    msLevel_ = 0;
    scanNumber_ = 0;
    totalIonCurrent_ = 0.0;
    precursorScanNumber_ = 0;
    precursorMz_ = 0.0;
    precursorCharge_ = 0;
    precursorAbundance_ = 0.0;
}

Spectrum Spectrum::subset(const double beginMz, 
  const double endMz) const
{
    const_iterator first = std::lower_bound(begin(), end(), beginMz, 
      LessThanMz<Spectrum::Element, double>());
    const_iterator last  = std::upper_bound(first, end(), endMz, 
      LessThanMz<double, Spectrum::Element>());
    Spectrum s(first, last);
    s.setRetentionTime(rt_);
    s.setMsLevel(msLevel_);
    s.setScanNumber(scanNumber_);
    s.setTotalIonCurrent(s.getTotalAbundance());
    s.setPrecursorScanNumber(precursorScanNumber_);
    s.setPrecursorMz(precursorMz_);
    s.setPrecursorCharge(precursorCharge_);
    s.setPrecursorAbundance(precursorAbundance_);
    return s;
}

void Spectrum::shiftTo(const double to)
{
    if (!empty()) {
        double diff = to - c_[0].mz;
        transform(c_.begin(), c_.end(), c_.begin(), 
          Spectrum::ShiftMz(diff));
    }
}

void Spectrum::shiftBy(const double by)
{
    transform(c_.begin(), c_.end(), c_.begin(), Spectrum::ShiftMz(by));
}

void Spectrum::shiftMaxToMonoisotopicMass(void)
{
    if (!empty()) {
        iterator maxIdx = std::max_element(c_.begin(), c_.end(), 
          LessThanAbundance<Spectrum::Element, Spectrum::Element>());
        if (maxIdx != c_.begin()) {
            shiftBy(c_[0].mz - (maxIdx->mz));
        }
    }
}

Spectrum::const_iterator Spectrum::getMaxAbundancePeak() const
{
    return std::max_element(begin(), end(), Spectrum::LessThanAbundance<
      Spectrum::Element, Spectrum::Element>());
}

Spectrum::iterator Spectrum::getMaxAbundancePeak()
{
    return std::max_element(begin(), end(), Spectrum::LessThanAbundance<
      Spectrum::Element, Spectrum::Element>());
}

double Spectrum::getTotalAbundance()
{
    return std::accumulate(begin(), end(), 0.0, Spectrum::SumAbundance());
}

void Spectrum::merge(const Spectrum& other)
{
    // FIXME: check ms levels and do something smart w/ the member variables
    iterator s1 = begin();
    const_iterator s2 = other.begin();

    while (s1 != end() && s2 != other.end()) {
        if (s1->mz < s2->mz) {
            //s1 is behind
            ++s1;
        }
        else if (s1->mz > s2->mz) {
            //s2 is behind
            s1 = insert(s1, Element(s2->mz, s2->abundance));
            ++s2;
        }
        else { /*s1->mz == s2->mz */
            *s1 = Element(s1->mz, s1->abundance + s2->abundance);
            ++s1;
            ++s2;
        }
    }
    // if we terminated because we reached the end of
    // this->sepctrum_, add the remainder
    std::copy(s2, other.end(), std::back_inserter(c_));
}

Spectrum Spectrum::removeDuplicates(double tol)
{
    if (empty() || size() == 1) {
        return *this;
    }
    Spectrum unique;
    Spectrum::iterator i=this->begin();
    Spectrum::iterator j=i+1;
    while (i != this->end()){
        while (j->mz - i->mz < tol && j!=this->end())
            ++j;
        double sumab = 0;
        double wmz = 0;
        for (Spectrum::iterator temp = i; temp != j; ++temp) {
            sumab += temp->abundance;
            wmz += temp->mz * temp->abundance;
        }
        wmz /= sumab;
        Spectrum::Element el(wmz, sumab);
        unique.push_back(el);
        i=j;
        if (i != this->end()) {
            j=i+1;
        }            
    }
    return unique;
}

std::ostream& operator<<(std::ostream& os, Spectrum& p)
{
    for (Spectrum::iterator i = p.begin(); i != p.end(); ++i) {
        os << i->mz << " " << i->abundance << std::endl;
    }
    return(os);
}

std::istream& operator>>(std::istream& is, Spectrum& s)
{
    double mz, ab;
    if (is.good()) {
        while (is >> mz >> ab) {
            // only push_back if abundance is > 0
            if (ab > 0.0) {
                s.push_back(Spectrum::Element(mz, ab));
            }
        }
    }

    return is;
}

std::ostream& operator<<(std::ostream& os, Spectrum::Element& e) {
    os << e.mz << " " << e.abundance;
    return os;
}

} // namespace fe

} // namespace mstk

