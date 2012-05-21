/*
 * Spectrum.hpp
 *
 * Copyright (c) 2007-2012 Marc Kirchner
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
#ifndef __MSTK_INCLUDE_MSTK_FE_TYPES_SPECTRUM_HPP__
#define __MSTK_INCLUDE_MSTK_FE_TYPES_SPECTRUM_HPP__
#include <MSTK/config.hpp>

#include <MSTK/common/Collection.hpp>
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace mstk {

namespace fe {

/** A single entry in a mass spectrum: a pair of \f$mz\f$ and abundance
 * values.
 */
struct SpectrumElement
{
    typedef double first_type;
    typedef double second_type;

    struct MzAccessor
    {
        double operator()(const SpectrumElement& e) const
        {
            return e.mz;
        }
    };

    struct AbundanceAccessor
    {
        double operator()(const SpectrumElement& e) const
        {
            return e.abundance;
        }
    };

    SpectrumElement(const SpectrumElement& rhs) :
        mz(rhs.mz), abundance(rhs.abundance)
    {
    }

    SpectrumElement(const double m, const double a) :
        mz(m), abundance(a)
    {
    }
    SpectrumElement(const std::pair<double, double>& p) :
        mz(p.first), abundance(p.second)
    {
    }

    bool operator==(const SpectrumElement& p) const
    {
        return mz == p.mz && abundance == p.abundance;
    }

    double mz;
    double abundance;
};

/**
 * A sparse representation of a mass spectrum.
 */
class MSTK_EXPORT Spectrum : public Collection<SpectrumElement>
{
public:
    typedef SpectrumElement Element;
    // for use with SpectrumTraits<T>:
    typedef SpectrumElement Value;
    /** Compare two mstk::Spectrum::Elements based on their \f$mz\f$ values
     */
    template<class A, class B>
    struct LessThanMz : std::binary_function<A, B, bool>
    {
        bool operator()(const A& lhs, const B& rhs);
    };

    /** Compare two mstk::Spectrum::Element based on their abundance values
     */
    template<class A, class B>
    struct LessThanAbundance : std::binary_function<A, B, bool>
    {
        bool operator()(const A& lhs, const B& rhs);
    };

    /** Compare two mstk::Spectrum objects based on the \f$mz\f$ values
     *  at their respective maximal abundances.
     */
    template<class A, class B>
    struct LessThanMzAtMaxAbundance : std::binary_function<A, B, bool>
    {
        bool operator()(const A& lhs, const B& rhs);
    };

    /** Compare two mstk::Spectrum objects based on their retention time
     */
    template<class A, class B>
    struct LessThanRt : std::binary_function<A, B, bool>
    {
        bool operator()(const A& lhs, const B& rhs);
    };

    /** Compare two mstk::Spectrum objects based on their MS level
     */
    template<class A, class B>
    struct LessThanMsLevel : std::binary_function<A, B, bool>
    {
        bool operator()(const A& lhs, const B& rhs);
    };

    /** Check if two mstk::Spectrum objects have the same MS level
     */
    struct EqualMsLevel : std::unary_function<Spectrum, bool>
    {
        EqualMsLevel(unsigned int msLevel) :
            msLevel_(msLevel)
        {
        }
        bool operator()(const Spectrum& s) const
        {
            return msLevel_ == s.msLevel_;
        }
        unsigned int msLevel_;
    };

    /** Functor to accumulate the abundance over a Spectrum object.
     */
    struct SumAbundance : std::binary_function<double,
            Spectrum::value_type, double>
    {
        double operator()(const double& val,
            const Spectrum::value_type& e)
        {
            return val + e.abundance;
        }
    };

    /**
     * functor to add a constant \f$mz\f$ value to a mstk::Spectrum::Element
     */
    struct ShiftMz : std::unary_function<value_type, value_type>
    {
        ShiftMz(double val) :
            val_(val)
        {
        }
        value_type operator()(const value_type& p)
        {
            return (Element(p.mz + val_, p.abundance));
        }
        double val_;
    };

    /** Default Constructor. Constructs an empty @Spectrum.
     */
    Spectrum();
    Spectrum(const std::vector<double>& mz,
        const std::vector<double>& abundances);

    /** Range constructor.
     */
    Spectrum(iterator first, iterator last);

    /** Range constructor.
     */
    Spectrum(const_iterator first, const_iterator last);

    /** Destructor.
     */
    ~Spectrum();

    /** Assignment operator.
     */
    Spectrum& operator=(const Spectrum& rhs);

    /** Test if two @Spectrum object are equal.
     */
    bool operator==(const Spectrum& s) const;

    /* Clear the @Spectrum object. Deletes all elements in the object
     * and resets all member values to zero.
     */
    void clear();

    void setRetentionTime(const double rt);

    double getRetentionTime() const;

    void setMsLevel(const unsigned int l);

    unsigned int getMsLevel() const;

    void setScanNumber(const unsigned int scanNumber);

    unsigned int getScanNumber() const;

    void setTotalIonCurrent(const double totalIonCurrent);

    double getTotalIonCurrent() const;

    void setPrecursorScanNumber(const unsigned int psn);

    unsigned int getPrecursorScanNumber() const;

    void setPrecursorMz(const double pmz);

    double getPrecursorMz() const;

    void setPrecursorCharge(const int pz);

    int getPrecursorCharge() const;

    void setPrecursorAbundance(const double pab);

    double getPrecursorAbundance() const;

    /** Get an iterator pointing at the maximum abundance peak.
     *   @return Iterator to the maximum abundance element
     */
    iterator getMaxAbundancePeak();

    /** Get an iterator pointing at the maximum abundance peak.
     *   @return Iterator to the maximum abundance element
     */
    const_iterator getMaxAbundancePeak() const;

    /** Get the sum over all abundances. This does *not* triangulate.
     *   @return The accumulated abundance.
     */
    double getTotalAbundance();

    /** Merges two Spectrum objects.
     *   @param Spectrum instance to merge with.
     *   @pre Requires the current instance and the merge instance to be sorted.
     *
     * This merges two SparseSpectra, combining their abundances if
     * they exhibit identical masses.
     */
    void merge(const Spectrum& ss);

    /** Remove duplicate m/z entries with tolerance.
     *   Add up the abundances
     */
    Spectrum removeDuplicates(double tol = 0.0016);

    /** Splice a Spectrum into two spectra.
     *   @param first Iterator to first Spectrum::Element that is
     *          supposed to be moved
     *   @param last Iterator to the first SparseSepctrum::Element that
     *          is NOT to be moved (or end(), following STL semantics)
     *   @note The elements are copied (the underlying data structure
     *         is a vector, not a list)
     */
    template<typename Out>
    void splice(const iterator first, const iterator last, Out out);

    /** @return the Spectrum containing the m/z range [beginMz, endMz)
     *  @param beginMz The lower bound of mz range
     *  @param endMz The upper bound of mz range
     *  @return A new Spectrum object
     */
    Spectrum subset(const double beginMz, const double endMz) const;

    void shiftTo(const double to);
    void shiftBy(const double diff);
    void shiftMaxToMonoisotopicMass();

private:
    double rt_;
    unsigned int msLevel_;
    unsigned int scanNumber_;
    double totalIonCurrent_;
    unsigned int precursorScanNumber_;
    double precursorMz_;
    int precursorCharge_;
    double precursorAbundance_;
}; // class Spectrum

MSTK_EXPORT std::istream
        & operator>>(std::istream& is, Spectrum& s);
MSTK_EXPORT std::ostream
        & operator<<(std::ostream& os, Spectrum& s);
MSTK_EXPORT std::ostream& operator<<(std::ostream& os,
    Spectrum::Element& e);

///
/// inline functions
///

inline void Spectrum::setRetentionTime(const double rt)
{
    rt_ = rt;
}

inline double Spectrum::getRetentionTime() const
{
    return rt_;
}

inline void Spectrum::setMsLevel(const unsigned int l)
{
    msLevel_ = l;
}

inline unsigned int Spectrum::getMsLevel() const
{
    return msLevel_;
}

inline void Spectrum::setScanNumber(const unsigned int scanNumber)
{
    scanNumber_ = scanNumber;
}

inline unsigned int Spectrum::getScanNumber() const
{
    return scanNumber_;
}

inline void Spectrum::setTotalIonCurrent(const double totalIonCurrent)
{
    totalIonCurrent_ = totalIonCurrent;
}

inline double Spectrum::getTotalIonCurrent(void) const
{
    return totalIonCurrent_;
}

inline void Spectrum::setPrecursorScanNumber(const unsigned int psn)
{
    precursorScanNumber_ = psn;
}

inline unsigned int Spectrum::getPrecursorScanNumber() const
{
    return precursorScanNumber_;
}

inline void Spectrum::setPrecursorMz(const double pmz)
{
    precursorMz_ = pmz;
}

inline double Spectrum::getPrecursorMz() const
{
    return precursorMz_;
}

inline void Spectrum::setPrecursorCharge(const int pz)
{
    precursorCharge_ = pz;
}

inline int Spectrum::getPrecursorCharge() const
{
    return precursorCharge_;
}

inline void Spectrum::setPrecursorAbundance(const double pab)
{
    precursorAbundance_ = pab;
}

inline double Spectrum::getPrecursorAbundance() const
{
    return precursorAbundance_;
}

//
// template function definitions
//

template<>
struct Spectrum::LessThanAbundance<Spectrum::Element,
        Spectrum::Element> : std::binary_function<Element, Element, bool>
{
    bool operator()(const Element& lhs, const Spectrum::Element& rhs) const
    {
        return lhs.abundance < rhs.abundance;
    }
};

template<>
struct Spectrum::LessThanAbundance<Spectrum::Element, double> : std::binary_function<
        Element, double, bool>
{
    bool operator()(const Element& lhs, const double& rhs) const
    {
        return lhs.abundance < rhs;
    }
};

template<>
struct Spectrum::LessThanAbundance<double, Spectrum::Element> : std::binary_function<
        double, Element, bool>
{
    bool operator()(const double& lhs, const Element& rhs) const
    {
        return lhs < rhs.abundance;
    }
};

template<>
struct Spectrum::LessThanMz<Spectrum::Element,
        Spectrum::Element> : std::binary_function<Element, Element, bool>
{
    bool operator()(const Element& lhs, const Element& rhs) const
    {
        return lhs.mz < rhs.mz;
    }
};

template<>
struct Spectrum::LessThanMz<Spectrum::Element, double> : std::binary_function<
        Element, double, bool>
{
    bool operator()(const Element& lhs, const double& rhs) const
    {
        return lhs.mz < rhs;
    }
};

template<>
struct Spectrum::LessThanMz<double, Spectrum::Element> : std::binary_function<
        double, Element, bool>
{
    bool operator()(const double& lhs, const Element& rhs) const
    {
        return lhs < rhs.mz;
    }
};

template<>
struct Spectrum::LessThanMzAtMaxAbundance<Spectrum, Spectrum> : std::binary_function<
        Spectrum, Spectrum, bool>
{
    bool operator()(const Spectrum& lhs, const Spectrum& rhs) const
    {
        Spectrum::const_iterator li = std::max_element(lhs.begin(),
            lhs.end(), Spectrum::LessThanAbundance<
                    Spectrum::Element, Spectrum::Element>());
        Spectrum::const_iterator ri = std::max_element(rhs.begin(),
            rhs.end(), Spectrum::LessThanAbundance<
                    Spectrum::Element, Spectrum::Element>());
        return li->mz < ri->mz;
    }
};

template<>
struct Spectrum::LessThanMzAtMaxAbundance<double, Spectrum> : std::binary_function<
        double, Spectrum, bool>
{
    bool operator()(const double lhs, const Spectrum& rhs) const
    {
        Spectrum::const_iterator ri = std::max_element(rhs.begin(),
            rhs.end(), Spectrum::LessThanAbundance<
                    Spectrum::Element, Spectrum::Element>());
        return lhs < ri->mz;
    }
};

template<>
struct Spectrum::LessThanMzAtMaxAbundance<Spectrum, double> : std::binary_function<
        Spectrum, double, bool>
{
    bool operator()(const Spectrum& lhs, const double rhs) const
    {
        Spectrum::const_iterator li = std::max_element(lhs.begin(),
            lhs.end(), Spectrum::LessThanAbundance<
                    Spectrum::Element, Spectrum::Element>());
        return li->mz < rhs;
    }
};

template<>
struct Spectrum::LessThanMsLevel<Spectrum, Spectrum> : std::binary_function<
        Spectrum, Spectrum, bool>
{
    bool operator()(const Spectrum& lhs, const Spectrum& rhs) const
    {
        return lhs.getMsLevel() < rhs.getMsLevel();
    }
};

template<>
struct Spectrum::LessThanMsLevel<Spectrum, double> : std::binary_function<
        Spectrum, double, bool>
{
    bool operator()(const Spectrum& lhs, const double& rhs) const
    {
        return lhs.getMsLevel() < rhs;
    }
};

template<>
struct Spectrum::LessThanMsLevel<double, Spectrum> : std::binary_function<
        double, Spectrum, bool>
{
    bool operator()(const double& lhs, const Spectrum& rhs) const
    {
        return lhs < rhs.getMsLevel();
    }
};

template<>
struct Spectrum::LessThanRt<Spectrum, Spectrum> : std::binary_function<
        Spectrum, Spectrum, bool>
{
    bool operator()(const Spectrum& lhs, const Spectrum& rhs) const
    {
        return lhs.getRetentionTime() < rhs.getRetentionTime();
    }
};

template<>
struct Spectrum::LessThanRt<Spectrum, double> : std::binary_function<
        Spectrum, double, bool>
{
    bool operator()(const Spectrum& lhs, const double& rhs) const
    {
        return lhs.getRetentionTime() < rhs;
    }
};
template<>
struct Spectrum::LessThanRt<double, Spectrum> : std::binary_function<
        double, Spectrum, bool>
{
    bool operator()(const double& lhs, const Spectrum& rhs) const
    {
        return lhs < rhs.getRetentionTime();
    }
};

template<typename Out>
void Spectrum::splice(const iterator first, const iterator last, Out out)
{
    std::copy(first, last, out);
    c_.erase(first, last);
}

} /* namespace fe */

} /* namespace mstk */

#endif /*__MSTK_INCLUDE_MSTK_FE_TYPES_SPECTRUM_HPP__*/
