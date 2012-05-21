/*
 * Centroid.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_TYPES_CENTROID_HPP__
#define __MSTK_INCLUDE_MSTK_FE_TYPES_CENTROID_HPP__

#include <MSTK/config.hpp>
#include <MSTK/common/Types.hpp>
#include <MSTK/fe/types/Spectrum.hpp>
#include <iosfwd>

namespace mstk {

namespace fe {

/** Representation of an m/z centroid, including retention time and scan number.
 */
class Centroid
{
public:
    /** Functor to compare two centroids based on their retention time.
     */
    struct LessThanRt : public std::binary_function<Centroid, Centroid, bool>
    {
        /** Comparison operator.
         * @param[in] lhs Left-hand-side centroid object for the comparison.
         * @param[in] rhs Right-hand-side centroid object for the comparison.
         * @return True if the left hand side object has a retention time that is
         *         smaller than the retention time of the right hand side centroid
         *         object.
         */
        bool operator()(const Centroid& lhs, const Centroid& rhs);
    };

    /** Functor to compare two centroids based on their mass/charge ratio.
     */
    struct LessThanMz : public std::binary_function<Centroid, Centroid, bool>
    {
        /** Comparison operator.
         * @param[in] lhs Left-hand-side centroid object for the comparison.
         * @param[in] rhs Right-hand-side centroid object for the comparison.
         * @return True if the left hand side object has a mass/charge ratio that is
         *         smaller than the mass/charge ratio of the right hand side centroid
         *         object.
         */
        bool operator()(const Centroid& lhs, const Centroid& rhs);
    };

    /** Functor to compare two centroids based on their retention time.
     */
    struct LessThanScanNumber : public std::binary_function<Centroid,
            Centroid, bool>
    {
        /** Comparison operator.
         * @param[in] lhs Left-hand-side centroid object for the comparison.
         * @param[in] rhs Right-hand-side centroid object for the comparison.
         * @return True if the left hand side object has a scan number that is
         *         smaller than the scan number of the right hand side centroid
         *         object.
         */
        bool operator()(const Centroid& lhs, const Centroid& rhs);
    };

    /** Functor to compare two centroids based on their retention time.
     */
    struct LessThanAbundance : public std::binary_function<Centroid, Centroid,
            bool>
    {
        /** Comparison operator.
         * @param[in] lhs Left-hand-side centroid object for the comparison.
         * @param[in] rhs Right-hand-side centroid object for the comparison.
         * @return True if the abundance of the left hand centroid object is
         *         smaller than the abundance of the right hand side centroid
         *         object.
         */
        bool operator()(const Centroid& lhs, const Centroid& rhs);
    };

    /** Accessor functor to access the m/z value of a centroid.
     */
    struct MzAccessor
    {
        typedef double value_type;
        double& operator()(Centroid& c);
        double operator()(const Centroid& c) const;
    };

    /** Accessor functor to access the rt value of a centroid.
     */
    struct RtAccessor
    {
        typedef double value_type;
        double& operator()(Centroid& c);
        double operator()(const Centroid& c) const;
    };

    /** Accessor functor to access the abundance value of a centroid.
     */
    struct AbundanceAccessor
    {
        typedef double value_type;
        double& operator()(Centroid& c);
        double operator()(const Centroid& c) const;
    };

    /** Constructor.
     */
    Centroid();

    /** Copy constructor.
     */
    Centroid(const Centroid& rhs);

    /** Alternate constructor.
     */
    Centroid(Double retentionTime, Double mz, UnsignedInt sn, Double ab,
        Spectrum::const_iterator first,
        Spectrum::const_iterator last);

    /** Comparison operator.
     * @param[in] rhs The right hand side comparison object.
     */
    bool operator==(const Centroid& rhs) const;

    /** Desctructor.
     */
    ~Centroid();

    /** Get the retention time associated with the centroid.
     * @return The retentiion time value, in seconds.
     */
    Double getRetentionTime() const;

    /** Set the retention time for the centroid.
     * @param[in] rt The retention time, in seconds.
     * @throw mstk::PreconditionViolation if \c rt < 0.0.
     */
    void setRetentionTime(Double rt);

    /** Get the mass/charge value of the centroid.
     * @return The m/z value.
     */
    Double getMz() const;

    /** Set the mass/charge value of the centroid.
     * @param[in] mz The mass/charge value. Must be non-negative.
     */
    void setMz(const Double);

    /** Get the scan number associated with the centroid.
     * @return The scan number.
     * @throw mstk::PreconditionViolation if \c mz < 0.0.
     */
    UnsignedInt getScanNumber() const;

    /** Set the scan number of the centroid.
     * @param[in] sn The scan number.
     */
    void setScanNumber(const UnsignedInt sn);

    /** Get the centroid abundance.
     * @return The abundance of the centroid.
     */
    Double getAbundance() const;

    /** Set the centroid abundance.
     * @param[in] ab The centroid abundance. Must be non-negative.
     */
    void setAbundance(const Double ab);

    /** Get the raw data underlying the centroid.
     * @return A \c Spectrum object holding the raw data.
     * @throw mstk::PreconditionViolation if \c ab < 0.0.
     */
    const Spectrum& getRawData() const;

    /** Set the raw data underlying the centroid.
     * @param[in] s An \c mstk::Spectrum object.
     */
    void setRawData(const Spectrum& s);

private:
    /** The retention time.
     */
    Double rt_;

    /** The m/z position.
     */

    Double mz_;
    /** The scan number.
     */

    UnsignedInt sn_;
    /** The centroid abundance.
     */
    Double ab_;

    /** The raw data. May be empty.
     */
    Spectrum raw_;
};

/** Streams a text representation of the centroid object.
 * @param[inout] os The stream into which to write.
 * @param[in] c The centroid object that should be streamed.
 */
std::ostream& operator<<(std::ostream& os, const Centroid& c);

} // namespace fe   

} // namespace mstk

#endif

