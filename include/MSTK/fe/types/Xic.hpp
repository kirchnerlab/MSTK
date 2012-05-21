/*
 * Xic.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_TYPES_XIC_HPP__
#define __MSTK_INCLUDE_MSTK_FE_TYPES_XIC_HPP__

#include <MSTK/config.hpp>
#include <MSTK/common/Collection.hpp>
#include <MSTK/fe/types/Centroid.hpp>

#include <functional>

namespace mstk {

namespace fe {

/** EXtracted Ion Current (XIC) class.
 * @note Clients using this class need to be aware that changing
 *       the container contents necessarily invalidates the XIC
 *       statistics like m/z and RT means, variances, as well as 
 *       overall abundance. Consequently, clients that change the
 *       internal state of the vector need to make use of 
 *       \c recalculate() to force the statistics update
 *       before any \c getXYZ() request.
 */
class MSTK_EXPORT Xic : public Collection<Centroid>
{
public:
    /** Comparison functor for abundance-based comparison of XIC objects.
     */
    struct LessThanAbundance : public std::binary_function<Xic, Xic, bool>
    {
        bool operator()(const Xic& lhs, const Xic& rhs) const;
    };

    /** Comparison functor for retention time-based comparison of XIC objects.
     */
    struct LessThanRt : public std::binary_function<Xic, Xic, bool>
    {
        bool operator()(const Xic& lhs, const Xic& rhs) const;
    };

    /** Comparison functor for m/z-based comparison of XIC objects.
     */
    struct LessThanMz : public std::binary_function<Xic, Xic, bool>
    {
        bool operator()(const Xic& lhs, const Xic& rhs) const;
    };

    /** Accessor functor to access the m/z value of a centroid.
     */
    struct MzAccessor
    {
        typedef double value_type;
        double& operator()(Xic& c);
        double operator()(const Xic& c) const;
    };

    /** Accessor functor to access the rt value of a centroid.
     */
    struct RtAccessor
    {
        typedef double value_type;
        double& operator()(Xic& c);
        double operator()(const Xic& c) const;
    };

    /** Accessor functor to access the abundance value of a centroid.
     */
    struct AbundanceAccessor
    {
        typedef double value_type;
        double& operator()(Xic& c);
        double operator()(const Xic& c) const;
    };

    // other typedefs
    typedef Collection<Centroid>::const_iterator const_iterator;

    /** Constructor.
     */
    Xic();

    /** Construct XIC from XIC subset.
     * @param[in] first The first element of the XIC that should be used to
     *                  construct the sub-XIC.
     * @param[in] last The one-beyond-the-end element of the XIC that
     *                  should be used to construct the sub-XIC.
     */
    Xic(const_iterator first, const_iterator last);

    /** Comparison operator for equality.
     * @param[in] rhs The XIC object to compare with.
     * @return True if the current object and \c rhs are the same.
     */
    bool operator==(const Xic& rhs);

    /** Returns the overall abundance of an XIC. The abundance is calculated by
     *  triangulation along the retention time domain.
     *  @return The overall abundance of the XIC.
     */
    double getAbundance() const;

    /** Abundance-weighted mean mass/charge ratio over all high-accuracy
     *  centroid masses in the XIC.
     *  @return The mean mass/charge ratio of the XIC.
     */
    double getMz() const;

    /** Get the mass/charge standard deviation of the XIC.
     * @return The m/z standard deviation of the XIC.
     */
    double getMzTolerance() const;

    /** Abundance-weighted mean retention time of all centroids in the XIC.
     * @return The retention time of the XIC (in seconds).
     */
    double getRetentionTime() const;

    /** Get the retention time standard deviation of the XIC.
     * @return The retention time standard deviation, in seconds.
     */
    double getRetentionTimeTolerance() const;

    /** Split the current XIC into pieces. The methods follows the idea
     * of [Cox and Mann, 2008]: if two consecutive local maxima are separated
     * by a local minimum that is at least \c mindepth smaller than the smaller
     * of the two local maxima, then the XIC should be split.
     *
     * @param[out] xics The XICs that result from the splitting.
     * @param[in] mindepth The factor by which the local minimum must be smaller
     *                     than the lower of the two surrounding local maxima
     *                     to warrant a split.
     */
    void split(std::vector<Xic>& xics, double mindepth = 0.76);

    /** Recalculate all XIC sufficient statistics.
     */
    void recalculate();

    /** Calculate uncentered Pearson correlation w/ another Xic.
     * @param[in] rhs \c Xic object to correlated with.
     * @return The uncentered Pearson correlation.
     */
    double correlate(Xic& rhs);


private:
    struct RtAbAccessor;

    /** Merge centroids inside the XIC that have the same scan number.
     */
    void mergeDuplicates();

    /** Return a 3-point average smoothed version of the XIC.
     * @param[out] xic The smooted XIC. Any previous values are overwritten.
     */
    void getSmoothedXic(Xic& xic);

    double rt_, rtSigma_;
    double mz_, mzSigma_, abundance_;
};

// stream operator for Xics (not const due to cached getXXXX calls).
std::ostream& operator<<(std::ostream&, const Xic&);

} // namespace fe

} // namespace mstk

#endif

