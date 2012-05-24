/*
 * SpectrumAlgorithm.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_PSF_SPECTRUMALGORITHM_HPP__
#define __MSTK_INCLUDE_MSTK_PSF_SPECTRUMALGORITHM_HPP__
#include <MSTK/config.hpp>

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include <MSTK/common/Log.hpp>
#include <MSTK/common/Error.hpp>

namespace mstk {

namespace psf {

/**
 * Compare two elements with regard to a certain aspect.
 *
 * Typically, an element in a spectrum represents more than one value 
 * (such as m/z, intensity, time etc.) Use this functor to compare two
 * elements with regard to one of these values.
 */
template< typename Element, typename Extractor >
class MSTK_EXPORT LessByExtractor {
  public:
  LessByExtractor( const Extractor& e ) : extract_(e) {};
  bool operator()( const Element& lhs, const Element& rhs ) const {
    return extract_(lhs) < extract_(rhs);    
  };

  private:
  Extractor extract_;
};

template< typename Element, typename Extractor >
class MSTK_EXPORT MoreThanValue {
   public:
   MoreThanValue( const Extractor& e, typename Extractor::result_type val ) : extract_(e), val_(val) {};
   bool operator()( const Element& e ) const {
     return val_ < extract_(e);    
  };

  private:
  Extractor extract_;
  typename Extractor::result_type val_;
};



/**
 * Finds the first 'bump' in a sequence.
 *
 * A bump is a range in a sequence containing a (local) maximum and stricly decreasing values
 * to the left and to the right of the maximum.
 *
 * The smallest possible bump consists of only three elements: .'.
 * 
 * @param first Iterator pointing to the first element of the input sequence. Has to point to
 *              an existing element. Else the behaviour is undefined.
 * @param last Iterator pointing to one past the last element of the input sequence. Else the
 *             behaviour is undefined.
 * @param comp comp(*iter, *(iter+1)) is called to compare two elements. Think of it as a 'less than' operator <.
 * 
 * @return The first bump found in the input sequence: [pair.first, pair.second]. Returns a pair of 'last', if no
 *         bump is found.
 */
template<typename FwdIter, typename Compare>
MSTK_EXPORT std::pair<FwdIter, FwdIter> findBump(FwdIter first, FwdIter last, Compare comp);

 /**
 * Sample the full width at a fraction of the maximum.
 *
 * Goes through a (sub-)spectrum and measures the full width at a fraction of the maximum for
 * every spectral peak considered pure.
 *
 * A pure peak fulfills the requirements of a 'bump' and is at least as low as the fraction
 * of its maximum. The true peak maximum is estimated as the most abundant element of the 
 * bump.
 *
 * Note: You should see this really as a measurement in the physical sense, meaning even
 * in the case of an exactly calculatable width, this function may return a slightly 
 * different value due to rounding errors and similar effects.
 *
 * The distance (last - first) may not be negative, else the behaviour is undefined.
 * If (last - first) is zero, an empty vector is returned. 
 *
 * @param first The first Element of the sequence to sample from.
 * @param last One past the last Element of the sequence to sample from.
 * @param fraction Has to be between 0.0 and 1.0, borders included.
 * @param minimalPeakHeight Even negative values are allowed, albeit in general not meaningful.
 * @return Pairs of (mz | width at mz) in ascending order of mz. If no pairs found, the
 *         vector is empty.
 *
 * @throw psf::PreconditionViolation Parameter fraction is out of the required range.
 *
 * @see psf::findBump
 * @see psf::SpectralPeak::lowness
 */
template<typename FwdIter, typename MzExtractor, typename IntensityExtractor> 
MSTK_EXPORT
std::vector<std::pair<typename MzExtractor::result_type, typename MzExtractor::result_type> > 
measureFullWidths(const MzExtractor&, const IntensityExtractor&, FwdIter first, FwdIter last, double fraction, typename IntensityExtractor::result_type minimalPeakHeight = 0);

/**
 * Aggregates algorithms working on a spectral peak.
 *
 * A spectral peak is a single peak in a mass spectrum; in contrary to a monoisotopic peak,
 * which represents a whole isotope pattern.
 *
 * A spectral peak is represented as a sequence of spectrum elements. The
 * algorithms work with iterators on this sequence. The elements have to be in ascending
 * order of their m/z values.
 * There are no further requirements to constitute a sequence a spectral peak. For 
 * example, even a set of equiabundant elements can be seen as a spectral peak.
 *
 * @see psf::Peak
 */
namespace SpectralPeak
{
    /**
     * The height of a spectral peak.
     *
     * The heighest abundance in a sequence of spectral elements is detected.
     * This abundance is interpreted as the peak height.
     *
     * The minimal requirements of the sequence are:
     * @li At least one element.
     *
     * @param firstElement Points to the first element of the sequences of spectral elements
     *                     to be considered a peak.
     * @param lastElement Points to the last element of the sequences of spectral elements
     *                    to be considered a peak.
     * @return The peak height, in units of the element abundance.
     *
     * @throw psf::PreconditionViolation Minimal sequence requirements aren't met.
     */
    template< typename FwdIter, typename IntensityExtractor >
    MSTK_EXPORT
    typename IntensityExtractor::result_type
    height(const IntensityExtractor&, FwdIter firstElement, FwdIter lastElement);     

    /**
      * The lowness of a spectral peak.
      *
      * The highest element in a sequence of spectral elements is detected and the two
      * elements with the lowest abundance are searched; one on the left and one on the
      * right of the maximum (the maximum itself can be selected). The more abundant peak of these two is chosen.
      * One minus the ratio of this abundance to the maximum is called  the 'peak lowness'.
      *
      * An equiabundant sequence of spectral elements has a lowness of 0.0. In contrast,
      * a maximum flanked by two elements with almost zero abundance has a lowness of
      * almost 1.00.
      *
      * The minimal requirements of the sequence are:
      * @li At least one element. The lowness is then 0.0.
      *       
      * @param firstElement Points to the first element of the sequences of spectral elements
      *                     to be considered a peak.
      * @param lastElement Points to the last element of the sequences of spectral elements
      *                    to be considered a peak.
      * @return The lowness is between 0.0 and 1.0, borders included.
      */
    template< typename FwdIter, typename IntensityExtractor >
    MSTK_EXPORT 
    double lowness(const IntensityExtractor&, FwdIter firstElement, FwdIter lastElement);

    /**
      * The full width at a fraction of the maximum of a spectral peak.
      *
      * The most abundant element in a sequence of spectral elements is found. Coming from 
      * the left and right, the two elements which are nearest to the fraction of the
      * abundance maximum are search: abundance fraction <= element abundance. The two
      * elements on both sides are then linearly interpolated with their neighboring elements
      * just below the fraction of the maximum. (In case of fraction == element abundance, this is avoided.)
      * The distance in m/z dimension between the two interpolated elements
      * flanking the maximum is then returned as the full width at the fraction of the maximum.
      *
      * For example, setting the fraction to 0.5 just returns the full width at half maximum.
      *
      * The minimal requirements of the sequence are:
      * @li At least one abundance maximum. If there are multiple maxima, the first such element
      *     is chosen.
      * @li At least one element below the fraction of the maximum on both flanks of the maximum.
      *     The maximum itself may be chosen as the element above the fraction.
      *
      * @param firstElement Points to the first element of the sequences of spectral elements
      *                     to be considered a peak.
      * @param lastElement Points to the last elmement of the sequences of spectral elements
      *                    to be considered a peak.
      * @param fraction Has to be between 0.0 and 1.0, borders included.
      * @return The full width at the fraction of the maximum, in units of m/z.
      *
      * @throw psf::PreconditionViolation Parameter fraction not in the required range.
      * @throw psf::Starvation The sequence doesn't satisfy the minimal requirements.
      */
    template< typename FwdIter, typename MzExtractor, typename IntensityExtractor >
    MSTK_EXPORT
    typename MzExtractor::result_type 
    fullWidthAtFractionOfMaximum(const MzExtractor&, const IntensityExtractor&, FwdIter firstElement, FwdIter lastElement, const double fraction);

} /* namespace SpectralPeak */

/******************/
/* implementation */
/******************/

template<typename FwdIter, typename Compare>
std::pair<FwdIter, FwdIter> findBump(FwdIter first, FwdIter last, Compare comp) {
    // First, let's assume the bump starts right on the first element.
    FwdIter leftEdge = first;
    // Furthermore, we are starting out on the first element, too.
    FwdIter currentElement = first;

    // Since we are working with ForwardIterators, we cannot rewind and have to store the
    // next position separately
    const FwdIter onePastLastElement = last;
    FwdIter onePastCurrentElement = currentElement; std::advance(onePastCurrentElement, 1);   

    // state of our search
    bool onIncreasingSlope = false;
    bool foundBumpTop = false;

    // Now let's go through the sequence element by element.
    while(onePastCurrentElement != onePastLastElement) {
        mstk_invariant(std::distance(currentElement, onePastCurrentElement) == 1, "findBump(): Distance between current and next elment is not 1.");
        
        /* Current element is smaller than its right neighbor. */      
        if(comp(*currentElement, *(onePastCurrentElement))) {
            // We already passed the top of a bump and are finished!
            if(foundBumpTop) {
                // currentElement points to the right edge of the bump
                break;
            }
            else {
                // We are on the bottom of an increasing slope.
                if(!onIncreasingSlope) {
                    onIncreasingSlope = true;
                    // That's were bumps are starting!
                    leftEdge = currentElement;
                }
            }
        }
        
        /* Current element is bigger than its right neighbor. */
        else if(comp(*(onePastCurrentElement), *currentElement)) { 
            // Were we on an increasing slope?
            if(onIncreasingSlope) {
                // So we are on the top of the bump!
                foundBumpTop = true;
            }
            // We are on a decreasing slope.
            else {
                onIncreasingSlope = false;
            }
        }

        /* Current element and its neighbor are equal. */
        else {
            // Great, we are finished!
            if(foundBumpTop) {
                // currentElement points to the right edge of our bump.
                break;
            }
            
            // That's bad. We have to start our search again.
            else {            
                leftEdge = onePastCurrentElement;
                onIncreasingSlope = false;
            }
        }

        ++currentElement;
        ++onePastCurrentElement;
    }
    
    if(foundBumpTop) { 
        // we return [leftEdge, rightEdge]
        return std::make_pair(leftEdge, currentElement);
    }
    else {
        return std::make_pair(last, last);
    }
}

template<typename FwdIter, typename MzExtractor, typename IntensityExtractor> 
std::vector<std::pair<typename MzExtractor::result_type, typename MzExtractor::result_type> > 
measureFullWidths(const MzExtractor& get_mz, const IntensityExtractor& get_int, FwdIter first, FwdIter last, double fraction, typename IntensityExtractor::result_type minimalPeakHeight = 0) {
    typedef typename MzExtractor::result_type Mz;
    typedef typename IntensityExtractor::result_type Intensity;

    mstk_precondition(0. <= fraction && fraction <= 1., 
        "measureFullWidths(): Parameter fraction out of required range.");

    // The result
    std::vector<std::pair<Mz, Mz> > widths;

    // Check for empty spectrum or only one element
    if((last - first) < 1) {
        return widths;
    }

    // prerequisites    
    const double requiredLowness = 1. - fraction;
    std::pair<FwdIter, FwdIter> bump;
    Mz positionOfMaximum = 0;
    Intensity bumpHeight = 0;    
    Mz width = 0;
    LessByExtractor< typename IntensityExtractor::element_type, IntensityExtractor > comp(get_int);
    
    // go through all bumps in the spectrum */       
    while(first < last) {
        bump = findBump(first, last, comp);
        // No new bump found
        if (bump.first == last) {
            break;
        }
        
        // calc full width if bump is low enough and has a minimal height
        mstk_invariant(bump.first <= bump.second && bump.second < last, "Bump in illegal state.");
        bumpHeight = get_int(*std::max_element(bump.first, bump.second + 1, comp));       
        if(SpectralPeak::lowness(get_int, bump.first, bump.second) >= requiredLowness && bumpHeight >= minimalPeakHeight) {
            // we don't have to try for exceptions here, because a bump fulfills the preconditions
	  width = SpectralPeak::fullWidthAtFractionOfMaximum(get_mz, get_int, bump.first, bump.second, fraction);
            positionOfMaximum = get_mz(*std::max_element(bump.first, bump.second + 1, comp));
            MSTK_LOG(logDEBUG) << "measureFullWidths(): Measured peak (mz | width): (" << positionOfMaximum << " | " << width << ")";
            widths.push_back(std::make_pair(positionOfMaximum, width));
        }
        
        // last element of the bump may be the first of the next one
        first = bump.second;
    }
    
    return widths;
}

template< typename FwdIter, typename IntensityExtractor >
typename IntensityExtractor::result_type
SpectralPeak::height(const IntensityExtractor& get_int, FwdIter firstElement, FwdIter lastElement) {
    mstk_precondition(distance(firstElement, lastElement) >=0, "SpectralPeak::height(): Distance between first and last input element is not nonnegative.");
    
    // Compare elements by intensity
    LessByExtractor<typename IntensityExtractor::element_type, IntensityExtractor> comp(get_int);
    // find maximum intensity
    FwdIter maximum = max_element(firstElement, ++lastElement, comp);

    return get_int(*maximum);
}

template< typename FwdIter, typename IntensityExtractor >
double psf::SpectralPeak::lowness(const IntensityExtractor& get_int, FwdIter firstElement, FwdIter lastElement) {
    // comply to STL standards.
    FwdIter last = lastElement;
    ++last;    

    // Compare elements by intensity
    LessByExtractor<typename IntensityExtractor::element_type, IntensityExtractor> comp(get_int);

    // find maximum intensity
    FwdIter maximum = max_element(firstElement, last, comp);

    // find least abundant element right of the maximum
    FwdIter rightMinimum = min_element(maximum, last, comp);
    // and to the left (both times with the maximum included as possible minimum)
    FwdIter leftMinimum = min_element(firstElement, ++maximum, comp);
    --maximum; // STL required [first, last)

    // more abundant element of the two
    typename IntensityExtractor::element_type moreAbundantOne = std::max(*leftMinimum, *rightMinimum, comp);

    return 1. - (get_int(moreAbundantOne)/get_int(*maximum));
}

namespace {
/**
 * Find Element below the target abundance given the Element above or on the target abundance.
 *
 * We are working on the sequence [firstElement, above]. Every Element besides the 'above'
 * element has to be lower than the 'target' abundance. The above Element may actually
 * be exactly on or above the target abundance. If these preconditions are not met,
 * the behaviour is undefined.
 *
 * @param firstElement Input const_iterator to the first spectrum element of the sequence to be searched.
 * @param above Input const_iterator to the spectrum element on or above the target abundance. The distance 
 *              above - firstElement may not be negative. Else, the behaviour is undefined.
 * @param target The target abundance
 *
 * @return The element below the target abundance nearest (in mz dimension) to the above Element.
 *         Is the above Element exactly on the target, below is the same as the above Element.
 *
 * @throw psf::Starvation No Element below found.
 */      
template <typename const_InIter, typename IntensityExtractor>
const_InIter findElementBelowTargetAbundance(const IntensityExtractor&, const const_InIter firstElement, const const_InIter above, const typename IntensityExtractor::result_type target);

/**
 * Takes two elements and blends them together with a specific target intensity.
 *
 * The interpolation is done linearly in the mz dimension so that the result has the 
 * target abundance. The order of element1 and element2 is not important. If element1 and element2 have the same mz value,
 * then the mz value is directly return without interpolation.
 * If element1 and element2 differ in mz, they have to differ in abundance, too.
 *
 * @return The mz value of the interpolated element with an abundance of 'target'.
 *
 * @throw psf::InvariantViolation Element1 and element2 differ in mz but not in abundance.
 *                               No interpolation to a general target is then possible.
 */
template< typename MzExtractor, typename IntensityExtractor, typename element_type>
typename MzExtractor::result_type interpolateElements(const MzExtractor&, const IntensityExtractor&, const element_type& element1, const element_type& element2, const typename IntensityExtractor::result_type target);

} /* anonymous namespace */

template< typename FwdIter, typename MzExtractor, typename IntensityExtractor >
typename MzExtractor::result_type 
SpectralPeak::
fullWidthAtFractionOfMaximum(const MzExtractor& get_mz, const IntensityExtractor& get_int, FwdIter firstElement, FwdIter lastElement, const double fraction) {
    mstk_precondition(0. <= fraction && fraction <= 1., "fullWidthAtFractionOfMaximum(): Fraction parameter out of range.");

    /* Determine target intensity */
    // Compare elements by intensity
    LessByExtractor<typename IntensityExtractor::element_type, IntensityExtractor> comp(get_int);
    // find maximum intensity
    FwdIter maximum = max_element(firstElement, lastElement + 1, comp);
    MSTK_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): Spectral peak maximum detected at (mz, intensity): " << get_mz(*maximum) << " ," << get_int(*maximum); 
    // calc target intensity
    const typename IntensityExtractor::result_type target = get_int(*maximum) * fraction;
    MSTK_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): Fraction of maximal intensity is: " << target;
    
    // we need that further below for finding elements
    MoreThanValue<typename IntensityExtractor::element_type, IntensityExtractor> compScalar(get_int, target);

    /* find utter left element nearest above or on target */
    // target <= above == !(above < target) 
    FwdIter aboveOnLeft = find_if(firstElement, maximum + 1, compScalar);
    MSTK_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): aboveOnLeft detected at (mz, intensity): " << get_mz(*aboveOnLeft) << " ," << get_int(*aboveOnLeft); 
    // determine belowOnLeft
    FwdIter belowOnLeft = findElementBelowTargetAbundance(get_int, firstElement, aboveOnLeft, target);
    MSTK_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): belowOnLeft detected at (mz, intensity): " << get_mz(*belowOnLeft) << " ," << get_int(*belowOnLeft);

    /* find utter right element */
    // we now start searching from the right
    std::reverse_iterator<FwdIter> rlast(lastElement + 1);
    std::reverse_iterator<FwdIter> rmaximum(maximum);
    // target <= above == !(above < target) 
    std::reverse_iterator<FwdIter> aboveOnRight = find_if(rlast, rmaximum, compScalar);
    MSTK_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): aboveOnRight detected at (mz, intensity): " << get_mz(*aboveOnRight) << " ," << get_int(*aboveOnRight);     
    // determine belowOnRight
    std::reverse_iterator<FwdIter> belowOnRight = findElementBelowTargetAbundance(get_int, rlast, aboveOnRight, target);
    MSTK_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): belowOnRight detected at (mz, intensity): " << get_mz(*belowOnRight) << " ," << get_int(*belowOnRight); 
    
    /* interpolate below and above elements */
    typename MzExtractor::result_type leftInterpolated = interpolateElements(get_mz, get_int, *belowOnLeft, *aboveOnLeft, target);
    MSTK_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): leftInterpolated is: " << leftInterpolated;
    typename MzExtractor::result_type rightInterpolated = interpolateElements(get_mz, get_int, *belowOnRight, *aboveOnRight, target);
    MSTK_LOG(logDEBUG1) << "fullWidthAtFractionOfMaximum(): rightInterpolated is: " << rightInterpolated;

    return rightInterpolated - leftInterpolated;
}

namespace {
    template <typename const_InIter, typename IntensityExtractor>
    const_InIter findElementBelowTargetAbundance(const IntensityExtractor& get_int, const const_InIter firstElement, const const_InIter above, const typename IntensityExtractor::result_type target) {
        const_InIter below;

        // Check if the above Element coincides with the firstElement
        if (firstElement == above) {
            // Rule out special case, where abundance of above Element equals target abundance
            if(target < get_int(*above)) { 
                throw Starvation("fullWidthAtFractionOfMaximum(): No elements on the left below target abundance.");
            }
            // special case: target == above
            else {
                below = above;
                MSTK_LOG(logDEBUG2) << "findElementBelowTargetAbundance(): Target abundance equals abundance of element above. Setting below element equal to above element.";
            }        
        } else {
            // Choose nearest neighbor in mz dimension (guaranteed to be below target abundance by our imposed preconditions).
            below = above - 1;
        }

        return below;
    }

    template< typename MzExtractor, typename IntensityExtractor, typename element_type>
    typename MzExtractor::result_type interpolateElements(const MzExtractor& get_mz, const IntensityExtractor& get_int, const element_type& element1, const element_type& element2, const typename IntensityExtractor::result_type target) {
        if (get_mz(element1) == get_mz(element2)) {
            return get_mz(element2);
        } else {
            mstk_invariant(get_int(element1) != get_int(element2), "interpolateElements(): Illegal abundance state: below < target && target <= above && above == below."); 

            // abundance = slope * mz + shift      
            double slope = 1.0 * (get_int(element2) - get_int(element1)) / (get_mz(element2) - get_mz(element1));
            MSTK_LOG(logDEBUG2) << "interpolateElements(): slope of linear interpolation: " << slope;              

            // just take one of the two Elements to determine the shift
            double shift = get_int(element1) - slope * get_mz(element1);
            MSTK_LOG(logDEBUG2) << "interpolateElements(): shift of linear interpolation: " << shift; 

            // => mz = (abundance - shift)/slope
            return static_cast<typename MzExtractor::result_type>( (target - shift) / slope );
        }
    }
} /* anonymous namespace */

} /* namespace psf */

} /* namespace mstk */

#endif /*__SPECTRUMALGORITHM_HPP__*/
