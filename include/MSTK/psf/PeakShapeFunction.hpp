/*
 * PeakShapeFunction.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_PSF_PEAKSHAPEFUNCTION_HPP__
#define __MSTK_INCLUDE_MSTK_PSF_PEAKSHAPEFUNCTION_HPP__

#include <MSTK/config.hpp>
#include <cmath>
#include <string>
#include <MSTK/common/Error.hpp>
#include <MSTK/common/Log.hpp>
#include <MSTK/psf/PeakParameter.hpp>
#include <MSTK/psf/PeakShape.hpp>

/**
 * @page peakshapefunction Peakshape Functions
 *
 *
 *
 * @section introductionofpeakshapefunction What is a Peakshape Function?
 *
 * A mass spectrum usually consists of intensities drawn over the mass-over-charge values,
 * called 'mass channels'.
 * In theory, ions with a distinct mass-over-charge ratio should only contribute to one
 * mass channel producing a sharp peak in the mass spectrum. Due to imprecisions in the
 * measurement caused by random physical processes the actual peak is smeared out over a
 * few neighboring mass channels, resembling roughly a bell-shaped-curve in the intensity
 * domain.
 * Thus, a function describing the envelope of a single peak in a mass spectrum is called
 * a peakshape function.
 *
 *
 * @section typesofpeakshapefunctions Types of Peakshape Functions in psf
 *
 * @subsection boxpsf Box
 * The box peakshape function is the rectangular function with a fixed width. It is implemented and used in the
 * transmogrifier unit test.
 *
 * This function is too crude an approximation and shouldn't be used in serious data analysis.
 *
 *
 * @subsection gaussianpsf Gaussian
 * The gaussian peakshape function is @f$ e^{-\frac{(\Delta\frac{m}{z})^{2}}{2\cdot\sigma^{2}}} @f$.
 * The only parameter is the full width at half maximum, which is connected to @f$ \sigma @f$ via
 * @f$ \mathrm{FWHM}=2\sqrt{2\ln2}\cdot\sigma @f$.
 *
 * This peakshape function is independent of the absolute value of the mass channel. In a real mass spectrum
 * the width of the peakshape is dependent on the mass channel (the resolution is usually decreasing
 * with higher masses). So, this function should only be used in a relatively small mass interval, where one
 * can safely neglect this change in the peak width.
 *
 * @see psf::GaussianPeakShapeFunction
 *
 *
 * @subsection orbipsf Orbitrap
 * The Orbitrap peak shape is gaussian. The full width at half maximum of the gaussian depends on the
 * absolute value of the mass channel and is calculated as follows: @f$ \mathrm{FWHM} = a\cdot\mathrm{mass}\sqrt{\mathrm{mass}} @f$.
 *
 * @see psf::OrbitrapPeakShapeFunction
 *
 *
 * @subsection orbiboxpsf Orbitrap Box
 * The Orbitrap box peakshape function is the rectangular function. The width of the rectangle depends on the
 * absolute value of the mass channel: @f$ \mathrm{width} = \mathrm{mass} \cdot \mathrm{sigma\_ppm}@f$.
 *
 * @f$ \mathrm{sigma\_ppm} @f$ is the only parameter and can be choosen as the resolution of an Orbitrap in parts-per-million.
 * Since for a real Orbitrap the peak width depends on the mass channel in a non-linear way and the peak shape is
 * far away from a box, this is a very crude approximation. Nevertheless, it can be used to speed up the calculation on good data sets
 * (high signal-to-noise ratio, few overlapping isotope patterns etc.) and can perform passable on such data.
 *
 * @see psf::OrbiBoxPeakShapeFunction
 *
 *
 * @subsection tofpsf Time-of-Flight
 * The time-of-flight peakshape function is gaussian. The full width at half maximum depends on the absolute value
 * of the mass channel as follows: @f$ \mathrm{FWHM} = a \cdot \sqrt{\mathrm{mass}} + b @f$.
 * @f$ a @f$ and @f$ b @f$ are parameters. They can be determined with the help of a peakshape function
 * estimator (by fitting the function on a given spectrum).
 * This formula can be deduced from the physics of a time-of-flight mass spectrometer.
 * 
 * @see psf::TOFPeakShapeFunction and psf::TOFPeakShapeFunctionEstimator
 */

namespace mstk {

namespace psf {
/**
 * The peak shape functions, which are implementing the abstract PeakShapeFunction interface.
 *
 * The 'box' type is only used for unit testing. The other types are included in the core
 * library.
 */
enum MSTK_EXPORT PeakShapeFunctionTypes {box, gaussian, orbi, orbiBox, tof};

/**
 * Encapsulates the PeakShapeFunctionTypes and provides conversion functions.
 */
class MSTK_EXPORT PeakShapeFunctionType
{
    public:
        /**
         * You may use the constructor in an implict way: 'PeakShapeFunctionType type = gaussian;'.
         */
        PeakShapeFunctionType(PeakShapeFunctionTypes type) : type_(type) {};

        PeakShapeFunctionTypes toEnum();

        /**
         * @return Depending on the type: 'box', 'gaussian', 'orbi', 'orbiBox' or 'time-of-flight'.
         *         If the type is not known, 'unknown'.
         */
        std::string toString();

    protected:
        PeakShapeFunctionTypes type_;
};



/**
 * Implements a generic psf::PeakShapeFunction.
 *
 * This generic implementation of the psf::PeakShapeFunction interface hast to be parameterized
 * by three template parameters.
 *
 * Use 'typedef PeakShapeFunctionTemplate<...> MyPeakShapeFunction' to implement new peak shape
 * functions quickly. Several peak shape functions defined in this way, can be found in the headerfile 'PeakShapeFunctions.h'.
 *
 * At the moment, this template supports up to two general paramters a,b to parameterize the
 * peak shape function. The actual number depends on the PeakParameterT used. If the PeakParameterT doesn't
 * support a or b, you will get a compile time error trying to use them.
 *
 * Furthermore, some PeakParameterT support autocalibration and/or setting a desired value at a specific
 * coordinate in the spectrum (for example the FWHM at 400 Th). In case of no support, your code will
 * not compile after calling such a function. (This is a desirable behaviour.)
 *
 * @param PeakShapeT Describes the general form of a peak. Can be found in the headerfile 'PeakShape.h'
 * @param PeakParameterT The PeakShapeT has to be parameterized according to peak parameters depending
 *                       on the type of mass spectrometer to describe. Can be found in the
 *                       headerfile 'PeakParameter.h'
 * @param PeakShapeFunctionTypeT The proper name of the peak shape function. Can be found
 *                               in the headerfile 'PeakShapeFunction.h'.
 *
 * @see psf::OrbitrapPeakShapeFunction
 * @see psf::GaussianPeakShapeFunction
 */
template <typename PeakShapeT, typename PeakParameterT, PeakShapeFunctionTypes PeakShapeFunctionTypeT>
class MSTK_EXPORT PeakShapeFunctionTemplate
{
public:
    PeakShapeFunctionTemplate();
    explicit PeakShapeFunctionTemplate(const double a);
    PeakShapeFunctionTemplate(const double a, const double b);

    /**
     * operator()
     * @param referenceMass the m/z value at the center of the PSF
     * @param observedMass the m/z value of the mass for which the value of the PSF is desired
     * @return the value of the PSF at (observedMass-referenceMass)
     */
    double operator()(const double referenceMass, const double observedMass) const;

    /**
     * Return the width of the PSF support at a specific m/z value.
     *
     * The threshold is a relative distance measured from the center of the peak shape
     * and is symmetrical around the center. After the threshold, the peak shape function
     * is set to zero.
     * @param mz the m/z value of the PSF center
     * @return the width of the PSF support at the given m/z position
     */
    double getSupportThreshold(const double mz) const;

    /**
     * Returns the actual implementation type of the abstract PeakShapeFunction interface.
     */
    PeakShapeFunctionType getType();

    void setA(const double a);
    double getA() const;

    void setB(const double b);
    double getB() const;

    /**
     * Autocalibrate peak shape function parameters using regression.
     *
     * Use this function to learn from a sequence of spectrum elements.
     *
     * Note, that there is no internal error threshold or similar for the quality of the
     * calibration. The calibration is performed as long as it is possible in any way.
     * To achieve a good result, one should filter out the noise of the input spectrum and/
     * or use high quality data in the first place.
     *
     * The distance (last - first) may not be negative, else the behaviour is undefined.
     * It doesn't violate the preconditions, if (last - first) is zero or small. Nevertheless,
     * it increases the chance of a Starvation exception to happen.
     *
     * The elements in the spectrum have to be in ascending order of their mz value. Furthermore,
     * no elements with duplicate mz values are allowed. Else, the behaviour is undefined.
     *
     * @param first Points to the first Element of the spectrum.
     * @param last Points to one past the last Element of the spectrum.
     *
     * @throw psf::Starvation To few or bad data extracted from the input spectrum to make a
     *                       calibration possible.
     */
    template< typename FwdIter, typename MzExtractor, typename IntensityExtractor >
    void calibrateFor(const MzExtractor&, const IntensityExtractor&, FwdIter first, FwdIter last);
    
    // setMinimalPeakHeightForCalibration()
    /**
     * Use only peaks with a minimal intensity for autocalibration.
     *
     * The internally used peak parameters have to support this feature. Else, calling the
     * function would result in a compile time error.
     *
     * @param minimalHeight The minimal absolute intensity. Negative values are possible,
     *                      albeit not meaningful.
     */
    void setMinimalPeakHeightForCalibration(double minimalHeight);

    // getMinimalPeakHeightForCalibration()
    /**
     * Uses only peaks with a minimal intensity for autocalibration.
     *
     * @see PeakShapeFunctionTemplate::setMinimalPeakHeightForCalibration
     */
    double getMinimalPeakHeightForCalibration();

private:
    mutable PeakShapeT peakshape_;
    PeakParameterT peakparameter_;
};



////////////////////////////////////////////////////////////////////////////////
/* Some predefined PeakShapeFunctionTemplates for out-of-the-box application. */
////////////////////////////////////////////////////////////////////////////////

/**
* A peak shape function as it occurs in Orbitrap mass spectra.
*
* The Orbitrap peak shape function is parameterized via a linear sqrt model, which goes through the origin: @f$ f(x) = a\cdot x\sqrt{x}@f$ .
* You can set the model parameter a via the corresponding getter/setter method in the
* PeakShapeFunctionTemplate interface. Furthermore, you may autocalibrate the parameter calling the calibrateFor() method.
*
* This function is robust concerning autocalibration, because the autocalibration cannot set the parameter a so, that the function becomes invalid in some
* mz ranges. That's why we choose a PeakParameter, which is constrained in the origin.
* @see psf::RobustOrbitrapPeakShapeFunction 
*/
typedef PeakShapeFunctionTemplate<GaussianPeakShape, OrbitrapWithOriginFwhm, orbi> OrbitrapPeakShapeFunction;

/**
* A peak shape function as it occurs in centroided Orbitrap mass spectra.
*
* The function is similar to the @see psf::OrbitrapPeakShapeFunction, the only difference being that 
* the window shape is a box. Support threshold calculation etc are identical to the @see OrbitrapPeakShapeFunction.
* The Orbitrap peak shape function is parameterized via a linear sqrt model, which goes through the origin: @f$ f(x) = a\cdot x\sqrt{x}@f$ .
* You can set the model parameter a via the corresponding getter/setter method in the
* PeakShapeFunctionTemplate interface. Furthermore, you may autocalibrate the parameter calling the calibrateFor() method.
*
* This function is robust concerning autocalibration, because the autocalibration cannot set the parameter a so, that the function becomes invalid in some
* mz ranges. That's why we choose a PeakParameter, which is constrained in the origin.
* @see psf::RobustOrbitrapPeakShapeFunction 
*/
typedef PeakShapeFunctionTemplate<BoxPeakShape, OrbitrapWithOriginFwhm, orbi> OrbitrapBoxPeakShapeFunction;

/**
* A peak shape function with a gaussian shape static everywhere in a mass spectrum.
*
* You can set the full width at half maximum of the gaussian via the 'a' getter and setter 
* in the PeakShapeFunctionTemplate
* interface (the 'b' functions are not supported and produce compile time errors if used).
* Furthermore, you may autocalibrate the parameter calling the calibrateFor() method.
*/
typedef PeakShapeFunctionTemplate<GaussianPeakShape, ConstantFwhm, gaussian> GaussianPeakShapeFunction;




////////////////////
/* Implementation */
////////////////////

// operator()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
double 
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
operator()(const double referenceMass, const double observedMass) const {
    double supportThreshold = this->getSupportThreshold( referenceMass );        
    double massDifference = observedMass - referenceMass;

    if((-supportThreshold <= massDifference) && (massDifference <= supportThreshold)) {
        peakshape_.setFwhm(peakparameter_.at(referenceMass));
        return peakshape_.at(massDifference);
    }
    else {
        return 0.0;
    }
}

// getSupportThreshold()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
double 
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
getSupportThreshold(const double mz) const {
    peakshape_.setFwhm(peakparameter_.at(mz));
    return peakshape_.getSupportThreshold();        
}

// getType()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
inline
PeakShapeFunctionType
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
getType() {
    return PeakShapeFunctionTypeT;
}

// PeakShapeFunctionTemplate()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
PeakShapeFunctionTemplate() {
}

template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
PeakShapeFunctionTemplate(const double a) {
    this->setA(a);
}

template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
PeakShapeFunctionTemplate(const double a, const double b) {
    this->setA(a);
    this->setB(b);
}

// setA()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
inline
void
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
setA(const double a) {
    peakparameter_.setA(a);
}
// getA()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
inline
double
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
getA() const {
    return peakparameter_.getA();
}

// setB()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
inline
void
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
setB(const double b) {
    peakparameter_.setB(b);
}
// getB()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
inline 
double
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
getB() const {
    return peakparameter_.getB();
}

// calibrateFor()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
template< typename FwdIter, typename MzExtractor, typename IntensityExtractor >
void
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
calibrateFor(const MzExtractor& get_mz, const IntensityExtractor& get_int, FwdIter first, FwdIter last) {
  peakparameter_.learnFrom(get_mz, get_int, first, last);
}

// setMinimalPeakHeightForCalibration()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
void 
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
setMinimalPeakHeightForCalibration(const double minimalHeight) {
    peakparameter_.setMinimalPeakHeightToLearnFrom(minimalHeight);
}

// getMinimalPeakHeightForCalibration()
template <typename PeakShapeT, typename PeakParameterT, psf::PeakShapeFunctionTypes PeakShapeFunctionTypeT>
double
PeakShapeFunctionTemplate<PeakShapeT, PeakParameterT, PeakShapeFunctionTypeT>::
getMinimalPeakHeightForCalibration() {
    return peakparameter_.getMinimalPeakHeightToLearnFrom();
}

} /* namespace psf */

} /* namespace mstk */

#endif /*__PEAKSHAPEFUNCTION _HPP___*/

