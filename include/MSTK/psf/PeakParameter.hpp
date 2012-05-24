/*
 * PeakParameter.hpp
 *
 * Copyright (C) 2009-2011 Bernhard Kausler
 * Copyright (C) 2009-2011 Marc Kirchner
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
#ifndef __MSTK_INCLUDE_MSTK_PSF_PEAKPARAMETER_HPP__
#define __MSTK_INCLUDE_MSTK_PSF_PEAKPARAMETER_HPP__

#include <MSTK/config.hpp>
#include <utility>
#include <vector>

#include <MSTK/common/Error.hpp>
#include <MSTK/psf/Error.hpp>
#include <MSTK/common/Log.hpp>
#include <MSTK/psf/SpectrumAlgorithm.hpp>

#include <vigra/matrix.hxx>
#include <vigra/regression.hxx>

namespace mstk {
    
namespace psf {
/**
 * The slope (including the bias) of a multidimensional linear function.
 *
 * Example: For the function @f$ f(x_1,x_2) = ax_1 + bx_2 + c + 3.0 @f$, the generalized slope is
 * @f$ (a,b,c,3.0) @f$ with a bias of 3.0 .
 */
typedef std::vector<double> GeneralizedSlope;


// class ConstantModel
/**
 * @f$ f(x) = a @f$
 *
 * ConstantModel is a policy class. The model depends on one parameter a.
 * @see http://en.wikipedia.org/wiki/Policy_based_design
 * @see psf::PeakParameterModel
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSTK_EXPORT ConstantModel
{
public:
    /**
     * This model has one parameter.
     */
    unsigned int numberOfParameters();

    void setParameter(unsigned index, double value);
    double getParameter(unsigned index);

protected:
    /**
     * Value of the model at position x.
     */
    double at(const double x) const;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
    ~ConstantModel() {}

    /**
     * Equal to @f$ (1.,0.) @f$.
     */
    GeneralizedSlope slopeInParameterSpaceFor(double x) const;

public:
    ConstantModel() : a_(0.1) {}

    void setA(const double a);
    double getA() const;

private:
    static const unsigned numberOfParameters_ = 1;
    double a_;
};


// class LinearSqrtModel()
/**
 * @f$ f(x) = a\cdot x\sqrt{x} + b @f$
 *
 * LinearSqrtModel is a policy class. The model depends on two parameter a and b.
 * @see http://en.wikipedia.org/wiki/Policy_based_design
 * @see psf::PeakParameterModel
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSTK_EXPORT LinearSqrtModel
{
public:
    /**
     * This model has two parameters.
     */
    unsigned int numberOfParameters();

    void setParameter(unsigned index, double value);
    double getParameter(unsigned index);

protected:
    /**
     * Value of the model at position x.
     */
    double at(const double x) const;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
    ~LinearSqrtModel() {}

    /**
     * Equal to @f$ (x\sqrt{x},1.,0.) @f$.
     */
    GeneralizedSlope slopeInParameterSpaceFor(double x) const;

public:
    LinearSqrtModel() : a_(0.1), b_(0.1) {}

    void setA(const double a);
    double getA() const;

    void setB(const double b);
    double getB() const;

private:
    static const unsigned numberOfParameters_ = 2;
    double a_, b_;
};


// class LinearSqrtOriginModel()
/**
 * @f$ f(x) = a\cdot x\sqrt{x}@f$
 *
 * LinearSqrtOriginModel is a policy class. The model depends on one parameter a.
 * @see http://en.wikipedia.org/wiki/Policy_based_design
 * @see psf::PeakParameterModel
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSTK_EXPORT LinearSqrtOriginModel
{
public:
    /**
     * This model has two parameters.
     */
    unsigned int numberOfParameters();

    void setParameter(unsigned index, double value);
    double getParameter(unsigned index);

protected:
    /**
     * Value of the model at position x.
     */
    double at(const double x) const;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
    ~LinearSqrtOriginModel() {}

    /**
     * Equal to @f$ (x\sqrt{x},0.) @f$.
     */
    GeneralizedSlope slopeInParameterSpaceFor(double x) const;

public:
    LinearSqrtOriginModel() : a_(0.1) {}

    void setA(const double a);
    double getA() const;

private:
    static const unsigned numberOfParameters_ = 1;
    double a_;
};


// class SqrtModel
/**
 * @f$ f(x) = a\cdot \sqrt{x} + b @f$
 *
 * SqrtModel is a policy class. The model depends on two parameter a and b.
 * @see http://en.wikipedia.org/wiki/Policy_based_design
 * @see psf::PeakParameterModel
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSTK_EXPORT SqrtModel
{
public:
    /**
     * This model has two parameters.
     */
    unsigned int numberOfParameters();

    void setParameter(unsigned index, double value);
    double getParameter(unsigned index);

protected:
    /**
     * Value of the model at position x.
     */
    double at(const double x) const;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
    ~SqrtModel() {}

    /**
     * Equal to @f$ (\sqrt{x},1.,0.) @f$.
     */    
    GeneralizedSlope slopeInParameterSpaceFor(double x) const;

public:
    SqrtModel() : a_(0.1), b_(0.1) {}

    void setA(const double a);
    double getA() const;

    void setB(const double b);
    double getB() const;

private:
    static const unsigned numberOfParameters_ = 2;
    double a_, b_;
};


// class QuadraticModel
/**
 * @f$ f(x) = a\cdot x^2 + b @f$
 *
 * QuadraticModel is a policy class. The model depends on two parameter a and b.
 * @see http://en.wikipedia.org/wiki/Policy_based_design
 * @see psf::PeakParameterModel
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class MSTK_EXPORT QuadraticModel
{
public:
    /**
     * This model has two parameters.
     */
    unsigned int numberOfParameters();

    void setParameter(unsigned index, double value);
    double getParameter(unsigned index);

protected:
    /**
     * Value of the model at position x.
     */
    double at(const double x) const;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
    ~QuadraticModel() {}

    /**
     * Equal to @f$ (x^2,1.,0.) @f$.
     */ 
    GeneralizedSlope slopeInParameterSpaceFor(double x) const;

public:
    QuadraticModel() : a_(0.1), b_(0.1) {}

    void setA(const double a);
    double getA() const;

    void setB(const double b);
    double getB() const;

private:
    static const unsigned numberOfParameters_ = 2;
    double a_, b_;
};



// class ParameterModel
/**
 * The interface for the peak parameter policy to be used in psf::PeakParameterFwhm.
 *
 * @attention This interface exists only for the purpose of documentation. Don't inherit
 *      from it, but simply implement it.
 * @see psf::PeakParameterFwhm
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
class ParameterModel
{
// mandatory interface
protected:
     /**
     * Value of the model at position x.
     */
     virtual double at(double x) const = 0;

    /**
     * Protected non-virtual destructor.
     *
     * This is a so called 'policy class'. The host class is inheriting from it to add its
     * functionality. But it is not intended to be used as a base class. In case the user
     * utilizes the policy as a baseclass anyway, we have to prevent calling 'delete' on it
     * (deleting a baseclass with a non-virtual destructor leads to undefined behaviour).
     * So we declare the destructor 'protected'.
     * @see http://en.wikipedia.org/wiki/Policy_based_design
     */
     virtual ~ParameterModel() {}

// optional interface
public:
     /**
      * The number of parameters in the model.
      *
      * May be any positive number including zero.
      * This is equivalent to the dimension of parameter space.
      * @attention: Don't implement, if the model has no parameters. This will cause
      * desirable compile time errors, if one tries to learn a model without parameters.
      */
     virtual unsigned int numberOfParameters() = 0;

     /**
      * 0 <= index < this->numberOfParameters() .
      * 
      * @attention: Don't implement, if the model has no parameters. This will cause
      * desirable compile time errors, if one tries to learn a model without parameters.
      * @throw psf::PreconditionViolation Parameter index is out of range or the number of
      *                                  parameters is zero.
      */
     virtual void setParameter(unsigned index, double value) = 0;
     /**
      * 0 <= index < this->numberOfParameters() .
      * 
      * @attention: Don't implement, if the model has no parameters. This will cause
      * desirable compile time errors, if one tries to learn a model without parameters.
      * @throw psf::PreconditionViolation Parameter index is out of range or the number of
      *                                  parameters is zero.
      */
     virtual double getParameter(unsigned index);
     

protected:
    /**
     * The slope of the linear function representing the model in parameter space.
     *
     * In parameter space the coordinates @f$ \vec{x} @f$ and
     * regular parameters @f$ \vec{p} @f$ of the model are switching their roles. Some models may have
     * a linear representation in this space: @f$ model(\vec{p}) = \vec{m}\cdot\vec{p} @f$.
     * @f$ \vec{m} @f$ is the generalized slope including the bias. In a n-dimensional
     * parameter space,  the  generalized slope is (n+1)-dimensional. So, the parameter vector 
     * @f$ \vec{p} @f$  has n+1 elements with the last element always set to unity.
     * Example: For the model @f$ f(x_1,x_2) = ax_1 + bx_2 + c @f$ with parameters a,b and c,
     * the generalized slope in parameter space is
     * @f$ (x_1,x_2,1,0) @f$ and @f$ \vec{p}= (a,b,c,1)^T @f$.
     * @attention This function is part of the optional interface, since not every model may have
     * a linear representation in parameter space.
     *
     * @see psf::PeakParameterModel::GeneralizedSlope
     * 
     * @param x The coordinate x now playing the role of a parameter.
     * @return The generalized slope, a mulitdimensional vector including the bias.
     */
    virtual GeneralizedSlope slopeInParameterSpaceFor(double x) const = 0;
};



// class PeakParameterFwhm
/**
 * 'Full width at half maximum' peak shape parameter.
 * 
 * The FWHM can be used to paramterize a peak shape.
 * 
 * Usually, the FWHM depends on the mass over charge ratio in a mass spectrum. This 
 * dependency can be linear, quadratic or any other function. Therefore, the
 * PeakParameterFwhm is a template, which has to be instantiated with a concrete
 * ParameterModel. Note, that the model should yield non-negative values for non-negative
 * parameters and mz values (There is no such thing as a negative FWHM). Else, you will
 * probably get very unexpected behaviour. 
 *
 * @see psf::ConstantModel, 
 *      psf::LinearSqrtModel,
        psf::LinearSqrtOriginModel
 *      psf::SqrtModel,
 *      psf::QuadraticModel
 * @see psf::PeakShape
 * @see psf::GaussianPeakShape
 *
 * @author Bernhard Kausler <bernhard.kausler@iwr.uni-heidelberg.de>
 */
template <typename ParameterModel>
class MSTK_EXPORT PeakParameterFwhm : public ParameterModel 
{
public:
    PeakParameterFwhm() : minimalPeakHeightToLearnFrom_(0) {}

    /**
     * The FWHM at a specific mass channel.
     *
     * @param mz Mass channel; has to be positive.
     * @throw psf::PreconditionViolation Parameter mz is not positive.
     * @throw psf::PostconditionViolation The computed fwhm is negative or zero. This may
     *      be caused by an invalid ParameterModel.
     */
    double at(const double mz) const {
        mstk_precondition(mz > 0, "PeakParameterFwhm::at(): Parameter mz has to be positive.");
        double fwhm = this->ParameterModel::at(mz);
        mstk_postcondition(fwhm > 0, "PeakParameterFwhm::at(): Model returned negative or zero fwhm.");
        return fwhm;
    }

    /**
     * Calibrates the internal model for a specific mass spectrum.
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
     * @param first Points to the first Element of the sequence.
     * @param last Points to one past the last Element of the sequence.
     *
     * @throw mstk::psf::Starvation Insufficient data to carry out a calibration.
     */
    template< typename FwdIter, typename MzExtractor, typename IntensityExtractor >
    void learnFrom(const MzExtractor&, const IntensityExtractor&, FwdIter first, FwdIter last);

    // setMinimalPeakHeightToLearnFrom()
    /**
     * Only use peaks with a minimum absolute intensity to learn from.
     *
     * @param minimalHeight May be negative, albeit that's not meaningful.
     */
    void setMinimalPeakHeightToLearnFrom(double minimalHeight);
    
    // getMinimalPeakHeightToLearnFrom()
    /**
     * Only peaks with this absolute minimal intensity are used for learning.
     */
    double getMinimalPeakHeightToLearnFrom();

private:
    static const double fractionOfMaximum_;

    double minimalPeakHeightToLearnFrom_; 

    /**
     * Fit the parameter model to measured mz-width pairs.
     *
     * Internally, a linear regression using least squares is done. The PeakParameterModel
     * has to support the optional slopeInParameterSpaceFor() function for this regression
     * to work. Else, calling this function will not compile.
     *
     * If new regression methods should be added, move the implementation of the learn_
     * function to policy classes.
     *
     * @throw psf::PreconditionViolation Parameter pairs is an empty vector.
     * @throw psf::InvariantViolation Numerical regression algorithm failed.
     */
    template< typename MzExtractor >
    void learn_(const std::vector<std::pair<typename MzExtractor::result_type, typename MzExtractor::result_type> >& pairs);
};

/**
 * Fwhm as it occurs in an Orbitrap mass spectrum.
 */
typedef PeakParameterFwhm<LinearSqrtModel> OrbitrapFwhm;
/**
 * Fwhm as it occurs in an Orbitrap mass spectrum. Is zero at zero Dalton.
 */
typedef PeakParameterFwhm<LinearSqrtOriginModel> OrbitrapWithOriginFwhm;

/**
 * Fwhm as it occurs in a FT-ICR mass spectrum.
 */
typedef PeakParameterFwhm<QuadraticModel> FtIcrFwhm;

/**
 * Fwhm as it occurs in a TOF mass spectrum.
 *
 * The specific Time-of-Flight mass analyzer should measure time internally (not velocity
 * or energy), for this peak parameter to be applicable.
 */
typedef PeakParameterFwhm<SqrtModel> TofFwhm; 

/**
 * A FWHM independent of the mass channel.
 */
typedef PeakParameterFwhm<ConstantModel> ConstantFwhm;



////////////////////
/* implementation */
////////////////////

template <typename ParameterModel>
const double PeakParameterFwhm<ParameterModel>::fractionOfMaximum_ = 0.5;  

// learnFrom()
template <typename ParameterModel>
template< typename FwdIter, typename MzExtractor, typename IntensityExtractor >
void PeakParameterFwhm<ParameterModel>::learnFrom(const MzExtractor& get_mz, const IntensityExtractor& get_int, FwdIter first, FwdIter last) {
    using namespace std;
    typedef std::vector<std::pair<typename MzExtractor::result_type, typename MzExtractor::result_type> > MzWidthPairs_;
 
    // sample some FWHMs from the spectrum
    MzWidthPairs_ pairs = measureFullWidths(get_mz, get_int, first, last, fractionOfMaximum_, getMinimalPeakHeightToLearnFrom());        

    if(pairs.empty()) {
        throw psf::Starvation("PeakParameterFwhm::learnFrom(): No (Mz | FWHM) could be measured in input spectrum to learn from.");
    }

    // fit the PeakParameterModel to the measured data
    try {
        learn_<MzExtractor>(pairs);
    } catch(const mstk::InvariantViolation& e) {
		MSTK_UNUSED(e);
        MSTK_LOG(logWARNING) << "PeakParameterFwhm::learnFrom(): Numerical regression failed.";
        throw psf::Starvation("PeakParameterFwhm::learnFrom(): Regression of the parameter model for the measured (Mz | FWHM) pairs failed.");
    }
    
    MSTK_LOG(logINFO) << "Learned peak parameter FWHM from spectrum. FWHM at 400 Th is now " << at(400)  << " Th. This corresponds to a resolution of " << 400./at(400) << ".";
}

template <typename ParameterModel>
void PeakParameterFwhm<ParameterModel>::setMinimalPeakHeightToLearnFrom(const double minimalHeight) {
    minimalPeakHeightToLearnFrom_ = minimalHeight;
}

template <typename ParameterModel>
double PeakParameterFwhm<ParameterModel>::getMinimalPeakHeightToLearnFrom() {
    return minimalPeakHeightToLearnFrom_;
}

// learn_()
template <typename ParameterModel>
template< typename MzExtractor >
void PeakParameterFwhm<ParameterModel>::learn_(const std::vector<std::pair<typename MzExtractor::result_type, typename MzExtractor::result_type> >& pairs) {
    using namespace vigra;
    typedef std::vector<std::pair<typename MzExtractor::result_type, typename MzExtractor::result_type> > MzWidthPairs_;
    mstk_precondition(pairs.empty() == false, "PeakParameterFwhm::learn_(): Called with empty input vector. This is not supposed to happen. A bug in the code preceding the call of learn_ probably caused it.");

    // We now fit the PeakParameterModel to the spectrum using linear regression.    
    // We minimize the residue |A*x - b|^2.
    // b is a column vector with all measured widths.
    // x is a column vector with the model parameters, for example (a, b)^T.
    // A is a matrix such that A*x resembles the original model in every row for the
    // corresponding mz value.
    // x gets optimized.

    /* Construct A and b */

    // Models with no parameter should not implement the numberOfParameters() function.
    // Calling it here while trying to learn this non-existing parameters would cause a
    // desirable compile time error. Nevertheless, we have to guard against wrongly implemented
    // ParameterModels.
    mstk_invariant(this->ParameterModel::numberOfParameters() > 0, "PeakParameterFwhm::learn_(): Number of model parameters is not greater than zero.");
    
    // A: #rows is number of measured pairs; #columns is dimension of parameter space
    linalg::Matrix<double> A(pairs.size(), this->ParameterModel::numberOfParameters());
    // b: column vector with as many elements as measured pairs
    linalg::Matrix<double> b(pairs.size(), 1);
    GeneralizedSlope slope;

    // Calc GeneralizedSlope for every measured pair and store it as rows of A
    // Store the measured width in b
    for(typename MzWidthPairs_::size_type pairIndex = 0; pairIndex < pairs.size(); ++pairIndex) {
        slope = this->ParameterModel::slopeInParameterSpaceFor(pairs[pairIndex].first);

        // copy the slope into a row of A
        // (we ignore the bias, because it can't be optimized.)
        mstk_invariant((slope.size()-1) == static_cast<GeneralizedSlope::size_type>(A.columnCount()), "PeakParameterFwhm::learn_(): Generalized slope has different dimension than the space, it is living in.");        
        for(linalg::Matrix<double>::difference_type_1 column = 0; column < A.columnCount(); ++column) {      
            A(static_cast<linalg::Matrix<double>::difference_type_1>(pairIndex), column) = slope.at(static_cast<GeneralizedSlope::size_type>(column));   
        }       

        b(static_cast<linalg::Matrix<double>::difference_type_1>(pairIndex), 0) = pairs[pairIndex].second;    
    }

    /* do least squares */
    // result: the optimized parameters
    // Note, that we don't include the bias.
    linalg::Matrix<double> x(this->ParameterModel::numberOfParameters(),1);
    
    // We have to enforce a positive FWHM for positive mz values, so we use a non-negative
    // least squares with x guaranteed to be non-negative. The model then has to yield
    // positive values too, of course. But that is in the responsibility of the caller. 
    linalg::nonnegativeLeastSquares(A, b, x);

    /* set the fitted parameters */
    for(linalg::Matrix<double>::difference_type_1 index = 0; index < x.rowCount(); ++index) {
        MSTK_LOG(logDEBUG2) << "PeakParameterFwhm::learn_(): Parameter " << index << " found: " <<  x(index, 0);
        
        mstk_invariant(index >= 0, "PeakParameterFwhm::learn_(): index may not be negative before calling setParameter(index,0).");
        this->ParameterModel::setParameter(index, x(static_cast<unsigned>(index), 0));    
    }
}

} /* namespace psf */

} /* namespace mstk */

#endif /*__PEAKPARAMETER _HPP___*/
