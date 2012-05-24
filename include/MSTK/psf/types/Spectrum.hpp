/*
 * Spectrum.hpp
 *
 * Copyright (C) 2011 Bernhard Kausler
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
#ifndef __MSTK_INCLUDE_MSTK_PSF_SPECTRUM_HPP__
#define __MSTK_INCLUDE_MSTK_PSF_SPECTRUM_HPP__

#include <fstream>
#include <istream>
#include <string>

/**
 * @page mstk_spectrum MSTK/psf: spectrum sample implementation
 *
 * @section spectrumdatatye Spectrum data type
 *
 * MSTK/psf can work with any datatype representing a mass spectrum as long as you provide
 * mz and intensity extractors, that adhere to the expected interface.
 * We provide a sample implementation of a spectrum. Feel free to use your own.
 *
 * @section extractorinterface Extractor interface
 *
 * struct MyExtractor {
 *  typedef ... result_type; // numeric type used for the value (either mz or intensity)
 *  typedef ... element_type; // type of entries in a spectrum
 *
 *  result_type operator()( const element_type& e ) const;
 * };
 */

namespace mstk {

namespace psf {

/**
* a single entry in a mass spectrum, characterized by \f$mz\f$ and intensity
* values
*/
struct SpectrumElement {
    typedef double mz_type;
    typedef double intensity_type;

    SpectrumElement(const double m, const double i) : mz(m), intensity(i) {}

    double mz;
    double intensity;
};

struct MzExtractor {
  typedef SpectrumElement::mz_type result_type;
  typedef SpectrumElement element_type;

  result_type operator()( const SpectrumElement& e ) const {
    return e.mz;
  };
};

struct IntensityExtractor {
  typedef SpectrumElement::intensity_type result_type;
  typedef SpectrumElement element_type;

  result_type operator()( const SpectrumElement& e ) const {
    return e.intensity;
  };
};



/**
 * A mass spectrum is a sequence of SpectrumElements ordered by mz.
 */ 
typedef std::vector<SpectrumElement> Spectrum;

std::istream& operator>>(std::istream& is, Spectrum& s) {
    double mz, intensity;
    if (is.good()) {
        while (is >> mz >> intensity) {
	    if (intensity > 0) {
		s.push_back(SpectrumElement(mz, intensity));
	    }
        }
    }
    return is;
}

void loadSpectrumElements(Spectrum& s, const std::string& filename){
    std::ifstream ifs(filename.c_str());
    if (ifs.good()) {
        ifs >> s;
    }
}

} /* namespace psf */

} /* namespace mstk */

#endif /*__SPECTRUM _HPP___*/
