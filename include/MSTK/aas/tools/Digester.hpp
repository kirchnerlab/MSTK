/*
 * Digester.hpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2009,2010,2011,2012 Marc Kirchner
 * Copyright (c) 2010 Nathan Huesken
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_DIGESTER_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_DIGESTER_HPP__

#include "MSTK/aas/AminoAcidSequence.hpp"
#include "MSTK/common/Types.hpp"

#include <boost/regex.hpp>

#include <vector>

namespace mstk {
namespace aas {
namespace tools {

/** @addtogroup mstk_aas
 * @{
 */

/**General digests class.
 * Digests peptides (AminoAcidSequences) using regular expressions.
 *
 * The regular expression are expected to be build in a way that the cleavage
 * appeares after the first submatch.
 */
class Digester
{
public:

    /**Convenience typedef of a list of amino acid sequences.
     */
    typedef std::vector<aas::AminoAcidSequence> AminoAcidSequences;

    /** List of enzymes for which the regular expression is stored in R.
     */
    enum EnzymEnum
    {
        ARG_C_PROTEINASE = 0,
        ASP_N_ENDOPETIDASE,
        CHYMOTRYPSIN,
        LYSC,
        PEPSIN_A,
        TRYPSIN
    };

    /** Constructor.
     *  The regular expression MUST be in the format.
     *  (beforeCleavage1)(afterCleavage1)|(beforeCleavage2)(afterCleavage2)| ...
     *  The parentheses are mandatory and the expression in the
     *  parentheses may not contain any parentheses resulting in more sub
     *  expressions.
     *  This ensures that the cleavage point is found by matching
     *  subexpressions.
     *  @param[in] re Regular expression identifying the positions of cleavages.
     */

    Digester(const mstk::String& re);

    /** Constructor taking a list of exceptions.
     *  The the same limitations applying to the regular expression apply
     *  to the exceptions. The exceptions mark points in the
     *  peptide where no cleavage may occur by the enzym.
     *  @param[in] re Regular expression identifying the positions of cleavages.
     *  @param[in] exceptions Regular expression identifying exceptions for cleavages.
     */
    Digester(const mstk::String& re, const mstk::String& exceptions);

    /** Constructor building a pre-defined enzym.
     *  Build a digester after the rules of a hardcoded enzym.
     *  @param[in] enzym Element of the EnzymElement list identifying the enzym.
     */
    Digester(const EnzymEnum enzym);

    /**
     * Digests the peptide \a seq; \a missedCleavages
     * indicates the maximum number of cleavages that are omitted.
     * @param[in] seq The sequence to digest.
     * @param[out] frags Vector that will be filled with the resulting sequences
     * @param[in] missedCleavages maximum number of missed cleavages that should
     *            be included in the result.
     * @return possible fragments are stored in \a frags.
     */
    void operator()(const aas::AminoAcidSequence& seq,
        AminoAcidSequences& frags, mstk::UnsignedInt missedCleavages = 0) const;

private:

    /**Set of regular expressions for the build in digesters.
     */
    static mstk::String R_[];

    /**Regular expression of the digester.
     */
    boost::regex regular_expression_;

    /**Regular expressions marking the exceptions to the cleavage points given
     * by \a regular_expression.
     */
    boost::regex exceptions_;

    /**Indicates, if the exceptions should be used.
     */
    mstk::Bool exceptions_enabled_;
};

/** @\ */

} // namespace tools
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_DIGESTER_HPP__ */
