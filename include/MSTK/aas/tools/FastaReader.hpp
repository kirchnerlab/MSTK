/*
 * FastaReader.hpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2010,2011,2012 Marc Kirchner
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_FASTAREADER_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_FASTAREADER_HPP__

#include "MSTK/aas/AminoAcidSequence.hpp"
#include "MSTK/aas/tools/Digester.hpp"
#include "MSTK/common/Types.hpp"

#include <boost/regex.hpp>

namespace mstk {
namespace aas {
namespace tools {

/** @addtogroup mstk_aas
 * @{
 */

/**General fasta-file reader.
 *
 */
class FastaReader
{
public:

    /**Convenience typedef for a list of amino acid sequences.
     */
    typedef std::vector<AminoAcidSequence> AminoAcidSequences;

    /**Constructor.
     *
     * @param[in] filename fasta file
     * @param[in] digester Digester used to digest the sequences
     * @param[in] fixedModifications Fixed modifications which are applied after the digester
     */
    FastaReader(const mstk::String& filename, const Digester& digester =
            Digester(""),
        const aas::AminoAcidSequence::ModificationList& fixedModifications =
                aas::AminoAcidSequence::ModificationList());

    /**Destructor.
     */
    ~FastaReader();

    /**Reads the fasta file, digests all sequences and applies fixed modifications.
     * @param[out] aminoAcidSequences Resulting amino acid sequences
     */
    void read(AminoAcidSequences& aminoAcidSequences) const;

private:

    /**Convenience typedef for a fasta entry list.
     */
    typedef std::vector<std::pair<mstk::String, mstk::String> > DescSeq;

    /**Parses the fasta file.
     *
     * @param[out] ds A list containing the header information and the sequence of
     */
    void parse(DescSeq& ds) const;

    /**Digests the given list of sequences.
     *
     * @param[in] ds List of header information (first) and the sequence (second)
     * @param[out] aminoAcidSequences List of resulting amino acid sequences
     */
    void digest(const DescSeq& ds,
        AminoAcidSequences& aminoAcidSequences) const;

    /**Applies fixed modifcations to all amino acid sequences.
     *
     * @param[in,out] aminoAcidSequences List of amino acid sequences
     */
    void modify(AminoAcidSequences& aminoAcidSequences) const;

    /**Fasta file name.
     */
    mstk::String filename_;

    /**Digester which is used to digest the amino acid sequences.
     */
    Digester digester_;

    /**List of modifications which are applied as fixed modifications.
     */
    aas::AminoAcidSequence::ModificationList fixedModifications_;
};

/** @\ */

} // namespace tools
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_FASTAREADER_HPP__ */
