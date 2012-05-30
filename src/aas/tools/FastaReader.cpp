/*
 * FastaReader.cpp
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

#include "MSTK/aas/tools/FastaReader.hpp"
#include "MSTK/common/Error.hpp"

#include <fstream>

namespace mstk {
namespace aas {
namespace tools {

FastaReader::FastaReader(const String& filename, const Digester& digester,
    const aas::AminoAcidSequence::ModificationList& fixedModifications) :
        filename_(filename), digester_(digester)
{
    typedef aas::AminoAcidSequence::ModificationList::const_iterator VSSI;
    try {
        // build modifications vectors: fixed mods
        for (VSSI modit = fixedModifications.begin();
                modit != fixedModifications.end(); ++modit) {
            fixedModifications_.push_back(*modit);
        }
    } catch (mstk::RuntimeError& e) {
        throw;
    }
}

FastaReader::~FastaReader()
{
}

void FastaReader::read(AminoAcidSequences& aminoAcidSequences) const
{
    DescSeq ds;
    // read the fasta file
    parse(ds);
    // digest the protein sequences and apply modifications
    digest(ds, aminoAcidSequences);
    // apply modifications
    modify(aminoAcidSequences);
}

void FastaReader::parse(DescSeq& ds) const
{
    // simple state machine to parse FASTA format files
    enum State
    {
        STATE_START, STATE_DESC, STATE_SEQ
    };
    State state = STATE_START;
    // open file
    std::ifstream ifs(filename_.c_str());
    if (!ifs) {
        mstk_fail("Cound not open " + filename_ + ".");
    }
    String desc, seq;
    while (!ifs.eof()) {
        String line;
        std::getline(ifs, line);
        if (line.empty()) {
            // ignore whitespace
            continue;
        }
        if (line[0] == '>' || line[0] == ';') {
            switch (state) {
                case STATE_START:
                    desc = line;
                    state = STATE_DESC;
                    break;
                case STATE_DESC:
                    // ignore
                    break;
                case STATE_SEQ:
                    // we finished a sequence. Store the sequence and
                    // the description.
                    ds.push_back(std::make_pair(desc, seq));
                    seq.clear();
                    state = STATE_DESC;
                    desc = line;
                    break;
            }
        } else {
            // we have a sequence
            switch (state) {
                case STATE_START:
                    // format error -- entries need to start with a description
                    mstk_fail("Syntax error in FASTA file.");
                    break;
                case STATE_DESC:
                    // start a new sequence
                    seq = line;
                    state = STATE_SEQ;
                    break;
                case STATE_SEQ:
                    // attach line to existing sequence
                    seq += line;
                    break;
            }
        }
    }
    // we allow empty files; hence the final state is either STATE_START
    // or STATE_SEQ. Everything else is an error.
    switch (state) {
        case STATE_DESC:
            mstk_fail("Syntax error in FASTA file.");
            break;
        case STATE_SEQ:
            // push the remaining sequence
            ds.push_back(std::make_pair(desc, seq));
            break;
        case STATE_START:
            // nothing to do
            break;
    }
}

void FastaReader::digest(const DescSeq& ds,
    AminoAcidSequences& aminoAcidSequences) const
{
    typedef Digester::AminoAcidSequences VAAS;
    VAAS peptides;
    for (DescSeq::const_iterator i = ds.begin(); i != ds.end(); ++i) {
        peptides.clear();
        digester_(i->second, peptides);
        for (VAAS::const_iterator pi = peptides.begin(); pi != peptides.end();
                ++pi) {
            aminoAcidSequences.push_back(*pi);
        }
    }
}

void FastaReader::modify(AminoAcidSequences& aminoAcidSequences) const
{
    for (AminoAcidSequences::iterator it = aminoAcidSequences.begin();
            it != aminoAcidSequences.end(); ++it) {
        it->applyFixedModifications(fixedModifications_);
    }
}

} // namespace tools
} // namespace aas
} // namespace mstk
