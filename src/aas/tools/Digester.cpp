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

#include "MSTK/aas/tools/Digester.hpp"
#include "MSTK/common/Error.hpp"

#include <iterator>
#include <set>

using namespace mstk;

namespace mstk {
namespace aas {
namespace tools {

String Digester::R_[] = { "(C)([^P])", // ARG_C_PROTEINASE
        "(.)(B|D)", // ASP_N_ENDOPETIDASE
        "(F|L|W|Y)([^P])", // CHYMOTRYPSIN
        "(K)([^P])", // LYSC
        "(F|L)(.)", // PEPSIN_A
        "(R|K)([^P])" // TRYPSIN
        };

Digester::Digester(const String& re) :
        regular_expression_(re), exceptions_(""), exceptions_enabled_(false)
{
}

Digester::Digester(const String& re, const String& exceptions) :
        regular_expression_(re), exceptions_(exceptions), exceptions_enabled_(
            true)
{
}
Digester::Digester(const EnzymEnum enzym) :
        regular_expression_(R_[enzym]), exceptions_(""), exceptions_enabled_(
            false)
{
}

void Digester::operator()(const aas::AminoAcidSequence& seq,
    AminoAcidSequences& frags, UnsignedInt missedCleavages) const
{
    if (regular_expression_.size() == 0) {
        frags.push_back(seq);
        return;
    }
    using namespace boost;

    String str = seq.toUnmodifiedSequenceString();
    String::const_iterator start = str.begin();
    String::const_iterator end = str.end();
    smatch m;

    // Build exception list
    std::set<String::const_iterator> cleavage_exceptions;
    if (exceptions_enabled_) {
        while (start != str.end() && regex_search(start, end, m, exceptions_)) {
            // Find the matched submatch
            end = str.end();
            for (UnsignedInt subMatch = 1; subMatch < m.size();
                    subMatch += 2) {
                if (m[subMatch].matched) {
                    end = m[subMatch].second;
                    break;
                }
            }
            // is this an allowed match?
            if (end == str.end()) {
                mstk_fail(
                    "Digester with invalid regular expression for exceptions found.");
            }
            // Add to unallowed
            cleavage_exceptions.insert(end);
            start = end;
            end = str.end();
        }
    }

    // Reset variables
    start = str.begin();
    end = str.end();
    AminoAcidSequence::const_iterator aa_start = seq.begin();
    AminoAcidSequence::const_iterator aa_end = seq.begin();
    String::const_iterator copy_start = str.begin();

    while ((start != str.end())
            && regex_search(start, end, m, regular_expression_)) {
        // Find the matched submatch
        end = str.end();
        for (UnsignedInt subMatch = 1; subMatch < m.size(); subMatch += 2) {
            if (m[subMatch].matched) {
                end = m[subMatch].second;
                break;
            }
        }
        //Is this an allowed mathc?
        if (end == str.end()) {
            mstk_fail("Digester with invalid regular expression found.");
        }
        // Check if this is allowed
        if (cleavage_exceptions.count(end)) {
            // Not allowed! update everything but copy_start
            start = end;
            end = str.end();
            continue;
        }
        // Translate the string iterators
        aa_start = seq.begin()
                + std::distance((String::const_iterator) str.begin(),
                    copy_start) + 1;
        aa_end = seq.begin()
                + std::distance((String::const_iterator) str.begin(), end)
                + 1;
        // If end is just before the C-term, set it to the C-term
        if (aa_end + 1 == seq.end()) {
            aa_end = seq.end();
        }
        // If start is just after the N-term, set it to the N-term
        if (aa_start == seq.begin() + 1) {
            aa_start = seq.begin();
        }

        frags.push_back(AminoAcidSequence(aa_start, aa_end));
        start = end;
        copy_start = end;
        end = str.end();
    }
    // Set new testing start to last end
    aa_start = aa_end;
    // something left?
    if (aa_start != seq.end() && aa_start != seq.end()) {
        frags.push_back(AminoAcidSequence(aa_start, seq.end()));
    }

    // Deal with missed cleavages
    if (missedCleavages > 0) {
        std::vector < AminoAcidSequence > missedCleavageSequences;
        for (std::vector<AminoAcidSequence>::iterator it = frags.begin();
                it != frags.end(); ++it) {
            for (unsigned int i = 1; i <= missedCleavages; ++i) {
                AminoAcidSequence s("");
                for (unsigned int j = 0; j <= i; ++j) {
                    if ((it + j) == frags.end()) {
                        break;
                    }
                    s.append(*(it + j));
                }
                missedCleavageSequences.push_back(s);
            }
        }
        std::copy(missedCleavageSequences.begin(),
            missedCleavageSequences.end(), std::back_inserter(frags));
    }

    // Sort
    std::sort(frags.begin(), frags.end(),
        AminoAcidSequence::LessThanSequenceUnmodified());
    // Remove duplicates
    std::vector<AminoAcidSequence>::iterator new_end = std::unique(
        frags.begin(), frags.end(),
        AminoAcidSequence::EqualToSequenceUnmodified());
    frags.erase(new_end, frags.end());
}

} // namespace tools
} // namespace aas
} // namespace mstk
