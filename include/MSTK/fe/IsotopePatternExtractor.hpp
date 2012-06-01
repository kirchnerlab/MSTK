/*
 * IsotopePatternExtractor.hpp
 *
 * Copyright (C) 2012 Marc Kirchner
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
#ifndef __MSTK_INCLUDE_MSTK_FE_ISOTOPEPATTERNEXTRACTOR_HPP__
#define __MSTK_INCLUDE_MSTK_FE_ISOTOPEPATTERNEXTRACTOR_HPP__

#include <MSTK/config.hpp>
#include <MSTK/common/Types.hpp>
#include <MSTK/common/Log.hpp>

#include <fbi/fbi.h>
#include <fbi/connectedcomponents.h>

namespace mstk {

namespace fe {

template<class Correlator, class Splitter>
class IsotopePatternExtractor : public Correlator, public Splitter
{
public:
    template<class XicContainer, class IsotopePatternContainer,
            class XicBoxGenerators>
    Size operator()(const XicContainer& xics,
        const XicBoxGenerators& boxGenerators,
        const typename Correlator::ThresholdType& correlationThreshold,
        const UnsignedInt minCardinality,
        IsotopePatternContainer& isotopePatterns);
};

} // namespace fe

} // namespace mstk

//
// template implementation
//
#include <MSTK/common/CardinalityLessThan.hpp>
#include <MSTK/common/Error.hpp>
#include <MSTK/common/Log.hpp>
#include <MSTK/fe/QuickCharge.hpp>
#include <MSTK/fe/XicTraits.hpp>
#include <MSTK/fe/IsotopePatternTraits.hpp>

namespace mstk {

namespace fe {

template<class Correlator, class Splitter>
template<class XicContainer, class IsotopePatternContainer,
        class XicBoxGenerators>
Size IsotopePatternExtractor<Correlator, Splitter>::operator()(
    const XicContainer& xics, const XicBoxGenerators& boxGenerators,
    const typename Correlator::ThresholdType& correlationThreshold,
    const UnsignedInt minCardinality, IsotopePatternContainer& isotopePatterns)
{
    // precondition
    mstk_precondition(!boxGenerators.empty(),
        "Require at least one search box generator.");
    MSTK_LOG(logDEBUG) << "Extracting isotope patterns from " << xics.size() << "XICs.";
    typedef typename XicContainer::value_type XicType;
    typedef typename fbi::SetA<XicType, 0, 1> SetA;
    // carry out the box intersection
    typedef typename SetA::ResultType AdjList;
    AdjList adjList = SetA::intersect(xics,
        boxGenerators[0], boxGenerators);

    // Correlation filter. Keep adjacency list entries only if
    // the correlation between XICs along rt is above a user-defined
    // threshold.
    size_t nXics = xics.size();
    AdjList filteredAdjList(nXics);
    for (size_t i = 0; i < nXics; ++i) {
        typedef typename AdjList::value_type::const_iterator SCI;
        for (SCI j = adjList[i].begin(); j != adjList[i].end(); ++j) {
            // The adjacency list models an undirected graph and the
            // pearson correlation is a symmetric measure; hence, avoid
            // a double calculation, only run the test for one of (i, *j)
            // or (*j, i) and insert the result into both places in the
            // adjacency list.
            // Also, make sure that every vertex has an edge to itself.
            if (*j <= i) {
                if (i == *j
                        || this->correlate(xics[i].begin(), xics[i].end(),
                            xics[*j].begin(), xics[*j].end())
                                >= correlationThreshold) {
                    filteredAdjList[i].push_back(*j);
                    filteredAdjList[*j].push_back(i);
                }
            } else {
                break;
            }
        }
    }
    // get rid of old adjacency list
    adjList.clear();
    // find the connected components in the filtered adjancency list
    typedef typename fbi::SetA<XicType, 0, 1>::IntType LabelType;
    std::vector < LabelType > labels;
    size_t nComponents = findConnectedComponents(filteredAdjList, labels);
    MSTK_LOG(logDEBUG) << "Found " << nComponents << " primary isotope patterns.";

    // Collect the XICs into pseudo-isotope patterns.
    typedef std::vector<XicType> XicSet;
    typedef std::vector<XicSet> XicSets;
    XicSets ips;
    ips.clear(); // make sure it is empty
    ips.resize(nComponents); // construct in a single allocation
    XicSets(ips).swap(ips); // trim excess memory
    // now fill in data
    for (Size i = 0; i < xics.size(); ++i) {
        // push the i-th centroid into the label[i]-th XIC
        XicSet& ip = ips[labels[i] - 1]; // labels start at 1
        ip.push_back(xics[i]);
    }
    MSTK_LOG(logDEBUG) << "Extracted " << ips.size() << " primary isotope patterns.";
    // filter the pre-isotope patterns (based on size requirements etc)
    typedef CardinalityLessThan<XicSet> InsufficientNumberOfXics;
    InsufficientNumberOfXics insufficientNumberOfXics(minCardinality);
    ips.erase(std::remove_if(ips.begin(), ips.end(), insufficientNumberOfXics),
        ips.end());
    MSTK_LOG(logDEBUG) << "Found " << ips.size() << " secondary isotope patterns.";
    // now make sure all pips are sorted by m/z inside the patterns
    // and determine the charges present in the isotope pattern
    //QuickCharge qc;
    typedef typename XicSets::iterator XSI;
    for (XSI i = ips.begin(); i != ips.end(); ++i) {
        typedef typename XicTraits<XicType>::LessThanMz LessThanMz;
        std::sort(i->begin(), i->end(), LessThanMz());
        //std::vector<int> charges;
        //qc(i->begin(), i->end(), std::back_inserter(charges));
        //i->setCharges(std::set<int>(charges.begin(), charges.end()));
    }
    // split the isotope patterns using the Splitter policy class
    this->split(ips);
    //
    // make sure that the size requirements still hold
    ips.erase(
        std::remove_if(ips.begin(), ips.end(), insufficientNumberOfXics),
        ips.end());

    // convert the results into real XICs
    MSTK_LOG(logDEBUG) << "Transforming " << ips.size() << " isotope patterns.";
    typedef typename IsotopePatternTraits<
            typename IsotopePatternContainer::value_type>::Creator Creator;
    isotopePatterns.clear();
    std::transform(ips.begin(), ips.end(), std::back_inserter(isotopePatterns),
        Creator());
    return isotopePatterns.size();
}

} // namespace fe

} // namespace mstk

#endif
