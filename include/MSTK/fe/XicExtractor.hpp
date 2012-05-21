/*
 * XicExtractor.hpp
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
#ifndef __MSTK_INCLUDE_MSTK_FE_XICEXTRACTOR_HPP__
#define __MSTK_INCLUDE_MSTK_FE_XICEXTRACTOR_HPP__

#include <MSTK/config.hpp>
#include <MSTK/common/Types.hpp>
#include <MSTK/common/Log.hpp>
#include <MSTK/fe/CentroidTraits.hpp>
#include <MSTK/fe/XicTraits.hpp>

#include <fbi/fbi.h>
#include <fbi/connectedcomponents.h>

namespace mstk {

namespace fe {

template<class Disambiguator, class Smoother, class Splitter>
class XicExtractor : public Disambiguator, public Smoother
{
public:
    template<class CentroidContainer, class XicContainer,
            class CentroidBoxGenerator>
    Size operator()(const CentroidContainer& centroids,
        const CentroidBoxGenerator& boxGenerator,
        const UnsignedInt minCardinality,
        const typename Splitter::ThresholdType splitThreshold,
        XicContainer& xics);
};

} // namespace fe

} // namespace mstk

//
// template implementation
//
#include <MSTK/common/CardinalityLessThan.hpp>

namespace mstk {

namespace fe {

template<typename Disambiguator, typename Smoother, typename Splitter>
template<typename CentroidContainer, typename XicContainer,
        typename CentroidBoxGenerator>
Size XicExtractor<Disambiguator, Smoother, Splitter>::operator ()(
    const CentroidContainer& centroids,
    const CentroidBoxGenerator& boxGenerator, const UnsignedInt minCardinality,
    const typename Splitter::ThresholdType splitThreshold, XicContainer& xics)
{
    MSTK_LOG(logDEBUG) << "XicExtractor::operator(): extracting XICs from "
            << centroids.size() << " centroids.";
    typedef typename CentroidContainer::value_type CentroidType;

    // calculate the adjacency list
    auto results = fbi::SetA<CentroidType, 0, 1>::intersect(centroids,
        boxGenerator, boxGenerator);

    // get the connected components
    typedef typename fbi::SetA<CentroidType, 0, 1>::IntType LabelType;
    std::vector < LabelType > labels;
    size_t nComponents = findConnectedComponents(results, labels);
    MSTK_LOG(logDEBUG) << "XicExtractor::operator(): Found " << nComponents
            << " primary XICs.";

    // Collect the centroids into pseudo-XICs.
    // We are using an internal vector of vectors because this allows
    // us to make the interface more generic, e.g. there is no guarantee that
    // other XIC implementations actually store their underlying centroids.
    typedef std::vector<CentroidType> CentroidSet; // a pseudo-XIC
    typedef std::vector<CentroidSet> CentroidSets; // a set of pseudo-XICs
    CentroidSets centroidSets;
    centroidSets.reserve(nComponents); // make sure we only allocate once
    centroidSets.resize(nComponents); // construct sets
    for (size_t i = 0; i < centroids.size(); ++i) {
        // push the i-th centroid into the label[i]-th centroid set
        CentroidSet& cs = centroidSets[labels[i] - 1]; // labels start at 1
        cs.push_back(centroids[i]);
    }

    // sort and get rid of any ambiguity in the sets
    typedef typename CentroidTraits<CentroidType>::LessThanRt Less;
    for (typename CentroidSets::iterator i = centroidSets.begin();
            i != centroidSets.end(); ++i) {
        std::sort(i->begin(), i->end(), Less());
        MSTK_LOG(logDEBUG) << "Have ambiguous XIC with " << i->size() << " entries.";
        i->erase(this->disambiguate(i->begin(), i->end()), i->end());
        MSTK_LOG(logDEBUG) << "Have disambiguated XIC with " << i->size() << " entries.";
    }

    // require a minimum cardinality for the centroid sets
    typedef CardinalityLessThan<std::vector<CentroidType> > InsufficientNumberOfCentroids;
    InsufficientNumberOfCentroids insufficientNumberOfCentroids(
        minCardinality);
    centroidSets.erase(
        std::remove_if(centroidSets.begin(), centroidSets.end(),
            insufficientNumberOfCentroids), centroidSets.end());
    MSTK_LOG(logDEBUG) << "XicExtractor::operator(): Found "
            << centroidSets.size() << " secondary XICs.";

    // Now walk through all XICs and attempt to split them.
    // Whenever a XIC is split, we collect the new parts in a splits vector
    // that is later merged with the existing XICs. We are working with two
    // containers because we don't know how many splits there are going to be
    // and need to make sure that there will not be any reallocations that
    // would invalidate the iterators.
    CentroidSets splits;
    typedef typename CentroidSets::iterator XI;
    for (XI i = centroidSets.begin(); i != centroidSets.end(); ++i) {
        CentroidSet smoothCopy(*i);
        this->smooth(smoothCopy.begin(), smoothCopy.end());
        Splitter splitter;
        splitter.split(i->begin(), i->end(), smoothCopy.begin(),
            smoothCopy.end(), splitThreshold);
        // store the first subXic in place
        *i = CentroidSet(splitter.begin()->first, splitter.begin()->second);
        // store the others
        if (splitter.size() > 1) {
            typename Splitter::const_iterator k = splitter.begin();
            std::advance(k, 1);
            while (k < splitter.end()) {
                splits.push_back(CentroidSet(k->first, k->second));
                ++k;
            }
        }
    }
    // join the lists
    centroidSets.insert(centroidSets.end(), splits.begin(), splits.end());
    splits.clear();
    MSTK_LOG(logDEBUG) << "XicExtractor::operator(): Found "
            << centroidSets.size() << " ternary XICs.";
    // once again, get rid of all XICs with an insufficient number of centroids
    centroidSets.erase(
        std::remove_if(centroidSets.begin(), centroidSets.end(),
            insufficientNumberOfCentroids), centroidSets.end());
    MSTK_LOG(logDEBUG) << "XicExtractor::operator(): Found "
            << centroidSets.size() << " quaternary XICs.";

    // convert the results into real XICs
    typedef typename XicTraits<typename XicContainer::value_type>::Creator Creator;
    xics.clear();
    std::transform(centroidSets.begin(), centroidSets.end(),
        std::back_inserter(xics), Creator());
    return xics.size();
}

} // namespace fe

} // namespace mstk

#endif
