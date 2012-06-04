/*
 * Stoichiometry.hpp
 *
 * Copyright (c) 2011,2012 Mathias Wilhelm
 * Copyright (c) 2009, 2010, 2011, 2012 Marc Kirchner
 * Copyright (c) 2010 Nathan Hueksen
 * Copyright (c) 2008 Thorben Kroeger
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

#ifndef __MSTK_INCLUDE_MSTK_AAS_STOICHIOMETRY_HPP__
#define __MSTK_INCLUDE_MSTK_AAS_STOICHIOMETRY_HPP__

#include "MSTK/aas/StoichiometryConfig.hpp"
#include "MSTK/aas/Element.hpp"
#include "MSTK/common/Types.hpp"

#include <map>
#include <vector>
#include <iostream>

namespace mstk {
namespace aas {
namespace stoichiometries {

/** @addtogroup mstk_aas
 * @{
 */

/**Representation of a stoichiometry.
 *
 * This class ensures that no element with count 0 exists in the mapping.
 */
class Stoichiometry
{

public:

    /** Convenience typedef of the data type.
     */
    typedef std::map<aas::elements::Element, mstk::Double> DataType;
    /** Convenience typedef of the data type const_iterator.
     */
    typedef DataType::const_iterator const_iterator;
    /** Convenience typedef of the data type iterator.
     */
    typedef DataType::iterator iterator;
    /** Convenience typedef of the data type value_type.
     */
    typedef DataType::value_type value_type;

    /** Creates an empty stoichiometry.
     * The annotation id is set to 0.
     */
    Stoichiometry();

    /**Sets the annotation id.
     * @param[in] id The annotation id
     */
    void setAnnotationId(const mstk::Int& id);

    /**Returns the annotation id.
     * @returns The annotation id
     */
    mstk::Int getAnnotationId(void) const;

    /**Clears the internal map.
     */
    void clear();

    /**Returns a const iterator referring to the first element in the stoichiometry.
     * @returns A const iterator to the beginning of the stoichiometry
     */
    const_iterator begin() const;

    /**Returns an iterator referring to the first element in the stoichiometry.
     * @returns An iterator to the beginning of the stoichiometry
     */
    iterator begin();

    /**Returns a const iterator referring to the past-the-end element in the stoichiometry.
     * @returns A const iterator to the element past the end of the stoichiometry
     */
    const_iterator end() const;

    /**Returns an iterator referring to the past-the-end element in the stoichiometry.
     * @returns An iterator to the element past the end of the stoichiometry
     */
    iterator end();

    /**Returns the number of elements in the stoichiometry.
     * @returns The number of elements that conform the stoichiometry content.
     */
    mstk::Size size() const;

    /**Returns whether the map container is empty, i.e. whether its size is 0.
     * @returns true if the vector size is 0, false otherwise.
     */
    mstk::Bool empty() const;

    /**Sets the count of an element.
     * @param[in] element The element
     * @param[in] count The count of the element
     */
    void set(const aas::elements::Element& element, const mstk::Double& count);

    /**Adds a certain amount of to an element.
     * In case the result is 0 the element is delted from the stoichiometry.
     * @param[in] element The element
     * @param[in] count The count of the element
     */
    void add(const aas::elements::Element& element, const mstk::Double& count);

    /**Returns the count of the element with the specified elementId.
     * @param[in] element Element id
     * @returns 0.0 if the elementId is not present in the stoichiometry, the
     * count otherwise
     */
    mstk::Double get(const aas::elements::Element& element) const;

    /**Checks whether all entries in the stoichiometry are non negative.
     * @returns True if all elements in the stoichiometry are non negative,
     * false otherwise
     */
    mstk::Bool nonNegative() const;

    /**Returns a human readable string containing the chemical formula of the
     * stoichiometry i.e. H(36)C(20)N(7)O(12)P(1)S(3)
     * @returns A string representing the stoichiometry, without the isotopic
     * data from the elements
     */
    mstk::String toString() const;

    /**Applies a stoichiometry configuration to this stoichiometry.
     * @param[in] config Stoichiomtery configuration
     */
    void applyStoichiometryConfiguration(const StoichiometryConfig& config);

    /**Recalculates the stoichiometry with the given stoichiometry configurtions
     * and returns a copy the result.
     * @param[in] config Stoichiometry config which is used to recalculate the
     * stoichiometry
     * @returns A copy of the recalculates stoichiometry
     */
    Stoichiometry recalculatesWithConfiguration(
        const StoichiometryConfig& config) const;

    /**Adds the content of the argument to the stoichiometry and returns a new
     * stoichiometry.
     *
     * @param[in] s Stoichiometry object
     * @returns A copy of *this + s
     */
    Stoichiometry operator+(const Stoichiometry& s);

    /**Subtracts the content of the argument from the stoichiometry and retruns
     * a new stoichiometry object.
     * @param[in] s Stoichiometry object
     * @returns A copy of *this - s
     */
    Stoichiometry operator-(const Stoichiometry& s);

    /**Adds a copy of the argument to the stoichiometry.
     * The new stoichiometry content is the content existing in the
     * stoichiometry object before the call plus the content of the argument.
     * @param[in] s Stoichiometry object
     * @return *this
     */
    Stoichiometry& operator+=(const Stoichiometry& s);

    /**Subtracts a copy of the argument to the stoichiometry.
     * The new stoichiometry content is the content existing in the
     * stoichiometry object before the call minus the content of the argument.
     * @param[in] s Stoichiometry object
     * @return *this
     */
    Stoichiometry& operator-=(const Stoichiometry& s);

    /**Sets a copy of the argument as the new content for the stoichiometry object.
     * The previous content is dropped.
     * @param[in] s Stoichiometry to copy
     * @returns *this
     */
    Stoichiometry& operator=(const Stoichiometry& s);

    /**Compares the stoichiometry against another.
     * @param[in] s Stoichiometry object to compare *this with
     * @returns true if both residues are the same, false otherwise
     */
    bool operator==(const Stoichiometry &s) const;

    /**Compares the stoichiometry against another, with opposite result of
     * Stoichiometry::operator==.
     * @param[in] s Stoichiometry object to compare *this with
     * @returns true if the stoichiometries are different, false otherwise.
     */
    bool operator!=(const Stoichiometry &s) const;

private:

    /** Annotation id.
     */
    mstk::Int annotationId_;
    /** Internal mapping of aas::elements::Element to the amount.
     */
    DataType counts_;

};
// class Stoichiometry

std::ostream& operator<<(std::ostream& o, const aas::stoichiometries::Stoichiometry& s);
std::ostream& operator<<(std::ostream& o,
    const std::vector<aas::stoichiometries::Stoichiometry>& s);

inline void Stoichiometry::setAnnotationId(const mstk::Int& id)
{
    annotationId_ = id;
}

inline mstk::Int Stoichiometry::getAnnotationId(void) const
{
    return annotationId_;
}

/** @\ */

} // namespace stoichiometries
} // namespace aas
} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_AAS_STOICHIOMETRY_HPP__ */
