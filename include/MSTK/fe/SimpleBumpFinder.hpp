/*
 * SimpleBumpFinder.hpp
 *
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
#include <MSTK/config.hpp>
#include <MSTK/common/Error.hpp>
#include <MSTK/common/Log.hpp>
#include <MSTK/fe/utils.hpp>
#include <iterator>
#include <utility>

#ifndef __MSTK_INCLUDE_MSTK_FE_SIMPLEBUMPFINDER_HPP__
#define __MSTK_INCLUDE_MSTK_FE_SIMPLEBUMPFINDER_HPP__

namespace mstk {

namespace fe {

/** Simple bump finder.
 * This implements a simple finite state machine that detects bumps in
 * a sequence of values using templated iterators and a templated comparison
 * operator. The class is intended as a base class in a policy-based
 * design [Alexandrescu, 2002]. It implements a single, non-virtual
 * function \c findBump that takes a range of values and returns iterators to
 * the left and right boundaries of the first bump in the range.
 * See the tests for usage examples.
 */
class SimpleBumpFinder
{
protected:
    /** Virtual destructor.
     */
    virtual ~SimpleBumpFinder() = 0;

    /** Finds the first bump in a range.
     * @param first Iterator, pointing at the beginning of the range.
     * @param last Iterator, pointing one beyond the last element of the range.
     * @return A pair of iterators, where the first and second element point at
     *         the beginning and (one beyond) the end of the first bump in the
     *         range [first, last).
     */
    template<typename InputIterator>
    std::pair<InputIterator, InputIterator> findBump(InputIterator first,
        InputIterator last);
private:
    /** Compare of there is a positive slope between \c left and \c right.
     * @return True if \c comp(*left, *right) evaluates to true.
     */
    template<typename InputIterator, typename Compare>
    bool slopeUp(InputIterator left, InputIterator right, Compare comp);

    /** Compare of there is a negative slope between \c left and \c right.
     * @return True if \c comp(*right, *left) evaluates to true.
     */
    template<typename InputIterator, typename Compare>
    bool slopeDown(InputIterator left, InputIterator right, Compare comp);

    enum State
    {
        STATE_START, STATE_RAMP_UP, STATE_RAMP_DOWN, STATE_BUMP, STATE_STOP
    };
};

} // namespace fe

} // namespace mstk

//
// template implementation
//

namespace mstk {

namespace fe {

template<typename InputIterator, typename Compare>
bool SimpleBumpFinder::slopeUp(InputIterator left, InputIterator right,
    Compare comp)
{
    return comp(*left, *right);
}

template<typename InputIterator, typename Compare>
bool SimpleBumpFinder::slopeDown(InputIterator left, InputIterator right,
    Compare comp)
{
    return comp(*right, *left);
}

template<typename InputIterator>
std::pair<InputIterator, InputIterator> SimpleBumpFinder::findBump(
    InputIterator first, InputIterator last)
{
    MSTK_LOG(logDEBUG2)
            << "findBump(): spectrum size = " << std::distance(first, last);
    // empty spectrum
    if (first == last) {
        return std::make_pair(first, last);
    }

    // set up iterators, only require prefix ++
    InputIterator left = first, right = first;
    ++right;

    // single entry
    if (right == last) {
        return std::make_pair(left, right);
    }

    // abundance lessThan comparison operator
    typedef typename std::iterator_traits<InputIterator>::value_type ValueType;
    LessByAbundance<ValueType> comp;

    // FSM
    State state = STATE_START;
    while (state != STATE_STOP) {
        // in the state machine, always check the end condition, i.e.
        // right == last first, otherwise an end() iterator could be
        // dereferenced.
        switch (state) {
            case STATE_START:
                if (right == last) {
                    MSTK_LOG(logDEBUG3)
                            << "STATE_START: got end";
                    state = STATE_STOP;
                } else {
                    if (slopeDown(left, right, comp)) {
                        MSTK_LOG(logDEBUG3)
                                << "STATE_START: got D";
                        state = STATE_BUMP;
                    } else {
                        MSTK_LOG(logDEBUG3)
                                << "STATE_START: got U or E";
                        state = STATE_RAMP_UP;
                    }
                }
                break;

            case STATE_BUMP:
                // note that peak pos is always at *left
                if (right == last || slopeUp(left, right, comp)) {
                    MSTK_LOG(logDEBUG3)
                            << "STATE_BUMP: got end or U";
                    state = STATE_STOP;
                } else {
                    // slopeDown, slopeEqual
                    MSTK_LOG(logDEBUG3)
                            << "STATE_BUMP: got D or E";
                    state = STATE_RAMP_DOWN;
                }
                break;

            case STATE_RAMP_UP:
                if (right == last) {
                    MSTK_LOG(logDEBUG3)
                            << "STATE_RAMP_UP: got end";
                    state = STATE_STOP;
                } else {
                    if (slopeDown(left, right, comp)) {
                        MSTK_LOG(logDEBUG3)
                                << "STATE_RAMP_UP: got D";
                        state = STATE_BUMP;
                    } else {
                        MSTK_LOG(logDEBUG3)
                                << "STATE_RAMP_UP: got U or E";
                    }
                }
                break;

            case STATE_RAMP_DOWN:
                if (right == last || slopeUp(left, right, comp)) {
                    MSTK_LOG(logDEBUG3)
                            << "STATE_RAMP_UP: got end or U";
                    state = STATE_STOP;
                } else {
                    MSTK_LOG(logDEBUG3)
                            << "STATE_RAMP_UP: got D or E";
                }
                break;

            default:
                mstk_assert(false, "Invalid FSM state.");
                break;
        }
        // Don't advance iterators in last step.
        // (We could also advance and then use 'left', but it's better to
        // make sure that our iterators never go beyond 'last'.
        if (state != STATE_STOP) {
            ++left;
            ++right;
        }
    }
    MSTK_LOG(logDEBUG3)
            << "STATE_STOP";
    return std::make_pair(first, right);
}

} // namespace fe

} // namespace mstk

#endif /* __MSTK_INCLUDE_MSTK_FE_SIMPLEBUMPFINDER_HPP__ */
