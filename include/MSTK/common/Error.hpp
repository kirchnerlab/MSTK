/*
 * Error.hpp
 *
 * Copyright (C) 2009-2012 Marc Kirchner
 * Copyright (C) 2009 Bernhard X. Kausler
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
#ifndef __MSTK_INCLUDE_MSTK_COMMON_ERROR_HPP__
#define __MSTK_INCLUDE_MSTK_COMMON_ERROR_HPP__

#include <exception>
#include <MSTK/common/Types.hpp>

#ifdef MSTK_DEBUG
#include <cassert>
#include <MSTK/common/Log.hpp>
#endif

namespace mstk {

/** @addtogroup mstk_common
 * @{
 */

#ifdef MSTK_DEBUG
#define mstk_assert(condition, message) {       \
    if (!(condition)) {                         \
        MSTK_LOG(mstk::logERROR)                \
            << "Assertion failed: "             \
            << #condition << ". " << message;   \
        assert(condition && (message));         \
    }                                           \
}

#else

/** Assertion macro that supports the MSTK logging facility.
 * Use this instead of the usual \c cassert-based call to \c assert().
 * Once, @c message has been logged into the MSTK logging facility,
 * the macro calls the standard C \c assert() function and terminates the
 * program.
 *
 * The macro is only available for non-release builds in which the
 * MSTK_DEBUG preprocessor directive is defined. If \c MSTK_DEBUG is not
 * defined, all calls are transformed into NOPs. Hence, never put any
 * program logic into these calls, ie. never use
 * \code
 *      ...
 *      mstk_assert((result = compute()) >= 0, "Result must be non-negative.")
 *      process(result);
 *      ...
 * \endcode
 * because the call to \c compute() will never be carried out if \c MSTK_DEBUG
 * is not defined. Instead write
 * \code
 *      ...
 *      result = compute();
 *      mstk_assert(result >= 0, "Result must be non-negative.")
 *      process(result);
 *      ...
 * \endcode
 *
 * @param condition The assertion condition that should be checked.
 *                  If the condition evaluates to \c false, then the
 *                  assertion will fire.
 * @param message A message describing the asserted condition. The message
 *                will be logged to the MSTK logging facility.
 */
#define mstk_assert(condition, message)

#endif

/** Base class for all MSTK exceptions.
 */
class Exception : public std::exception
{
public:
    /** Constructor.
     *  @param message String that describes the error condition.
     *                 Contents of the string will be copied into an
     *                 internal buffer; the exception will not take
     *                 ownership of the memory \c message points to.
     */
    explicit Exception(const Char* message);

    /** Constructor.
     *  @param message String describing the error condition.
     */
    explicit Exception(const String& message);

    /** Virtual destructor.
     */
    virtual ~Exception() throw ();

    /** C-string-style access to the error description.
     *  @return A pointer to C-string-style representation of the
     *          error description.
     */
    virtual const Char * what() const throw ();

protected:
    /** Stores the string that describes the error condition.
     * The implementation leaves the variable as 'protected' to
     * simplify reimplementations of 'what()' in derived classes.
     */
    String what_;
};

/**
 *   Base class for all logic errors.
 *
 *   Use LogicError for defects that, in principle, could be detected by static
 *   analysis.  A prominent example are data access attempts outside expected
 *   parameter ranges.
 */
class LogicError : public mstk::Exception
{
public:
    explicit LogicError(const Char * message);
    explicit LogicError(const String& message);
    virtual ~LogicError() throw ();
};

/**
 *   Base class for all runtime errors.
 *
 *   Use RuntimeError for defects, which could only happen or be detected during
 *   runtime.  Runtime errors are caused by not acquirable system resources (like
 *   memory, file handles, network etc.), race conditions of different threads or
 *   processes and other unforseeable failures.
 */
class RuntimeError : public mstk::Exception
{
public:
    explicit RuntimeError(const Char* message);
    explicit RuntimeError(const String& message);
    virtual ~RuntimeError() throw ();
};

class PreconditionViolation : public LogicError
{
public:
    explicit PreconditionViolation(const Char* message);
    explicit PreconditionViolation(const String& message);
    virtual ~PreconditionViolation() throw ();
};

class PostconditionViolation : public LogicError
{
public:
    explicit PostconditionViolation(const Char* message);
    explicit PostconditionViolation(const String& message);
    virtual ~PostconditionViolation() throw ();
};

class InvariantViolation : public LogicError
{
public:
    explicit InvariantViolation(const Char* message);
    explicit InvariantViolation(const String& message);
    virtual ~InvariantViolation() throw ();
};

//////////////////////////////////////////////////////
// some helper functions to throw exceptions easily //
// dont't use them directly, use the macros below   //
//////////////////////////////////////////////////////
inline
void throw_invariant_error(bool predicate, const Char* message)
{
    if (!predicate)
        throw mstk::InvariantViolation(message);
}

inline
void throw_precondition_error(bool predicate, const Char* message)
{
    if (!predicate)
        throw mstk::PreconditionViolation(message);
}

inline
void throw_postcondition_error(bool predicate, const Char* message)
{
    if (!predicate)
        throw mstk::PostconditionViolation(message);
}

inline
void throw_invariant_error(bool predicate, const String& message)
{
    if (!predicate)
        throw mstk::InvariantViolation(message);
}

inline
void throw_precondition_error(bool predicate, const String& message)
{
    if (!predicate)
        throw mstk::PreconditionViolation(message);
}

inline
void throw_postcondition_error(bool predicate, const String& message)
{
    if (!predicate)
        throw mstk::PostconditionViolation(message);
}

/*
 * Convenience macros to write quick throw statements.
 */

/**
 * Throws a mstk::PreconditionViolation, if the PREDICATE is false.
 */
#define mstk_precondition(PREDICATE, MESSAGE) mstk::throw_precondition_error((PREDICATE), MESSAGE)

/**
 * Throws a mstk::PostconditionViolation, if the PREDICATE is false.
 */
#define mstk_postcondition(PREDICATE, MESSAGE) mstk::throw_postcondition_error((PREDICATE), MESSAGE)

/**
 * Throws a mstk::InvariantViolation, if the PREDICATE is false.
 */
#define mstk_invariant(PREDICATE, MESSAGE) mstk::throw_invariant_error((PREDICATE), MESSAGE)

/**
 * Throws a RuntimeError.
 */
#define mstk_fail(MESSAGE) throw mstk::RuntimeError(MESSAGE)

/** @} */

} // namespace mstk

#endif // __MSTK_INCLUDE_MSTK_ERROR_HPP__
