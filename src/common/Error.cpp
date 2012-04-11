/*
 * Error.cpp
 *
 * Copyright (C) 2009-2012 Marc Kirchner
 * Copyright (C) 2009 Bernhard X. Kausler
 * 
 * This file is part of the Mass Spectrometry Toolkit (MSTK).
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include <MSTK/common/Error.hpp>
#include <MSTK/common/Log.hpp>
#include <MSTK/common/Types.hpp>

namespace mstk {

Exception::Exception(const Char* message)
{
    what_ = message;
}

Exception::Exception(const String& message)
{
    what_ = message;
}

Exception::~Exception() throw ()
{
}

const Char * Exception::what() const throw ()
{
    return what_.c_str();
}

LogicError::LogicError(const Char * message) :
    Exception(message)
{
}
LogicError::LogicError(const String& message) :
    Exception(message)
{
}

LogicError::~LogicError() throw ()
{
}

RuntimeError::RuntimeError(const Char* message) :
    Exception(message)
{
}
RuntimeError::RuntimeError(const String& message) :
    Exception(message)
{
}

RuntimeError::~RuntimeError() throw ()
{
}

PreconditionViolation::PreconditionViolation(const Char* message) :
    LogicError(message)
{
}
PreconditionViolation::PreconditionViolation(const String& message) :
    LogicError(message)
{
}
PreconditionViolation::~PreconditionViolation() throw ()
{
}

PostconditionViolation::PostconditionViolation(const Char* message) :
    LogicError(message)
{
}
PostconditionViolation::PostconditionViolation(const String& message) :
    LogicError(message)
{
}

PostconditionViolation::~PostconditionViolation() throw ()
{
}

InvariantViolation::InvariantViolation(const Char* message) :
    LogicError(message)
{
}
InvariantViolation::InvariantViolation(const String& message) :
    LogicError(message)
{
}
InvariantViolation::~InvariantViolation() throw ()
{
}

} // namespace mstk
