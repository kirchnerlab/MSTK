/* 
 * Error-test.cpp
 *
 * Copyright (c) 2012 Marc Kirchner
 * Copyright (c) 2009 Bernhard Kausler 
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

#include <cstring>
#include <iostream>
#include <ostream>

#include "unittest.hxx"
#include "MSTK/common/Error.hpp"

using namespace mstk;

struct ErrorTestSuite : vigra::test_suite {
    ErrorTestSuite() : vigra::test_suite("Error reporting") {
        add( testCase(&ErrorTestSuite::testExceptions));
        add( testCase(&ErrorTestSuite::testHelperFunctions));
        add( testCase(&ErrorTestSuite::testErrorMacros));
    }

    void testExceptions() {
        Exception e1("");
        Exception e2("ex123");
        LogicError le1("");
        LogicError le2("le123");
        RuntimeError re1("");
        RuntimeError re2("re123");
        PreconditionViolation pv1("");
        PreconditionViolation pv2("pv123");
        PreconditionViolation pv3(std::string("pv123"));
        PostconditionViolation pov1("");
        PostconditionViolation pov2("pov123");
        InvariantViolation iv1("");
        InvariantViolation iv2("iv123");

        should(std::strcmp(e1.what(), "") == 0);
        should(std::strcmp(e2.what(), "ex123") == 0);
        should(std::strcmp(le1.what(), "") == 0);
        should(std::strcmp(le2.what(), "le123") == 0);
        should(std::strcmp(re1.what(), "") == 0);
        should(std::strcmp(re2.what(), "re123") == 0);
        should(std::strcmp(pv1.what(), "") == 0);
        should(std::strcmp(pv2.what(), "pv123") == 0);
        should(std::strcmp(pv3.what(), "pv123") == 0);
        should(std::strcmp(pov1.what(), "") == 0);
        should(std::strcmp(pov2.what(), "pov123") == 0);
        should(std::strcmp(iv1.what(), "") == 0);
        should(std::strcmp(iv2.what(), "iv123") == 0);
    }

    void testHelperFunctions() {
        // these shouldn't throw anything
        try {
            throw_invariant_error(true, "");
            throw_precondition_error(true, "");
            throw_postcondition_error(true, "");
        }
        catch (...) {
            failTest("Helper function throws unknown exception");
        }

        bool thrown = false;
        try {
            throw_invariant_error(false, "");
        }
        catch (const InvariantViolation &e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        if (!thrown) failTest("InvariantViolation not thrown");

        thrown = false;
        try {
            throw_precondition_error(false, "");
        }
        catch (const PreconditionViolation &e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        if (!thrown) failTest("PreconditionViolation not thrown");

        thrown = false;
        try {
            throw_postcondition_error(false, "");
        }
        catch (const PostconditionViolation &e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        if (!thrown) failTest("PostconditionViolation not thrown");
    }

    void testErrorMacros() {
        // shouldn't throw exceptions
        try {
            mstk_precondition(true, "");
            mstk_postcondition(true, "");
            mstk_invariant(true, "");
        }
        catch (...) {
            failTest("Macro throws unknown exception");
        }

        bool thrown = false;
        try {
            mstk_precondition(false, "");
        }
        catch (const PreconditionViolation &e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        if (!thrown) failTest("mstk_precondition not throwing");

        thrown = false;
        try {
            mstk_postcondition(false, "");
        }
        catch (const PostconditionViolation &e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        if (!thrown) failTest("mstk_postcondition not throwing");

        thrown = false;
        try {
            mstk_invariant(false, "");
        }
        catch (const InvariantViolation &e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        if (!thrown) failTest("mstk_invariant not throwing");

        thrown = false;
        try {
            mstk_fail("");
        }
        catch (const RuntimeError &e) {
			MSTK_UNUSED(e);
            thrown = true;
        }
        if (!thrown) failTest("mstk_fail not throwing");
    }
};

int main()
{
    ErrorTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


