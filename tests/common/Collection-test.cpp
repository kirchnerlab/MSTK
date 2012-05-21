/*
 * Collection-test.cpp
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
#include "unittest.hxx"
#include <iostream>
#include <vector>
// expose the class
#define private public
#define protected public
#include <MSTK/common/Collection.hpp>
#undef protected
#undef private


struct CollectionTestSuite : vigra::test_suite {
    CollectionTestSuite() : vigra::test_suite("Collection") {
        add( testCase(&CollectionTestSuite::testConstructors));
        add( testCase(&CollectionTestSuite::testOperators));
        add( testCase(&CollectionTestSuite::testAssignment));
        add( testCase(&CollectionTestSuite::testStackOperations));
        add( testCase(&CollectionTestSuite::testListOperations));
        add( testCase(&CollectionTestSuite::testIterators));
        add( testCase(&CollectionTestSuite::testElementAccess));
        add( testCase(&CollectionTestSuite::testSizeAndCapacity));
        add( testCase(&CollectionTestSuite::testOther));
    }

    void testConstructors() {
        // default construction
        mstk::Collection<int> c;
        
        // construction from vector
        std::vector<int, std::allocator<int> > vector;
        mstk::Collection<int, std::allocator<int> > cFromVector = mstk::Collection<int, std::allocator<int> >(vector);
        should(vector == cFromVector.c_);

        // construction with preinit
        mstk::Collection<int> cPreinit(5);
        should(cPreinit.c_.size() == 5);

        // construction from sequence
        int myints[] = {1,2,3,4};
        mstk::Collection<int> cFromSeq (myints, myints + sizeof(myints) / sizeof(int) );

        should(cFromSeq.c_.size() == 4);
        int i = 1;
        for(mstk::Collection<int>::iterator it = cFromSeq.begin(); it < cFromSeq.end(); ++it) {           
            should(*it == i);
            ++i;         
        }

        // copy constructor        
        mstk::Collection<int> c_copy = c; // This wouldn't compile in case of an explicit copy constructor (which we don't want!)       
        should(c == c_copy);
    }

    
    void testOperators() {
        // operator=
        mstk::Collection<int> c1(2);
        mstk::Collection<int> c2(3);
        should(!(c1 == c2));
        c1 = c2;
        should(c1 == c2);

        // operator==
        mstk::Collection<int> c3(2);
        mstk::Collection<int> c4(3);
        should(!(c3 == c4));
        c3.c_ = c4.c_;
        should(c3 == c4);

        // operator<
        mstk::Collection<int> c5;
        mstk::Collection<int> c6;
        should((c5.c_ < c6.c_) == (c5 < c6));

        // swap operator
        mstk::Collection<int> c7(2);
        mstk::Collection<int> c7_copy = c7;        
        mstk::Collection<int> c8(3);
        mstk::Collection<int> c8_copy = c8;

        should(!((c7_copy == c8) && (c8_copy == c7)));        
        mstk::swap(c7, c8);
        should((c7_copy == c8) && (c8_copy == c7));
    }

    
    void testAssignment() {
        // assign sequence
        mstk::Collection<int> c(10);
        should(c.size() == 10);

        int myints[] = {1,2,3,4};
        c.assign(myints, myints + sizeof(myints) / sizeof(int) );

        should(c.size() == 4);
        int i = 1;
        for(mstk::Collection<int>::iterator it = c.begin(); it < c.end(); ++it) {           
            should(*it == i);
            ++i;       
        }

        // assign values
        mstk::Collection<int> c2(10);
        should(c2.size() == 10);

        // cast to size_type to prevent confusion with the other insert signature
        c2.assign((mstk::Collection<int>::size_type)4, 2317);
        should(c2.size() == 4);

        for(mstk::Collection<int>::iterator it = c2.begin(); it < c2.end(); ++it) {
            should(*it == 2317);
        }
    }

    
    void testStackOperations() {
        // push_back
        mstk::Collection<int> c;
        should(c.size() == 0);
        
        c.push_back(17);
        should(c.size() == 1);
        should(c[0] == 17);
    
        // pop_back
        int myints[] = {1,2,3,4};
        c.assign(myints, myints + sizeof(myints) / sizeof(int) );
        should(c[3] == 4);
        
        c.pop_back();
        should(c.size() == 3);
        should(c[2] == 3);
    }


    void testListOperations() {
        mstk::Collection<int> c;
        int myints[] = {46,243,45};
        c.assign(myints, myints + sizeof(myints) / sizeof(int) );

        // insert
        should(!(c[1] == 23)); 
        mstk::Collection<int>::iterator it_insert = c.insert(c.begin()+1, 23);
        should(c.size() == 4);
        should(c[1] == 23);
        should(*it_insert == 23);
        should((c.begin() + 1) == it_insert);

        // multiple insert
        // c is now {46,23,243,45}
        should(c.size() == 4);
        // avoid confusion with the templated insert signature by explicitly casting to size_type
        c.insert(c.begin()+1, (mstk::Collection<int>::size_type)2, 123); // 2x 123 at position
        should(c.size() == 6);
        should(c[1] == 123 && c[2] == 123);

        // c is now {46,123,123,23,243,45}
        // templated insert
        int toBeInserted[] = {823, 329, 198};
        c.insert(c.begin()+1, toBeInserted, toBeInserted + 3);
        should(c.size() == 9);
        should(c[1] == 823 && c[2] == 329 && c[3] == 198);

        // c is now {46,823,329,198,123,123,23,243,45}
        // erase
        mstk::Collection<int>::iterator it_erase = c.erase(c.begin()+1);
        should(c.size() == 8);
        should(*it_erase == 329);
        should(c.begin()+1 == it_erase);        

        // c is now {46,329,198,123,123,23,243,45}
        // erase sequence
        mstk::Collection<int>::iterator it_erase_seq = c.erase(c.begin()+1,c.begin()+4);
        should(c.size() == 5);
        should(*it_erase_seq == 123);
        should(c.begin()+1 == it_erase_seq);

        // clear
        should(c.size() == 5);
        should(!(c.empty()));
        c.clear();
        should(c.size() == 0);
        should(c.empty());
    }

    
    void testIterators() {
        mstk::Collection<int> c(5);

        // begin, rbegin, end and rend
        should(c.begin() == c.c_.begin());
        should(c.rbegin() == c.c_.rbegin());
        should(c.end() == c.c_.end());
        should(c.rend() == c.c_.rend());
        should(c.begin() == c.c_.begin());

        // const begin and end
        const mstk::Collection<int> c_const(5);
        should(c_const.begin() == c_const.c_.begin());
        should(c_const.end() == c_const.c_.end());
    }


    void testElementAccess() {
        mstk::Collection<int> c;
        int myints[] = {46,243,45};
        c.assign(myints, myints + sizeof(myints) / sizeof(int) );

        // operator[]
        should(c[0] == 46 && c[1] == 243 && c[2] == 45);
        c[1] = 127;
        should(*(c.begin() + 1) ==  127);
        
        // at
        should(c.at(0) == 46 && c.at(1) == 127 && c.at(2) == 45);
        c.at(1) = 259;
        should(*(c.begin() + 1) ==  259);

        
        int myints_for_const[] = {46,243,45};
        const mstk::Collection<int> c_const(myints_for_const, myints_for_const + sizeof(myints_for_const) / sizeof(int) );

        // const operator[]
        should(c_const[0] == 46 && c_const[1] == 243 && c_const[2] == 45);

        // const at
        should(c_const.at(0) == 46 && c_const.at(1) == 243 && c_const.at(2) == 45);
    }

    
    void testSizeAndCapacity() {
        // size
        int myints[] = {26,45,12};
        mstk::Collection<int> c(myints, myints + sizeof(myints) / sizeof(int));
        should(c.size() == 3);
        c.push_back(32);
        should(c.size() == 4);
        c.pop_back();
        should(c.size() == 3);
        c.pop_back();
        should(c.size() == 2);
        c.clear();
        should(c.size() == 0);
        should(c.empty());
        
        // maxSize
        mstk::Collection<int> c_maxSize(4);
        should(c_maxSize.max_size() == c_maxSize.c_.max_size());

        // empty
        mstk::Collection<int> c_empty(4);
        should(!(c_empty.empty()));
        c_empty.clear();
        should(c_empty.empty());

        // resize
        mstk::Collection<int> c_resize(5);
        should(c_resize.size() == 5);
        c_resize.resize(10, 13);
        should(c_resize.size() == 10);
        should(c_resize.at(9) == 13);
        c_resize.resize(2);
        should(c_resize.size() == 2);
        should(c_resize.at(1) == 0);

        // capacity
        mstk::Collection<int> c_capacity(5);
        should(c_capacity.capacity() == c_capacity.c_.capacity());

        // reserve
        mstk::Collection<int> c_reserve;
        c_reserve.reserve((mstk::Collection<int>::size_type)1000000);
        should(c_reserve.capacity() >= 1000000);
    }

    
    void testOther() {
        // swap
        mstk::Collection<int> c1(2);
        mstk::Collection<int> c1_copy = c1;        
        mstk::Collection<int> c2(3);
        mstk::Collection<int> c2_copy = c2;

        should(!((c1_copy == c2) && (c2_copy == c1)));        
        c1.swap(c2);
        should((c1_copy == c2) && (c2_copy == c1));


        // get_allocator
        mstk::Collection<int> c_alloc;
        should(c_alloc.get_allocator() == c_alloc.c_.get_allocator());
    }
};

int main()
{
    CollectionTestSuite test;
    int success = test.run();
    std::cout << test.report() << std::endl;
    return success;
}


