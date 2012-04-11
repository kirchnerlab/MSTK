#!/usr/bin/env python

import sys, os

if len(sys.argv) != 3:
    print 'Usage: ./generatetest.py <classname> <templatename>'
    sys.exit('ERROR: Wrong parameters.')

classname = sys.argv[1]
templatename = sys.argv[2]
print classname
testfilename = classname+'-test.cpp'
print testfilename
testfile = open(testfilename, 'w')

template = open(templatename)
template = template.read()
template = template.replace('CLASSNAME', classname)
template = template.replace('FILENAME', testfilename)

testfile.write(template)

cmakefile = open('CMakeLists.txt', 'r')
cmake = cmakefile.read()
cmakefile.close()

srcs = 'SET(SRCS_'+classname.upper()+' '+testfilename+')'
print srcs
cmake = cmake.replace('#### Sources', '#### Sources\n'+srcs)


test = 'ADD_MSTK_TEST("'+classname+'" test_'+classname.lower()+' ${SRCS_'+classname.upper()+'})'
print test
cmake = cmake.replace('#### Unit tests', '#### Unit tests\n'+test)

cmakefile = open('CMakeLists.txt', 'w')
cmakefile.write(cmake)
cmakefile.close()
