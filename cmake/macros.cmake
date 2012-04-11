################################################################
# logging_level_to_define()
# Gives the define corresponding to a global logging level.
#
# Copyright (c) 2009 Bernhard Kausler
# For example: The logging level 'INFO' corresponds to the define 'logINFO'.
#
# Parameters:
#   LOGGING_LEVEL   STRING  One of the logging levels:
#                           NO_LOGGING, ERROR, WARNING, INFO, DEBUG, DEBUG1, (...), DEBUG4
#                           input parameter
#
#   DEFINE          STRING  output parameter
################################################################
MACRO(LOGGING_LEVEL_TO_DEFINE LOGGING_LEVEL DEFINE)
    IF(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
        SET(${DEFINE} "mstk::logNO_LOGGING")
    ELSEIF(${LOGGING_LEVEL} STREQUAL "ERROR")
        SET(${DEFINE} "mstk::logERROR")
    ELSEIF(${LOGGING_LEVEL} STREQUAL "WARNING")
        SET(${DEFINE} "mstk::logWARNING")
    ELSEIF(${LOGGING_LEVEL} STREQUAL "INFO")
        SET(${DEFINE} "mstk::logINFO")
    ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG")
        SET(${DEFINE} "mstk::logDEBUG")
    ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG1")
        SET(${DEFINE} "mstk::logDEBUG1")
    ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG2")
        SET(${DEFINE} "mstk::logDEBUG2")
    ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG3")
        SET(${DEFINE} "mstk::logDEBUG3")
    ELSEIF(${LOGGING_LEVEL} STREQUAL "DEBUG4")
        SET(${DEFINE} "mstk::logDEBUG4")
    ELSE(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
        MESSAGE(SEND_ERROR "Unknown LOGGING_LEVEL: ${LOGGING_LEVEL}. Default to INFO.")
        SET(${DEFINE} "mstk::logINFO")
    ENDIF(${LOGGING_LEVEL} STREQUAL "NO_LOGGING")
ENDMACRO(LOGGING_LEVEL_TO_DEFINE)
    


################################################################
# MACRO_ENSURE_OUT_OF_SOURCE_BUILD(<errorMessage>)
#
# Copyright (c) 2006, Alexander Neundorf, <neundorf@kde.org>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
################################################################
macro (MACRO_ENSURE_OUT_OF_SOURCE_BUILD _errorMessage)
   string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" _insource)
   if (_insource)
     message(SEND_ERROR "${_errorMessage}")
     message(FATAL_ERROR "Remove the file CMakeCache.txt in ${CMAKE_SOURCE_DIR} first.")
   endif (_insource)
endmacro (MACRO_ENSURE_OUT_OF_SOURCE_BUILD)

################################################################
# Add an MSTK test.
# 
# Parameters:
#
#   lib  -- Name of the MSTK lib that the test belongs to,
#           e.g. 'common', or 'fe'
#   classname -- Name of the test; generally the name of the class
#           that is being tested.
#   src  -- Source file 
#
# The macro will generate targets called 'lib_test_name' and
# 'lib_memtest_name'.
#
################################################################
MACRO(ADD_MSTK_TEST lib classname src)
    # derive the test file name
    SET(testName "${lib}_${classname}_test") 
    SET(testNameExe "${testName}_exe")

    # build the test
    ADD_EXECUTABLE(${testNameExe} ${src})
    # and link
    TARGET_LINK_LIBRARIES(${testNameExe} ${TEST_LIBS})

    # add test to global list of unit tests
    MESSAGE(STATUS "Adding test for ${lib}/${classname}: ${testName}.")
    ADD_TEST(${testName} ${testNameExe})
    LIST(APPEND MSTK_TEST_NAMES ${testName})
    LIST(APPEND "MSTK_${lib}_TEST_NAMES" ${testName})
    # Add target for the test
    ADD_CUSTOM_TARGET(${testName} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${testNameExe})

    # if we have valgrind, also add memory test
    IF (HAVE_VALGRIND)
        SET(memtestName "${lib}_${classname}_memtest") 
        LIST(APPEND MSTK_MEMTEST_NAMES ${memtestName})
        LIST(APPEND "MSTK_${lib}_MEMTEST_NAMES" ${memtestName})
        # Add a target for the memory test
        #ADD_CUSTOM_TARGET(${memtestName} COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/memcheck.py ${CMAKE_CURRENT_BINARY_DIR}/${exe})
        MESSAGE(STATUS "Adding memory test for ${lib}/${classname}: ${memtestName}.")
        ADD_TEST(${memtestName}
            ${CMAKE_SOURCE_DIR}/tests/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/${testNameExe} ${CMAKE_BINARY_DIR})
        ADD_CUSTOM_TARGET(${memtestName}
            ${CMAKE_SOURCE_DIR}/tests/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/${testNameExe} ${CMAKE_BINARY_DIR})
    ENDIF(HAVE_VALGRIND)
ENDMACRO(ADD_MSTK_TEST lib classname src)
