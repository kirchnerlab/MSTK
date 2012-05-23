PROJECT(MSTK)
#############################################################################
# cmake setup
#############################################################################
# require 2.8
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
# set macro search path
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/macros.cmake)
# include module path
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

#############################################################################
# versioning
#############################################################################
INCLUDE(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC MSTK_VERSION)

#############################################################################
# data diectories
#############################################################################
SET(DATA_INSTALL_DIR ${CMAKE_BINARY_DIR/data})

#############################################################################
# user flags
#############################################################################
# Which libs should be built into MSTK? Add new libs using a semicolon
# as a separator: "common;fe;ngq" etc.
SET(MSTK_COMPONENTS "common;fe;ipaca" CACHE STRING "Which libraries should be built?")
# build type: Debug, RelWithDebInfo, and Release
SET(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "MSTK build type.")
# enable testing: this will also attempt to build memory tests
OPTION(ENABLE_TESTING "Enable regression/unit tests." ON)
# gcov support
OPTION(ENABLE_COVERAGE "Enable GCov coverage analysis (defines a 'coverage' target and enforces static build of MSTK)" OFF)
IF(ENABLE_COVERAGE AND NOT ENABLE_TESTING)
    MESSAGE(FATAL_ERROR "Coverage analysis requires the unit tests.")
ENDIF(ENABLE_COVERAGE AND NOT ENABLE_TESTING)
IF(ENABLE_COVERAGE AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    MESSAGE(FATAL_ERROR "Can not enable coverage analysis in release mode.")
ENDIF(ENABLE_COVERAGE AND(CMAKE_BUILD_TYPE STREQUAL "Release"))

#############################################################################
# build type
#############################################################################
# build type-dependent compiler switches
IF(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    ADD_DEFINITIONS(-DMSTK_DEBUG)
ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
IF(MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS -DEXP_STL)
ELSE(MSVC)
    ADD_DEFINITIONS(-Wall)
    SET(CMAKE_CXX_FLAGS_RELENABLEDEBINFO "-O3 -g -std=c++0x")
    SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -std=c++0x")
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -std=c++0x")
    SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -Wall -ggdb3 -std=c++0x")
    IF(ENABLE_COVERAGE)
        IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
            SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
        ELSE(CMAKE_BUILD_TYPE STREQUAL "Debug")
            MESSAGE(FATAL_ERROR "Coverage analysis requires a DEBUG build")
        ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    ENDIF(ENABLE_COVERAGE)
ENDIF(MSVC)
# build static or shared libraries
OPTION(BUILD_SHARED_LIBS "Build MSTK and helper libraries as shared libs (default=OFF)." OFF)

#############################################################################
# global logging level
#############################################################################
SET(LOGGING_LEVEL "INFO" CACHE STRING "Choose a global logging level: NO_LOGGING, ERROR, WARNING, INFO, DEBUG, DEBUG1, ..., DEBUG4")
LOGGING_LEVEL_TO_DEFINE(LOGGING_LEVEL LOG_DEFINE)
ADD_DEFINITIONS(-D FILELOG_MAX_LEVEL=${LOG_DEFINE})

#############################################################################
# require out-of-source build
#############################################################################
MACRO_ENSURE_OUT_OF_SOURCE_BUILD("MSTK requires out-of-source builds.")

#############################################################################
# config files
#############################################################################
CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/MSTK/config.hpp.cmake 
    ${MSTK_BINARY_DIR}/include/MSTK/config.hpp
)
CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/config.cpp.cmake 
    ${MSTK_BINARY_DIR}/src/config.cpp
)

#############################################################################
# header locations
#############################################################################
INCLUDE_DIRECTORIES(
    ${MSTK_BINARY_DIR}
    ${MSTK_BINARY_DIR}/include
    ${MSTK_SOURCE_DIR}
    ${MSTK_SOURCE_DIR}/src
    ${MSTK_SOURCE_DIR}/include
)

#############################################################################
# build: libraries
#############################################################################
# linker paths
#LINK_DIRECTORIES(
#    # add include dirs of libs we depend on here
#)

# require boost
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_USE_MULTITHREADED OFF)
#FIND_PACKAGE(Boost COMPONENTS filesystem program_options regex system REQUIRED)
FIND_PACKAGE(Boost REQUIRED)
#IF (NOT ( Boost_MAJOR_VERSION EQUAL 1 AND Boost_MINOR_VERSION EQUAL 47))
#    MESSAGE(SEND_ERROR "Found Boost ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION} but need 1.47.0")
#ENDIF (NOT(Boost_MAJOR_VERSION EQUAL 1 AND Boost_MINOR_VERSION EQUAL 47))
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

#############################################################################
# build: subdirectories
#############################################################################
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(lib)

#############################################################################
# custom targets
#############################################################################
# add uninstall at some point
#ADD_CUSTOM_TARGET(uninstall
#    "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
#)

#############################################################################
# regression tests
#############################################################################
IF (ENABLE_TESTING)
    INCLUDE(CTest)
    ENABLE_TESTING()
    ADD_SUBDIRECTORY(tests)
ENDIF(ENABLE_TESTING)

#############################################################################
# documentation
#############################################################################
CONFIGURE_FILE(${MSTK_SOURCE_DIR}/doc/Doxyfile.cmake
    ${MSTK_BINARY_DIR}/doc/Doxyfile
    @ONLY IMMEDIATE
)
FIND_PACKAGE(Doxygen)
# add w/ the ALL option to make sure the docs are built for the install
ADD_CUSTOM_TARGET(doc ALL
    COMMAND ${DOXYGEN} ${MSTK_BINARY_DIR}/doc/Doxyfile)
INSTALL(DIRECTORY ${MSTK_BINARY_DIR}/doc/html
    DESTINATION doc
    COMPONENT doc
)

#############################################################################
# summary
#############################################################################
MESSAGE(STATUS)
MESSAGE(STATUS "-----------------------------------------")
MESSAGE(STATUS "System version: ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_VERSION}")
MESSAGE(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
MESSAGE(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "MSTK Version: ${MSTK_VERSION}")
MESSAGE(STATUS "Global logging level: ${LOGGING_LEVEL}")
MESSAGE(STATUS "Regression tests: ${ENABLE_TESTING}")
MESSAGE(STATUS "Coverage analysis: ${ENABLE_COVERAGE}")
MESSAGE(STATUS "Boost version: ${Boost_VERSION}=${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost.SUBMINOR_VERSION}")
MESSAGE(STATUS "Boost include dir: ${Boost_INCLUDE_DIRS}")
MESSAGE(STATUS "Boost library dir:  ${Boost_LIBRARY_DIRS}")
MESSAGE(STATUS "-----------------------------------------")
MESSAGE(STATUS)
