##############################################################################
# required libraries
##############################################################################
FIND_PACKAGE(LIBFBI REQUIRED) 
INCLUDE(${LIBFBI_USE_FILE}) 
IF (NOT ( LIBFBI_VERSION_MAJOR EQUAL 1 AND LIBFBI_VERSION_MINOR EQUAL 3))
    MESSAGE(SEND_ERROR "Found LIBFBI ${LIBFBI_VERSION_MAJOR}.${LIBFBI_VERSION_MINOR} but need 1.3")
ENDIF (NOT(LIBFBI_VERSION_MAJOR EQUAL 1 AND LIBFBI_VERSION_MINOR EQUAL 3))

##############################################################################
# build
##############################################################################
# We add all mstk-fe sources here because they should end up in a single
# library and this is easier that way.
SET(SRCS
    CentroidWeightedMeanDisambiguator.cpp
    GaussianMeanAccumulator.cpp
    RunningMeanSmoother.cpp
    SimpleBumpFinder.cpp
    SumAbundanceAccumulator.cpp
    UncenteredCorrelation.cpp
    types/Centroid.cpp
    types/IsotopePattern.cpp
    types/Spectrum.cpp
    types/Xic.cpp
    #types/XicFbiTraits.cpp
)

ADD_LIBRARY(mstk-fe ${SRCS})

##############################################################################
# installation
##############################################################################
# All install headers/targets are collected here as well.
FILE(GLOB mstkFeHeaders 
    "${MSTK_SOURCE_DIR}/include/MSTK/fe/*.hpp"
    "${MSTK_BINARY_DIR}/include/MSTK/fe/*.hpp")
INSTALL(FILES ${mstkFeHeaders} DESTINATION include/MSTK/fe)
FILE(GLOB mstkFeTypesHeaders 
    "${MSTK_SOURCE_DIR}/include/MSTK/fe/types/*.hpp"
    "${MSTK_BINARY_DIR}/include/MSTK/fe/types/*.hpp")
INSTALL(FILES ${mstkFeTypesHeaders} DESTINATION include/MSTK/fe/types)
INSTALL(TARGETS mstk-fe
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

