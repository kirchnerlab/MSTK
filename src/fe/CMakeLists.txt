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

