##############################################################################
# variables
##############################################################################
SET(VIGRA_INCLUDE_DIR "/usr/local/include" CACHE STRING "Path to VIGRA headers?")
INCLUDE_DIRECTORIES(${VIGRA_INCLUDE_DIR})

##############################################################################
# build
##############################################################################
SET(SRCS 
    BoxPeakShape.cpp
    ConstantModel.cpp
    GaussianPeakShape.cpp
    LinearSqrtModel.cpp
    LorentzianPeakShape.cpp
    PeakShapeFunction.cpp
    QuadraticModel.cpp
    SqrtModel.cpp
)

ADD_LIBRARY(mstk-psf ${SRCS})

##############################################################################
# installation
##############################################################################
# All install headers/targets are collected here as well.
FILE(GLOB mstkPsfHeaders 
    "${MSTK_SOURCE_DIR}/include/MSTK/psf/*.hpp"
    "${MSTK_BINARY_DIR}/include/MSTK/psf/*.hpp")
INSTALL(FILES ${mstkPsfHeaders} DESTINATION include/MSTK/psf)
FILE(GLOB mstkPsfTypesHeaders 
    "${MSTK_SOURCE_DIR}/include/MSTK/psf/types/*.hpp"
    "${MSTK_BINARY_DIR}/include/MSTK/psf/types/*.hpp")
INSTALL(FILES ${mstkPsfTypesHeaders} DESTINATION include/MSTK/psf/types)
INSTALL(TARGETS mstk-psf
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

