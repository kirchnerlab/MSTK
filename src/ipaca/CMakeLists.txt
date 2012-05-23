##############################################################################
# build
##############################################################################
SET(SRCS
    Mercury7Impl.cpp
    Spectrum.cpp
    Stoichiometry.cpp
    Traits.cpp
)

ADD_LIBRARY(mstk-ipaca ${SRCS})

##############################################################################
# installation
##############################################################################
# All install headers/targets are collected here as well.
FILE(GLOB mstkIpacaHeaders 
    "${MSTK_SOURCE_DIR}/include/MSTK/ipaca/*.hpp"
    "${MSTK_BINARY_DIR}/include/MSTK/ipaca/*.hpp")
INSTALL(FILES ${mstkIpacaHeaders} DESTINATION include/MSTK/ipaca)
INSTALL(TARGETS mstk-ipaca
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

