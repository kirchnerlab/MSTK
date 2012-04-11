SET(SRCS
    Error.cpp
)

ADD_LIBRARY(mstk-common ${SRCS})

##############################################################################
# installation
##############################################################################
FILE(GLOB headers 
    "${MSTK_SOURCE_DIR}/include/MSTK/common/*.hpp"
    "${MSTK_BINARY_DIR}/include/MSTK/common/*.hpp")
INSTALL(FILES ${headers} DESTINATION include/MSTK/common)
INSTALL(TARGETS mstk-common
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
