FOREACH(mstkComponent ${MSTK_COMPONENTS})
    ADD_SUBDIRECTORY(${mstkComponent})
ENDFOREACH(mstkComponent)

SET(SRCS
    ${CMAKE_CURRENT_BINARY_DIR}/config.cpp
)

IF(NOT ENABLE_COVERAGE)
    ADD_LIBRARY(mstk ${SRCS})
ELSE(NOT ENABLE_COVERAGE)
    ADD_LIBRARY(mstk STATIC ${SRCS})
ENDIF(NOT ENABLE_COVERAGE)

ADD_DEPENDENCIES(mstk mstk-common)
TARGET_LINK_LIBRARIES(mstk mstk-common)

##############################################################################
# installation
##############################################################################
FILE(GLOB headers 
    "${MSTK_SOURCE_DIR}/include/MSTK/*.hpp"
    "${MSTK_BINARY_DIR}/include/MSTK/*.hpp")
INSTALL(FILES ${headers} DESTINATION include/MSTK)
INSTALL(TARGETS mstk
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
)
