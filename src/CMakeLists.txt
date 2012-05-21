FOREACH(mstkComponent ${MSTK_COMPONENTS})
    ADD_SUBDIRECTORY(${mstkComponent})
    SET(MSTK_DEPENDENCIES ${MSTK_DEPENDENCIES} mstk-${mstkComponent})
    SET(MSTK_LIBRARIES ${MSTK_LIBRARIES} mstk-${mstkComponent})
ENDFOREACH(mstkComponent)

SET(SRCS
    ${CMAKE_CURRENT_BINARY_DIR}/config.cpp
)

IF(NOT ENABLE_COVERAGE)
    ADD_LIBRARY(mstk ${SRCS})
ELSE(NOT ENABLE_COVERAGE)
    ADD_LIBRARY(mstk STATIC ${SRCS})
ENDIF(NOT ENABLE_COVERAGE)

ADD_DEPENDENCIES(mstk ${MSTK_DEPENDENCIES})
TARGET_LINK_LIBRARIES(mstk ${MSTK_LIBRARIES})

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
