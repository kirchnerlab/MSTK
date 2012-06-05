
FOREACH(mstkComponent ${MSTK_COMPONENTS})
    IF(IS_DIRECTORY "${MSTK_SOURCE_DIR}/examples/${mstkComponent}")
        ADD_SUBDIRECTORY(${mstkComponent})
    ELSE(IS_DIRECTORY "${MSTK_SOURCE_DIR}/examples/${mstkComponent}")
        MESSAGE(STATUS "No examples for MSTK/${mstkComponent}.")
    ENDIF(IS_DIRECTORY "${MSTK_SOURCE_DIR}/examples/${mstkComponent}")
ENDFOREACH(mstkComponent)

