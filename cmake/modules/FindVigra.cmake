# Copyright : ETH Zurich
# License   : three-clause BSD license
# Authors   : Witold Wolski
# for full text refer to files: LICENSE, AUTHORS and COPYRIGHT

FIND_PATH(VIGRA_INCLUDE_DIR vigra/accessor.hxx)
FIND_LIBRARY(VIGRA_LIBRARY NAMES vigraimpex)

IF (VIGRA_INCLUDE_DIR AND VIGRA_LIBRARY)
   MESSAGE(STATUS "XXXXXXXXXXXXX ${VIGRA_LIBRARY} XXXXXXXXXXXXXXX")
   SET(VIGRA_FOUND TRUE)
ENDIF (VIGRA_INCLUDE_DIR AND VIGRA_LIBRARY)

IF (VIGRA_FOUND)
   # show which CppUnit was found only if not quiet
   IF (NOT VIGRA_FIND_QUIETLY)
      MESSAGE(STATUS "Found Pwiz: ${VIGRA_LIBRARY}")
   ENDIF (NOT VIGRA_FIND_QUIETLY)
ELSE (VIGRA_FOUND)
   # fatal error if CppUnit is required but not found
   IF (VIGRA_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find Pwiz")
   ENDIF (VIGRA_FIND_REQUIRED)
ENDIF (VIGRA_FOUND)

