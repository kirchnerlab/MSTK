# require boost
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREAD OFF)
FIND_PACKAGE(Boost COMPONENTS regex REQUIRED)
# TODO check minimal version we need for boost regex
IF (NOT ( Boost_MAJOR_VERSION EQUAL 1 AND Boost_MINOR_VERSION GREATER 32))
	MESSAGE(SEND_ERROR "libaas requires Boost >= 1.40.0")
ENDIF (NOT(Boost_MAJOR_VERSION EQUAL 1 AND Boost_MINOR_VERSION GREATER 32))
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

#require threads
FIND_PACKAGE(Threads REQUIRED)

SET(SRCS
	AminoAcid.cpp
	RawAminoAcid.cpp
	RawAminoAcidImpl.cpp
	AminoAcidSequence.cpp
	Element.cpp
	ElementImpl.cpp
	Isotope.cpp
	Modification.cpp
	RawModification.cpp
	RawModificationImpl.cpp
	Residue.cpp
	Specificity.cpp
	Stoichiometry.cpp
	StoichiometryConfig.cpp
	StoichiometryConfigImpl.cpp
	tools/Digester.cpp
	tools/FastaReader.cpp
)

ADD_LIBRARY(mstk-aas ${SRCS})
TARGET_LINK_LIBRARIES(mstk-aas
	${Boost_REGEX_LIBRARY}
	${CMAKE_THREAD_LIBS_INIT}
)

##############################################################################
# installation
##############################################################################
FILE(GLOB headers 
    "${MSTK_SOURCE_DIR}/include/MSTK/aas/*.hpp"
    "${MSTK_SOURCE_DIR}/include/MSTK/aas/tools/*.hpp"
    "${MSTK_SOURCE_DIR}/include/MSTK/aas/adapter/*.hpp"
    "${MSTK_BINARY_DIR}/include/MSTK/aas/*.hpp")
INSTALL(FILES ${headers} DESTINATION include/MSTK/aas)
INSTALL(TARGETS mstk-aas
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
