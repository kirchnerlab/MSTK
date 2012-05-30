
MACRO(LIST_CONTAINS var value)
	SET(${var})
	FOREACH (value2 ${ARGN})
		IF (${value} STREQUAL ${value2})
			SET(${var} TRUE)
		ENDIF (${value} STREQUAL ${value2})
	ENDFOREACH (value2)
ENDMACRO(LIST_CONTAINS)

SET(TEST_LIBS mstk-aas mstk-common)

#### Sources
SET(SRCS_FASTAREADER FastaReader-test.cpp)
SET(SRCS_DIGESTER Digester-test.cpp)
SET(SRCS_STOICHIOMETRYCONFIG StoichiometryConfig-test.cpp)
SET(SRCS_RAWAMINOACID RawAminoAcid-test.cpp)
SET(SRCS_RAWMODIFICATION RawModification-test.cpp)
SET(SRCS_AMINOACIDSEQUENCE AminoAcidSequence-test.cpp)
SET(SRCS_ISOTOPE Isotope-test.cpp)
SET(SRCS_SPECIFICITY Specificity-test.cpp)
SET(SRCS_RESIDUE Residue-test.cpp)
SET(SRCS_STOICHIOMETRY Stoichiometry-test.cpp)
SET(SRCS_MODIFICATION Modification-test.cpp)
SET(SRCS_AMINOACID AminoAcid-test.cpp)
SET(SRCS_ELEMENT Element-test.cpp)
SET(SRCS_COMBINATION Combination-test.cpp)

#########  List of tests
ADD_MSTK_TEST("aas" "FastaReader" ${SRCS_FASTAREADER})
ADD_MSTK_TEST("aas" "Digester" ${SRCS_DIGESTER})
ADD_MSTK_TEST("aas" "StoichiometryConfig" ${SRCS_STOICHIOMETRYCONFIG})
ADD_MSTK_TEST("aas" "RawAminoAcid" ${SRCS_RAWAMINOACID})
ADD_MSTK_TEST("aas" "RawModification" ${SRCS_RAWMODIFICATION})
ADD_MSTK_TEST("aas" "AminoAcidSequence" ${SRCS_AMINOACIDSEQUENCE})
ADD_MSTK_TEST("aas" "Isotope" ${SRCS_ISOTOPE})
ADD_MSTK_TEST("aas" "Specificity" ${SRCS_SPECIFICITY})
ADD_MSTK_TEST("aas" "Residue" ${SRCS_RESIDUE})
ADD_MSTK_TEST("aas" "Stoichiometry" ${SRCS_STOICHIOMETRY})
ADD_MSTK_TEST("aas" "Modification" ${SRCS_MODIFICATION})
ADD_MSTK_TEST("aas" "AminoAcid" ${SRCS_AMINOACID})
ADD_MSTK_TEST("aas" "Element" ${SRCS_ELEMENT})

######### iff ipaca is in MSTK_COMPONENTS
LIST_CONTAINS(ipacaInMSTK ipaca ${MSTK_COMPONENTS})
IF (ipacaInMSTK)
    MESSAGE(STATUS "Adding aas<->ipaca test (Combination)")
    SET(TEST_LIBS mstk-aas mstk-common mstk-ipaca)
    ADD_MSTK_TEST("aas" "Combination" ${SRCS_COMBINATION})
ENDIF (ipacaInMSTK)

MESSAGE(STATUS "Tests for 'aas': ${MSTK_aas_TEST_NAMES}")
MESSAGE(STATUS "Memory tests for 'aas': ${MSTK_aas_MEMTEST_NAMES}")

ADD_CUSTOM_TARGET(aas_test 
    DEPENDS ${MSTK_aas_TEST_NAMES} ${MSTK_aas_MEMTEST_NAMES}
)

