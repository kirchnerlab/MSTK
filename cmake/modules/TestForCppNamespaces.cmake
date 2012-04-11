MACRO(TestForCppNamespaces)
    MESSAGE(STATUS "Testing for C++ namespace support")
    TRY_COMPILE(HAVE_CPP_NAMESPACES
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/cmake/modules/TestForCppNamespaces.cpp
    )
    IF (HAVE_CPP_NAMESPACES)
        MESSAGE(STATUS "Testing for C++ namespace support - available")
    ELSE (HAVE_CPP_NAMESPACES)
        MESSAGE(STATUS "Testing for C++ namespace support - not available")
    ENDIF (HAVE_CPP_NAMESPACES)
ENDMACRO(TestForCppNamespaces)
