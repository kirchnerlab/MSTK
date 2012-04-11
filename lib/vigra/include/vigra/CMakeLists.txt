FILE(GLOB headers
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${headers} DESTINATION include/MSTK/lib/vigra)
