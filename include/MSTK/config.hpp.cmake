/*
 * config.hpp
 *
 * Copyright (c) 2010 Marc Kirchner <mail@marc-kirchner.de>
 *
 */

#ifndef __MSTK_CONFIG_HPP__
#define __MSTK_CONFIG_HPP__

#ifdef _MSC_VER
	#include "winsock2.h"
	#include "vigra/windows.h"
#endif

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#cmakedefine HAVE_UNIX_ISNAN
#cmakedefine HAVE_UNIX_ISINF
#cmakedefine ENABLE_TESTING

#ifdef _WIN32
	#define MSTK_EXPORT __declspec( dllexport )
	/* Disable a template related MSVC warning.
	   See: http://www.unknownroad.com/rtfm/VisualStudio/warningC4251.html */
	#pragma warning( disable: 4251 )
#else
	#define MSTK_EXPORT
#endif

/* Use this to avoid unsued variable warnings */
#define MSTK_UNUSED(x) (void)x;

/* Direcory with all the data */
#define MSTK_DATA_DIR "@DATA_INSTALL_DIR@"

// version info
namespace mstk {

extern const char MSTK_VERSION[];

} // namespace mstk

#endif

