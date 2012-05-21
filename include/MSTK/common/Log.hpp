/*
 * Log.hpp
 *
 * Copyright (C) 2012 Marc Kirchner
 * Copyright (C) 2009 Bernhard X. Kausler
 * 
 * This file is part of the Mass Spectrometry Toolkit (MSTK).
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef __MSTK_INCLUDE_MSTK_LOG_HPP__
#define __MSTK_INCLUDE_MSTK_LOG_HPP__

#include <ctime>
#include <cstdio>

#include <string>
#include <sstream>
#include <iostream>

/**
 *
 * @page common_log Logging
 *
 * @section common_log_usage Logging: Usage
 * @c MSTK provides severity-based multi-level logging to @c stderr via the @c MSTK_LOG macro:
 * @code
 * #include <MSTK/common/Log.hpp>
 * ...
 *     std::string username = getUser();
 *     MSTK_LOG(logINFO) << "Hello" << username << "std::endl is appended automatically";
 *     MSTK_LOG(logDEBUG) << "username has length " << username.size();
 * ...
 * @endcode
 *
 * You can set a logging level:
 * @code
 * #undef MSTK_LOG_MAX_LEVEL
 * #define MSTK_LOG_MAX_LEVEL mstk::logWARNING
 * @endcode
 *
 * In the simplest case, this can be used to define a global logging level; in
 * more complicated setups, this enables fine-grained logging configuration
 * across different compilation units, and even within a single file.
 *
 * The @c MSTK_LOG macro guanrantees that the respective logging code will only
 * be compiled into the final binary if the specified logging level is above
 * the global logging level. Consequently, clients should not use the @c Log API
 * directly but should instead make use of @c MSTK_LOG.
 *
 * Valid logging levels, sorted by increasing level of detail, are
 * @code
 * logERROR, logWARNING, logINFO, logDEBUG, logDEBUG1, logDEBUG2, logDEBUG3, logDEBUG4
 * @endcode
 *
 * @section common_log_disablelogging Logging: Turning off logging
 * MSTK also supports turning logging off completely. For this case, a special
 * logging level is provided. Using
 * @code
 * #undef MSTK_LOG_MAX_LEVEL
 * #define MSTK_LOG_MAX_LEVEL mstk::logNO_LOGGING
 * @endcode
 * will turn off all logging.
 *
 * @note Don't use logNO_LOGGING as a paramter to the logging macro.
 */

namespace mstk {

/** @addtogroup mstk_common
 * @{
 */

// nowTime()
/**
 * The current time as a string.
 *
 * The exact presentation of the string is locale dependent, but should be the following
 * in most cases:
 * @code
 * 15:50:57.979
 * @endcode
 * The precision is up to milliseconds.
 */
inline std::string nowTime();

// LogLevel
/**
 * The logging levels, which can be used.
 *
 * Assign one of these logging levels to each message.
 * Don't use the highest level 'logNO_LOGGING'. It is needed to turn off logging and would
 * default to logINFO if used.
 */
enum LogLevel
{
    logNO_LOGGING,
    logERROR,
    logWARNING,
    logINFO,
    logDEBUG,
    logDEBUG1,
    logDEBUG2,
    logDEBUG3,
    logDEBUG4
};

// Output2FILE
/**
 * Redirector of the logging stream to a file handle.
 *
 * At the moment, the file handle is hardcoded to 'stderr'.
 *
 * Use it in conjunction with the Log<T> class: Log<Output2FILE>.
 *
 */
class Output2FILE
{
public:
    // getRedirect()
    /**
     * The file handle to which the logging stream is redirected.
     */
    static FILE*& getRedirect();

    // output()
    /**
     * Writes message to file handle.
     *
     * This function is used in the Log<T> and mandatory for every Redirector.
     */
    static void output(const std::string& msg);
};

// getRedirect()
inline FILE*& mstk::Output2FILE::getRedirect()
{
    static FILE* pStream = stderr;
    return pStream;
}

// output()
inline void mstk::Output2FILE::output(const std::string& msg)
{
    FILE* pStream = getRedirect();
    if (!pStream) {
        return;
    }
    fprintf(pStream, "%s", msg.c_str());
    fflush(pStream);
}

//Log<T>
/**
 * A thread-safe logging tool
 *
 * The Log<T> passes its internal logging stream to a Redirector T.
 *
 * You have to provide a Redirector with the following interface:
 * @code
 * void T::output(const std::string& msg)
 * @endcode
 *
 */
template<typename T>
class Log
{
public:
    Log();
    virtual ~Log();

    // get()
    /**
     * The internal logging stream.
     *
     * Calling this function writes something like the following into the internal logging
     * stream (depending on the chosen logging level and the current time):
     * @code
     * '- 16:17:23.714 WARNING: '
     * @endcode
     * Afterwards, it returns a reference to that internal stream.
     */
    std::ostringstream& get(LogLevel level = logINFO);

    // getReportingLevel
    /**
     * Returns the deepest logging level available.
     *
     * This function can be used as a safeguard in the logging macros to defend against illegal
     * user defined global logging levels.
     *
     * @return The deepest mstk::LogLevel available.
     */
    static LogLevel& getReportingLevel();

    // toString()
    /**
     * Converts enumerated logging level to a string.
     *
     * For example: logINFO to "INFO".
     */
    static std::string toString(LogLevel level);

    // fromString()
    /**
     * Converts from string to a enumerated logging level
     *
     * For example: "INFO" to logINFO.
     */
    static LogLevel fromString(const std::string& level);

protected:
    // os_
    /**
     * Internal string stream used for logging messages.
     *
     * In the destructor, this stream is written to the Redirector T.
     */
    std::ostringstream os_;

private:
    // Declare copy constructor etc. as private, since we don't want them to be used.
    Log(const Log&);
    Log& operator =(const Log&);
};

// Log()
template<typename T>
mstk::Log<T>::Log()
{
}

// get()
template<typename T>
std::ostringstream& mstk::Log<T>::get(LogLevel level)
{
    // check for valid logging level
    LogLevel ll = level;
    if (ll <= logNO_LOGGING || ll > getReportingLevel()) {
        mstk::Log<T>().get(mstk::logWARNING)
                << "Log<T>::get(): Invalid logging level '" << ll
                << "'. Using INFO level as default.";
        ll = logINFO;
    }

    // print standard logging preambel to logging stream
    os_ << "- " << mstk::nowTime();
    os_ << " " << toString(ll) << ": ";
    os_ << std::string(ll > mstk::logDEBUG ? ll - mstk::logDEBUG : 0, '\t');

    return os_;
}

// ~Log()
template<typename T>
mstk::Log<T>::~Log()
{
    os_ << std::endl;
    T::output(os_.str());
}

// getReportingLevel()
template<typename T>
mstk::LogLevel& mstk::Log<T>::getReportingLevel()
{
    static mstk::LogLevel reportingLevel = mstk::logDEBUG4;
    return reportingLevel;
}

// toString()
template<typename T>
std::string mstk::Log<T>::toString(LogLevel level)
{
    if (level > getReportingLevel() || level < logNO_LOGGING) {
        mstk::Log<T>().get(mstk::logWARNING)
                << "Log<T>::toString(): Unknown logging level '" << level
                << "'. Using INFO level as default.";
        return "INFO";
    }

    static const char* const buffer[] = { "NO_LOGGING", "ERROR", "WARNING",
                                          "INFO", "DEBUG", "DEBUG1", "DEBUG2",
                                          "DEBUG3", "DEBUG4" };
    return buffer[level];
}

// fromString()
template<typename T>
mstk::LogLevel mstk::Log<T>::fromString(const std::string& level)
{
    if (level == "DEBUG4")
        return mstk::logDEBUG4;
    if (level == "DEBUG3")
        return mstk::logDEBUG3;
    if (level == "DEBUG2")
        return mstk::logDEBUG2;
    if (level == "DEBUG1")
        return mstk::logDEBUG1;
    if (level == "DEBUG")
        return mstk::logDEBUG;
    if (level == "INFO")
        return mstk::logINFO;
    if (level == "WARNING")
        return mstk::logWARNING;
    if (level == "ERROR")
        return mstk::logERROR;
    if (level == "NO_LOGGING")
        return mstk::logNO_LOGGING;

    // else
    mstk::Log<T>().get(mstk::logWARNING)
            << "Log<T>::fromString(): Unknown logging level '" << level
            << "'. Using INFO level as default.";
    return mstk::logINFO;
}

// FILELOG_DECLSPEC
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#   if defined (BUILDING_FILELOG_DLL)
#       define FILELOG_DECLSPEC   __declspec (dllexport)
#   elif defined (USING_FILELOG_DLL)
#       define FILELOG_DECLSPEC   __declspec (dllimport)
#   else
#       define FILELOG_DECLSPEC
#   endif 
#else
#   define FILELOG_DECLSPEC
#endif // _WIN32

// FILELog
/**
 * An instance of LOG<T> which is writing to a FILE.
 */
class FILELOG_DECLSPEC FILELog : public Log<Output2FILE>
{
};
//typedef Log<Output2FILE> FILELog;


// FILELOG_MAX_LEVEL
#ifndef FILELOG_MAX_LEVEL
/**
 * The deepest logging level to be compiled into the code.
 *
 * Every logging message depper than that level will not be compiled into the code.
 */
#define FILELOG_MAX_LEVEL mstk::logDEBUG4
#endif

// MSTK_LOG()
/**
 * Logs to a file handle.
 *
 * This macro checks, if the logging level should be compiled. After that, it creates an
 * anonymous instance of FILELog and writes to its logging stream. Afterwards, the
 * anonymous object is destroyed and the logging stream flushed out to the FILE.
 *
 * Use it like this:
 * @code
 * MSTK_LOG(logINFO) << "some logging" << 1224 << "no endl, will be appended automatically";
 * @endcode
 */
#define MSTK_LOG(level) \
    if (level > FILELOG_MAX_LEVEL) ;\
    else if (level > mstk::FILELog::getReportingLevel() || !mstk::Output2FILE::getRedirect()) ; \
    else mstk::FILELog().get(level)

// nowTime()
// We have to do the following yaketiyak, because the standard <ctime> is not thread safe.
// (It is using static internal buffers in some functions like ctime() .)
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
} // Temporarily close ms namespace to include the external Windows headers.

// winsocks2.h has always to be included BEFORE windows.h
// We don't use winsocks2 here, but it may be used in a file including this header.
#include <winsock2.hpp>
#include <windows.hpp>

// Reopen the ms namespace.
namespace mstk {
inline std::string nowTime()
{
    const int MAX_LEN = 200;
    char buffer[MAX_LEN];

    // get time
    if (GetTimeFormatA(
                    LOCALE_USER_DEFAULT, // locale
                    0, // time format flags
                    0, // optional ptr to systemtime structure
                    "HH':'mm':'ss", // format
                    buffer, // ptr to output buffer
                    MAX_LEN) // size of output buffer
            == 0) {
        return "Error in nowTime()";
    }

    // format time according to our format: "hh:mm:ss.ms"
    static DWORD first = GetTickCount();
    char result[100] = {0};
    std::sprintf(result, "%s.%03ld", buffer, (long)(GetTickCount() - first) % 1000);

    return result;
}

#else
} // Temporarily close ms namespace to inclue header files.
#include <sys/time.h>

// Reopen namespace mstk.
namespace mstk {
inline std::string nowTime()
{
    // get time
    time_t t;
    t = time(NULL);
    if (t == static_cast<std::time_t> (-1)) {
        return "Error_in_nowTime().time";
    }

    // convert time to local time
    tm r = { 0 };
    if (localtime_r(&t, &r) == NULL) {
        return "Error_in_nowTime().localtime_r";
    }

    // convert localtime to a string
    char buffer[101];
    if (strftime(buffer, sizeof(buffer), "%X", &r) == 0) {
        return "Error_in_nowTime().strftime";
    }

    // format the string according to our format: "hh:mm:ss.ms"
    struct timeval tv;
    gettimeofday(&tv, 0);
    char result[101] = { 0 };
    std::sprintf(result, "%s.%03ld", buffer, (long) tv.tv_usec / 1000);

    return result;
}

#endif //WIN32

/** @} */

} /* namespace mstk */

#endif /* __MSTK_INCLUDE_MSTK_LOG_HPP__ */
