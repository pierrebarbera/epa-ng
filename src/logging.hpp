#ifndef GENESIS_UTILS_CORE_LOGGING_H_
#define GENESIS_UTILS_CORE_LOGGING_H_

/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2017 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

/**
 * @brief Provides easy and fast logging functionality.

 *
 * For more information on the logging, see Logging class.
 *
 * @file
 * @ingroup utils
 */

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace genesis {
namespace utils {

// =============================================================================
//     Macro definitions
// =============================================================================

// TODO add TIME1..4 (or TMR maybe) for indented timing logs
// TODO change macros for timing to be out of usual log levels
// TODO define a function that logs the detail column headers (difficult because
// of length difference for file name log detail)
// TODO offer csv as output format
// TODO offer remote streams

#ifndef LOG_LEVEL_MAX
    /**
     * @brief Static maximal logging level.
     *
     * Everything above this level will be pruned by the compiler.
     */
#    define LOG_LEVEL_MAX genesis::utils::Logging::kDebug4
#endif

// TODO make DEBUG a special macro with proper usage makefile etc,
// also add maybe stuff like RELEASE TEST etc, prepend ENV_ or so!
// TODO do an undefine first

// override this setting when debugging is turned on (eg by makefile)
#ifdef DEBUG
#    define LOG_LEVEL_MAX genesis::utils::Logging::kDebug4
#endif

// try to find a macro that expands to the current function name
#ifdef __cplusplus
#    define GENESIS_FUNC __PRETTY_FUNCTION__
#else
#    if defined __STDC_VERSION__
#        define GENESIS_FUNC __func__
#    else
#        define GENESIS_FUNC ((const char *) 0)
#    endif
#endif

// define the actual log macro, so that the compiler can prune calls
// to a log with levels above the static max log level
#define GENESIS_LOG(level) \
    if (level > LOG_LEVEL_MAX) ; \
    else if (level > genesis::utils::Logging::max_level()) ; \
    else genesis::utils::Logging().get(__FILE__, __LINE__, GENESIS_FUNC, level)

// define a similar log macro, this time changing the details of the message
#define GENESIS_LOG_DETAILS(level, ...) \
    if (level > LOG_LEVEL_MAX) ; \
    else if (level > genesis::utils::Logging::max_level()) ; \
    else genesis::utils::Logging().get(__FILE__, __LINE__, GENESIS_FUNC, level, {__VA_ARGS__})

// define standard logging types as macro shortcuts:

/** @brief Log an error. See genesis::utils::LoggingLevel. */
#define LOG_ERR  GENESIS_LOG(genesis::utils::Logging::kError)

/** @brief Log a warning. See genesis::utils::LoggingLevel. */
#define LOG_WARN GENESIS_LOG(genesis::utils::Logging::kWarning)

/** @brief Log an info message. See genesis::utils::LoggingLevel. */
#define LOG_INFO  GENESIS_LOG(genesis::utils::Logging::kInfo)

/** @brief Log a debug message. See genesis::utils::LoggingLevel. */
#define LOG_DBG  GENESIS_LOG(genesis::utils::Logging::kDebug)

/** @brief Log a debug message. See genesis::utils::LoggingLevel. */
#define LOG_DBG1 GENESIS_LOG(genesis::utils::Logging::kDebug1)

/** @brief Log a debug message. See genesis::utils::LoggingLevel. */
#define LOG_DBG2 GENESIS_LOG(genesis::utils::Logging::kDebug2)

/** @brief Log a debug message. See genesis::utils::LoggingLevel. */
#define LOG_DBG3 GENESIS_LOG(genesis::utils::Logging::kDebug3)

/** @brief Log a debug message. See genesis::utils::LoggingLevel. */
#define LOG_DBG4 GENESIS_LOG(genesis::utils::Logging::kDebug4)

// define special log shortcuts: the list of bools represent
// the members of struct LogDetails and indicate which parts shall be included.

/**
 * @brief %Logging of a message that is always displayed.
 *
 * It does not include any details with its message stream (thus, is
 * independent of Logging::details) and is always logged (independent of the max
 * levels). This is for example used to log the program header and footer on startup and
 * termination.
 */
#define LOG_BOLD GENESIS_LOG_DETAILS(genesis::utils::Logging::kNone, \
    false, false, false, false, false, false, false, false, false)

/**
 * @brief %Logging of a message with timing information.
 *
 * It includes the run time difference to the last log message in seconds
 * as its only detail (independent of Logging::details). This is particularly
 * useful for timing and profiling code sections. Its level is Logging::kDebug, so
 * that in can be easily turned of for production code.
 */
#define LOG_TIME GENESIS_LOG_DETAILS(genesis::utils::Logging::kDebug, \
    false, false, false, false , true,  false, false, false, false)

/**
 * @brief This macro sets the logging level to a certain value for this scope and sets
 * it back to the original value on exiting the scope.
 *
 * This is useful to disable debug messages in a function while still beging able to activate them
 * when needed. The macro creates class instance of LoggingScopeLevel with a long and hopefully
 * unique name.
 */
#define LOG_SCOPE_LEVEL(level) \
    genesis::utils::LoggingScopeLevel genesis_logging_scope_level_temp_object(level);

/**
 * @brief Hack function to make sure that the value arugment in #LOG_PROG is only evaluated once.
 *
 * Without this function, #LOG_PROG would include two appearances of its variable `value`, which
 * means that a statement like
 *
 *     LOG_PROG(++i, n) << "of progress.";
 *
 * would lead to a double evaluation of the increment statement `++i`. That is not intended, thus
 * we need this hack function.
 */
inline long logging_progress_value (long value = -1)
{
    static long v = -1;
    if (value > -1) {
        v = value;
    }
    return v;
}

/**
 * @brief %Logging of a progress message.
 *
 * This special logging mechanism provides an easy way to report progress while running long
 * operations. It uses LoggingLevel::kProgress for its messages.
 *
 * It takes two arguments, the first one being the number of the current iteration of the loop,
 * the second the total number of iterations, both are casted to long.
 *
 * The integer value of Logging::report_percentage can be used to control how often progress is
 * reported in terms of a percentage counter. Default is 5, meaning every 5% progress is reported.
 *
 * The trailing message is optional -- when left out, simply the percentage number will be logged.
 * However, the statement must always be ended by a semicolon.
 *
 * The default use case looks like this:
 *
 *     Logging::report_percentage(2);
 *     for (int i = 0; i < n; ++i) {
 *         // long loop from 0 to n
 *         // do stuff...
 *         LOG_PROG(i, n) << "of the loop finished.";
 *     }
 *
 * This will report the progress of the loop every 2 percent with a message like this:
 *
 *     0% of the loop finished.
 *     2% of the loop finished.
 *     4% of the loop finished.
 *     ...
 *
 * It is also possible to increment the iteration counter inline:
 *
 *     LOG_PROG(++i, n) << "of stuff finished.";
 *
 * It automatically checks whether it is time to report, so that the user does not have to take care
 * of this.
 *
 * Caveat: For small number of iterations (<100) it is not entirely accurate. Because
 * of the involved integer modulo calculations, it will trigger too many messages. However, as this
 * type of logging is usually used for loops with many iterations, this should rarely be an issue.
 *
 * There is a slight overhead of ~60ms per 1mio invocations because of the needed calculations.
 */
#define LOG_PROG( value, quantity ) \
    if( genesis::utils::Logging::kProgress > LOG_LEVEL_MAX ) ; \
    else if( genesis::utils::Logging::kProgress > genesis::utils::Logging::max_level() ) ; \
    else if( (long) genesis::utils::logging_progress_value(value) % \
        (((long) (quantity) * genesis::utils::Logging::report_percentage() / 100) > 0 ? \
        ((long) (quantity) * genesis::utils::Logging::report_percentage() / 100) : 1) != 0 ) ; \
    else genesis::utils::Logging().get(__FILE__, __LINE__, GENESIS_FUNC, genesis::utils::Logging::kProgress) << \
    (int) round(100.0 * (double) genesis::utils::logging_progress_value() / ((quantity) > 0 ? (quantity) : 1)) << "% "

// =============================================================================
//     Logging Details
// =============================================================================

/**
 * @brief POD stuct containing the settings for which information is included
 * with each logging message.
 *
 * The details are activated via accessing the static variable of the
 * Logging class:
 *
 *     Logging::details.level = true;
 *
 * All active details are prepended to the actual log message and
 * separated by spaces (execpt file and line, they are separated by a
 * colon). Their order is fixed.
 *
 * A message with all details activates looks like this
 *
 *     0003 2014-10-28 11:40:47 0.001859 0.000103 src/main/main.cc:28 (int main (int argc, char* argv[])) INFO Hello World.
 *
 * It was the third message being logged in this run of the program, at
 * a date and time, 0.001859 sec after the program started and 0.000103
 * sec after the last log message. It was called from main.cc line 28 in function main
 * and has LoggingLevel Info.
 */
typedef struct {
    /** @brief Include a counter of how many messages have been logged so far. */
    bool count;

    /** @brief Include the current date. */
    bool date;

    /** @brief Include the current time. */
    bool time;

    /** @brief Include the current run time of the program in sec. */
    bool runtime;

    /**
     * @brief Include the run time difference to the last log message
     * in sec.
     *
     * Useful for timing and profiling code sections. Is 0.0 at the first
     * log message.
     */
    bool rundiff;

    /** @brief Include the filename where the log message was generated. */
    bool file;

    /**
     * @brief Include the line of the file where the log message was
     * generated.
     */
    bool line;

    /**
     * @brief Include the function name where the log message was generated.
     *
     * This might not be available on certain compilers. In this case, it is empty.
     */
    bool function;

    /** @brief Include the level (e.g. Info, Debug) of the message. */
    bool level;
} LoggingDetails;

// =============================================================================
//     Logging
// =============================================================================

/**
 * @brief %Logging class with easy and fast usage.
 *
 * The basic usage of this class is to invoke the macros for the different types
 * of log messages and send a stream of messages to them:
 *
 *     LOG_DBG << "you are here";
 *     LOG_ERR << "there was an error: " << 42;
 *
 * The provided types of macros are: #LOG_ERR, #LOG_WARN, #LOG_INFO ,
 * #LOG_DBG, #LOG_DBG1, #LOG_DBG2, #LOG_DBG3, #LOG_DBG4 for all levels of
 * logging explained in LoggingLevel.
 *
 * The details that are logged with each message can be changed via accessing
 * Logging::details -- see LoggingDetails for more on that.
 *
 * In order to use this class, at least one output stream has to be added first
 * by invoking either AddOutputStream or AddOutputFile.
 *
 * The depths of logging can be changed in order to reduce the amount of written
 * messages. First, at compile time the macro constant LOG_LEVEL_MAX can be set
 * (e.g. in the Makefile) to the highest level that shall be logged. All log
 * invoakations with a higher level are never called, and moreover will be
 * pruned from the code by most compilers completetly. This makes the log class
 * fast -- instead of deleting all log invokations by hand, it is sufficient to
 * set the constant to a low level.
 * Second, the depths of logging can be changed dynamically at run time by
 * setting Logging::max_level to the desired value. Of course, this value
 * cannot be higher than LOG_LEVEL_MAX, because those logs are already pruned
 * by the compiler.
 *
 * There are also three more special log types that create a different output than
 * the previously mentioned types: #LOG_BOLD, #LOG_TIME and #LOG_PROG. See their
 * respective documentation for more information.
 *
 * The inner working of this class is as follows: Upon invokation via one of the
 * macros, an instance is created that stays alive only for the rest of the
 * source code line. In this line, the log message is inserted to the buffer
 * stream. After the line finishes, the object is automatically destroyed, so
 * that its destructor is called. In the destructor, the message is composed and
 * logged to all given output streams.
 *
 * Caveat: Because the macro contains conditions depending on the log level,
 * do not use code in a log line that changes the program state.
 * Therefore, instead of
 *
 *     LOG_INFO  << "this is iteration " << ++i;
 *
 * use
 *
 *     ++i;
 *     LOG_INFO  << "this is iteration " << i;
 *
 * because the former will not work when the max log level is below info level.
 */
class Logging
{
public:

    // -------------------------------------------------------------------
    //     Class Functions
    // -------------------------------------------------------------------

    /**
     * @brief Levels of severity used for logging.
     *
     * The levels are in ascending order and are used both to signal what
     * kind of message is being logged and to provide a threshold for less
     * important messages that can be filtered out, for example debug
     * messages in the production build of the program. Because some of the
     * filtering is already done at compile time, log messages with a level
     * higher than #LOG_LEVEL_MAX do not produce any overhead. See also
     * Logging class for more on this.
     */
    enum LoggingLevel {
        /** @brief Special messages that are always logged, e.g. program header. */
        kNone = 0,

        /** @brief Errors, usually non-recoverable. See #LOG_ERR. */
        kError,

        /** @brief Warnings if somthing went wront, but program can continue. See #LOG_WARN. */
        kWarning,

        /** @brief Infos, for example when a file was written. See #LOG_INFO . */
        kInfo,

        /** @brief Progess, used in long executing functions. See #LOG_PROG. */
        kProgress,

        /** @brief Basic debugging message. See #LOG_DBG. */
        kDebug,

        /** @brief Debugging message with indent level 1 (e.g. for loops). See #LOG_DBG1. */
        kDebug1,

        /** @brief Debugging message with indent level 2. See #LOG_DBG2. */
        kDebug2,

        /** @brief Debugging message with indent level 3. See #LOG_DBG3. */
        kDebug3,

        /** @brief Debugging message with indent level 4. See #LOG_DBG4. */
        kDebug4
    };

    Logging() {};
    ~Logging();

    // Logging is kind of singleton, its instances are only provided via the
    // get functions. do not allow other instances by blocking the copy
    // constructors and copy assignment
    Logging (const Logging&) {};
    Logging& operator = (const Logging&) = delete;

    // getter for the singleton instance of log, is called by the macros
    std::ostringstream& get (
        const std::string& file, const int line, const std::string& function,
        const LoggingLevel level
    );
    std::ostringstream& get (
        const std::string& file, const int line, const std::string& function,
        const LoggingLevel level, const LoggingDetails dets
    );

    // -------------------------------------------------------------------
    //     Logging Settings
    // -------------------------------------------------------------------

    // TODO make logging accept one file, and a setter for cout and for cerr
    // TODO allow different levels to be logged to different streams?!

    // methods to handle the output streams to write the log messages to
    static void log_to_stdout (const bool b = true);
    static void log_to_stream (std::ostream& os);
    static void log_to_file   (const std::string& fn);

    /**
     * @brief Settings for which information is included with each log message.
     * See LoggingDetails for usage.
     */
    static LoggingDetails details;

    /** @brief Get the highest log level that is reported. */
    static inline LoggingLevel max_level ()
    {
        return max_level_;
    }
    static void max_level (const LoggingLevel level);

    /** @brief Get the current percentage for reporting #LOG_PROG messages. */
    static inline int report_percentage ()
    {
        return report_percentage_;
    }
    static void report_percentage (const int percentage);

    // return a string representation for a log level
    static std::string level_to_string (const LoggingLevel level);

    /** @brief Indention string for Debug Levels 1-4. */
    static std::string debug_indent;

    // -------------------------------------------------------------------
    //     Static Wrappers
    // -------------------------------------------------------------------

    static inline void log_error (const std::string& msg)
    {
        LOG_ERR << msg;
    }

    static inline void log_warning (const std::string& msg)
    {
        LOG_WARN << msg;
    }

    static inline void log_info (const std::string& msg)
    {
        LOG_INFO  << msg;
    }

    static inline void log_debug (const std::string& msg)
    {
        LOG_DBG << msg;
    }

    static inline void log_debug_1 (const std::string& msg)
    {
        LOG_DBG1 << msg;
    }

    static inline void log_debug_2 (const std::string& msg)
    {
        LOG_DBG2 << msg;
    }

    static inline void log_debug_3 (const std::string& msg)
    {
        LOG_DBG3 << msg;
    }

    static inline void log_debug_4 (const std::string& msg)
    {
        LOG_DBG4 << msg;
    }

    // -------------------------------------------------------------------
    //     Internal Data Members
    // -------------------------------------------------------------------

protected:
    // storage for information needed during one invocation of a log
    std::ostringstream buff_;
    std::string        file_;
    int                line_;
    std::string        function_;
    LoggingLevel       level_;
    LoggingDetails     details_;

    // dynamic log level limit
    static LoggingLevel max_level_;

    // how often to report progress messages
    static int report_percentage_;

    // how many log calls were made so far
    static long    count_;

    // when was the last call to logging (used for time measurements)
    static clock_t last_clock_;

    // array of streams that are used for output
    static std::vector<std::ostream*> ostreams_;
};

// =============================================================================
//     Logging Scope Level
// =============================================================================

/**
 * @brief Class that sets the Logging Level to a value von construction and set it back on
 * destruction. This is used by the log scope level macro.
 */
class LoggingScopeLevel
{
public:
    LoggingScopeLevel(Logging::LoggingLevel scope_level)
    {
        previous_level = Logging::max_level();
        Logging::max_level(scope_level);
    }

    ~LoggingScopeLevel()
    {
        Logging::max_level(previous_level);
    }

private:
    Logging::LoggingLevel previous_level;
};

// =============================================================================
//     Logging Progress
// =============================================================================

/*
 * This was a test to make LOG_PROG work with incrementing counters like
 *
 *     LOG_PROG(++progress, total) << "of something";
 *
 * but it is slow and dirty. If used, Logging().get(...) has to be expanded by a
 * std::string init parameters, that is used to initialize buff_.str(init).
 *
 * Also, the definition of LOG_PROG that was used for this is:
 *
 *     #define LOG_PROG(value, quantity) \
 *         if (genesis::utils::Logging::kProgress > LOG_LEVEL_MAX) ; \
 *         else if (genesis::utils::Logging::kProgress > genesis::utils::Logging::max_level()) ; \
 *         else LOG_PROG_FUNC(value, quantity, __FILE__, __LINE__, GENESIS_FUNC)
 * %
inline std::ostringstream& LOG_PROG_FUNC (const int value, const int quantity, const std::string& file, const int line, const std::string& function)
{
    static std::ostringstream trash;
    int base = quantity * Logging::report_percentage() / 100;
    if (value % (base > 0 ? base : 1) != 0) {
        trash.str("");
        return trash;
    } else {
        std::string init = std::to_string((int) round(100.0 * (double) value / quantity)) + "% ";
        return Logging().get(file, line, function, Logging::kProgress, Logging::details, init);
    }
}
*/

} // namespace utils
} // namespace genesis

#endif // include guard
