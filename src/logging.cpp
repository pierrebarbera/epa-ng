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
 * @brief Implementation of Logging functions.
 *
 * @file
 * @ingroup utils
 */

#include "logging.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>

#ifdef GENESIS_PTHREADS
#    include <mutex>
#endif

// #include "genesis/utils/io/output_stream.hpp"
// #include "genesis/utils/text/string.hpp"
// #include "genesis/utils/tools/date_time.hpp"

namespace genesis {
namespace utils {

#include <cstdio>
#include <ctime>

// to avoid includes
static std::string current_date()
{
    time_t now = time(0);
    tm*    ltm = localtime(&now);

    char out[12];
    sprintf (out, "%u-%02u-%02u",
        ltm->tm_year + 1900, ltm->tm_mon + 1, ltm->tm_mday
    );
    return out;
}

static std::string current_time()
{
    time_t now = time(0);
    tm*    ltm = localtime(&now);

    char out[10];
    sprintf (out, "%02u:%02u:%02u", ltm->tm_hour, ltm->tm_min, ltm->tm_sec);
    return out;
}

static std::string replace_all (
    std::string const& text, std::string const& search, std::string const& replace
) {
    std::string tmp = text;
    for (size_t pos = 0; ; pos += replace.length()) {
        pos = tmp.find(search, pos);

        if (pos == std::string::npos)
            break;

        tmp.erase(pos, search.length());
        tmp.insert(pos, replace);
    }
    return tmp;
}

static std::string trim_right (
    std::string const& s,
    std::string const& delimiters = " \f\n\r\t\v"
) {
    auto const pos = s.find_last_not_of(delimiters);
    if( std::string::npos == pos ) {
        return "";
    } else {
        return s.substr( 0, pos + 1 );
    }
}

// =============================================================================
//     Settings
// =============================================================================

#ifdef GENESIS_PTHREADS
    static std::mutex log_mutex;
#endif

// TODO use different init for log details depending on DEBUG

// init static members
LoggingDetails Logging::details = {
    false, // count
    false, // date
    false, // time
    false, // runtime
    false, // rundiff
    false, // file
    false, // line
    false, // function
    true   // level
};
Logging::LoggingLevel      Logging::max_level_  = kDebug4;
long                       Logging::count_      = 0;
clock_t                    Logging::last_clock_ = 0;
std::vector<std::ostream*> Logging::ostreams_;
int                        Logging::report_percentage_ = 5;
std::string                Logging::debug_indent       = "    ";

/**
 * @brief Set the highest log level that is reported.
 *
 * Invocations of log with higher levels will create no output.
 * It creates a warning if the set level is higher than the static compile time
 * level set by #LOG_LEVEL_MAX.
 */
void Logging::max_level (const LoggingLevel level)
{
    if (level > LOG_LEVEL_MAX) {
        LOG_WARN << "Logging max level set to " << level << ", but compile "
                 << "time max level is " << LOG_LEVEL_MAX << ", so that "
                 << "everything above that will not be logged.";
    }
    max_level_ = level;
}

/**
 * @brief set the percentage for reporting #LOG_PROG messages.
 */
void Logging::report_percentage (const int percentage)
{
    if (percentage <= 0) {
        LOG_WARN << "Logging report percentage less than 1% not possible.";
        report_percentage_ = 1;
        return;
    }
    if (percentage > 100) {
        LOG_WARN << "Logging report percentage greater than 100% not meaningful.";
        report_percentage_ = 100;
        return;
    }
    report_percentage_ = percentage;
}

/**
 * @brief Return a string representation of a log level.
 */
std::string Logging::level_to_string(const LoggingLevel level)
{
    static const char* const buffer[] = {
        "NONE", "ERR ", "WARN", "INFO", "PROG", "DBG ", "DBG1", "DBG2", "DBG3", "DBG4"
    };
    return buffer[level];
}

/**
 * @brief Add stdout as output stream to which log messages are written.
 */
void Logging::log_to_stdout (const bool b)
{
    if (b) {
        // check whether stdout was already added.
        for (std::ostream* os : ostreams_) {
            if (os == &std::cout) {
                return;
            }
        }

        // if not, add it as output stream.
        ostreams_.push_back (&std::cout);
    } else {
        ostreams_.erase(
            std::remove(std::begin(ostreams_), std::end(ostreams_), &std::cout),
            std::end(ostreams_)
        );
    }
}

/**
 * @brief Add an output stream to which log messages are written.
 */
void Logging::log_to_stream (std::ostream& os)
{
    ostreams_.push_back (&os);
}

/**
 * @brief Add an output file to which log messages are written.
 *
 * This creates a stream to the file.
 */
void Logging::log_to_file( std::string const& filename )
{
    // TODO the log file stream is never deleted. this is not a big leak,
    // as commonly only one file is used for logging, but still is a smell.
    // use uniqute stored statically in the log class instead!

    std::ofstream* file_stream = new std::ofstream();
    // utils::file_output_stream( filename, *file_stream );
    file_stream->open(filename, std::ios_base::out);
    ostreams_.push_back( file_stream );
    if( file_stream->fail() ) {
        throw std::runtime_error( "Cannot open file '" + filename + "'." );
    }
}

// =============================================================================
//     Destructor (does the actual work)
// =============================================================================

/**
 * @brief Destructor that is invoked at the end of each log line and does the actual
 * output.
 */
Logging::~Logging()
{
    // build the details for the log message into a buffer
    clock_t now_clock = clock();
    std::ostringstream det_buff;
    det_buff.str("");
    if (details_.count) {
        det_buff.fill('0');
        det_buff.width(4);
        det_buff << count_ << " ";
    }
    if (details_.date) {
        det_buff << current_date() << " ";
    }
    if (details_.time) {
        det_buff << current_time() << " ";
    }
    if (details_.runtime) {
        det_buff << std::fixed
                 << std::setprecision(6)
                 << double(now_clock) / CLOCKS_PER_SEC
                 << " ";
    }
    if (details_.rundiff) {
        double val = 0.0;
        if (last_clock_ > 0) {
            val = (double) (now_clock - last_clock_) / CLOCKS_PER_SEC;
        }
        det_buff << std::fixed
                 << std::setprecision(6)
                 << val
                 << " ";
        last_clock_ = now_clock;
    }
    if (details_.file) {
        det_buff << file_ << (details_.line ? "" : " ");
    }
    if (details_.line) {
        det_buff << ":" << line_ << " ";
    }
    if (details_.function) {
        det_buff << "(" << function_ << ") ";
    }
    if (details_.level) {
        det_buff << level_to_string(level_) << " ";
    }

    // add spaces for nested debug levels
    if (level_ > kDebug) {
        for (int i = 0; i < level_ - kDebug; i++) {
            det_buff << debug_indent;
        }
    }

    // make multi line log messages align to the length of the detail header,
    // and trim trailing whitespace, as we only want one newline at the end
    std::string msg = det_buff.str();
    if (msg.length() > 0) {
        msg += utils::replace_all(
            buff_.str(), "\n", "\n" + std::string(msg.length(), ' ')
        );
    } else {
        msg += buff_.str();
    }
    msg = utils::trim_right(msg);

    // output the message to every stream, thread safe!
#   ifdef GENESIS_PTHREADS
    log_mutex.lock();
#   endif
    for (std::ostream* out : ostreams_) {
        (*out) << msg << std::endl << std::flush;
    }
#   ifdef GENESIS_PTHREADS
    log_mutex.unlock();
#   endif

    // inc log message counter
    count_++;
}

// =============================================================================
//     Singleton accessors
// =============================================================================

/**
 * @brief Getter for the singleton instance of log, is called by the standard macros.
 *
 * It returns the string stream buffer used to capture the log messages.
 */
std::ostringstream& Logging::get(
    const std::string& file, const int line, const std::string& function,
    const LoggingLevel level
)
{
    return get(file, line, function, level, details);
}

/**
 * @brief Getter for the singleton instance of log, is called by special macros
 * that change the details of the log message.
 *
 * It stores some relevant information and returns the string stream buffer
 * used to capture the log messages.
 */
std::ostringstream& Logging::get(
    const std::string& file, const int line, const std::string& function,
    const LoggingLevel level, const LoggingDetails dets
)
{
    // save the information given when called from the macros
    file_     = file;
    line_     = line;
    function_ = function;
    level_    = level;
    details_  = dets;
    buff_.str("");
    return buff_;
}

} // namespace utils
} // namespace genesis
