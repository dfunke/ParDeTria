#pragma once

#include <boost/noncopyable.hpp>
#include <vector>
#include <iostream>
#include <sstream>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>

// fix bug in TBB where emplaced_back is only defined when LIBCPP is defined
#ifndef _LIBCPP_VERSION
#define _LIBCPP_VERSION 1
#define _UNDEF_LIBCPP 1
#endif

#include <tbb/concurrent_vector.h>

#ifdef _UNDEF_LIBCPP
#undef _LIBCPP_VERSION
#undef _UNDEF_LIBCPP
#endif

#include "ccolor.h"

class Logger : private boost::noncopyable {

#define LOGGER Logger::getInstance()
#define LOG Logger::getInstance().addLogEntry(Logger::Verbosity::NORMAL)
#define VLOG Logger::getInstance().addLogEntry(Logger::Verbosity::VERBOSE)
#define PLOG Logger::getInstance().addLogEntry(Logger::Verbosity::PROLIX)
#define LLOG std::cout << Logger::getInstance().indent()

#define CONT Logger::getInstance().continueLogEntry()
#define LCONT std::cout

#define INDENT Logger::getInstance().incIndent();
#define DEDENT Logger::getInstance().decIndent();

#define IS_VERBOSE                                                             \
  Logger::abs(Logger::getInstance().getLogLevel()) >= Logger::Verbosity::VERBOSE
#define IS_PROLIX                                                              \
  Logger::abs(Logger::getInstance().getLogLevel()) >= Logger::Verbosity::PROLIX

public:
  enum class Verbosity : short {
    LIVEPROLIX = -3,
    LIVEVERBOSE = -2,
    LIVE = -1,
    SILENT = 0,
    NORMAL = 1,
    VERBOSE = 2,
    PROLIX = 3
  };

  static Verbosity abs(const Verbosity &v) {
    return static_cast<Verbosity>(std::abs(static_cast<int>(v)));
  }

public:
  typedef std::pair<Verbosity, std::unique_ptr<std::stringstream>> LogEntry;
  typedef tbb::concurrent_vector<LogEntry> LogEntries;

  static Logger &getInstance();

  const Verbosity &getLogLevel() const { return logLevel; }
  void setLogLevel(const Verbosity &l) { logLevel = l; }

  std::ostream &addLogEntry(Verbosity level);
  std::ostream &continueLogEntry();

  std::ostream &printLog(std::ostream &out) const;

  std::ostream &printLog(std::ostream &out, Verbosity level) const;

  std::string indent() const;
  void incIndent() { ++indentLevel; }
  void decIndent() { --indentLevel; }

  template <typename InputIt>
  void logContainer(const InputIt &first, const InputIt &last, Verbosity level,
                    const std::string &prefix = "",
                    const std::string &sep = " ") {
    auto &stream = addLogEntry(level);

    if (prefix != "")
      stream << prefix << sep;

    for (auto it = first; it != last; ++it) {
      if (it != first)
        stream << sep;
      stream << *it;
    }
    stream << std::endl;
  }

  template <typename Container>
  void logContainer(const Container &c, Verbosity level,
                    const std::string &prefix = "",
                    const std::string &sep = " ") {

    logContainer(c.begin(), c.end(), level, prefix, sep);
  }

private:
  Logger()
      : logLevel(Logger::Verbosity::NORMAL), indentLevel(0),
        nullStream((boost::iostreams::null_sink())){};

  LogEntries logEntries;

  Verbosity logLevel;
  uint indentLevel;

  boost::iostreams::stream<boost::iostreams::null_sink> nullStream;

  static thread_local Verbosity contVerbosity; // continue log level for LIVE
  static thread_local LogEntries::iterator
      contIt; // continue stream for non-live
};

std::ostream &operator<<(std::ostream &s, const Logger &log);
