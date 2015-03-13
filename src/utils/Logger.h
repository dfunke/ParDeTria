#pragma once

#include <boost/noncopyable.hpp>
#include <vector>
#include <iostream>
#include <sstream>

#include "utils/ASSERT.h"

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

#define IS_NORMAL                                                              \
  Logger::abs(Logger::getInstance().getLogLevel()) >= Logger::Verbosity::NORMAL
#define IS_VERBOSE                                                             \
  Logger::abs(Logger::getInstance().getLogLevel()) >= Logger::Verbosity::VERBOSE
#define IS_PROLIX                                                              \
  Logger::abs(Logger::getInstance().getLogLevel()) >= Logger::Verbosity::PROLIX

#define LOGGER Logger::getInstance()
#define LOG(msg)                                                               \
  if (IS_NORMAL)                                                               \
  Logger::getInstance().addLogEntry(Logger::Verbosity::NORMAL) << msg
#define VLOG(msg)                                                              \
  if (IS_VERBOSE)                                                              \
  Logger::getInstance().addLogEntry(Logger::Verbosity::VERBOSE) << msg

#ifndef NDEBUG
#define PLOG(msg)                                                              \
  if (IS_PROLIX)                                                               \
  Logger::getInstance().addLogEntry(Logger::Verbosity::PROLIX) << msg
#else // NDEBUG
#define PLOG(msg) ((void)(0))
#endif

#define INDENT Logger::getInstance().incIndent();
#define DEDENT Logger::getInstance().decIndent();

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

  std::ostream &printLog(std::ostream &out) const;

  std::ostream &printLog(std::ostream &out, Verbosity level) const;

  std::string indent() const;
  void incIndent() { ++indentLevel; }
  void decIndent() {
    ASSERT(indentLevel > 0);
    --indentLevel;
  }
  void setIndent(const uint i) { indentLevel = i; }
  uint getIndent() { return indentLevel; }

  template <typename InputIt>
  void logContainer(const InputIt &first, const InputIt &last, Verbosity level,
                    const std::string &prefix = "",
                    const std::string &sep = " ") {

    if (Logger::abs(Logger::getInstance().getLogLevel()) >= level) {
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
  }

  template <typename Container>
  void logContainer(const Container &c, Verbosity level,
                    const std::string &prefix = "",
                    const std::string &sep = " ") {

    logContainer(c.begin(), c.end(), level, prefix, sep);
  }

private:
  Logger() : logLevel(Logger::Verbosity::NORMAL){};

  LogEntries logEntries;

  Verbosity logLevel;
  static thread_local uint indentLevel;
};

std::ostream &operator<<(std::ostream &s, const Logger &log);
