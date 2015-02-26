/*
 * Logger.cpp
 *
 *  Created on: May 15, 2013
 *      Author: dfunke
 */

#include "Logger.h"

Logger &Logger::getInstance() {
  static Logger instance; // Guaranteed to be destroyed and thread safe
  // Instantiated on first use.
  return instance;
}

thread_local uint Logger::indentLevel = 0;

thread_local Logger::Verbosity
    Logger::contVerbosity; // continue log level for LIVE
thread_local Logger::LogEntries::iterator
    Logger::contIt; // continue stream for non-live

std::ostream &Logger::addLogEntry(Verbosity level) {

  if (logLevel <= Logger::Verbosity::LIVE) {
    contVerbosity = level; // thread local
    if (level <= abs(logLevel)) {
      return std::cout << indent();
    } else {
      return nullStream;
    }
  }

  if (logLevel == Logger::Verbosity::SILENT)
    return nullStream;

  contIt = logEntries.emplace_back(
      level, std::make_unique<std::stringstream>()); // thread local

  return *(contIt->second) << indent();
}

std::ostream &Logger::continueLogEntry() {

  if (logLevel <= Logger::Verbosity::LIVE) {
    if (contVerbosity <= abs(logLevel)) {
      return std::cout;
    } else {
      return nullStream;
    }
  }

  if (logLevel == Logger::Verbosity::SILENT)
    return nullStream;

  return *(contIt->second);
}

std::ostream &Logger::printLog(std::ostream &out) const {
  return printLog(out, logLevel);
}

std::ostream &Logger::printLog(std::ostream &out, Verbosity level) const {

  // shortcut SILENT processing
  if (level == Logger::Verbosity::SILENT)
    return out;

  for (const LogEntry &t : logEntries) {
    if (t.first <= level) {
      out << t.second->str();
    }
  }

  return out;
}

std::string Logger::indent() const {
  std::stringstream ss;
  for (uint i = 0; i < indentLevel; ++i)
    ss << "\t";
  return ss.str();
}

std::ostream &operator<<(std::ostream &s, const Logger &log) {

  return log.printLog(s);
}
