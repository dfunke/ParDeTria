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

std::ostream &Logger::addLogEntry(Verbosity level) {

    if (logLevel <= Logger::Verbosity::LIVE) {
        return std::cout << indent();
    }

    auto it = logEntries.emplace_back(
            level, std::make_unique<std::stringstream>()); // thread local

    return *(it->second) << indent();
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
