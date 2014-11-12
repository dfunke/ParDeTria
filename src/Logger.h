#pragma once

#include <boost/noncopyable.hpp>
#include <vector>
#include <iostream>
#include <sstream>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/null.hpp>

#include <tbb/concurrent_vector.h>

class Logger : private boost::noncopyable {

#define LOG  Logger::getInstance().addLogEntry(Logger::Verbosity::NORMAL)
#define VLOG Logger::getInstance().addLogEntry(Logger::Verbosity::VERBOSE)
#define PLOG Logger::getInstance().addLogEntry(Logger::Verbosity::PROLIX)
#define LLOG std::cout << Logger::getInstance().indent()

#define CONT  Logger::getInstance().continueLogEntry()
#define LCONT std::cout

#define INDENT Logger::getInstance().incIndent();
#define DEDENT Logger::getInstance().decIndent();

#define IS_VERBOSE abs(Logger::getInstance().getLogLevel()) >= Logger::Verbosity::VERBOSE
#define IS_PROLIX  abs(Logger::getInstance().getLogLevel()) >= Logger::Verbosity::PROLIX

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

	 static Verbosity abs(const Verbosity & v){
		 return static_cast<Verbosity>(std::abs(static_cast<int>(v)));
	 }

public:

	typedef std::pair<Verbosity, std::stringstream*> tLogEntry;

	static Logger& getInstance();

	const Verbosity & getLogLevel() const { return logLevel; }
	void setLogLevel(const Verbosity & l) { logLevel = l; }

	std::ostream & addLogEntry(Verbosity level);
	std::ostream & continueLogEntry();

	std::ostream& printLog(std::ostream & out);

	std::ostream& printLog(std::ostream & out, Verbosity level);

	std::string indent() const;
	void incIndent() { ++indentLevel; }
	void decIndent() { --indentLevel; }

	~Logger();

private:
	Logger() : logLevel(Logger::Verbosity::NORMAL),
			   indentLevel(0),
			   nullStream( ( boost::iostreams::null_sink() ) ) { };

	tbb::concurrent_vector<tLogEntry> logEntries;

	Verbosity logLevel;
	uint indentLevel;

	boost::iostreams::stream< boost::iostreams::null_sink > nullStream;

	static thread_local Verbosity contVerbosity; //continue log level for LIVE
	static thread_local tbb::concurrent_vector<tLogEntry>::iterator contIt; //continue stream for non-live
};

std::ostream & operator<<(std::ostream& s, const Logger & log);
