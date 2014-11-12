/*
 * Logger.cpp
 *
 *  Created on: May 15, 2013
 *      Author: dfunke
 */

#include "Logger.h"

Logger& Logger::getInstance()
{
	static Logger    instance; // Guaranteed to be destroyed and thread safe
	// Instantiated on first use.
	return instance;
}

thread_local Logger::Verbosity Logger::contVerbosity; //continue log level for LIVE
thread_local tbb::concurrent_vector<Logger::tLogEntry>::iterator Logger::contIt; //continue stream for non-live

std::ostream & Logger::addLogEntry(Verbosity level) {

	if(logLevel <= Logger::Verbosity::LIVE){
		if(level <= abs(logLevel)){
			contVerbosity = level; //thread local
			return std::cout << indent();
		} else {
			return nullStream;
		}
	}

	if(logLevel == Logger::Verbosity::SILENT)
		return nullStream;

	std::stringstream* s = new std::stringstream();

	tLogEntry entry = std::make_pair(level, s);
	contIt = logEntries.push_back(entry); //thread local

	return *(contIt->second) << indent();
}

std::ostream & Logger::continueLogEntry() {

	if(logLevel <= Logger::Verbosity::LIVE){
		if(contVerbosity <= abs(logLevel)){
			return std::cout;
		} else {
			return nullStream;
		}
	}

	if(logLevel == Logger::Verbosity::SILENT)
		return nullStream;

	return *(contIt->second);
}

std::ostream & Logger::printLog(std::ostream & out){
	return printLog(out, logLevel);
}

std::ostream & Logger::printLog(std::ostream & out, Verbosity level){

	//shortcut SILENT processing
	if(level == Logger::Verbosity::SILENT)
		return out;

	for(tLogEntry t : logEntries){
		if( t.first <= level){
			out << t.second->str();
		}
	}

	return out;
}

Logger::~Logger(){
	for(tLogEntry t : logEntries){
		delete t.second;
	}
}

std::string Logger::indent() const {
	std::stringstream ss;
	for(uint i = 0; i < indentLevel; ++i)
		ss << "\t";
	return ss.str();
}

std::ostream & operator<<(std::ostream& s, const Logger & log){

	return Logger::getInstance().printLog(s);
}



