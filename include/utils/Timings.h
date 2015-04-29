/*
 * Timings.h
 *
 *  Created on: Oct 28, 2014
 *      Author: dfunke
 */

#pragma once

#include <chrono>
#include <string>
#include <vector>
#include <sstream>
#include <string>
#include <map>

#include <boost/noncopyable.hpp>

#include "utils/TBB_Containers.h"

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::nanoseconds tDuration;

class ExperimentRun {

public:
    typedef ulong tMeasurement;
    typedef tbb::concurrent_vector<tMeasurement> tMeasurements;

public:

    ExperimentRun();

    ExperimentRun(const std::string &run, std::string _traitSep = "   ", std::string _innerSep = ": ");

    void addTrait(const std::string &name, const std::string &value) {
        m_traits.emplace(name, value);
    }

    template<typename T>
    void addTrait(const std::string &name, const T &value) {
        std::stringstream ss;
        ss << value;
        addTrait(name, ss.str());
    }

    std::string getTrait(const std::string &name) const {
        return m_traits.at(name);
    }

    template<typename T>
    T getTrait(const std::string &name) const {
        std::istringstream ss(m_traits.at(name));
        T o;
        ss >> o;
        return o;
    }

    void updateTrait(const std::string &name, const std::string &value) {
        m_traits[name] = value;
    }

    template<typename T>
    void updateTrait(const std::string &name, const T &value) {
        std::stringstream ss;
        ss << value;
        updateTrait(name, ss.str());
    }

    const auto &traits() const {
        return m_traits;
    }

    void addMeasurement(const std::string &name, const tMeasurement &value) {
        m_measurements[name].emplace_back(value);
    }

    const auto &measurements() const {
        return m_measurements;
    }

    const auto &measurements(const std::string &name) const {
        return m_measurements.at(name);
    }

    tMeasurement avgMeasurement(const std::string &name) const;

    void incCounter(const std::string &name, const tMeasurement &value = 1) {
        m_counters[name] += value;
    }

    const auto &counters() const {
        return m_counters;
    }

    const auto &counter(const std::string &name) const {
        return m_counters.at(name);
    }

    std::string str(std::string _traitSep = "   ", std::string _innerSep = ": ") const;

    bool operator==(const ExperimentRun &o) const;

private:
    std::map<std::string, std::string> m_traits;
    tbb::concurrent_unordered_map<std::string, tMeasurements> m_measurements;
    tbb::concurrent_unordered_map<std::string, tMeasurement> m_counters;

private:
    static std::vector<std::string> c_ignored_fields;

};

class ExperimentRunAccessor : private boost::noncopyable {

#define PROFILER ExperimentRunAccessor::getInstance()

#ifdef ENABLE_PROFILING

#define PROFILER_INC(name) ExperimentRunAccessor::getInstance().incCounter(name)
#define PROFILER_ADD(name, value) ExperimentRunAccessor::getInstance().incCounter(name, value)

#else // ENABLE_PROFILING

#define PROFILER_INC(name) ((void)(0))
#define PROFILER_ADD(name, value) ((void)(0))

#endif

public:
    static ExperimentRunAccessor &getInstance();

    void addMeasurement(const std::string &name, const ExperimentRun::tMeasurement &value) {
        if(m_run)
            m_run->addMeasurement(name, value);
    }

    void incCounter(const std::string &name, const ExperimentRun::tMeasurement &value = 1) {
        if(m_run)
            m_run->incCounter(name, value);
    }

    void setRun(ExperimentRun *run){
        m_run = run;
    }

private:
    ExperimentRunAccessor() : m_run(nullptr) { }

    ExperimentRun * m_run;

};

std::ostream &operator<<(std::ostream &o, const ExperimentRun &p);