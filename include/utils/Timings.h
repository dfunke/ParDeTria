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

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::nanoseconds tDuration;

class ExperimentRun {

public:
    typedef ulong tMeasurement;

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

    void addMeasurement(const std::string &name, const tMeasurement &value) {
        m_measurements[name].emplace_back(value);
    }

    const auto &traits() const {
        return m_traits;
    }

    const auto &measurements() const {
        return m_measurements;
    }

    const auto &measurements(const std::string &name) const {
        return m_measurements.at(name);
    }

    tMeasurement avgMeasurement(const std::string &name) const;

    std::string str(std::string _traitSep = "   ", std::string _innerSep = ": ") const;

    bool operator==(const ExperimentRun &o) const;

private:
    std::map<std::string, std::string> m_traits;
    std::map<std::string, std::vector<tMeasurement>> m_measurements;

private:
    static std::vector<std::string> c_ignored_fields;

};

std::ostream &operator<<(std::ostream &o, const ExperimentRun &p);