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

class ExperimentRun
{

public:

  ExperimentRun();

  ExperimentRun(const std::string & run, std::string _traitSep = "   ", std::string _innerSep = ": ");

  void addTrait(const std::string & name, const std::string & value){
    m_traits.emplace(name, value);
  }

  template<typename T>
  void addTrait(const std::string & name, const T & value){
    std::stringstream ss;
    ss << value;
    addTrait(name, ss.str());
  }

  std::string getTrait(const std::string & name) const {
    return m_traits.at(name);
  }

  template<typename T>
  T getTrait(const std::string & name) const {
    std::istringstream ss(m_traits.at(name));
    T o;
    ss >> o;
    return o;
  }

  void updateTrait(const std::string & name, const std::string & value){
    m_traits[name] = value;
  }

  template<typename T>
  void updateTrait(const std::string & name, const T & value){
    std::stringstream ss;
    ss << value;
    updateTrait(name, ss.str());
  }

  void addTime(const tDuration & t){
    m_times.emplace_back(t);
  }

  void addMemory(const std::size_t & m){
    m_mem.emplace_back(m);
  }

  const auto & traits() const {
    return m_traits;
  }

  const auto & times() const {
    return m_times;
  }

  const auto & mem() const {
    return m_mem;
  }

  tDuration avgTime() const;
  std::size_t avgMem() const;

  std::string str(std::string _traitSep = "   ", std::string _innerSep = ": ") const;

  bool operator ==(const ExperimentRun & o) const;

private:
  std::map<std::string, std::string> m_traits;
  std::vector<tDuration> m_times;
  std::vector<std::size_t> m_mem;

private:
  static std::vector<std::string> c_ignored_fields;

};

std::ostream &operator<<(std::ostream &o, const ExperimentRun &p);