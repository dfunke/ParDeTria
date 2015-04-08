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
#include <map>

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::nanoseconds tDuration;

class ExperimentRun
{

public:

  ExperimentRun();

  void addTrait(const std::string & name, const std::string & value){
    m_traits.emplace(name, value);
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

private:
  std::map<std::string, std::string> m_traits;
  std::vector<tDuration> m_times;
  std::vector<std::size_t> m_mem;

};
