/*
 * Timings.h
 *
 *  Created on: Oct 28, 2014
 *      Author: dfunke
 */

#pragma once

#include <chrono>
#include <vector>
#include <cmath>

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::nanoseconds tDuration;
typedef std::vector<tDuration> tDurations;

struct tMeas {
  double avg, std;
};

typedef std::vector<tMeas> tMeass;

tMeas stats(const tDurations &meas);
