/*
 * Random.h
 *
 *  Created on: Oct 28, 2014
 *      Author: dfunke
 */

#pragma once

#include <random>
#include <sstream>
#include <functional>

//#############################################################################

const uint START_SEED = 1986;
typedef std::mt19937 tGenerator;
extern tGenerator startGen;

/* Prototype
 * std::uniform_int_distribution<uint> distribution(0,100);
 * auto dice = std::bind ( distribution, generator );
 */
