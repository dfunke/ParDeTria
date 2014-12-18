/*
 * utils.h
 *
 *  Created on: May 24, 2013
 *      Author: dfunke
 */

#pragma once

#include <sstream>
#include <fstream>

class CSV {
public:
  static long double nsToMs(long double ns) { return ns * 1E-6; }

  static long double sToNs(long double ns) { return ns * 1E9; }

  template <class... Types> static std::string csv(const Types &... args) {
    std::stringstream ss;
    _csv(ss, args...);
    std::string s = ss.str();

    // ignore the last separator in the string
    return s.substr(0, s.length() - 1);
  }

private:
  static void _csv(__attribute((unused)) std::stringstream &ss) {}

  template <typename T, class... Types>
  static void _csv(std::stringstream &ss, const T &value,
                   const Types &... args) {
    ss << value << SEP;
    _csv(ss, args...);
  }

public:
  static constexpr char SEP = ' ';
};
