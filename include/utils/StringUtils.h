#pragma once

#include <iostream>
#include <regex>
#include <string>

bool is_int(const std::string & s){
    return !s.empty() &&
           ((s.size() == 1 && std::isdigit(s[0]))
           || std::regex_match(s, std::regex("[(-|+)|][0-9]+")));
}

bool is_float(const std::string & s){
    return !s.empty() && (
           (s.size() == 1 && std::isdigit(s[0]))
           || std::regex_match(s, std::regex("[(-|+)|](([1-9][0-9]*\\.?[0-9]*)|(\\.[0-9]+))([Ee][+-]?[0-9]+)?")));
}