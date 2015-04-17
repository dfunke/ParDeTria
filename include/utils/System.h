#pragma once

#include <string>

#ifdef __cplusplus
extern "C" {
#endif

size_t getPeakRSS();
size_t getCurrentRSS();

#ifdef __cplusplus
}
#endif

std::string getHostname();
std::string getDatetime();
