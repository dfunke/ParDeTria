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

std::string print_stacktrace(const std::string & file, const int line, const uint max_frames = 63);

#define ALLOCATE(size) \
try { \
    _allocate(size); \
} catch (std::bad_alloc &e) { \
    std::string s = e.what() \
                    + std::string(" in ") \
                    + typeid(*this).name() \
                    + " for size " + std::to_string(size) + "\n" \
                    + "current/peak RSS: " + std::to_string(getCurrentRSS() / 1e6) + "/" + std::to_string(getPeakRSS() / 1e6) + " MB\n" \
                    + print_stacktrace(__FILE__, __LINE__); \
    std::cerr << s << std::endl; \
    throw e; \
}

