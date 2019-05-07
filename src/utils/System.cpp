//
// Created by dfunke on 4/9/15.
//

#include <unistd.h>
#include "System.h"
#include "Random.h"

#include <sstream>

tGenerator startGen(START_SEED);


std::string getHostname() {
    char hostname[128];
    gethostname(hostname, sizeof(hostname));
    return std::string(hostname);
}

std::string getDatetime() {
    char datetime[64];
    time_t tnow = time(NULL);
    strftime(datetime, sizeof(datetime), "%Y-%m-%d %H:%M:%S", localtime(&tnow));

    return std::string(datetime);
}

// stacktrace.h (c) 2008, Timo Bingmann from http://idlebox.net/
// published under the WTFPL v2.0

#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>
#include <cxxabi.h>

/** Print a demangled stack backtrace of the caller function to FILE* out. */
std::string print_stacktrace(const std::string & file, const int line, const uint max_frames)
{

    std::stringstream ss;

    ss << "stack trace from (" << file <<":" << line << "):\n";

    // storage array for stack trace address data
    void** addrlist = new void*[max_frames+1];

    // retrieve current stack addresses
    int addrlen = backtrace(addrlist, max_frames+1);

    if (addrlen == 0) {
        ss << "\t<empty, possibly corrupt>\n";
        return ss.str();
    }

    // resolve addresses into strings containing "filename(function+address)",
    // this array must be free()-ed
    char** symbollist = backtrace_symbols(addrlist, addrlen);

    // allocate string which will be filled with the demangled function name
    size_t funcnamesize = 256;
    char* funcname = (char*)malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for (int i = 1; i < addrlen; i++)
    {
        char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

        // find parentheses and +address offset surrounding the mangled name:
        // ./module(function+0x15c) [0x8048a6d]
        for (char *p = symbollist[i]; *p; ++p)
        {
            if (*p == '(')
                begin_name = p;
            else if (*p == '+')
                begin_offset = p;
            else if (*p == ')' && begin_offset) {
                end_offset = p;
                break;
            }
        }

        if (begin_name && begin_offset && end_offset
            && begin_name < begin_offset)
        {
            *begin_name++ = '\0';
            *begin_offset++ = '\0';
            *end_offset = '\0';

            // mangled name is now in [begin_name, begin_offset) and caller
            // offset in [begin_offset, end_offset). now apply
            // __cxa_demangle():

            int status;
            char* ret = abi::__cxa_demangle(begin_name,
                                            funcname, &funcnamesize, &status);
            if (status == 0) {
                funcname = ret; // use possibly realloc()-ed string
                ss << "\t" << symbollist[i] << " : " << funcname << "+" << begin_offset << "\n";
            }
            else {
                // demangling failed. Output function name as a C function with
                // no arguments.
                ss << "\t" << symbollist[i] << " : " << begin_name << "()+" << begin_offset << "\n";
            }
        }
        else
        {
            // couldn't parse the line? print the whole line.
            ss << "\t " << symbollist[i] << "\n";
        }
    }

    free(funcname);
    free(symbollist);
    delete[] addrlist;

    return ss.str();
}

