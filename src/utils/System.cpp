//
// Created by dfunke on 4/9/15.
//

#include <unistd.h>
#include "System.h"


std::string getHostname() {
    char hostname[128];
    gethostname(hostname, sizeof(hostname));
    return std::string(hostname);
}