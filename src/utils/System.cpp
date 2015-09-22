//
// Created by dfunke on 4/9/15.
//

#include <unistd.h>
#include "System.h"
#include "Random.h"

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