#include "Timings.h"

#include <time.h>

#include "System.h"

#include "version.h"

ExperimentRun::ExperimentRun() {

    // output datetime, hostname and git commit

    char datetime[64];
    time_t tnow = time(NULL);
    strftime(datetime,sizeof(datetime),"%Y-%m-%d %H:%M:%S", localtime(&tnow));
    addTrait("datetime", std::string(datetime));

    addTrait("host", getHostname());

    addTrait("git-rev", std::string(GIT_COMMIT));

}

std::string ExperimentRun::str(std::string _sep) const {

    std::stringstream ss;

    std::string sep = "";
    for(const auto & trait : m_traits) {
        ss << sep << trait.first << ": " << trait.second;
        sep = _sep;
    }

    return ss.str();

}

tDuration ExperimentRun::avgTime() const {

    tDuration avg(0);

    for(const auto & t : m_times)
        avg += t;

    return avg / m_times.size();
}