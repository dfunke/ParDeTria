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

ExperimentRun::ExperimentRun(const std::string &run, std::string _traitSep, std::string _innerSep)
        : ExperimentRun()
{

    // iterate through string in search for traits
    std::string traitName, traitValue;
    std::size_t innerSep, traitSep, currPos = 0;

    while(currPos < run.length()){
        innerSep = run.find(_innerSep, currPos);
        traitSep = run.find(_traitSep, currPos);

        //we have the last trait in the string
        traitSep = traitSep != std::string::npos ? traitSep : run.length() + 1;

        traitName = run.substr(currPos, (innerSep - currPos));
        traitValue = run.substr(innerSep + _innerSep.length(), (traitSep-(innerSep + _innerSep.length())));

        if(   traitName != "datetime"
           && traitName != "host"
           && traitName != "git-rev")

            addTrait(traitName, traitValue);

        currPos = traitSep + _traitSep.length();
    }
}

std::string ExperimentRun::str(std::string _traitSep, std::string _innerSep) const {

    std::stringstream ss;

    std::string sep = "";
    for(const auto & trait : m_traits) {
        ss << sep << trait.first << _innerSep << trait.second;
        sep = _traitSep;
    }

    return ss.str();

}

tDuration ExperimentRun::avgTime() const {

    tDuration avg(0);

    for(const auto & t : m_times)
        avg += t;

    return avg / m_times.size();
}

std::size_t ExperimentRun::avgMem() const {

    std::size_t avg(0);

    for(const auto & t : m_mem)
        avg += t;

    return avg / m_mem.size();
}

bool ExperimentRun::operator==(const ExperimentRun &o) const {

    if(m_traits.size() != o.m_traits.size())
        return false;

    bool eq = true;
    for(const auto & t : m_traits) {

        if (t.first == "datetime"
            && t.first == "host"
            && t.first == "git-rev")

            //skip datetime host and git version
            continue;

        if (o.m_traits.count(t.first) == 0)
            // the other object does not contain this trait
            return false;
        else if (t.second != o.m_traits.at(t.first))
            eq = false;
    }


    return eq;

}

std::ostream &operator<<(std::ostream &o, const ExperimentRun &p){
    return o << p.str();
}