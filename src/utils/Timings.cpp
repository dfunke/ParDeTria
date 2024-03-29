#include "Timings.h"

#include <algorithm>

#include "System.h"
#include "version.h"

#define IGNORED_FIELDS {"git-rev", "start-time", "end-time", "run-number", "host"}
std::vector<std::string> ExperimentRun::c_ignored_fields(IGNORED_FIELDS);

ExperimentRun::ExperimentRun() {

    // hostname and git commit

    addTrait("host", getHostname());
    addTrait("git-rev", std::string(GIT_COMMIT));
}

ExperimentRun::ExperimentRun(const std::string &run, std::string _traitSep, std::string _innerSep)
        : ExperimentRun() {

    // iterate through string in search for traits
    std::string traitName, traitValue;
    std::size_t innerSep, traitSep, currPos = 0;

    while (currPos < run.length()) {
        innerSep = run.find(_innerSep, currPos);
        traitSep = run.find(_traitSep, currPos);

        //we have the last trait in the string
        traitSep = traitSep != std::string::npos ? traitSep : run.length() + 1;

        traitName = run.substr(currPos, (innerSep - currPos));
        traitValue = run.substr(innerSep + _innerSep.length(), (traitSep - (innerSep + _innerSep.length())));

        if (std::find(c_ignored_fields.begin(), c_ignored_fields.end(), traitName) == c_ignored_fields.end())
            //trait not in ignored fields
            addTrait(traitName, traitValue);

        currPos = traitSep + _traitSep.length();
    }
}

std::string ExperimentRun::str(std::string _traitSep, std::string _innerSep, bool printMeas) const {

    std::stringstream ss;

    std::string sep = "";
    for (const auto &trait : m_traits) {
        ss << sep << trait.first << _innerSep << trait.second;
        sep = _traitSep;
    }

    if(printMeas){
        for (const auto &meas : m_measurements){
            ss << sep << meas.first << _innerSep << avgMeasurement(meas.first);
        }
    }

    return ss.str();

}

ExperimentRun::tMeasurement ExperimentRun::avgMeasurement(const std::string &name) const {

    tMeasurement avg = 0;

    for (const auto &m : m_measurements.at(name))
        avg += m;

    if (m_measurements.at(name).size() > 0)
        avg /= m_measurements.at(name).size();

    return avg;
}

bool ExperimentRun::operator==(const ExperimentRun &o) const {

    if (m_traits.size() != o.m_traits.size())
        return false;

    bool eq = true;
    for (const auto &t : m_traits) {

        if (std::find(c_ignored_fields.begin(), c_ignored_fields.end(), t.first) != c_ignored_fields.end())
            //trait is ignored field
            continue;

        if (o.m_traits.count(t.first) == 0)
            // the other object does not contain this trait
            return false;
        else if (t.second != o.m_traits.at(t.first))
            eq = false;
    }


    return eq;

}

std::ostream &operator<<(std::ostream &o, const ExperimentRun &p) {
    return o << p.str();
}

ExperimentRunAccessor &ExperimentRunAccessor::getInstance() {
    static ExperimentRunAccessor instance; // Guaranteed to be destroyed and thread safe
    // Instantiated on first use.
    return instance;
}