//
// Created by dfunke on 5/7/15.
//
#pragma once

#include <boost/noncopyable.hpp>
#include "utils/Timings.h"
#include "utils/Logger.h"


class TaskTimer : private boost::noncopyable {

public:

    TaskTimer(const std::string & name, int rank, bool autostart = true)
            : m_name(name), m_rank(rank), m_started(false) {

        if(autostart)
            start();
    }

    ~TaskTimer(){
        end();
    }

    void start() {
        m_started = true;
        m_startTime = Clock::now();
    }

    void end() {
        if(m_started) {
            auto end = Clock::now();
            tDuration duration = std::chrono::duration_cast<tDuration>(end - m_startTime);
            m_started = false;
            LOG("TIMING " << m_rank << " " << m_name << " " << m_startTime.time_since_epoch().count() << " " << end.time_since_epoch().count() << " " << duration.count() << std::endl);
        }
    }

private:
    std::string m_name;
    int m_rank;
    bool m_started;
    Clock::time_point m_startTime;

};

#ifdef ENABLE_TASK_TIMER

#define TIME_TASK(task) TaskTimer __taskTMR_##task(#task, this->world.rank())
#define TIME_PREP_TASK(task) TaskTimer __taskTMR_##task(#task, this->world.rank(), false)
#define TIME_START_TASK(task) __taskTMR_##task.start()
#define TIME_END_TASK(task) __taskTMR_##task.end()

#else // ENABLE_TASK_TIMER

#define TIME_TASK(task) ((void)(0))
#define TIME_PREP_TASK(task) ((void)(0))
#define TIME_START_TASK(task) ((void)(0))
#define TIME_END_TASK(task) ((void)(0))

#endif
