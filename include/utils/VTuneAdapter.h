//
// Created by dfunke on 5/7/15.
//
#pragma once

#ifdef ENABLE_VTUNE

#include <boost/noncopyable.hpp>
#include <ittnotify.h>
#include <string>
#include <unordered_map>


class VTuneDomain : private boost::noncopyable {

public:
    static VTuneDomain MAIN;

public:
    VTuneDomain(const std::string & domain){
        m_name = domain;
        m_domain = __itt_domain_create(m_name.c_str());
    }

    __itt_domain * handle(){
        return m_domain;
    }

    const std::string & name() const {
        return m_name;
    }

private:
    __itt_domain * m_domain;
    std::string m_name;

};

class VTuneStringHandle : private boost::noncopyable {

public:
    VTuneStringHandle(const std::string & string){
        m_string = string;
        m_handle = __itt_string_handle_create(m_string.c_str());
    }

    __itt_string_handle * handle(){
        return m_handle;
    }

    const std::string & str() const {
        return m_string;
    }

private:
    __itt_string_handle * m_handle;
    std::string m_string;

};

class VTuneTask : private boost::noncopyable {

public:

    VTuneTask(VTuneDomain & domain, VTuneStringHandle & name, bool autostart = true)
            : m_domain(domain), m_name(name), m_running(false) {

        if(autostart)
            start();
    }

    ~VTuneTask(){
        end();
    }

    void start() {
        __itt_task_begin(m_domain.handle(), __itt_null, __itt_null, m_name.handle());
        m_running = true;
    }

    void end() {
        if(m_running) {
            __itt_task_end(m_domain.handle());
            m_running = false;
        }
    }

private:
    VTuneDomain & m_domain;
    VTuneStringHandle & m_name;
    bool m_running;

};

#define VTUNE_TASK(task) static VTuneStringHandle __vtuneStr_##task(#task); \
                         VTuneTask __vtuneTsk_##task(VTuneDomain::MAIN, \
                                                     __vtuneStr_##task)
#define VTUNE_PREP_TASK(task) static VTuneStringHandle __vtuneStr_##task(#task); \
                              VTuneTask __vtuneTsk_##task(VTuneDomain::MAIN, \
                                                          __vtuneStr_##task, false)
#define VTUNE_START_TASK(task) __vtuneTsk_##task.start()
#define VTUNE_END_TASK(task) __vtuneTsk_##task.end()

#else // ENABLE_VTUNE

#define VTUNE_TASK(task) ((void)(0))
#define VTUNE_PREP_TASK(task) ((void)(0))
#define VTUNE_START_TASK(task) ((void)(0))
#define VTUNE_END_TASK(task) ((void)(0))

#endif
