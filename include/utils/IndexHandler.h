#pragma once

#include <atomic>
#include <queue>
#include <tbb/concurrent_queue.h>

namespace _detail {
    class IndexHandler {

    public:
        IndexHandler() : m_id(0) { }

        IndexHandler(const tIdType i) : m_id(i) { }

        tIdType getId() const { //needs to be const because of CGAL::Compact_container
            if(m_freeIds.empty()){
                return m_id++;
            } else {
                auto ret = m_freeIds.front();
                m_freeIds.pop();
                return ret;
            }
        }

        void releaseId(const tIdType & id) const { //needs to be const because of CGAL::Compact_container
            if(m_enabled)
                m_freeIds.push(id);
        }

        tIdType maxId() const {
            return m_id;
        }

        void disable() {
            m_enabled = false;
        }

    private:
        mutable tIdType m_id;
        mutable std::queue<tIdType> m_freeIds;
        bool m_enabled = true;

    };

    class Concurrent_IndexHandler {

    public:
        Concurrent_IndexHandler() : m_id(0) { }

        Concurrent_IndexHandler(const tIdType i) : m_id(i) { }

        tIdType getId() const { //needs to be const because of CGAL::Compact_container
            tIdType id;
            if(m_freeIds.try_pop(id))
                return id;
            else
                return m_id++;
        }

        void releaseId(const tIdType & id) const { //needs to be const because of CGAL::Compact_container
            if(m_enabled)
                m_freeIds.push(id);
        }

        tIdType maxId() const {
            return m_id.load();
        }

        void disable() {
            m_enabled = false;
        }

    private:
        mutable std::atomic<tIdType> m_id;
        mutable tbb::concurrent_queue<tIdType> m_freeIds;
        bool m_enabled = true;

    };
}