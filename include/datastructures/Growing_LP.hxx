#pragma once

#include <memory>
#include <mutex>

#include "LP_Set.hxx"

template<class HT>
class GrowingHashTableHandle;

template<class HT>
class ConstGrowingHashTableHandle;


template<class HT>
class GrowingHashTable {
    friend class GrowingHashTableHandle<HT>;

    friend class ConstGrowingHashTableHandle<HT>;

public:
    using Handle = GrowingHashTableHandle<HT>;

    template<typename... Args>
    GrowingHashTable(Args... args)
            : g_epoch_r(0), g_epoch_w(0), helperThreads(0) {
        g_table_r = std::make_shared<HT>(args...);
        g_table_w = g_table_r;
    }

    GrowingHashTable(GrowingHashTable<HT> &&other)
            : g_epoch_r(other.g_epoch_r.load()),
              g_table_r(other.g_table_r),
              g_epoch_w(other.g_epoch_w.load()),
              g_table_w(other.g_table_w),
              helperThreads(0) { }

    ~GrowingHashTable() = default;

    GrowingHashTable<HT> &operator=(GrowingHashTable<HT> &&other) {
        g_epoch_r.store(other.g_epoch_r.load());
        g_table_r = other.g_table_r;
        g_epoch_w.store(other.g_epoch_w.load());
        g_table_w = other.g_table_w;
        helperThreads.store(0);

        return *this;
    }

    uint get_epoch() { return g_epoch_r.load(); }

    GrowingHashTableHandle<HT> handle();

    ConstGrowingHashTableHandle<HT> handle() const;

private:
    using HashPtr = std::shared_ptr<HT>;

    std::atomic_uint g_epoch_r;
    HashPtr g_table_r;
    std::atomic_uint g_epoch_w;
    HashPtr g_table_w;

    size_t size() const {
        return g_table_r->size();
    }

    std::atomic_uint helperThreads;
    std::mutex grow_mutex;

};

template<class HT>
class ConstGrowingHashTableHandle {
    friend class GrowingHashTable<HT>;

public:

    ConstGrowingHashTableHandle() = delete;

    ConstGrowingHashTableHandle(const GrowingHashTable<HT> &parent_) : parent(parent_) {
        getTable();
    }

    template<typename... Args>
    bool contains(Args... args) const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->contains(args...);
    }

    template<typename... Args>
    auto get(Args... args) const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->get(args...);
    }

    template<typename... Args>
    std::size_t count(Args... args) const {
        return contains(args...);
    }

    bool empty() const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->empty();
    };

    bool empty(const std::size_t idx) const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->empty(idx);
    };

    auto at(const std::size_t idx) const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->at(idx);
    }

    auto capacity() const {
        uint t_epoch = parent.g_epoch_w.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->capacity();
    }

    auto size() const {
        uint t_epoch = parent.g_epoch_w.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->size();
    }

    auto begin() const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->begin();
    }

    auto end() const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->end();
    }

    auto range() const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->range();
    }

private:
    void getTable() const {
        l_epoch = parent.g_epoch_r.load();
        l_table = parent.g_table_r;

        if (l_table->m_version != l_epoch) {
            getTable();
        }
    }

private:
    using HashPtr = std::shared_ptr<HT>;

    const GrowingHashTable<HT> &parent;
    mutable uint l_epoch;
    mutable HashPtr l_table;

    const static uint c_BlockSize = 256;
};

template<class HT>
class GrowingHashTableHandle {
    friend class GrowingHashTable<HT>;

public:

    GrowingHashTableHandle() = delete;

    GrowingHashTableHandle(GrowingHashTable<HT> &parent_) : parent(parent_) {
        getTable();
    }

    template<typename... Args>
    InsertReturn insert(Args... args) {
        uint t_epoch = parent.g_epoch_w.load();
        if (t_epoch > l_epoch) {
            getTable();
            if (t_epoch > l_epoch) {
                grow();
            }
        }
        auto result = l_table->insert(args...);
        if (result == InsertReturn::State::Full) {
            grow();
            return insert(args...);
        } else
            return result;
    }

    template<typename... Args>
    bool contains(Args... args) const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->contains(args...);
    }

    template<typename... Args>
    auto get(Args... args) const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->get(args...);
    }

    template<typename... Args>
    std::size_t count(Args... args) const {
        return contains(args...);
    }

    bool empty() const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->empty();
    };

    bool empty(const std::size_t idx) const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->empty(idx);
    };

    auto at(const std::size_t idx) const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->at(idx);
    }

    auto capacity() const {
        uint t_epoch = parent.g_epoch_w.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->capacity();
    }

    auto size() const {
        uint t_epoch = parent.g_epoch_w.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->size();
    }

    template<class Source, class Filter>
    void unsafe_merge(Source &&other, const Filter &filter) {
        uint t_epoch = parent.g_epoch_w.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->unsafe_merge(std::move(other), filter);
    }

    auto begin() const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->begin();
    }

    auto end() const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->end();
    }

    auto range() const {
        uint t_epoch = parent.g_epoch_r.load();
        if (t_epoch > l_epoch) {
            getTable();
        }
        return l_table->range();
    }

private:
    void getTable() const {
        l_epoch = parent.g_epoch_r.load();
        l_table = parent.g_table_r;

        if (l_table->m_version != l_epoch) {
            getTable();
        }
    }

    void grow() {
        //generate table
        //CAS table into position

        //std::cout << "Growing " << typeid(HT).name() <<std::endl;

        HashPtr w_table;
        { //should be atomic
            std::lock_guard<std::mutex> lock(parent.grow_mutex);
            if (parent.g_table_w->m_version == l_table->m_version) {
                auto temp = std::make_shared<HT>(l_table->capacity() << 1,
                                                 l_table->m_version + 1);
                parent.g_table_w = temp;
            }
            w_table = parent.g_table_w;
        }

        //CAS epoch_w up
        auto t_epoch = l_table->m_version;
        parent.g_epoch_w.compare_exchange_strong(t_epoch, w_table->m_version);

        //nhelper ++
        parent.helperThreads.fetch_add(1);

        //get block + while block legal migrate and get new block
        auto begin = l_table->m_currentCopyBlock.fetch_add(c_BlockSize);
        while (begin < l_table->capacity()) {

            auto end = std::min(begin + c_BlockSize, l_table->capacity());
            for (auto i = begin; i < end; ++i) {
                l_table->_migrate(i, *w_table);
            }

            begin = l_table->m_currentCopyBlock.fetch_add(c_BlockSize);
        }

        //nhelper --
        parent.helperThreads.fetch_sub(1);

        //wait for other helpers
        while (parent.helperThreads.load(std::memory_order_acquire))
        {};

        //CAS table into R-Position
        {
            std::lock_guard<std::mutex> lock(parent.grow_mutex);
            if (parent.g_table_r->m_version < w_table->m_version) {
                parent.g_table_r = w_table;
            }
        }

        //CAS epoch_r up
        t_epoch = l_table->m_version;
        parent.g_epoch_r.compare_exchange_strong(t_epoch, w_table->m_version);
        getTable();
    }

private:
    using HashPtr = std::shared_ptr<HT>;

    GrowingHashTable<HT> &parent;
    mutable uint l_epoch;
    mutable HashPtr l_table;

    const static uint c_BlockSize = 256;
};

template<class HT>
GrowingHashTableHandle<HT> GrowingHashTable<HT>::handle() {
    return GrowingHashTableHandle<HT>(*this);
}

template<class HT>
ConstGrowingHashTableHandle<HT> GrowingHashTable<HT>::handle() const {
    return ConstGrowingHashTableHandle<HT>(*this);
}