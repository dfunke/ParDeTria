#pragma once

#include <vector>
#include <memory>
#include <atomic>

#include "utils/TBB_Containers.h"

template<typename T, std::size_t BS = 512>
class BlockedArray {

private:
    typedef std::unique_ptr<T[]> block_ptr;

public:
    BlockedArray(const std::size_t &initial_capacity = BS)
            : m_nBlocks(0) {

        reserve(initial_capacity);
    }

    BlockedArray(BlockedArray && other){
        m_nBlocks = other.m_nBlocks;
        m_blocks = std::move(other.m_blocks);
    }

    void reserve(const std::size_t &capacity) {

        while (m_nBlocks * BS < capacity) {
            m_blocks.emplace_back(new T[BS]);
            ++m_nBlocks;
        }
    }

    const T& operator[](const std::size_t idx) const {
        return m_blocks[idx / BS][idx % BS];
    }

    T& operator[](const std::size_t idx) {
        return m_blocks[idx / BS][idx % BS];
    }

    const T& at(const std::size_t idx) const {
        if(idx > m_nBlocks * BS)
            throw std::out_of_range("index: " + std::to_string(idx) + " size: " + std::to_string(m_nBlocks * BS));
        return m_blocks[idx / BS][idx % BS];
    }

    T& at(const std::size_t idx) {
        if(idx > m_nBlocks * BS)
            throw std::out_of_range("index: " + std::to_string(idx) + " size: " + std::to_string(m_nBlocks * BS));
        return m_blocks[idx / BS][idx % BS];
    }

private:
    std::size_t m_nBlocks;
    std::vector<block_ptr> m_blocks;

};

template<typename T, std::size_t BS = 512>
class Concurrent_BlockedArray {

private:
    typedef std::unique_ptr<T[]> block_ptr;

public:
    Concurrent_BlockedArray(const std::size_t &initial_capacity = BS)
            : m_nBlocks(0) {

        reserve(initial_capacity);
    }

    Concurrent_BlockedArray(Concurrent_BlockedArray && other){
        m_nBlocks.store(other.m_nBlocks.load());
        m_blocks = std::move(other.m_blocks);
    }

    void reserve(const std::size_t &capacity) {

        while (m_nBlocks.load() * BS < capacity) {
            m_blocks.emplace_back(new T[BS]);
            ++m_nBlocks;
        }
    }

    const T& operator[](const std::size_t idx) const {
        return m_blocks[idx / BS][idx % BS];
    }

    T& operator[](const std::size_t idx) {
        return m_blocks[idx / BS][idx % BS];
    }

    const T& at(const std::size_t idx) const {
        if(idx > m_nBlocks * BS)
            throw std::out_of_range("index: " + std::to_string(idx) + " size: " + std::to_string(m_nBlocks * BS));
        return m_blocks[idx / BS][idx % BS];
    }

    T& at(const std::size_t idx) {
        if(idx > m_nBlocks * BS)
            throw std::out_of_range("index: " + std::to_string(idx) + " size: " + std::to_string(m_nBlocks * BS));
        return m_blocks[idx / BS][idx % BS];
    }

private:
    std::atomic<std::size_t> m_nBlocks;
    tbb::concurrent_vector<block_ptr> m_blocks;

};