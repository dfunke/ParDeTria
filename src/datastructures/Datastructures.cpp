#include "datastructures/LP_Set.hxx"
#include "datastructures/Bit_Set.hxx"

//conversion constructors

LP_Set::LP_Set(Concurrent_LP_Set &&other)
        : m_arraySize(other.m_arraySize),
          m_items(other.m_items),
          m_array(reinterpret_cast<uint *>(other.m_array.release())),
          m_hasher(std::move(other.m_hasher)) { }

Concurrent_LP_Set::Concurrent_LP_Set(LP_Set &&other)
        : m_arraySize(other.m_arraySize),
          m_items(other.m_items),
          m_array(reinterpret_cast<std::atomic<uint> *>(other.m_array.release())),
          m_version(0),
          m_hasher(std::move(other.m_hasher)),
          m_currentCopyBlock(0) { }

Bit_Set::Bit_Set(Concurrent_Bit_Set &&other)
        : m_arraySize(other.m_arraySize),
          m_capacity(other.m_capacity),
          m_array(reinterpret_cast<uint *>(other.m_array.release())),
          m_ones(other.m_ones.load()) { }

Concurrent_Bit_Set::Concurrent_Bit_Set(Bit_Set &&other)
        : m_arraySize(other.m_arraySize),
          m_capacity(other.m_capacity),
          m_array(reinterpret_cast<std::atomic<uint> *>(other.m_array.release())),
          m_ones(other.m_ones) { }