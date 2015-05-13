#include <atomic>
#include <functional>
#include <cstring>

#include "Misc.h"
#include "ASSERT.h"

class LP_Set {

  public:
    typedef uint tKeyType;

  public:

    LP_Set(std::size_t size)
    : m_items(0),
      m_array(nullptr)
    {
      // Initialize cells
      m_arraySize = nextPow2(size);
      m_array = new std::atomic<tKeyType>[m_arraySize];
      std::fill(m_array, m_array + m_arraySize, 0);
    }

    ~LP_Set(){
      if(m_array)
        delete[] m_array;
    }

    bool insert(const tKeyType & key){
      ASSERT(key != 0);

      tKeyType startIdx = m_hasher(key) & (m_arraySize - 1);
      for (tKeyType idx = startIdx; idx != ((startIdx - 1) & (m_arraySize - 1)); ++idx) {
        idx &= m_arraySize - 1;
        ASSERT(idx < m_arraySize);

        // Load the key that was there.
        tKeyType probedKey = m_array[idx].load(std::memory_order_relaxed);

        if (probedKey == key)
          return false; // the key is already in the set, return false;
        else {
          // The entry was either free, or contains another key.
          if (probedKey != 0)
            continue; // Usually, it contains another key. Keep probing.
          // The entry was free. Now let's try to take it using a CAS.
          tKeyType prevKey = 0;
          bool cas = m_array[idx].compare_exchange_strong(prevKey, key, std::memory_order_relaxed);
          if (cas){
            ++m_items;
            return true; // we just added the key to the set
          }
          else if(prevKey == key)
            return false; // the key was already added by another thread
          else
            continue; // another thread inserted a different key in this position
        }
      }

      return false;
    }

    bool contains(const tKeyType & key){
      ASSERT(key != 0);
      for (tKeyType idx = m_hasher(key);; idx++) {
        idx &= m_arraySize - 1;
        tKeyType probedKey = m_array[idx].load(std::memory_order_relaxed);
        if (probedKey == key)
          return true;;
        if (probedKey == 0)
          return false;
      }

    }

    auto capacity() const { return m_arraySize; }
    auto size() const { return m_items.load(); }

    void merge(LP_Set && other){

      std::size_t oldSize = m_arraySize;
      m_arraySize = nextPow2(capacity() + other.capacity());

      std::atomic<tKeyType> * oldArray = m_array;

      m_array = new std::atomic<tKeyType>[m_arraySize];
      std::fill(m_array, m_array + m_arraySize, 0);

      for(std::size_t i = 0; i < oldSize; ++i){
        m_array[2*i + 0].store(oldArray[i].load(std::memory_order_relaxed), std::memory_order_relaxed);
      }

      for(std::size_t i = 0; i < other.capacity(); ++i){
        m_array[2*i + 1].store(other.m_array[i].load(std::memory_order_relaxed), std::memory_order_relaxed);
      }

      delete[] oldArray;
    }


  private:
    std::size_t m_arraySize;
    std::atomic<std::size_t> m_items;
    std::atomic<tKeyType> * m_array;
    std::hash<tKeyType> m_hasher;

};
