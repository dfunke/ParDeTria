#include <atomic>
#include <mutex>
#include <boost/noncopyable.hpp>

class SpinMutex : private boost::noncopyable {

public:
  void lock() {
    while (m_flag.test_and_set(std::memory_order_acquire)) // acquire lock
      ;                                                    // spin
  }

  bool try_lock() { return m_flag.test_and_set(std::memory_order_acquire); }

  void unlock() {
    m_flag.clear(std::memory_order_release); // release lock
  }

private:
  std::atomic_flag m_flag = ATOMIC_FLAG_INIT;
};