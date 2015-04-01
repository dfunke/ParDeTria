#include <boost/noncopyable.hpp>
#include <iostream>
#include <atomic>
#include <mutex>

//  progress_display  --------------------------------------------------------//
//  based on boost::progress_display but thread-safe

//  progress_display displays an appropriate indication of
//  progress at an appropriate place in an appropriate form.

// NOTE: (Jan 12, 2001) Tried to change unsigned long to boost::uintmax_t, but
// found some compilers couldn't handle the required conversion to double.
// Reverted to unsigned long until the compilers catch up.

class ProgressDisplay : private boost::noncopyable {
public:
  explicit ProgressDisplay(unsigned long expected_count_,
                           std::ostream &os = std::cout,
                           const std::string &s1 = "\n", // leading strings
                           const std::string &s2 = "",
                           const std::string &s3 = "")
      // os is hint; implementation may ignore, particularly in embedded systems
      : noncopyable(),
        m_os(os),
        m_s1(s1),
        m_s2(s2),
        m_s3(s3) {
    restart(expected_count_);
  }

  void restart(unsigned long expected_count_)
  //  Effects: display appropriate scale
  //  Postconditions: count()==0, expected_count()==expected_count_
  {
    _count = _next_tic_count = _tic = 0;
    _expected_count = expected_count_;

    m_os << m_s1 << "0%   10   20   30   40   50   60   70   80   90   100%\n"
         << m_s2 << "|----|----|----|----|----|----|----|----|----|----|"
         << std::endl // endl implies flush, which ensures display
         << m_s3;
    if (!_expected_count)
      _expected_count = 1; // prevent divide by zero
  }                        // restart

  unsigned long operator+=(unsigned long increment)
  //  Effects: Display appropriate progress tic if needed.
  //  Postconditions: count()== original count() + increment
  //  Returns: count().
  {
    while (_mtx.test_and_set(std::memory_order_acquire)) // acquire lock
      ;                                                  // spin
    if ((_count += increment) >= _next_tic_count) {
      display_tic();
    }

    _mtx.clear(std::memory_order_release); // release lock

    return _count;
  }

  unsigned long operator++() { return operator+=(1); }
  unsigned long count() const { return _count; }
  unsigned long expected_count() const { return _expected_count; }

private:
  std::ostream &m_os;     // may not be present in all imps
  const std::string m_s1; // string is more general, safer than
  const std::string m_s2; //  const char *, and efficiency or size are
  const std::string m_s3; //  not issues

  unsigned long _count, _expected_count, _next_tic_count;
  unsigned int _tic;
  std::atomic_flag _mtx = ATOMIC_FLAG_INIT;

  void display_tic() {
    // use of floating point ensures that both large and small counts
    // work correctly.  static_cast<>() is also used several places
    // to suppress spurious compiler warnings.
    unsigned int tics_needed = static_cast<unsigned int>(
        (static_cast<double>(_count) / _expected_count) * 50.0);
    do {
      m_os << '*' << std::flush;
    } while (++_tic < tics_needed);
    _next_tic_count =
        static_cast<unsigned long>((_tic / 50.0) * _expected_count);
    if (_count == _expected_count) {
      if (_tic < 51)
        m_os << '*';
      m_os << std::endl;
    }
  } // display_tic
};