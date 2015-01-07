#include "Timings.h"

tMeas stats(const tDurations &meas) {

  // assert(meas.size() > 0);

  size_t N = 0;
  tMeas result;
  result.avg = 0, result.std = 0;
  double Mprev = 0;

  for (const auto &x : meas) {
    ++N;
    Mprev = result.avg;
    result.avg += (x.count() - Mprev) / N;
    result.std += (x.count() - Mprev) * (x.count() - result.avg);
  }

  if (N > 1)
    result.std = std::sqrt(result.std / (N - 1));

  return result;
}
