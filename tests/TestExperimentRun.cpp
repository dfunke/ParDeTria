#include "utils/Timings.h"
#include <gtest/gtest.h>

TEST(ExperimentRun, StoreLoadString) {

  ExperimentRun run1;

  run1.addTrait("alg", 'c');
  run1.addTrait("basecase", 500);

  std::string s = run1.str();

  ExperimentRun run2(s);

  EXPECT_EQ(run1, run2);
}
