#include "utils/Timings.h"
#include "utils/StringUtils.h"
#include <gtest/gtest.h>

TEST(ExperimentRun, StoreLoadString) {

  ExperimentRun run1;

  run1.addTrait("alg", 'c');
  run1.addTrait("basecase", 500);

  std::string s = run1.str();

  ExperimentRun run2(s);

  EXPECT_EQ(run1, run2);
}

TEST(ExperimentRun, NumberRecognition) {
  EXPECT_TRUE(is_int("124"));
  EXPECT_TRUE(is_int("1"));
  EXPECT_TRUE(is_int("0"));
  EXPECT_TRUE(is_int("+230"));
  EXPECT_TRUE(is_int("-230"));

  EXPECT_FALSE(is_int(""));
  EXPECT_FALSE(is_int("24.3"));
  EXPECT_FALSE(is_int("24a"));

  EXPECT_TRUE(is_float("0"));
  EXPECT_TRUE(is_float("1"));
  EXPECT_TRUE(is_float("0.0"));
  EXPECT_TRUE(is_float("1.0"));
  EXPECT_TRUE(is_float("234"));
  EXPECT_TRUE(is_float("+234"));
  EXPECT_TRUE(is_float("-234"));

  EXPECT_TRUE(is_float("234.3"));
  EXPECT_TRUE(is_float("+234.2"));
  EXPECT_TRUE(is_float("-234.3"));

  EXPECT_TRUE(is_float("1.4e3"));
  EXPECT_TRUE(is_float("3.5E3"));

  EXPECT_TRUE(is_float(".23432"));
  EXPECT_TRUE(is_float("-.234"));

  EXPECT_FALSE(is_float("234a"));
  EXPECT_FALSE(is_float("234.3.2"));
  EXPECT_FALSE(is_float("234..234"));
  EXPECT_FALSE(is_float(""));

}