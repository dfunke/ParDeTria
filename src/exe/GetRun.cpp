// own
#include "utils/System.h"
#include "utils/ASSERT.h"
#include "utils/DBConnection.h"

DBConnection db("db_" + getHostname() + ".dat",
#ifndef ENABLE_PROFILING
                "benchmarks"
#else
                "profiling"
#endif
);

int main() {

    return db.getMaximum<uint>("run-number") + 1;
}
