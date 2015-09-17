//
// Created by dfunke on 5/26/15.
//

#include "datastructures/LP_Set.hxx"
#include "datastructures/LP_Map.hxx"
#include "datastructures/LP_MultiMap.hxx"
#include "datastructures/BlockedArray.hxx"

#include "utils/TBB_Containers.h"
#include <tbb/task_scheduler_init.h>

#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <csignal>

#include "utils/Timings.h"

const uint N = 1e7;
const uint R = 10;

typedef std::vector<uint> tVals;

template<typename T>
tDuration timeOps(const tVals & vals){

    T set(vals.size());
    tbb::blocked_range<std::size_t> range(std::size_t(0), vals.size(), 1e4);

    auto t1 = Clock::now();
    tbb::parallel_for(range, [&set, &vals] (const auto & r){
        for(auto i = r.begin(); i != r.end(); ++i){
            set.insert(vals[i]);
        }
    });
    auto t2 = Clock::now();

    auto t3 = Clock::now();
    tbb::parallel_for(range, [&set, &vals] (const auto & r){
        for(uint k = 0; k < R; ++k) {
            for (auto i = r.begin(); i != r.end(); ++i) {
                set.count(vals[i]);
            }
        }
    });
    auto t4 = Clock::now();

    return (t2 - t1) + (t4 - t3);
}

int main() {

    tVals vals;
    std::unordered_set<uint> cmp;

    auto distKey = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
    auto diceKey = std::bind(distKey, startGen);

    while(cmp.size() < N){
        uint val = diceKey();

        if(cmp.insert(val).second)
            vals.push_back(val);
    }

    std::ofstream f("benchmark_set.csv");
    f << "threads tbb own" << std::endl;

    for(int t = 1; t <= tbb::task_scheduler_init::default_num_threads(); ++t){
        tbb::task_scheduler_init init(t);

        std::cout << "Testing with " << t << " threads: " << std::flush;

        tDuration tTBB(0);
        tDuration tOwn(0);
        for(uint r = 0; r < R; ++r){
            tTBB += timeOps<tbb::concurrent_unordered_set<uint>>(vals);
            tOwn += timeOps<Concurrent_LP_Set<uint>>(vals);
            std::cout << "." << std::flush;
        }
        tTBB /= R; tOwn /= R;

        std::cout << std::endl;

        f << t << " " << tTBB.count() << " " << tOwn.count() << std::endl;
    }

}