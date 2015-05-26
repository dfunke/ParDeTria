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

const uint N = 1e6;
const uint R = 10;

struct tOp {
    char op;
    uint val;
    bool res;
};

typedef std::vector<tOp> tOps;

template<typename T>
tDuration timeOps(const tOps & ops){

    T set(ops.size());

    auto t1 = Clock::now();
    tbb::parallel_for(tbb::blocked_range<std::size_t>(std::size_t(0), ops.size()), [&set, &ops] (const auto & r){
        for(auto i = r.begin(); i != r.end(); ++i){
            const tOp & op = ops[i];

            switch(op.op){
                case 0:
                    set.insert(op.val);
                    break;
                case 1:
                    set.count(op.val);
                    break;
            }
        }
    });
    auto t2 = Clock::now();

    return t2 - t1;
}

int main() {

    tOps ops;
    std::unordered_set<uint> cmp;

    auto distKey = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
    auto diceKey = std::bind(distKey, startGen);

    auto distOp = std::uniform_int_distribution<char>(0, 1);
    auto diceOp = std::bind(distOp, startGen);

    for (uint i = 0; i < N; ++i) {
        tOp op;
        op.op = diceOp();
        op.val = diceKey();

        switch(op.op){
            case 0:
                op.res = cmp.insert(op.val).second; break;
            case 1:
                op.res = cmp.count(op.val); break; //only meaningful in sequential case
        }

        ops.push_back(op);
    }

    std::ofstream f("benchmark_set.csv");
    f << "threads tbb own" << std::endl;

    for(int t = 1; t <= tbb::task_scheduler_init::default_num_threads(); ++t){
        tbb::task_scheduler_init init(t);

        std::cout << "Testing with " << t << " threads: " << std::flush;

        tDuration tTBB(0);
        tDuration tOwn(0);
        for(uint r = 0; r < R; ++r){
            tTBB += timeOps<tbb::concurrent_unordered_set<uint>>(ops);
            tOwn += timeOps<Concurrent_LP_Set>(ops);
            std::cout << "." << std::flush;
        }
        tTBB /= R; tOwn /= R;

        std::cout << std::endl;

        f << t << " " << tTBB.count() << " " << tOwn.count() << std::endl;
    }

}