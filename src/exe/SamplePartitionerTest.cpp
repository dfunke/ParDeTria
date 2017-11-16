
#include "utils/Generator.h"
#include "utils/MakeIds.h"
#include "LoadBalancedDCTriangulator.h"
#include "load_balancing/Partitioner.h"
#include "load_balancing/sample_partitioner/BoundsDistancePointAssigningSamplePartitioner.h"
#include "load_balancing/OldPartitionerPartitioner.h"
#include "load_balancing/PartitionTreeDump.h"

int main()
{
    namespace lb = LoadBalancing;
    constexpr size_t d = 2;
    //const size_t numPoints = 100;
    tGenerator gen(START_SEED);
    
    const dBox<d, double> boundingBox({0, 0}, {1, 1});
    
    dPoints<d, double> pointCloud;
    pointCloud.reserveUpToIdx(4*4 + 1);
    pointCloud[1] = dPoint<d, double>{{0.5, 0.5}};
    pointCloud[1+4] = dPoint<d, double>{{0.4, 0.4}};
    pointCloud[1+8] = dPoint<d, double>{{0.4, 0.6}};
    pointCloud[1+12] = dPoint<d, double>{{0.6, 0.4}};
    
    pointCloud[2] = dPoint<d, double>{{0.5, -0.5}};
    pointCloud[2+4] = dPoint<d, double>{{0.4, -0.4}};
    pointCloud[2+8] = dPoint<d, double>{{0.4, -0.6}};
    pointCloud[2+12] = dPoint<d, double>{{0.6, -0.4}};
    
    pointCloud[3] = dPoint<d, double>{{-0.5, 0.5}};
    pointCloud[3+4] = dPoint<d, double>{{-0.4, 0.4}};
    pointCloud[3+8] = dPoint<d, double>{{-0.4, 0.6}};
    pointCloud[3+12] = dPoint<d, double>{{-0.6, 0.4}};
    
    pointCloud[4] = dPoint<d, double>{{-0.5, -0.5}};
    pointCloud[4+4] = dPoint<d, double>{{-0.4, -0.4}};
    pointCloud[4+8] = dPoint<d, double>{{-0.4, -0.6}};
    pointCloud[4+12] = dPoint<d, double>{{-0.6, -0.4}};
    auto ids = makePointIds(pointCloud);

    std::cout << pointCloud.size() << "\n";
    for(const auto& point : pointCloud) {
        std::cout << point << "\n";
    }

    auto partitioner = std::make_unique<lb::BoundsDistancePointAssigningSamplePartitioner<d, double>>(8, rand(), 2);
    auto partitionTree = partitioner->partition(boundingBox, pointCloud, ids);
    lb::dumpPartition(partitionTree);
    
    lb::DCTriangulator<2, double> triangulator(boundingBox, pointCloud, std::move(partitioner), 1, false, false, true);
    auto simplices = triangulator.triangulate();
}
