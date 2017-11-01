
#include "utils/Generator.h"
#include "Partitioner.h"
#include "LoadBalancedDCTriangulator.h"
#include "load_balancing/Partitioner.h"
#include "load_balancing/SimplePartitioner.h"
#include "load_balancing/OldPartitionerPartitioner.h"

int main()
{
    namespace lb = LoadBalancing;
    constexpr size_t d = 2;
    const size_t numPoints = 100;
    const dBox<d, double> boundingBox({0, 0}, {1, 1});

    tGenerator gen(START_SEED);

    BubblePointGenerator<d, double> generator;
    dPoints<d, double> pointCloud = generator.generate(numPoints, boundingBox, gen);

    std::cout << pointCloud.size() << "\n";
    for(const auto& point : pointCloud) {
        std::cout << point << "\n";
    }

    auto oldPartitioner = std::make_unique<dWayPartitioner<d, double>>();
    auto partitioner = std::make_unique<lb::OldPartitionerPartitioner<d, double>>(std::move(oldPartitioner), 10, 10);
    //auto partitioner = std::make_unique<lb::SimplePartitioner<d, double>>();
    lb::DCTriangulator<2, double> triangulator(boundingBox, pointCloud, std::move(partitioner), 1, false, false, true);
    auto simplices = triangulator.triangulate();
}
