#include <string>
#include <tbb/task_scheduler_init.h>
#include "Painter.h"
#include "utils/MakeIds.h"
#include "load_balancing/CommandLineInterface.h"

using Precision = double;
constexpr uint D = 2;

template <typename Precision, typename RandomIt>
RandomIt paintPartitionTree(const dPoints<2, Precision>& points,
                        const lb::PartitionTree<2, Precision>& tree,
                        Painter<2, Precision>& painter,
                        RandomIt first, RandomIt last, RandomIt current) {
    if(current == last)
        current = first;
    if(tree.isLeaf()) {
        painter.setColor(*current++);
        for(auto id : std::get<Point_Ids>(tree.attachment)) {
            painter.draw(points[id]);
        }
    } else {
        const auto& children = std::get<typename lb::PartitionTree<2, Precision>::ChildContainer>(tree.attachment);
        for(const auto& child : children) {
            current = paintPartitionTree(points, child, painter, first, last, current);
        }
    }
    return current;
}

template <typename Precision, typename RandomIt>
void paintPartitionTree(const dPoints<2, Precision>& points,
                        const lb::PartitionTree<2, Precision>& tree,
                        Painter<2, Precision>& painter,
                        RandomIt first, RandomIt last) {
    paintPartitionTree(points, tree, painter, first, last, first);
}

int main(int argc, char *argv[]) {
    uint threads = tbb::task_scheduler_init::default_num_threads();
    //std::string pointFile;    
    std::string partitionerName;

    auto cCommandLine = defaultOptions<Precision>();
    cCommandLine.add_options()("out", po::value<std::string>());
    cCommandLine.add_options()("threads", po::value(&threads), "specify number of threads");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << cCommandLine << std::endl;
        return EXIT_SUCCESS;
    }

    LOGGER.setLogLevel(static_cast<Logger::Verbosity>(0));
    
    if (!(vm.count("n") || vm.count("points"))) {
        std::cout << "Please specify number of points or point file" << std::endl;
        return EXIT_FAILURE;
    }
    
    if(!vm.count("partitioner")) {
        std::cout << "Please specify a partitioner" << std::endl;
        return EXIT_FAILURE;
    }
    
    if(!vm.count("distribution")) {
        std::cout << "Please specify a distribution" << std::endl;
        return EXIT_FAILURE;
    }
    
    auto filename = vm["out"].as<std::string>();
    
    dBox<D, Precision> bounds;
    bounds.low.fill(0);
    bounds.high.fill(100);

    dPoints<D, Precision> points;
    /*if (vm.count("points")) {
        points = loadObject<dPoints<D, Precision>>(pointFile);
    } else {*/
        auto pg = createGenerator<D, Precision>(vm);
        auto N = vm["n"].as<tIdType>();
        points = pg->generate(N, bounds, startGen);
    //}
    
    auto partitioner = createPartitioner<D, Precision>(vm, threads, startGen);
    assert(partitioner);

    tbb::task_scheduler_init init(threads);
    
    auto pointIds = makePointIds(points);
    auto partitioning = partitioner->partition(bounds, points, pointIds);
    
    Painter<2, double> pointsPainter(bounds);
    for(const auto& point : points) {
        pointsPainter.draw(point);
    }
    pointsPainter.save(filename + "_points");
    
    std::vector<tRGB> colors = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 1.0},
        {1.0, 0.0, 1.0},
    };
    
    Painter<2, double> partitionPainter(bounds);
    paintPartitionTree(points, partitioning, partitionPainter, colors.begin(), colors.end());
    partitionPainter.save(filename + "_partitioning");
}
