#include <string>
#include <tuple>
#include <tbb/task_scheduler_init.h>
#include "Painter.h"
#include "utils/MakeIds.h"
#include "load_balancing/CommandLineInterface.h"

using Precision = double;
constexpr uint D = 2;

template <typename Precision>
struct PartitionTreePainter
{
    PartitionTreePainter(Painter<2, Precision>& painter, const dPoints<2, Precision>& points)
     : painter(&painter), points(&points)
    {}
    
    void operator()(const lb::PartitionTree<2, Precision>& tree) {
        size_t partitions = 0;
        paintPartitionTree(tree, colors.begin(), partitions);
        
        painter->setColor(tRGB(0, 0, 0));
        for(const auto label : labels)
            painter->drawText(std::get<0>(label), std::get<1>(label));
        labels.clear();
    }
    

    std::vector<tRGB>::iterator paintPartitionTree(const lb::PartitionTree<2, Precision>& tree,
                                                   std::vector<tRGB>::iterator current,
                                                   size_t& currentNumPartitions) {
        if(current == colors.end())
            current = colors.begin();
        
        if(tree.isLeaf()) {
            auto color = *current++;
            painter->setColor(color);
            dVector<2, Precision> center;
            size_t currentNumPoints = 0;
            for(auto id : std::get<Point_Ids>(tree.attachment)) {
                painter->draw((*points)[id]);
                center = center * ((Precision)currentNumPoints/(1 + currentNumPoints))
                    + (*points)[id].coords * ((Precision)1/(1 + currentNumPoints));
            }
            
            painter->drawBox(tree.bounds);
            
            labels.push_back({std::to_string(currentNumPartitions++), center});
        } else {
            const auto& children = std::get<typename lb::PartitionTree<2, Precision>::ChildContainer>(tree.attachment);
            for(const auto& child : children) {
                current = paintPartitionTree(child, current, currentNumPartitions);
            }
        }
        return current;
    }
    
private:
    Painter<2, Precision>* painter;
    const dPoints<2, Precision>* points;
    std::vector<std::tuple<std::string, dVector<2, Precision>>> labels;
    
    std::vector<tRGB> colors = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 1.0},
        {1.0, 0.0, 1.0},
    };
};

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
    
    Painter<2, double> partitionPainter(bounds);
    PartitionTreePainter<double> paintTree(partitionPainter, points);
    paintTree(partitioning);
    
    const auto* sp = dynamic_cast<const lb::SamplePartitioner<2, double>*>(partitioner.get());
    if(sp){
        partitionPainter.setColor(tRGB(0, 0, 0));
        const auto& samplePoints = sp->sampling().points;
        for(const auto& point : samplePoints) {
            partitionPainter.draw(point);
        }
    }
    
    partitionPainter.save(filename + "_partitioning");
}
