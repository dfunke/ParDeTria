#include <string>
#include <tuple>
#include <tbb/task_scheduler_init.h>
#include "Painter.h"
#include "utils/MakeIds.h"
#include "load_balancing/CommandLineInterface.h"
#include "load_balancing/GridIntersectionChecker.h"

using Precision = double;
constexpr uint D = 2;
    
const std::vector<tRGB> colors = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
    {1.0, 1.0, 0.0},
    {0.0, 1.0, 1.0},
    {1.0, 0.0, 1.0},
};

template <typename Precision>
struct PartitionTreePainter
{
    PartitionTreePainter(Painter<2, Precision>& painter, const dPoints<2, Precision>& points)
     : painter(&painter), points(&points)
    {}
    
    void operator()(const lb::PartitionTree<2, Precision>& tree) {
        size_t partitions = 0;
        paintPartitionTree(tree, colors.cbegin(), partitions);
        
        painter->setColor(tRGB(0, 0, 0));
        for(const auto label : labels)
            painter->drawText(std::get<0>(label), std::get<1>(label));
        labels.clear();
    }
    

    std::vector<tRGB>::const_iterator paintPartitionTree(const lb::PartitionTree<2, Precision>& tree,
                                                   std::vector<tRGB>::const_iterator current,
                                                   size_t& currentNumPartitions) {
        if(current == colors.cend())
            current = colors.cbegin();
        
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
            
			auto gridChecker = dynamic_cast<lb::GridIntersectionChecker<D, Precision>*>(tree.intersectionChecker.get());
			if(!gridChecker) {
	            painter->drawBox(tree.intersectionChecker->bounds());
			} else {
				for(auto it = gridChecker->cellsBegin(); it != gridChecker->cellsEnd(); ++it) {
					auto hasNeighbour = [it, gridChecker] (const dIndex<D, Precision>& direction) -> bool {
						auto cell = *it;
						std::transform(cell.begin(), cell.end(), direction.begin(), cell.begin(), [] (auto l, auto r) {
						               return l + r;
									   });
						return std::find(gridChecker->cellsBegin(), gridChecker->cellsEnd(), cell) != gridChecker->cellsEnd();
					};
					auto bounds = gridChecker->grid().boundsOf(*it);
					auto bottomLeft = bounds.low;
					auto topRight = bounds.high;
					auto bottomRight = bottomLeft;
					bottomRight[0] = topRight[0];
					auto topLeft = bottomLeft;
					topLeft[1] = topRight[1];
					if(hasNeighbour({0, 1})) {
						painter->drawLine(topLeft, topRight);
					}
					if(hasNeighbour({1, 0})) {
						painter->drawLine(topRight, bottomRight);
					}
					if(hasNeighbour({0, -1})) {
						painter->drawLine(bottomLeft, bottomRight);
					}
					if(hasNeighbour({-1, 0})) {
						painter->drawLine(topLeft, bottomLeft);
					}
				}
			}
            
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
    
    Painter<2, double> samplePainter(bounds);

    Painter<2, double> partitionPainter(bounds);
    PartitionTreePainter<double> paintTree(partitionPainter, points);
    paintTree(partitioning);
    
    const auto* sp = dynamic_cast<const lb::SamplePartitioner<2, double>*>(partitioner.get());
    if(sp){
	const auto& sampling = sp->sampling();
	const auto& graph = sampling.graph;
	for(size_t i = 0; i < graph.nodeRecords.size() - 1; ++i) {
	    for(int k = graph.nodeRecords[i]; k < graph.nodeRecords[i + 1]; ++k) {
		size_t j = graph.adjacency[k];
	    	samplePainter.setColor(tRGB(0, 0, 0));
		samplePainter.drawLine(sampling.points[i], sampling.points[j]);
	    }
	    samplePainter.setColor(colors[sampling.partition[i] % colors.size()]);
	    samplePainter.draw(sampling.points[i]);
	}

        partitionPainter.setColor(tRGB(0, 0, 0));
        const auto& samplePoints = sampling.points;
        for(const auto& point : samplePoints) {
	    partitionPainter.draw(point);
        }

	const auto* cdpasp = dynamic_cast<const lb::CenterDistancePointAssigningSamplePartitioner<2, double>*>(sp);
	if(cdpasp) {
		partitionPainter.setColor(tRGB(0.5, 0.5, 0.5));
		for(const auto& point : cdpasp->partitionCenters()) {
		    partitionPainter.draw(point);
		}
	}
    }
    
    samplePainter.save(filename + "_sampling");
    partitionPainter.save(filename + "_partitioning");
}
