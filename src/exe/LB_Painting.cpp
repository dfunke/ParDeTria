#include <string>
#include <tuple>
#include <tbb/task_scheduler_init.h>
#include "Painter.h"
#include "utils/MakeIds.h"
#include "load_balancing/CommandLineInterface.h"
#include "load_balancing/PartitionTreePainter.h"
#include "load_balancing/monitors/DrawingMonitor.h"

using Precision = double;
constexpr uint D = 2;
    

template <typename Precision, typename RandomIt>
void paintPartitionTree(const dPoints<2, Precision>& points,
                        const lb::PartitionTree<2, Precision>& tree,
                        Painter<2, Precision>& painter,
                        RandomIt first, RandomIt last) {
    paintPartitionTree(points, tree, painter, first, last, first);
}

template <typename Precision>
void execute(const po::variables_map& vm, uint threads, const std::string& out,
             dPoints<D, Precision> /*points*/) {
	startGen = tGenerator(START_SEED);
		
    dBox<D, Precision> bounds;
    bounds.low.fill(0);
    bounds.high.fill(100);

    dPoints<D, Precision> points;
	auto pg = createGenerator<D, Precision>(vm);
	auto N = vm["n"].as<tIdType>();
	points = pg->generate(N, bounds, startGen);

    auto partitioner = createPartitioner<D, Precision>(vm, threads, startGen);
    assert(partitioner);

	std::vector<std::pair<std::string, Painter<2, double>>> painters;
	std::mutex mutex;
	lb::DrawingMonitor<Precision> monitor(points, colors.begin(), colors.end(), mutex,
	               [&painters, &bounds, &out](const auto& name) -> auto& {
	               painters.emplace_back(out + "_triang_" + name, Painter<2, double>(bounds));
	               return painters.back().second;
	               });
   
	lb::DCTriangulator<D, Precision, decltype(monitor)> triangulator(bounds, points,
																	  std::move(partitioner), 100,
																	  false, false, true,
																	  std::move(monitor));
	triangulator.triangulate();
	for(auto& painter : painters) {
		painter.second.save(painter.first);
	}
}

int main(int argc, char *argv[]) {
    uint threads = tbb::task_scheduler_init::default_num_threads();
    //std::string pointFile;
    std::string partitionerName;

    auto cCommandLine = defaultOptions<Precision>();
    cCommandLine.add_options()("out", po::value<std::string>());
    cCommandLine.add_options()("threads", po::value(&threads), "specify number of threads");
    cCommandLine.add_options()("mode", po::value<std::string>());
    cCommandLine.add_options()("no-triang", "just partition, do not triangulate");

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
        /*const auto& samplePoints = sampling.points;
        for(const auto& point : samplePoints) {
	    partitionPainter.draw(point);
        }*/

	/*const auto* cdpasp = dynamic_cast<const lb::CenterDistancePointAssigningSamplePartitioner<2, double>*>(sp);
	if(cdpasp) {
		partitionPainter.setColor(tRGB(0.5, 0.5, 0.5));
		for(const auto& point : cdpasp->partitionCenters()) {
		    partitionPainter.draw(point);
		}
	}*/
    }
    
	if(vm.count("mode") == 0) {
		pointsPainter.save(filename + "_points");
		samplePainter.save(filename + "_sampling");
		partitionPainter.save(filename + "_partitioning");
	} else {
		auto mode = vm["mode"].as<std::string>();
		if("points" == mode)
			pointsPainter.save(filename + "_points");
		else if("sampling" == mode)
			samplePainter.save(filename + "_sampling");
		else if("partitioning" == mode)
			partitionPainter.save(filename + "_partitioning");
	}
     
	if(vm.count("no-triang") == 0)
	   execute<double>(vm, threads, filename, points);
}
