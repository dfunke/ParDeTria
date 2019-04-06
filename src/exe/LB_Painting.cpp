#include <string>
#include <tuple>
#include <tbb/task_scheduler_init.h>
#include "Painter.h"
#include "utils/MakeIds.h"
#include "load_balancing/CommandLineInterface.h"
#include "load_balancing/PartitionTreePainter.h"
#include "load_balancing/monitors/DrawingMonitor.h"
#include "load_balancing/monitors/SampleDrawingMonitor.h"
#include "load_balancing/monitors/UnpluggedMonitor.h"

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
bool execute(const po::variables_map& vm, uint threads, const std::string& out,
             dPoints<D, Precision> /*points*/, bool verify) {
	startGen = tGenerator(vm.count("seed") ? vm["seed"].as<uint>() : START_SEED);
		
    dBox<D, Precision> bounds;
    bounds.low.fill(0);
    bounds.high.fill(100);

    dPoints<D, Precision> points;
	auto pg = createGenerator<D, Precision>(vm);
	auto N = vm["n"].as<tIdType>();
	points = pg->generate(N, bounds, startGen);

	std::vector<std::pair<std::string, Painter<2, double>>> painters;
	std::mutex mutex;
	lb::DrawingMonitor<Precision> monitor(points, colors.begin(), colors.end(), mutex,
	               [&painters, &bounds, &out](const auto& name) -> auto& {
	               painters.emplace_back(out + "_triang_" + name, Painter<2, double>(bounds));
	               return painters.back().second;
	               });

    auto partitioner = createPartitioner<D, Precision, decltype(monitor)>(vm, threads, startGen);
    assert(partitioner);

	using Triangulator = lb::DCTriangulator<D, Precision, decltype(monitor)>;
	Triangulator triangulator(bounds, points, threads,
	                          std::move(partitioner), 100,
	                          false,
	                          (partitioner->info().find("bipartitioner") == std::string::npos),
	                          true,
	                          std::move(monitor));

	auto dt = triangulator.triangulate();
	for(auto& painter : painters) {
		painter.second.save(painter.first);
	}


	if(verify) {
		CGALTriangulator<D, Precision, false> cgal(bounds, points);
		auto realDT = cgal.triangulate();

		auto vr = dt.verify(points);
		auto ccr = dt.crossCheck(realDT);

        return vr.valid && ccr.valid;
	} else {
		return true;
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
    cCommandLine.add_options()("print-tree", "just partition, do not triangulate");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cCommandLine), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << cCommandLine << std::endl;
        return EXIT_SUCCESS;
    }

    LOGGER.setLogLevel(static_cast<Logger::Verbosity>(0));

	if(vm.count("seed")) {
		startGen.seed(vm["seed"].as<uint>());
	}
    
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
    
	std::vector<std::pair<std::string, Painter<2, double>>> painters;
	std::mutex mutex;
	lb::SampleDrawingMonitor<Precision> monitor(points, colors.begin(), colors.end(), mutex,
	               [&painters, &bounds, &out=filename](const auto& name) -> auto& {
	               painters.emplace_back(out + "_sampling_" + name, Painter<2, double>(bounds));
	               return painters.back().second;
	               });

    auto partitioner = createPartitioner<D, Precision, decltype(monitor)>(vm, threads, startGen);
    assert(partitioner);

    tbb::task_scheduler_init init(threads);
    
    auto pointIds = makePointIds(points);
    auto partitioning = partitioner->partition(bounds, points, pointIds, monitor);

	if(vm.count("print-tree") > 0)
		printPartitionTree(partitioning, std::cout);
    
    Painter<2, double> pointsPainter(bounds);
    for(const auto& point : points) {
        pointsPainter.draw(point);
    }
   
	lb::DCTriangulator<D, Precision, decltype(monitor)> triangulator(bounds, points, threads,
																	  std::move(partitioner), 100,
																	  false, false, true,
																	  std::move(monitor));

    Painter<2, double> partitionPainter(bounds);
    PartitionTreePainter<double> paintTree(partitionPainter, points);
    paintTree(partitioning);
    
	if(vm.count("mode") == 0) {
		pointsPainter.save(filename + "_points");
		partitionPainter.save(filename + "_partitioning");
		for(auto& painter : painters) {
			painter.second.save(painter.first);
		}
	} else {
		auto mode = vm["mode"].as<std::string>();
		if("points" == mode)
			pointsPainter.save(filename + "_points");
		else if("sampling" == mode)
			for(auto& painter : painters) {
				painter.second.save(painter.first);
			}
		else if("partitioning" == mode)
			partitionPainter.save(filename + "_partitioning");
	}
     
	int result = 0;
	if(vm.count("no-triang") == 0)
	   result = execute<double>(vm, threads, filename, points, vm.count("validate") == 1) ? 0 : -1;

	return result;
}
