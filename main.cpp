#include <pthread.h>
#include <time.h>
#include "graph.hpp"
#include "parser.hpp"
#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>
#include "vec2d.hpp"
#include "global.hpp"
#include "order.hpp"

using namespace std;
using namespace boost::timer;
using namespace boost::program_options;

#define TEST
#undef TEST

std::string getFileExtension(const std::string& fileName) {
    size_t dotPos = fileName.find_last_of(".");
    if (dotPos != std::string::npos && dotPos < fileName.length() - 1) {
        return fileName.substr(dotPos + 1);
    }
    return ""; 
}

//////////////////////////////////////////
//main function
//////////////////////////////////////////
int main(int argc, char** argv)
{
    options_description desc{"Options"};
    desc.add_options()
        ("data,d", value<string>()->default_value(""), "input_graph")
        ("in_feat,i", value<string>()->default_value(""), "input_feat")
        ("output,o", value<string>()->default_value(""), "output_graph")
        ("out_reorder,r", value<string>()->default_value(""), "output_vtx_map")
        ("size,s", value<int>()->default_value(512), "part_size")
        ("algorithm,a", value<int>()->default_value(0), "reorder_ID")
        ("kv,k", value<int>()->default_value(16), "K value");


    variables_map vm;
    try
    {
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);
    } catch (error& e) {
        cerr << "ERROR: " << e.what() << '\n' << '\n' << desc << '\n';
        return 1;
    }

    string data_file = vm["data"].as<string>();
    string write_file = vm["output"].as<string>();
    string reorder_file = vm["out_reorder"].as<string>();
    string feat_file = vm["in_feat"].as<string>();

    int size = vm["size"].as<int>();
    int reorder = vm["algorithm"].as<int>();
    int feat = 15;
    int kv = vm["kv"].as<int>();
    params::num_threads = omp_get_max_threads() / 2;

    if (data_file.empty()) {
        cout << desc << '\n';
        exit(1);
    } 
    bool if_out_reorder = false;
    if(reorder_file != "")
        if_out_reorder = true;

    Algo algo = static_cast<Algo>(reorder);
    // graph object
    Graph graph;

    /**************************************************
     Compute the preprocessing time
     *************************************************/
    // cpu_timer timer;
    cpu_times times;

    graph.cluster_num = kv;
    if(parseGraph(data_file, graph)) {
        // times = timer.elapsed();
        // cout << times.wall/(1e9) << "s: parse done for " << data_file << '\n';
        graph.printGraph();
    }
    
    int part_size = size * 1024 / sizeof(float);
    //int part_size = size;
    int NUM_THREADS = params::num_threads;
    int numVerticesPerBin= (graph.num_vertex/(NUM_THREADS*4));
    params::partition_size = (numVerticesPerBin < part_size) ? numVerticesPerBin : part_size;
    //params::partition_size = part_size;
    int pow2=1;
    while(pow2<=params::partition_size)
        pow2*=2;
    pow2/=2;
    if(pow2==0) params::partition_size=4;
    else params::partition_size = pow2;

    
    //////////////////////////////////////////
    // read csr file
    //////////////////////////////////////////
  //  writeEdgelist("edgelist.txt", graph);
    params::num_partitions = (graph.num_vertex-1)/params::partition_size + 1;
    printf("[Graph Info]\n");
    std::cout << "num vertex: " << graph.num_vertex << ", num edges: " << graph.num_edges << ", num partitions: " << params::num_partitions << '\n';
    //////////////////////////////////////////
    // output Degree array
    //////////////////////////////////////////
    graph.computeOutDegree();
    // times = timer.elapsed();
    // cout << times.wall/(1e9) << "s: outdegree is computed "  << '\n';


    if(feat_file != "" && !feat_file.empty()){
        graph.in_feat = feat_file;
        graph.initAttributeFile(feat_file, feat);
    }
    // TODO: add a logic to handle missing feature file
    Orderer orderer(&graph);
    
    ///////////////////////////////
    // re-order the graph 
    //////////////////////////////
    std::cout << "==========" << std::endl;
    cpu_timer reorder_timer;
    float order_time = 0.0;

    orderer.reorder(algo);
    cout << "reordering latency = "  << reorder_timer.elapsed().wall/(1e9) << '\n';

    orderer.getNewGraph(algo);
    // cout << times.wall/(1e9) << "s: a new graph is constructed, total time is: "  << reorder_timer.elapsed().wall/(1e9) << '\n';
    
    if(if_out_reorder){
        writeReorder(reorder_file, orderer, graph.num_edges, algo);
    }
    
    if(!write_file.empty()){
        writeGraph(write_file, graph);
    }

}
