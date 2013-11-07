//
// Author: Christian Schulz <christian.schulz@kit.edu>
// 


#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string.h> 

#include <argtable2.h>
#include <regex.h>

#include "../lib/data_structure/graph_access.h"
#include "../lib/io/graph_io.h"
#include "../lib/tools/timer.h"
#include "../lib/tools/quality_metrics.h"
#include "../lib/tools/macros_assertions.h"
#include "../lib/tools/random_functions.h"
#include "../lib/partition/partition_config.h"
#include "../lib/partition/graph_partitioner.h"
#include "../lib/partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"

#include "parse_parameters.h"
#include <math.h>

using namespace std;
int main(int argn, char **argv) {

        PartitionConfig partition_config;
        string graph_filename;
        bool is_graph_weighted = false;
        bool suppress_output = false;
        bool recursive = false;
       
        int ret_code = parse_parameters(argn, argv, partition_config, graph_filename, is_graph_weighted, suppress_output, recursive); 

        if(ret_code) {
                return 0;
        }

        streambuf* backup = cout.rdbuf();
        ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
               cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.LogDump(stdout);

        graph_access G;     
        timer t;
        if(!is_graph_weighted) {
                graph_io::readGraphUnweighted(G, graph_filename);
        } else {
                graph_io::readGraphWeighted(G, graph_filename);
        }
        cout << "io time: " << t.elapsed()  << endl;
       
        G.set_partition_count(partition_config.k); 
 
        NodeWeight largest_graph_weight = 0;
        forall_nodes(G, node) {
                largest_graph_weight += G.getNodeWeight(node);
        } endfor
        
        double epsilon = 0;
        if( partition_config.kaffpa_perfectly_balance ) {
                epsilon                              = (partition_config.imbalance+1)/100.0;
        } else {
                epsilon                              = (partition_config.imbalance)/100.0;
        }
        partition_config.upper_bound_partition      = (1+epsilon)*ceil(largest_graph_weight/(double)partition_config.k);
        cout <<  "upper bound " <<  partition_config.upper_bound_partition  << endl;
        partition_config.largest_graph_weight       = largest_graph_weight;
        partition_config.graph_allready_partitioned = false;
        partition_config.kway_adaptive_limits_beta  = log(largest_graph_weight);

       
        //double epsilon                              = partition_config.imbalance/100.0;
        //partition_config.upper_bound_partition      = ceil((1+epsilon)*G.number_of_nodes()/(double)partition_config.k);
        //partition_config.largest_graph_weight       = G.number_of_nodes();
        //partition_config.graph_allready_partitioned = false;
        //partition_config.kway_adaptive_limits_beta  = log(G.number_of_nodes());

        if(partition_config.input_partition != "") {
                std::cout <<  "reading input partition" << std::endl;
                graph_io::readPartition(G, partition_config.input_partition);
                partition_config.graph_allready_partitioned  = true;
                partition_config.no_new_initial_partitioning = true;
        }

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << endl;
        // ***************************** perform partitioning ***************************************       
        graph_partitioner partitioner;
        quality_metrics qm;
        
        unsigned tries = 10;
        double sum_time = 0;
        double sum_weight = 0;
        EdgeWeight cut = std::numeric_limits<EdgeWeight>::max();
        for( unsigned i = 0; i < tries; i++) {
                                partition_config.graph_allready_partitioned = false;
                                forall_nodes(G, node) {
                                        G.setPartitionIndex(node,0);
                                } endfor
                                
                                t.restart();
                                partitioner.perform_partitioning(partition_config, G);
                                sum_time += t.elapsed();
                                EdgeWeight cur_cut = qm.edge_cut(G);
                                sum_weight += cur_cut;

                                if( cur_cut < cut ) {
                                        cut = cur_cut;
                                }
        }
        // ******************************* done partitioning *****************************************       
        ofs.close();
        cout.rdbuf(backup);
        cout <<  "time spent for partitioning " << t.elapsed()  << endl;
       
        // output some information about the partition that we have computed 
        cout << "cut \t\t"   << qm.edge_cut(G)       << endl;
        cout << "finalobjective  "   << qm.edge_cut(G)       << endl;
        cout << "avgobjective  "   << setprecision(20)   <<  sum_weight/(double)tries      << endl;
        cout << "avgtime  "   << setprecision(20)<< sum_time/(double)tries      << endl;
        cout << "minimumcut "   << setprecision(20)<< cut      << endl;
        cout << "bnd \t\t"   << qm.boundary_nodes(G) << endl;
        cout << "balance \t" << qm.balance(G)        << endl;
        cout << "finalbalance \t" << qm.balance(G)        << endl;
        cout << "max_comm_vol \t" << qm.max_communication_volume(G)        << endl;

        // write the partition to the disc 
        string partition("tmppartition");
        stringstream noparts;
        noparts << "tmppartition" << partition_config.k;
        graph_io::writePartition(G, noparts.str());
        
}
