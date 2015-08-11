//
// Author: Christian Schulz <christian.schulz@kit.edu>
// 

#include <fstream>
#include "initial_node_separator.h"
#include "graph_partitioner.h"
#include "tools/quality_metrics.h"
#include "tools/random_functions.h"
#include "partition/uncoarsening/separator/vertex_separator_algorithm.h"
#include "partition/uncoarsening/refinement/node_separators/fm_ns_local_search.h"


initial_node_separator::initial_node_separator() {
                
}

initial_node_separator::~initial_node_separator() {
                
}

NodeWeight initial_node_separator::single_run( const PartitionConfig & config, graph_access & G) {
        std::cout <<  "IP run"  << std::endl;

        std::streambuf* backup = std::cout.rdbuf();
        std::ofstream ofs;
        ofs.open("/dev/null");
        std::cout.rdbuf(ofs.rdbuf()); 

        graph_partitioner partitioner;
        PartitionConfig partition_config         = config;
        partition_config.mode_node_separators    = false;
        partition_config.global_cycle_iterations = 1;
        partition_config.repetitions             = 1;
        partition_config.edge_rating             = config.sep_edge_rating_during_ip;


        //if( rnd == 0 ) partition_config.edge_rating = SEPARATOR_MULTX;
        //if( rnd == 1 ) partition_config.edge_rating = SEPARATOR_ADDX;
        //if( rnd == 2 ) partition_config.edge_rating = SEPARATOR_MAX;
        //if( rnd == 3 ) partition_config.edge_rating = SEPARATOR_LOG;

        //computing a partition
        partitioner.perform_partitioning(partition_config, G);

        complete_boundary boundary(&G);
        boundary.build();

        ofs.close();
        std::cout.rdbuf(backup);

        vertex_separator_algorithm vsa; std::vector<NodeID> separator;
        //create a very simple separator from that partition
        if( partition_config.sep_full_boundary_ip ) {
                vsa.compute_vertex_separator_simpler(partition_config, G, boundary, separator);
        } else {
                vsa.compute_vertex_separator_simple(partition_config, G, boundary, separator);
        }
 
	//fm_ns_local_search fmnsls;
	//std::vector< NodeWeight > block_weights(3,0);
	//forall_nodes(G, node) {
		//if( G.getPartitionIndex(node) == 0) {
			//block_weights[0] += G.getNodeWeight(node);
		//} else if( G.getPartitionIndex(node) == 1 ) {
			//block_weights[1] += G.getNodeWeight(node);
		//} else {
			//block_weights[2] += G.getNodeWeight(node);
		//}
	//} endfor

	//if(block_weights[0] > block_weights[1]) {
		//fmnsls.perform_refinement(partition_config, G, true, 1);
	//} else {
		//fmnsls.perform_refinement(partition_config, G, true, 0);
	//}
	//fmnsls.perform_refinement(partition_config, G);
	//separator.clear();
	//forall_nodes(G, node) {
		//if(G.getPartitionIndex(node) == 2) {
			//separator.push_back(node);
		//}
	//} endfor

        std::vector<NodeID> output_separator;
        //improve the separator using flow based techniques
        vsa.improve_vertex_separator(partition_config, G, separator, output_separator);

        quality_metrics qm;
        return qm.separator_weight(G);

}

void initial_node_separator::compute_node_separator( const PartitionConfig & config, graph_access & G) {
        if(config.graph_allready_partitioned) return;

        std::vector< NodeID > best_separator(G.number_of_nodes(),0);
        NodeWeight best_separator_size = std::numeric_limits< NodeWeight >::max();

        for( int i  = 0; i  < config.max_initial_ns_tries; i++) {
                NodeWeight cur_separator_size = single_run(config, G);
                if(cur_separator_size < best_separator_size) {
                        forall_nodes(G, node) {
                                best_separator[node] = G.getPartitionIndex(node);
                        } endfor
                        best_separator_size = cur_separator_size;
                
                        std::cout <<  "improved initial separator size " <<  cur_separator_size  << std::endl;
                }
        }

        forall_nodes(G, node) {
                G.setPartitionIndex(node, best_separator[node]);
        } endfor
        
}
