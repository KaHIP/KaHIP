//
// Author: Christian Schulz <christian.schulz@kit.edu>
// 

#include "initial_node_separator.h"
#include "graph_partitioner.h"
#include "tools/quality_metrics.h"
#include "partition/uncoarsening/separator/vertex_separator_algorithm.h"


initial_node_separator::initial_node_separator() {
                
}

initial_node_separator::~initial_node_separator() {
                
}

void initial_node_separator::compute_node_separator( const PartitionConfig & config, graph_access & G) {
        graph_partitioner partitioner;
        quality_metrics qm;

        PartitionConfig partition_config = config;
        std::cout <<  "initially computing a node separator"  << std::endl;
        partition_config.mode_node_separators = false;
        partitioner.perform_partitioning(partition_config, G);

        complete_boundary boundary(&G);
        boundary.build();

        std::cout <<  "now transforming "  << std::endl;
        vertex_separator_algorithm vsa;
        std::vector<NodeID> separator;
        vsa.compute_vertex_separator_simple(partition_config, G, boundary, separator);
        vsa.is_vertex_separator(G, separator);
        std::cout <<  "separator size " <<  separator.size()  << std::endl;
        for( NodeID v : separator) {
                std::cout <<  "separator " <<  v  << std::endl;
        }
        std::cout <<  "now improving"  << std::endl;
        
        std::vector<NodeID> output_separator;
        vsa.improve_vertex_separator(partition_config, G, separator, output_separator);
        vsa.is_vertex_separator(G, output_separator);
        std::cout <<  "here"  << std::endl;

}
