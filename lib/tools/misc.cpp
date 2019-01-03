/******************************************************************************
 * misc.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "misc.h"
#include "quality_metrics.h"

misc::misc() {

}

misc::~misc() {

}

void misc::balance_singletons(const PartitionConfig & config, graph_access & G) {
        quality_metrics qm;
        std::vector< NodeID > singletons;
        std::vector< NodeWeight > block_sizes(config.k,0);

        forall_nodes(G, node) {
                block_sizes[G.getPartitionIndex(node)] += G.getNodeWeight(node);

                if(G.getNodeDegree(node) == 0) {
                        singletons.push_back(node);
                }
        } endfor

        // use buckets?
        for( unsigned i = 0; i < singletons.size(); i++) {
                NodeWeight min = block_sizes[0];
                PartitionID p  = 0;
                for( unsigned j = 0; j < config.k; j++) {
                        if( block_sizes[j] < min ) {
                                min = block_sizes[j];
                                p   = j;
                        }
                }

                NodeID node = singletons[i];
                block_sizes[G.getPartitionIndex(node)] -= G.getNodeWeight(node);
                block_sizes[p] += G.getNodeWeight(node);
                G.setPartitionIndex(node, p);
        }
        std::cout <<  "log> balance after assigning singletons " <<  qm.balance(G)  << std::endl;
}
