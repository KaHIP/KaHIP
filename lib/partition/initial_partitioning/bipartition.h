/******************************************************************************
 * bipartition.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BIPARTITION_7I4IR31Y
#define BIPARTITION_7I4IR31Y

#include "initial_partitioner.h"

class bipartition : public initial_partitioner {
        public:
                bipartition();
                virtual ~bipartition();

                void initial_partition( const PartitionConfig & config, 
                                        const unsigned int seed, 
                                        graph_access & G, 
                                        int* partition_map); 

                void initial_partition( const PartitionConfig & config, 
                                const unsigned int seed,  
                                graph_access & G, 
                                int* xadj,
                                int* adjncy, 
                                int* vwgt, 
                                int* adjwgt,
                                int* partition_map); 

        private:
                void grow_regions_bfs(const PartitionConfig & config, graph_access & G);
                void grow_regions_fm(const PartitionConfig & config, graph_access & G);
                NodeID find_start_node( const PartitionConfig & config, graph_access & G);
                void post_fm(const PartitionConfig & config, graph_access & G);
                inline Gain compute_gain( graph_access & G, NodeID node, PartitionID targeting_partition);

};

inline Gain bipartition::compute_gain( graph_access & G, NodeID node, PartitionID targeting_partition) {
        //compute how connected is the target to the current partition
        Gain gain = 0;
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                if( G.getPartitionIndex(target) == targeting_partition) {
                        gain += G.getEdgeWeight(e); 
                } else {
                        gain -= G.getEdgeWeight(e); 
                }
        } endfor

        return gain;
}


#endif /* end of include guard: BIPARTITION_7I4IR31Y */
