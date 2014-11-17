/******************************************************************************
 * bipartition.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
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
