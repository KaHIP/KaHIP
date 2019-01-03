/******************************************************************************
 * quality_metrics.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>
#include <cmath>

#include "quality_metrics.h"
#include "data_structure/union_find.h"

#include <unordered_map>

quality_metrics::quality_metrics() {
}

quality_metrics::~quality_metrics () {
}

EdgeWeight quality_metrics::edge_cut(graph_access & G) {
        EdgeWeight edgeCut = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = G.getPartitionIndex(n);
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);

                        if (partitionIDSource != partitionIDTarget) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
        return edgeCut/2;
}

EdgeWeight quality_metrics::edge_cut(graph_access & G, int * partition_map) {
        EdgeWeight edgeCut = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = partition_map[n];
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = partition_map[targetNode];

                        if (partitionIDSource != partitionIDTarget) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
        return edgeCut/2;
}

EdgeWeight quality_metrics::edge_cut(graph_access & G, PartitionID lhs, PartitionID rhs) {
        EdgeWeight edgeCut = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = G.getPartitionIndex(n);
                if(partitionIDSource != lhs) continue;
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);

                        if(partitionIDTarget == rhs) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                } endfor 
        } endfor
        return edgeCut;
}

EdgeWeight quality_metrics::edge_cut_connected(graph_access & G, int * partition_map) {
        EdgeWeight edgeCut = 0;
        EdgeWeight sumEW = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = partition_map[n];
                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = partition_map[targetNode];

                        if (partitionIDSource != partitionIDTarget) {
                                edgeCut += G.getEdgeWeight(e);
                        }
                        sumEW+=G.getEdgeWeight(e);
                } endfor 
        } endfor
        union_find uf(G.number_of_nodes());
        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] == partition_map[target]) {
                                uf.Union(node, target); 
                        }
                } endfor
        } endfor

        std::unordered_map<NodeID, NodeID> size_right;
        forall_nodes(G, node) {
                size_right[uf.Find(node)] = 1;
        } endfor


        std::cout <<  "number of connected comp " <<  size_right.size()  << std::endl;
        if( size_right.size() == G.get_partition_count()) {
                return edgeCut/2;
        } else {
                return edgeCut/2+sumEW*size_right.size();
        }

}


EdgeWeight quality_metrics::max_communication_volume(graph_access & G, int * partition_map) {
    std::vector<EdgeWeight> block_volume(G.get_partition_count(),0);
    forall_nodes(G, node) {
        PartitionID block = partition_map[node];
        std::vector<bool> block_incident(G.get_partition_count(), false);
        block_incident[block] = true;

        int num_incident_blocks = 0;

        forall_out_edges(G, e, node) {
            NodeID target = G.getEdgeTarget(e);
            PartitionID target_block = partition_map[target];
            if(!block_incident[target_block]) {
                block_incident[target_block] = true;
                num_incident_blocks++;
            }
        } endfor
        block_volume[block] += num_incident_blocks;
    } endfor

    EdgeWeight max_comm_volume = *(std::max_element(block_volume.begin(), block_volume.end()));
    return max_comm_volume;
}

EdgeWeight quality_metrics::max_communication_volume(graph_access & G) {
    std::vector<EdgeWeight> block_volume(G.get_partition_count(),0);
    forall_nodes(G, node) {
        PartitionID block = G.getPartitionIndex(node);
        std::vector<bool> block_incident(G.get_partition_count(), false);
        block_incident[block] = true;
        int num_incident_blocks = 0;

        forall_out_edges(G, e, node) {
            NodeID target = G.getEdgeTarget(e);
            PartitionID target_block = G.getPartitionIndex(target);
            if(!block_incident[target_block]) {
                block_incident[target_block] = true;
                num_incident_blocks++;
            }
        } endfor
        block_volume[block] += num_incident_blocks;
    } endfor

    EdgeWeight max_comm_volume = *(std::max_element(block_volume.begin(), block_volume.end()));
    return max_comm_volume;
}


int quality_metrics::boundary_nodes(graph_access& G) {
        int no_of_boundary_nodes = 0;
        forall_nodes(G, n) { 
                PartitionID partitionIDSource = G.getPartitionIndex(n);

                forall_out_edges(G, e, n) {
                        NodeID targetNode = G.getEdgeTarget(e);
                        PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);

                        if (partitionIDSource != partitionIDTarget) {
                                no_of_boundary_nodes++;
                                break; 
                        }
                } endfor 
        }       endfor
        return no_of_boundary_nodes;
}


double quality_metrics::balance(graph_access& G) {
        std::vector<PartitionID> part_weights(G.get_partition_count(), 0);

        double overallWeight = 0;

        forall_nodes(G, n) {
                PartitionID curPartition = G.getPartitionIndex(n);
                part_weights[curPartition] += G.getNodeWeight(n);
                overallWeight += G.getNodeWeight(n);
        } endfor

        double balance_part_weight = ceil(overallWeight / (double)G.get_partition_count());
        double cur_max             = -1;

        forall_blocks(G, p) {
                double cur = part_weights[p];
                if (cur > cur_max) {
                        cur_max = cur;
                }
        } endfor

        double percentage = cur_max/balance_part_weight;
        return percentage;
}

EdgeWeight quality_metrics::objective(const PartitionConfig & config, graph_access & G, int* partition_map) {
        if(config.mh_optimize_communication_volume) {
                return max_communication_volume(G, partition_map);
        } else if(config.mh_penalty_for_unconnected) {
                return edge_cut_connected(G, partition_map);
        } else {
                return edge_cut(G, partition_map);
        }
}


