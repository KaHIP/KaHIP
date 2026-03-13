/******************************************************************************
 * connectivity_check.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef CONNECTIVITY_CHECK_7J3X9K2M
#define CONNECTIVITY_CHECK_7J3X9K2M

#include <queue>
#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"

// Check whether moving `node` out of `block` would disconnect the block.
// Returns true if the move would disconnect, false if it is safe.
//
// Correctness argument: if removing `node` splits the block into components
// C1, C2, ..., Ck, each component must contain at least one neighbor of `node`
// (because the block was connected before removal). Therefore, checking
// reachability among same-block neighbors is sufficient.
inline bool would_disconnect_block(graph_access & G, NodeID node, PartitionID block) {
        // Collect neighbors of `node` that are in `block`
        std::vector<NodeID> same_block_neighbors;
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                if(target != node && G.getPartitionIndex(target) == block) {
                        same_block_neighbors.push_back(target);
                }
        } endfor

        // 0 neighbors: isolated in block, safe to remove
        // 1 neighbor: leaf in block-induced subgraph, removing a leaf cannot disconnect
        if(same_block_neighbors.size() <= 1) {
                return false;
        }

        // BFS from first same-block neighbor, restricted to block \ {node}
        // Check if all other same-block neighbors are reachable
        std::vector<bool> visited(G.number_of_nodes(), false);
        visited[node] = true; // exclude the node being moved
        visited[same_block_neighbors[0]] = true;

        unsigned int found = 1; // count of same_block_neighbors reached (starting node counts)
        unsigned int need = same_block_neighbors.size();

        std::queue<NodeID> bfsqueue;
        bfsqueue.push(same_block_neighbors[0]);

        while(!bfsqueue.empty() && found < need) {
                NodeID source = bfsqueue.front();
                bfsqueue.pop();

                forall_out_edges(G, e, source) {
                        NodeID target = G.getEdgeTarget(e);
                        if(!visited[target] && G.getPartitionIndex(target) == block) {
                                visited[target] = true;
                                bfsqueue.push(target);

                                // Check if this is one of the same-block neighbors
                                for(unsigned int i = 1; i < same_block_neighbors.size(); i++) {
                                        if(target == same_block_neighbors[i]) {
                                                found++;
                                                break;
                                        }
                                }
                        }
                } endfor
        }

        return found < need;
}

#endif /* end of include guard: CONNECTIVITY_CHECK_7J3X9K2M */
