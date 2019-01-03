/******************************************************************************
 * boundary_bfs.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <queue>
#include <vector>

#include "boundary_bfs.h"
#include "random_functions.h"

boundary_bfs::boundary_bfs() {
                
}

boundary_bfs::~boundary_bfs() {
                
}

bool boundary_bfs::boundary_bfs_search(graph_access & G, 
                                       std::vector<NodeID> & start_nodes, 
                                       PartitionID partition, 
                                       NodeWeight upper_bound_no_nodes, 
                                       std::vector<NodeID> & reached_nodes,
                                       NodeWeight & stripe_weight,
                                       bool flow_tiebreaking) {

        std::queue<NodeID> node_queue;
	std::vector<int> deepth(G.number_of_nodes(), -1);
	int cur_deepth = 0;
       
        if(flow_tiebreaking) {
                random_functions::permutate_vector_good(start_nodes, false);
        } 
	/***************************
	 * Initialize the Queue
	 * *************************/
        NodeWeight accumulated_weight = 0;
	for(unsigned int i = 0; i < start_nodes.size(); i++) {
		node_queue.push(start_nodes[i]);
		ASSERT_TRUE(G.getPartitionIndex(start_nodes[i]) == partition);
		deepth[start_nodes[i]] = cur_deepth;
		reached_nodes.push_back(start_nodes[i]);
                accumulated_weight += G.getNodeWeight(start_nodes[i]); 
	}
	++cur_deepth;

	if(accumulated_weight >= upper_bound_no_nodes) {
                stripe_weight = accumulated_weight;
                return false;
        }
	/***************************
	 * Do the BFS
	 ***************************/
	while (!node_queue.empty()) {
		if(accumulated_weight >= upper_bound_no_nodes) break;
		NodeID n = node_queue.front();
		node_queue.pop();

		if (deepth[n] == cur_deepth) {
			cur_deepth++;
		}
		forall_out_edges(G,e,n) {
			NodeID t = G.getEdgeTarget(e);
			if(deepth[t] == -1 && G.getPartitionIndex(t) == partition 
                        && accumulated_weight + G.getNodeWeight(t) <= upper_bound_no_nodes) {
				deepth[t] = cur_deepth;
				node_queue.push(t);
				reached_nodes.push_back(t);
                                accumulated_weight += G.getNodeWeight(t);
			}
		} endfor
	}
        bool some_to_do = stripe_weight != accumulated_weight;
        stripe_weight   = accumulated_weight;
        return some_to_do;

}

