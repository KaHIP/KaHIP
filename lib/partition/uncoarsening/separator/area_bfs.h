/*
* Source of KaHIP -- Karlsruhe High Quality Partitioning.
* Christian Schulz <christian.schulz.phone@gmail.com>
*/ 

#ifndef AREA_BFS_OD13WIZM
#define AREA_BFS_OD13WIZM

#include <algorithm>
#include "partition_config.h"
#include "data_structure/graph_access.h"
#include "tools/random_functions.h"

class area_bfs {
	public:
		area_bfs();
		virtual ~area_bfs();

		void perform_bfs(const PartitionConfig & config, 
				graph_access & G, 
				std::vector< NodeID > & input_separator, 
				PartitionID block, 
				std::vector< NodeWeight > & block_weights,
				std::vector< NodeID > & reached_nodes) {


			// for correctness, in practice will almost never be called
			if( round == std::numeric_limits<int>::max()) {
				round = 0;
				for( int i = 0; m_deepth.size(); i++) {
					m_deepth[i] = 0;
				}
			}

			round++; std::queue<NodeID> node_queue;

			random_functions::permutate_vector_good(input_separator, false);

			/***************************
			 * Initialize the Queue
			 * *************************/
			for(unsigned int i = 0; i < input_separator.size(); i++) {
				node_queue.push(input_separator[i]);
				m_deepth[input_separator[i]] = round;
			}

			NodeWeight size_lhs = block_weights[0];
			NodeWeight size_rhs = block_weights[1];
			NodeWeight size_sep = block_weights[2];

			NodeWeight accumulated_weight = 0;
			NodeID upper_bound_no_nodes;

			if( block == 0 ) {
				upper_bound_no_nodes = std::max((int)(config.region_factor_node_separators*config.upper_bound_partition - size_rhs - size_sep), 0);
			} else {
				upper_bound_no_nodes = std::max((int)(config.region_factor_node_separators*config.upper_bound_partition - size_lhs - size_sep), 0);
			}
			upper_bound_no_nodes = std::min(upper_bound_no_nodes, block_weights[block]-1);

			/***************************
			 * Do the BFS
			 ***************************/
			while (!node_queue.empty()) {
				if(accumulated_weight >= upper_bound_no_nodes) break;
				NodeID n = node_queue.front();
				node_queue.pop();

				forall_out_edges(G,e,n) {
					NodeID target = G.getEdgeTarget(e);
					if(m_deepth[target] != round && G.getPartitionIndex(target) == block 
					   && accumulated_weight + G.getNodeWeight(target) <= upper_bound_no_nodes) {
						m_deepth[target] = round;
						node_queue.push(target);
						reached_nodes.push_back(target);
						accumulated_weight += G.getNodeWeight(target);
					}
				} endfor
			}
		}

		static std::vector<int> m_deepth;
		static int round;

};


#endif /* end of include guard: AREA_BFS_OD13WIZM */
