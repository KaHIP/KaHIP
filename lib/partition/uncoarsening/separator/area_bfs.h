//
// Author: Christian Schulz <christian.schulz@kit.edu>
// 

#ifndef AREA_BFS_OD13WIZM
#define AREA_BFS_OD13WIZM

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
                                 PartitionID block_weight,
                                 std::vector< NodeID > & reached_nodes) {

                        std::queue<NodeID> node_queue;
                        std::vector<int> deepth(G.number_of_nodes(), -1);
                        int cur_deepth = 0;

                        random_functions::permutate_vector_good(input_separator, false);
                        /***************************
                         * Initialize the Queue
                         * *************************/
                        for(unsigned int i = 0; i < input_separator.size(); i++) {
                                node_queue.push(input_separator[i]);
                                deepth[input_separator[i]] = cur_deepth;
                        }
                        ++cur_deepth;

                        NodeWeight size_lhs = 0;
                        NodeWeight size_rhs = 0;
                        NodeWeight size_sep = 0;

                        NodeWeight count_lhs = 0;
                        NodeWeight count_rhs = 0;
                        forall_nodes(G, node) {
                                if(G.getPartitionIndex(node) == 0) {
                                        size_lhs += G.getNodeWeight(node);
                                        count_lhs++;
                                } else if( G.getPartitionIndex(node) == 1) {
                                        size_rhs += G.getNodeWeight(node);
                                        count_rhs++;
                                } else if(G.getPartitionIndex(node) == 2) {
                                        size_sep += G.getNodeWeight(node);
                                }
                        } endfor
                        //std::cout <<  "size lhs " <<  size_lhs << std::endl;
                        //std::cout <<  "size rhs " <<  size_rhs << std::endl;
                        //std::cout <<  "size sep " <<  size_sep << std::endl;
                        //std::cout <<  "count lhs " <<  count_lhs << std::endl;
                        //std::cout <<  "count rhs " <<  count_rhs << std::endl;
                        NodeWeight accumulated_weight = 0;
                        NodeID upper_bound_no_nodes;

                        if( block == 0 ) {
                                upper_bound_no_nodes = std::max((int)(config.region_factor_node_separators*config.upper_bound_partition - size_rhs - size_sep), 0);
                        } else {
                                upper_bound_no_nodes = std::max((int)(config.region_factor_node_separators*config.upper_bound_partition - size_lhs - size_sep), 0);
                        }
                        upper_bound_no_nodes = std::min(upper_bound_no_nodes, block_weight-1);

                        //upper_bound_no_nodes *= config.region_factor_node_separators ;
                        //std::cout <<  "region facotr " <<  config.region_factor_node_separators  << std::endl;
                        //std::cout <<  "upper bound no nodes " <<  upper_bound_no_nodes  << std::endl;
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
                                        NodeID target = G.getEdgeTarget(e);
                                        if(deepth[target] == -1 && G.getPartitionIndex(target) == block 
                                                        && accumulated_weight + G.getNodeWeight(target) <= upper_bound_no_nodes) {
                                                deepth[target] = cur_deepth;
                                                node_queue.push(target);
                                                reached_nodes.push_back(target);
                                                accumulated_weight += G.getNodeWeight(target);
                                        }
                                } endfor
                        }
                        //std::cout <<  "accumulated weight is " <<  accumulated_weight  << std::endl;
                }

};


#endif /* end of include guard: AREA_BFS_OD13WIZM */
