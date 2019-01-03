/******************************************************************************
 * balance_configuration.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef BALANCE_CONFIGURATION_JCQB9ZGV
#define BALANCE_CONFIGURATION_JCQB9ZGV

#include "data_structure/graph_access.h"
#include "partition/partition_config.h"

class balance_configuration {
public:
        balance_configuration() {};
        virtual ~balance_configuration() {};

        void configurate_balance( PartitionConfig & partition_config, graph_access & G ) {
                NodeWeight largest_graph_weight = 0;
                forall_nodes(G, node) {
                        largest_graph_weight += G.getNodeWeight(node);
                } endfor

                NodeWeight edge_weights = 0;
                if(partition_config.balance_edges && partition_config.imbalance != 0) {
                        // balancing edges is disabled for the perfectly balanced case since this case requires uniform node weights
                        forall_nodes(G, node) {
                                NodeWeight weighted_degree = 0;
                                forall_out_edges(G, e, node) {
                                        weighted_degree += G.getEdgeWeight(e);
                                } endfor          

                                edge_weights += weighted_degree;
                                G.setNodeWeight(node, G.getNodeWeight(node) + weighted_degree);
                        } endfor
                        
                }

                double epsilon  = (partition_config.imbalance)/100.0;
                if(  partition_config.imbalance == 0 && !partition_config.kaffpaE) {
                        partition_config.upper_bound_partition    = (1+epsilon+0.01)*ceil(largest_graph_weight/(double)partition_config.k);
                        partition_config.kaffpa_perfectly_balance = true;
                } else {
                        NodeWeight load                           = largest_graph_weight + edge_weights;
                        partition_config.upper_bound_partition    = (1+epsilon)*ceil(load/(double)partition_config.k);
                }

                partition_config.largest_graph_weight       = largest_graph_weight;
                partition_config.graph_allready_partitioned = false;
                partition_config.kway_adaptive_limits_beta  = log(G.number_of_nodes());
                partition_config.work_load                  = largest_graph_weight + edge_weights;

                //std::cout <<  "block weight upper bound " <<  partition_config.upper_bound_partition  << std::endl;

        }

};


#endif /* end of include guard: BALANCE_CONFIGURATION_JCQB9ZGV */
