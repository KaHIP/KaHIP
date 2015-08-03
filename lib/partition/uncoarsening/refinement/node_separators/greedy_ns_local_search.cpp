//
// Author: Christian Schulz <christian.schulz@kit.edu>
// 

#include "greedy_ns_local_search.h"
#include "tools/quality_metrics.h"

greedy_ns_local_search::greedy_ns_local_search() {
                
}

greedy_ns_local_search::~greedy_ns_local_search() {
                
}

EdgeWeight greedy_ns_local_search::perform_refinement(const PartitionConfig & config, graph_access & G) {
        NodeWeight weight_blockA = 0;
        NodeWeight weight_blockB = 0;

        forall_nodes(G, node) {
                if( G.getPartitionIndex(node) == 0) {
                        weight_blockA += G.getNodeWeight(node);
                } else if( G.getPartitionIndex(node) == 1 ) {
                        weight_blockB += G.getNodeWeight(node);
                }
        } endfor
        
        forall_nodes(G, node) {
               if( G.getPartitionIndex(node) == 2) { 
                       // now this is a separator node
                       // check if we can move it to one of the blocks
                       Gain gainArhs = 0;
                       Gain gainBrhs = 0;

                       forall_out_edges(G, e, node) {
                               NodeID target = G.getEdgeTarget(e);
                               if( G.getPartitionIndex(target) == 0) {
                                        gainArhs += G.getNodeWeight(target);
                               } else if( G.getPartitionIndex(target) == 1 ) {
                                        gainBrhs += G.getNodeWeight(target);
                               }
                       } endfor

                       Gain gainToA = G.getNodeWeight(node) - gainBrhs;
                       Gain gainToB = G.getNodeWeight(node) - gainArhs;

                               //std::cout <<  "gainA " <<  gainToA <<  " gainB " <<  gainToB  << std::endl;
                       if( gainToA > 0 || gainToB > 0 ) {
                                // try to move it
                                quality_metrics qm;
                                std::cout <<  "size of separator before movement " <<  qm.separator_weight(G)  << std::endl;
                                Gain top_gain  = gainToA > gainToB ? gainToA : gainToB;
                                PartitionID top_block = top_gain == gainToA ? 0 : 1;
                                std::cout <<  "topblock " <<  top_block  << std::endl;
                                std::cout <<  "topgain " <<   top_gain  << std::endl;

                                if( top_block == 0 ) {
                                        if( weight_blockA + G.getNodeWeight(node) < config.upper_bound_partition) {
                                                G.setPartitionIndex(node, 0);
                                                forall_out_edges(G, e, node) {
                                                        NodeID target = G.getEdgeTarget(e);

                                                        if( G.getPartitionIndex( target ) == 1 ) {
                                                                G.setPartitionIndex(target,2);
                                                        }
                                                } endfor
                                                weight_blockA += G.getNodeWeight(node);
                                                weight_blockB -= gainBrhs;
                                        } else if( gainToB > 0 && weight_blockB + G.getNodeWeight(node) < config.upper_bound_partition ) {
                                                G.setPartitionIndex(node, 1);
                                                forall_out_edges(G, e, node) {
                                                        NodeID target = G.getEdgeTarget(e);

                                                        if( G.getPartitionIndex( target ) == 0 ) {
                                                                G.setPartitionIndex(target,2);
                                                        }
                                                } endfor

                                                weight_blockB += G.getNodeWeight(node);
                                                weight_blockA -= gainArhs;

                                        }
                                        std::cout <<  "moved a node !"   << std::endl;
                                }

                                if( top_block == 1 ) {
                                        if( weight_blockB + G.getNodeWeight(node) < config.upper_bound_partition) {
                                                G.setPartitionIndex(node, 1);
                                                forall_out_edges(G, e, node) {
                                                        NodeID target = G.getEdgeTarget(e);

                                                        if( G.getPartitionIndex( target ) == 0 ) {
                                                                G.setPartitionIndex(target,2);
                                                        }
                                                } endfor

                                                weight_blockB += G.getNodeWeight(node);
                                                weight_blockA -= gainArhs;


                                        } else if( gainToA > 0 && weight_blockA + G.getNodeWeight(node) < config.upper_bound_partition ) {
                                                G.setPartitionIndex(node, 0);
                                                forall_out_edges(G, e, node) {
                                                        NodeID target = G.getEdgeTarget(e);

                                                        if( G.getPartitionIndex( target ) == 1 ) {
                                                                G.setPartitionIndex(target,2);
                                                        }
                                                } endfor
                                                weight_blockA += G.getNodeWeight(node);
                                                weight_blockB -= gainBrhs;



                                        }

                                        std::cout <<  "moved a node !"   << std::endl;
                                }
                                std::cout <<  "size of separator after movement " <<  qm.separator_weight(G)  << std::endl;
                                //exit(0);

                       }
               
               }
        } endfor

        return 0;
}

