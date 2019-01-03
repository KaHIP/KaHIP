/******************************************************************************
 * graph_extractor.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPH_EXTRACTOR_PDUTVIEF
#define GRAPH_EXTRACTOR_PDUTVIEF

#include "data_structure/graph_access.h"
#include "definitions.h"

class graph_extractor {
        public:
                graph_extractor();
                virtual ~graph_extractor();

                void extract_block(graph_access & G, 
                                   graph_access & extracted_block, 
                                   PartitionID block, 
                                   std::vector<NodeID> & mapping);

                void extract_two_blocks(graph_access & G, 
                                        graph_access & extracted_block_lhs, 
                                        graph_access & extracted_block_rhs, 
                                        std::vector<NodeID> & mapping_lhs, 
                                        std::vector<NodeID> & mapping_rhs,
                                        NodeWeight & partition_weight_lhs,
                                        NodeWeight & partition_weight_rhs); 

               void extract_two_blocks_connected(graph_access & G, 
                                                 std::vector<NodeID> lhs_nodes,
                                                 std::vector<NodeID> rhs_nodes,
                                                 PartitionID lhs, 
                                                 PartitionID rhs,
                                                 graph_access & pair,
                                                 std::vector<NodeID> & mapping) ;


};


#endif /* end of include guard: GRAPH_EXTRACTOR_PDUTVIEF */
