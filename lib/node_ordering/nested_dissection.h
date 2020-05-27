/*
 * Author: Wolfgang Ost
 */

#ifndef NESTED_DISSECTION
#define NESTED_DISSECTION

#include <iosfwd>
#include <memory>
#include <unordered_map>
#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "node_ordering/reductions.h"
#include "partition/partition_config.h"

class nested_dissection {
public:
        
        nested_dissection(graph_access * const G);
        nested_dissection(graph_access * const G, int recursion_level);

        void perform_nested_dissection(PartitionConfig &config);

        // Obtain a reference to the ordering. If 'perform_nested_dissection' was not called,
        // the ordering is uninitialized and thus invalid
        const std::vector<NodeID>& ordering() const;

private:
        // pointer to the graph passed to perform_nested_dissection
        graph_access * const original_graph;

        // How often 'recurse_dissection' was called to get to this level
        int m_recursion_level;

        // computed elimination order
        // node x is eliminated in step m_label[x]
        std::vector<NodeID> m_label;

        // elimination order in the reduced graph
        std::vector<NodeID> m_reduced_label;

        std::vector<std::unique_ptr<Reduction>> m_reduction_stack;

        // Compute a separator of the graph G
        void compute_separator(PartitionConfig &config, graph_access &G);

        // Apply nested dissection to the subgraph of G induced by the partition with ID block
        // new labels start at order_begin, which is updated to the value past the new largest label
        void recurse_dissection(PartitionConfig &config, graph_access &G, PartitionID block, NodeID &order_begin);

};

#endif /* end of include guard: NESTED_DISSECTION */ 
