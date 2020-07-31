/*
 * Author: Wolfgang Ost
 */

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

#include "balance_configuration.h"
#include "node_ordering/min_degree_ordering.h"
#include "node_ordering/nested_dissection.h"
#include "node_ordering/reductions.h"
#include "partition/graph_partitioner.h"
#include "partition/uncoarsening/separator/area_bfs.h"
#include "partition/uncoarsening/separator/vertex_separator_algorithm.h"
#include "tools/graph_extractor.h"
#include "tools/macros_assertions.h"
#include "tools/quality_metrics.h"

nested_dissection::nested_dissection(graph_access * const G) :
        original_graph(G), m_recursion_level(0) {}

nested_dissection::nested_dissection(graph_access * const G, int recursion_level) :
        original_graph(G), m_recursion_level(recursion_level) {}


void nested_dissection::perform_nested_dissection(PartitionConfig &config) {
        if (original_graph->number_of_nodes() == 0) {
                return;
        }

        // 'reduced_graph' is a copy of 'original_graph', with reductions applied.
        // If no reductions were applied, 'reduced_graph' is empty, and we use 'original_graph' instead.
        graph_access reduced_graph;
        graph_access *active_graph;
        bool used_reductions = apply_reductions(config, *original_graph, m_reduction_stack, m_recursion_level);
        if (used_reductions) {
                active_graph = &m_reduction_stack.back()->get_reduced_graph();
        } else {
                active_graph = original_graph;
        }

        m_reduced_label.resize(active_graph->number_of_nodes());

        if (active_graph->number_of_nodes() > 0) {
                if (active_graph->number_of_nodes() < config.dissection_rec_limit) {
                        // Stop nested dissection and use the min degree algorithm instead
                        MinDegree(active_graph).perform_ordering(m_reduced_label);
                } else {
                        NodeID order_begin = 0;
                        // continue nested dissection
                        compute_separator(config, *active_graph);

                        // perform nested dissection on subgraphs
                        forall_blocks((*active_graph), p) {
                                if (p != active_graph->getSeparatorBlock()) {
                                        recurse_dissection(config, (*active_graph), p, order_begin);
                                }
                        } endfor
                        // Perform nested dissection on separator block
                        recurse_dissection(config, (*active_graph), active_graph->getSeparatorBlock(), order_begin);
                }
        }

        if (used_reductions) {
                // Map the ordering from the reduced graph to the original graph
                map_ordering(m_reduction_stack, m_reduced_label, m_label);
        } else {
                m_label = m_reduced_label;
        }
}

void nested_dissection::compute_separator(PartitionConfig &config, graph_access &G) {
        // set up the graph and config for computing a node separator
        config.k = 2;
        G.set_partition_count(config.k + 1);
        config.mode_node_separators = true;
        config.graph_allready_partitioned = false;
        balance_configuration bc;
        bc.configurate_balance(config, G);

        // compute the separator
        area_bfs::m_deepth.resize(G.number_of_nodes(), 0);
        forall_nodes(G, node) {
                area_bfs::m_deepth[node] = 0;
        } endfor

        graph_partitioner partitioner;
        partitioner.perform_partitioning(config, G);
}

void nested_dissection::recurse_dissection(PartitionConfig &config, graph_access &G, PartitionID block, NodeID &order_begin) {
        std::vector<NodeID> mapping;
        graph_extractor extractor;
        graph_access subgraph;
        extractor.extract_block(G, subgraph, block, mapping);
        nested_dissection dissection(&subgraph, m_recursion_level + 1);
        dissection.perform_nested_dissection(config);

        // Transfer labels from the subgraph to the reduced graph
        for (size_t i = 0; i < mapping.size(); ++i) {
                m_reduced_label[mapping[i]] = dissection.m_label[i] + order_begin;
        }
        order_begin += mapping.size();
}

const std::vector<NodeID>& nested_dissection::ordering() const {
        return m_label;
}
