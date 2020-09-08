/*
 * Author: Wolfgang Ost
 */

#ifndef H_ISOLATED_CLIQUE
#define H_ISOLATED_CLIQUE

#include <array>
#include <iostream>
#include <iomanip>
#include <memory>
#include <unordered_map>
#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "node_ordering/ordering_tools.h"
#include "partition/partition_config.h"

class reduction_stat_counter final {
private:
        typedef std::array<double, nested_dissection_reduction_type::num_types> percent_array;
        typedef std::array<int, nested_dissection_reduction_type::num_types> count_array;

        // Sum up to what percentage a graph was reduced by a reduction, stored per recursion level
        std::vector<percent_array> percent_sums;
        // Count how often a reduction has been applied, stored per recursion level
        std::vector<count_array> application_counts;

        int deg2_separator_count;

        inline reduction_stat_counter() {
                percent_sums.push_back(percent_array{});
                application_counts.push_back(count_array{});
        }

public:
        reduction_stat_counter(reduction_stat_counter const&) = delete;
        void operator=(reduction_stat_counter const&) = delete;

        inline static reduction_stat_counter &get_instance() {
                static reduction_stat_counter instance;
                return instance;
        }

        inline void count_reduction(nested_dissection_reduction_type type,
                                    int num_original_nodes,
                                    int num_reduced_nodes,
                                    int recursion_level = 0) {
                while (recursion_level >= (int)percent_sums.size()) {
                        percent_sums.push_back(percent_array{});
                        application_counts.push_back(count_array{});
                }
                // Only count, if there were nodes that could be removed
                if (num_original_nodes > 0) {
                        percent_sums[recursion_level][type] += (double)num_reduced_nodes/(double)num_original_nodes;
                        application_counts[recursion_level][type]++;
                }
        }

        inline void count_deg2_separators(int count) {
                deg2_separator_count += count;
        }

        inline void print_summary(std::ostream &stream) const {
                stream << "Eliminated " << deg2_separator_count << " degree-2 separators." << std::endl;
                stream << "type : average reduction, number of applications" << std::endl;
                for (int i = 0; i < nested_dissection_reduction_type::num_types; ++i) {
                        double percent = 0;
                        int count = 0;
                        for (size_t level = 0; level < percent_sums.size(); ++level) {
                                percent += percent_sums[level][i];
                                count += application_counts[level][i];
                        }
                        stream << i << ": " << std::setw(8) << percent/count
                                            << ", " << std::setw(3) << std::right << count << std::endl;
                }
        }

        inline void print_level(std::ostream &stream, int recursion_level) const {
                stream << std::setw(5) << recursion_level;
                for (int i = 0; i < nested_dissection_reduction_type::num_types; ++i) {
                        stream << " " << std::setw(8) << percent_sums[recursion_level][i]/application_counts[recursion_level][i];
                }
        }

        inline void print_histogram(std::ostream &stream) const {
                stream << "level";
                for (unsigned int i = 0; i < nested_dissection_reduction_type::num_types; ++i) {
                        stream << " " << std::setw(8) << i;
                }
                stream << std::endl;
                for (unsigned int i = 0; i < percent_sums.size(); ++i) {
                        print_level(stream, i);
                        stream << std::endl;
                }
        }

};

/**************/
/* REDUCTIONS */
/**************/

class Reduction {

public:
        inline Reduction(graph_access &graph_before) : graph_before(graph_before) {};

        virtual ~Reduction() {};

        // Apply the reduction to the graph passed to the constructor
        virtual void apply() = 0;

        // Map the ordering on the reduced graph ('reduced_label') to the original graph.
        // The mapped ordering is stored in 'new_label'
        virtual void map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const = 0;

        inline graph_access &get_reduced_graph() {
                return graph_after;
        }

        // Find the number of nodes for which the reachable set size in the reduced graph
        // is not the same as the reachable set size in the original graph
        virtual int test_node_degrees() = 0;

        // Print the mapping for this reduction to a stream
        virtual void print_mapping(std::ostream &stream) = 0;

protected:
        // Reference to unreduced graph
        graph_access &graph_before;
        // Reduced graph
        graph_access graph_after;

};

class EliminationReduction : public Reduction {

public:
        inline EliminationReduction(graph_access &graph_before) : Reduction(graph_before) {};

        // This test doesn't really make sense for eliminations, since we expect degrees to change
        // We still implement it, just for completeness.
        virtual inline int test_node_degrees() override {
                int count = 0;
                forall_nodes(graph_after, node) {
                        auto reach = compute_reachable_set_size(graph_after, node);
                        auto degree = compute_reachable_set_size(graph_before, mapping[node]);
                        if (reach != degree) {
                                count++;
                        }
                } endfor
                return count;
        }

        virtual inline void print_mapping(std::ostream &stream) override {
                for (const auto node: mapping) {
                        stream << node << std::endl;
                }
        }

protected:
        // Mapping from the reduced graph to the original graph
        std::vector<NodeID> mapping;
        // Eliminated nodes in elimination order
        std::vector<NodeID> label_first;

};

class ContractionReduction : public Reduction {

public:
        inline ContractionReduction(graph_access &graph_before) : Reduction(graph_before) {};

        virtual inline int test_node_degrees() override {
                int count = 0;
                forall_nodes(graph_after, node) {
                        auto reach = compute_reachable_set_size(graph_after, node);
                        for (const auto old_node: mapping[node]) {
                                auto degree = compute_reachable_set_size(graph_before, old_node);
                                if (reach != degree) {
                                        count++;
                                }
                        }
                } endfor
                return count;
        }

        virtual inline void print_mapping(std::ostream &stream) override {
                for (auto it = mapping.begin(); it != mapping.end(); ++it) {
                        for (const auto node: it->second) {
                                stream << node << " ";
                        }
                        stream << std::endl;
                }
        }

protected:
        // Mapping from the reduced graph to groups of nodes in the original graph
        std::unordered_map<NodeID, std::vector<NodeID>> mapping;

};

// Iteratively eliminates nodes whose neighborhood is a clique
class SimplicialNodeReduction : public EliminationReduction {

public:
        inline SimplicialNodeReduction(graph_access &graph_before, Count degree_limit) : EliminationReduction(graph_before),
                                                                                         degree_limit(degree_limit) {};

        void apply() override;

        void map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const override;

private:
        Count degree_limit;

};

// Contract indistinguishable nodes and twins
class IndistinguishableNodeReduction : public ContractionReduction {

public:
        IndistinguishableNodeReduction(graph_access &graph_before) : ContractionReduction(graph_before) {}
        
        void apply() override;

        void map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const override;

};

// Contract indistinguishable nodes and twins
class TwinReduction : public ContractionReduction {

public:
        TwinReduction(graph_access &graph_before) : ContractionReduction(graph_before) {}
        
        // Detect nodes that have a common neighborhood or neighborhood set.
        // Low-degree nodes are tested first. Stops, once more than stop_factor * graph_before.number_of_nodes() nodes have been tested. 
        // stop_factor == 1 -> find all indistinguishable nodes
        // stop_factor == 0.5 -> only test half of the nodes
        // stop_factor == 0 -> don't find any indistinguishable nodes
        void apply() override;

        void map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const override;

};

// Contract paths of nodes with degree 2
class PathCompression : public ContractionReduction {

public:
        inline PathCompression(graph_access &graph_before) : ContractionReduction(graph_before) {};

        void apply() override;

        void map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const override;

protected:
        // We reverse paths and parts of paths in the map function, so we make the mapping mutable for this reduction
        mutable std::unordered_map<NodeID, std::vector<NodeID>> mapping;

};

// Eliminate nodes of degree 2. This takes away a step of the min-degree algorithm.
// It is only useful if there are no simplicial nodes in the graph.
class Degree2Elimination : public EliminationReduction {

public:
        inline Degree2Elimination(graph_access &graph_before) : EliminationReduction(graph_before) {};

        void apply() override;

        void map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const override;

};

// Contract nodes of degree 3 that share at least one neighbor.
class TriangleContraction : public ContractionReduction {

public:
        inline TriangleContraction(graph_access &graph_before) : ContractionReduction(graph_before) {};

        // Contract nodes of degree 3 that share exactly one neighbor.
        void apply() override;

        void map(std::vector<NodeID> &reduced_label, std::vector<NodeID> &new_label) const override;

};

// Apply reductions to in_graph, based on the given configuration.
// The reductions applied to the graph are stored in 'reduction_stack'.
// Returns false if no reductions were applied, true otherwise.
// 'recursion_level' is an optional parameter for recording efficiency of reductions per level of recursion.
bool apply_reductions(const PartitionConfig &config,
                      graph_access &in_graph,
                      std::vector<std::unique_ptr<Reduction>> &reduction_stack,
                      int recursion_level = 0);

// Map an ordering on the reduced graph to the original graph.
// The reductions to "undo" are given in 'reduction_stack',
// the ordering on the reduced graph is given in 'reduced_label'.
// WARNING: The contents of 'reduced_label' are destroyed in the process.
// The resulting ordering is stored in 'final_label'.
void map_ordering(const std::vector<std::unique_ptr<Reduction>> &reduction_stack,
                  std::vector<NodeID> &reduced_label,
                  std::vector<NodeID> &final_label);

#endif /* end of include guard */
