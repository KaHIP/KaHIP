/*
 * Author: Wolfgang Ost
 */

#include <argtable3.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <regex.h>
#include <memory>
#include <vector>
#include <string>
#include "metis.h"

#include "data_structure/graph_access.h"
#include "node_ordering/reductions.h"
#include "io/graph_io.h"
#include "parse_parameters.h"
#include "partition/partition_config.h"
#include "tools/timer.h"

int main(int argn, char **argv) {
        PartitionConfig partition_config;
        std::string graph_filename;

        bool is_graph_weighted = false;
        bool suppress_output = false;
        bool recursive = false;

        int ret_code = parse_parameters(argn, argv,
                                        partition_config,
                                        graph_filename,
                                        is_graph_weighted,
                                        suppress_output, recursive);

        if (ret_code) {
                return 0;
        }

        // Backup stdout
        std::streambuf* backup = std::cout.rdbuf();
        if(suppress_output) {
                std::cout.rdbuf(nullptr);
        }

        graph_access input_graph;

        timer t;
        graph_io::readGraphWeighted(input_graph, graph_filename);
        std::cout << "io time: " << t.elapsed() << std::endl;

        // Make input_graph unweighted
        bool has_node_weights = false;
        forall_nodes(input_graph, node) {
                if (input_graph.getNodeWeight(node) != 1) {
                        has_node_weights = true;
                }
                input_graph.setNodeWeight(node, 1);
        } endfor
        if (has_node_weights) {
                std::cout << "There were nodes with weight != 1" << std::endl;
        }
        std::cout << "graph has " << input_graph.number_of_nodes() << " nodes and "
                                  << input_graph.number_of_edges() << " edges" << std::endl;

        timer full_nd_timer;
        full_nd_timer.restart();
        t.restart();
        // 'active_graph' is the graph to use after reductions have been applied.
        // If no reductions have been applied, 'active_graph' points to 'input_graph'.
        // Otherwise, it points to 'reduction_stack.back()->get_reduced_graph()'.
        graph_access *active_graph;
        std::vector<std::unique_ptr<Reduction>> reduction_stack;
        bool used_reductions = apply_reductions(partition_config, input_graph, reduction_stack);
        if (used_reductions) {
                active_graph = &reduction_stack.back()->get_reduced_graph();
        } else {
                active_graph = &input_graph;
        }
        //std::cout << "Time for reductions: " << t.elapsed() << std::endl;

        t.restart();
        idx_t num_nodes = active_graph->number_of_nodes();
        // convert the graph into metis style
        idx_t* xadj = new idx_t[num_nodes + 1];
        forall_nodes((*active_graph), node) {
                xadj[node] = (idx_t)active_graph->get_first_edge(node);
        } endfor
        xadj[num_nodes] = (idx_t)active_graph->number_of_edges();
        idx_t* adjncy = new idx_t[active_graph->number_of_edges()];
        forall_edges((*active_graph), edge) {
                adjncy[edge] = (idx_t)active_graph->getEdgeTarget(edge);
        } endfor

        idx_t* perm = new idx_t[active_graph->number_of_nodes()];
        idx_t* iperm = new idx_t[active_graph->number_of_nodes()];      // inverse ordering. This is the one we are interested in.
        idx_t* metis_options = new idx_t[METIS_NOPTIONS];

        // Perform nested dissection with Metis
        if (num_nodes > 0) {
                METIS_SetDefaultOptions(metis_options);
                METIS_NodeND(&num_nodes, xadj, adjncy, nullptr, metis_options, perm, iperm);
        }
        //std::cout << "Time spent in metis: " << t.elapsed() << std::endl;
        
        t.restart();
        // Place labels of reduced graph in a vector for uncontracting
        std::vector<NodeID> reduced_labels(active_graph->number_of_nodes(), 0);
        for (size_t i = 0; i < active_graph->number_of_nodes(); ++i) {
                reduced_labels[i] = iperm[i];
        }

        // Map ordering of reduced graph to input graph
        std::vector<NodeID> final_labels;
        if (used_reductions) {
                map_ordering(reduction_stack, reduced_labels, final_labels);
        } else {
                final_labels = reduced_labels;
        }
        //std::cout << "Time for mapping: " << t.elapsed() << std::endl;
        auto nd_time = full_nd_timer.elapsed();

        // Restore cout output stream
        std::cout.rdbuf(backup);

        // Write resut to file
        std::string filename;
        if (!partition_config.filename_output.compare("")) {
                filename = "tmpnodeordering";
        } else {
                filename = partition_config.filename_output;
        }
        std::ofstream result_fs;
        result_fs.open(filename);
        if (result_fs.is_open()) {
                print_ordering(result_fs, final_labels);
        } else {
                std::cout << "Failed to open file " << filename << std::endl;
        }

        std::cout << "time spent to compute node ordering " << nd_time << std::endl;

        //std::cout << "Reduction statistics:" << std::endl;
        //reduction_stat_counter::get_instance().print_summary(std::cout);

        delete[] xadj;
        delete[] adjncy;
        delete[] perm;
        delete[] iperm;
        delete[] metis_options;

        return 0;
}
