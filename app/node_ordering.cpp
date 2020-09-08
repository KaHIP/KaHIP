/*
 * Author: Wolfgang Ost
 */

#include <argtable3.h>
#include <iostream>
#include <regex.h>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "node_ordering/nested_dissection.h"
#include "node_ordering/ordering_tools.h"
#include "node_ordering/reductions.h"
#include "io/graph_io.h"
#include "parse_parameters.h"
#include "partition/partition_config.h"
#include "timer.h"
#include "tools/random_functions.h"

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

        graph_access G;

        timer t;
        graph_io::readGraphWeighted(G, graph_filename);
        std::cout << "io time: " << t.elapsed() << std::endl;

        // Make G unweighted
        bool has_node_weights = false;
        forall_nodes(G, node) {
                if (G.getNodeWeight(node) != 1) {
                        has_node_weights = true;
                }
                G.setNodeWeight(node, 1);
        } endfor
        if (has_node_weights) {
                std::cout << "There were nodes with weight != 1" << std::endl;
        }

        std::cout << "imbalance is set to " << partition_config.imbalance << "%" << std::endl;
        balance_configuration bc;
        bc.configurate_balance(partition_config, G);

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        std::cout << "graph has " << G.number_of_nodes() << " nodes and " << G.number_of_edges() << " edges" << std::endl;

        t.restart();
        nested_dissection dissection(&G);
        dissection.perform_nested_dissection(partition_config);

        // Restore cout output stream
        std::cout.rdbuf(backup);

        std::cout << "time spent to compute node ordering " << t.elapsed() << std::endl;

        std::cout << "Number of fill-edges: " << compute_fill(G, dissection.ordering()) << std::endl;

        std::string filename;
        if (!partition_config.filename_output.compare("")) {
                filename = "tmpnodeordering";
        } else {
                filename = partition_config.filename_output;
        }
        std::ofstream result_fs;
        result_fs.open(filename);
        if (result_fs.is_open()) {
                print_ordering(result_fs, dissection.ordering());
        } else {
                std::cout << "Failed to open file " << filename << std::endl;
        }

        //std::cout << "Reduction statistics:" << std::endl;
        //const auto &stats = reduction_stat_counter::get_instance();
        //stats.print_histogram(std::cout);
        //stats.print_summary(std::cout);

        return 0; 
}
