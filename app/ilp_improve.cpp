/******************************************************************************
 * ilp_improve.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/
#include <argtable3.h>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "ilp_improve/ilp_helpers.h"
#include "ilp_improve/ilp_improve.h"
#include "macros_assertions.h"
#include "mapping/mapping_algorithms.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"

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

        ilp_improve ilp;
        graph_access G;

        timer t;
        graph_io::readGraphWeighted(G, graph_filename);
        std::cout << "io time: " << t.elapsed() << std::endl;
        G.set_partition_count(partition_config.k);
        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        // partition graph
        quality_metrics qm;
        graph_io::readPartition(G, partition_config.input_partition);
        ilp_helpers ilp_h;
        auto pid1 = ilp_h.extractPartition(G);
        ilp_h.setFirstPartition(G, pid1);


        // output some information about the initial partition
        std::cout << "Info on start graph: " << std::endl;
        std::cout << "Total num nodes in uncoarsened graph: "
                << G.number_of_nodes() << std::endl;
        ilp.graphInfo(qm, G);

        EdgeWeight before_ec = qm.edge_cut(G);

        auto limit_nonzeroes = (size_t) partition_config.ilp_limit_nonzeroes;
        // compute BFS
        std::unordered_set<NodeID> reachable;
        ilp.computeBFS(G, reachable, partition_config, limit_nonzeroes);
        PartitionID num_blocks = ilp.computeBlocks(G, reachable, partition_config.k);
        PartitionID no_of_coarse_vertices = reachable.size() + num_blocks;

        EdgeWeight base_ec = before_ec;

        graph_access coarser;
        std::vector<bool> coarse_presets(coarser.number_of_nodes(), false);
        if (partition_config.ilp_mode == OptimizationMode::OVERLAP) {
                auto partitions = ilp.createPartitions(G, partition_config);
                size_t best_index;
                std::tie(no_of_coarse_vertices, base_ec, best_index) =
                        ilp.buildOverlapGraph(G, partition_config, pid1, partitions);

                G.set_partition_count(no_of_coarse_vertices);
                complete_boundary bnd(&G);
                bnd.build();
                bnd.getUnderlyingQuotientGraph(coarser);
                coarse_presets = ilp.findPresets(G, partition_config, partitions, best_index, coarser);

                for (uint32_t n = 0; n < G.number_of_nodes(); ++n) {
                        coarser.setPartitionIndex(G.getPartitionIndex(n), pid1[n]);
                }

                num_blocks = 0;
        } else {
                // condense graph
                G.set_partition_count(no_of_coarse_vertices);
                ilp.setPartitionIDs(G, reachable, num_blocks);
                complete_boundary bnd(&G);
                bnd.build();
                bnd.getUnderlyingQuotientGraph(coarser);
        }

        forall_nodes(G, node) {
                NodeID c_node = G.getPartitionIndex(node);
                coarser.setPartitionIndex(c_node, pid1[node]);
        } endfor


        //run 2 hours then use best solution
        double timelimit = partition_config.ilp_timeout;

        timer ilptime;
        ilptime.restart();

        // compute ILP on coarser graph
        ilp.computeIlp(coarser, partition_config, num_blocks,
                       timelimit, pid1, coarse_presets);

        std::cout <<  "ILP took " << ilptime.elapsed() << std::endl;

        // transfer ILP to graph
        ilp.transferIlp(G, coarser);
        G.set_partition_count(partition_config.k);

        std::stringstream filename;
        if(!partition_config.filename_output.compare("")) {
                filename << "tmppartition" << partition_config.k << ".impr";
        } else {
                filename << partition_config.filename_output;
        }

        graph_io::writePartition(G, filename.str());

        return 0;
}

