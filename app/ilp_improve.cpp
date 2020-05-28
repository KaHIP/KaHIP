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

    double imbalance_param;

    bool is_graph_weighted = false;
    bool suppress_output = false;
    bool recursive = false;

    int ret_code = parse_parameters(argn, argv,
                                    partition_config,
                                    graph_filename,
                                    is_graph_weighted,
                                    suppress_output, recursive);

    imbalance_param = partition_config.imbalance / 100 + 1;


    if (ret_code) {
        return 0;
    }

    ilp_improve ilp;


    //filename parsing
    size_t dir_pos = graph_filename.find_last_of('/');
    std::string dir = graph_filename.substr(0, dir_pos);
    std::string graph = graph_filename.substr(
            dir_pos + 1, graph_filename.length() - 7 - dir_pos);

    size_t parent_dir_pos = dir.find_last_of('/');
    std::string parent_dir = graph_filename.substr(0, parent_dir_pos);

    std::string partition_filename = parent_dir + "/partitions/"
                                     + std::to_string((int) partition_config.imbalance) + "/" + graph
                                     + "." + std::to_string(partition_config.k) + ".ptn";

    std::cout << partition_filename << std::endl;

    graph_access G;

    timer t;
    graph_io::readGraphWeighted(G, graph_filename);
    std::cout << "io time: " << t.elapsed() << std::endl;
    G.set_partition_count(partition_config.k);
    balance_configuration bc;
    bc.configurate_balance( partition_config, G);

    // partition graph
    quality_metrics qm;
    graph_io::readPartition(G, partition_filename);
    ilp_helpers ilp_h;
    auto pid1 = ilp_h.extractPartition(G);
    ilp_h.setFirstPartition(G, pid1);



    // output some information about the initial partition
    std::cout << "Info on start graph: " << std::endl;
    std::cout << "Total num nodes in uncoarsened graph: "
              << G.number_of_nodes() << std::endl;
    ilp.graphInfo(qm, G);

    EdgeWeight before_ec = qm.edge_cut(G);
    double before_balance = qm.balance(G);

    auto limit_nonzeroes = (size_t) partition_config.ilp_limit_nonzeroes;
    // compute BFS
    std::unordered_set<NodeID> reachable;
    size_t trees = ilp.computeBFS(G, reachable, partition_config, limit_nonzeroes);
    PartitionID num_blocks = ilp.computeBlocks(G, reachable, partition_config.k);
    PartitionID no_of_coarse_vertices = reachable.size() + num_blocks;

    EdgeWeight base_ec = before_ec;

    graph_access coarser;
    std::vector<bool> coarse_presets(coarser.number_of_nodes(), false);
    if (partition_config.ilp_mode == OptimizationMode::OVERLAP) {
        auto partitions = ilp.createPartitions(G, partition_config);
        size_t best_index;
        std::tie(no_of_coarse_vertices, base_ec, best_index) = ilp.buildOverlapGraph(G, partition_config, pid1,
                                                                                     partitions);

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
    double timelimit = 3600 * 2;

    EdgeWeight after_ec = before_ec;
    double after_balance = before_balance;

    size_t maxcut = 0, maxgain = 0;

    bool improved_solution = false;

    timer ilptime;

    // compute ILP on coarser graph
    ilp.computeIlp(coarser, partition_config, num_blocks,
                   timelimit, pid1, coarse_presets);

    double rawtime = ilptime.elapsed();

        // transfer ILP to graph
        ilp.transferIlp(G, coarser);
        G.set_partition_count(partition_config.k);

        after_ec = qm.edge_cut(G);
        after_balance = qm.balance(G);
        // write the partition to the disc
        std::stringstream Nfilename;


        if (partition_config.filename_output.empty()) {
            Nfilename << graph << "." << partition_config.k << ".ptn";
        } else {
            Nfilename << partition_config.filename_output;
        }

        improved_solution = (after_ec < before_ec ||
                             (after_ec == before_ec
                              && after_balance < before_balance));

        auto pid2 = ilp_h.extractPartition(G);

        std::tie(maxcut, maxgain) = ilp_h.comparePartitionsMoreInfo(G, pid2);

        if (improved_solution) {
            if (after_balance > imbalance_param) {
                std::cout << "BUG - PARTITIONS NOT BALANCED ENOUGH" << std::endl;
            }

            std::cout << "BETTER RESULT! WRITING!" << std::endl;
            graph_io::writePartition(G, parent_dir + "/better_partitions/"
                                        + std::to_string((int) partition_config.imbalance)
                                        + "/" + Nfilename.str() + "_c=" +
                                        std::to_string(after_ec) + "_m=" + std::to_string(after_balance));
        } else {
            std::cout << "Result is not better." << std::endl;
        }

        // output some information about the partition that we have computed
        std::cout << "New partition: " << std::endl;
        //ilp.graphInfo(qm, G);

        rawtime = ilptime.elapsed();

    std::string modestr;

    switch (partition_config.ilp_mode) {
        case OptimizationMode::GAIN :
            modestr = "gain";
            break;
        case OptimizationMode::BOUNDARY :
            modestr = "boundary";
            break;
        case OptimizationMode::TREES :
            modestr = "trees";
            break;
        case OptimizationMode::OVERLAP :
            modestr = "overlap";
    }


    std::ofstream logfile("log.txt", std::ios_base::app | std::ios_base::out);
    std::stringstream log;
    log << "RESULT"
        << " variant=start"
        << " graph=" << graph
        << " k=" << partition_config.k
        << " imbalance=" << (int) partition_config.imbalance
        << " imp_cut=" << (after_ec < before_ec)
        << " imp_sol=" << improved_solution
        << " time=" << t.elapsed()
        << " rawtime=" << rawtime
        << " cut=" << after_ec
        << " actual_imbalance=" << after_balance
        << " t_limit=" << (t.elapsed() >= timelimit)
        << " mode=" << modestr
        << " gain=" << partition_config.ilp_min_gain
        << " depth=" << partition_config.ilp_bfs_depth
        << " trees=" << trees
        << " maxcut=" << maxcut
        << " maxgain=" << maxgain
        << " nonzeroes=" << limit_nonzeroes
        << " overlap_runs=" << partition_config.ilp_overlap_runs
        << " vertices=" << no_of_coarse_vertices
        << " edges=" << (coarser.number_of_edges() / 2)
        << " improved_base=" << (after_ec < base_ec)
        << " base_cut=" << base_ec
        << " optimal=" << partition_config.ilp_optimality
        << std::endl;

    logfile << log.str();
    std::cout << log.str();


    return 0;
}

