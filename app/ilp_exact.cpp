/******************************************************************************
 * ilp_exact.cpp
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
#include "ilp_improve/ilp_exact.h"
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

        ilp_exact ilp;
        graph_access G;

        timer t;
        graph_io::readGraphWeighted(G, graph_filename);
        std::cout << "io time: " << t.elapsed() << std::endl;
        G.set_partition_count(partition_config.k);

        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        // compute BFS

        timer ilptime;
        ilptime.restart();

        // compute ILP on coarser graph
        ilp.computeIlp(G, partition_config);

        std::cout <<  "ILP took " << ilptime.elapsed() << std::endl;

        quality_metrics qm;
        // ******************************* done partitioning *****************************************       
        // output some information about the partition that we have computed 
        std::cout << "cut \t\t"         << qm.edge_cut(G)                 << std::endl;
        std::cout << "balance \t"       << qm.balance(G)                  << std::endl;

        std::stringstream filename;
        if(!partition_config.filename_output.compare("")) {
                filename << "tmppartition" << partition_config.k << ".exact";
        } else {
                filename << partition_config.filename_output;
        }

        graph_io::writePartition(G, filename.str());

        return 0;
}

