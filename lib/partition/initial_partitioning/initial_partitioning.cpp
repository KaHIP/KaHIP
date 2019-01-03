/******************************************************************************
 * initial_partitioning.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <math.h>

#include "bipartition.h"
#include "graph_partition_assertions.h"
#include "graph_partitioner.h"
#include "initial_partition_bipartition.h"
#include "initial_partitioning.h"
#include "initial_refinement/initial_refinement.h"
#include "initial_node_separator.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"

initial_partitioning::initial_partitioning() {

}

initial_partitioning::~initial_partitioning() {

}

void initial_partitioning::perform_initial_partitioning(const PartitionConfig & config, graph_hierarchy & hierarchy) {
        graph_access& G = *hierarchy.get_coarsest();
        if(config.mode_node_separators) {
                perform_initial_partitioning_separator(config, G);
        } else {
                perform_initial_partitioning(config, G);
        }
}


void initial_partitioning::perform_initial_partitioning(const PartitionConfig & config, graph_access &  G) {

        initial_partitioner* partition = NULL;
        switch(config.initial_partitioning_type) {
                case INITIAL_PARTITIONING_RECPARTITION:
                        partition = new initial_partition_bipartition();
                        break;
                case INITIAL_PARTITIONING_BIPARTITION:
                        partition = new bipartition();
                        break;


        }       

        quality_metrics qm;
        EdgeWeight best_cut;
        int* best_map = new int[G.number_of_nodes()];
        if(config.graph_allready_partitioned && !config.omit_given_partitioning) {
                best_cut = qm.edge_cut(G);
                forall_nodes(G, n) {
                        best_map[n] = G.getPartitionIndex(n); 
                } endfor
        } else {
                best_cut = std::numeric_limits<EdgeWeight>::max();
        }
        
        timer t;
        t.restart();
        int* partition_map  = new int[G.number_of_nodes()];
        unsigned reps_to_do = (unsigned) std::max((int)ceil(config.initial_partitioning_repetitions/(double)log2(config.k)),2);
         
        if(config.initial_partitioning_repetitions == 0) {
                reps_to_do = 1;
        }
        if(config.eco) {
                //bound the number of initial partitioning repetions
                reps_to_do = std::min((int)config.minipreps, (int)reps_to_do);
        }

        PRINT(std::cout << "no of initial partitioning repetitions = " << reps_to_do                     << std::endl;);
        PRINT(std::cout << "no of nodes for partition = "              << G.number_of_nodes()            << std::endl;);
        if(!((config.graph_allready_partitioned && config.no_new_initial_partitioning) || config.omit_given_partitioning)) {
                for(unsigned int rep = 0; rep < reps_to_do; rep++) {
                        unsigned seed = random_functions::nextInt(0, std::numeric_limits<int>::max()); 
                        PartitionConfig working_config = config;
                        working_config.combine = false;
                        partition->initial_partition(working_config, seed, G, partition_map);
                        
                        EdgeWeight cur_cut = qm.edge_cut(G, partition_map); 
                        if(cur_cut < best_cut) {
                                PRINT(std::cout << "log>" << "improved the current initial partitiong from " << best_cut 
                                                << " to " << cur_cut  << std::endl;)

                                forall_nodes(G, n) {
                                        best_map[n] = partition_map[n];
                                } endfor

                                best_cut = cur_cut; 
                                if(best_cut == 0) break;
                        }
                }

                forall_nodes(G, n) {
                        G.setPartitionIndex(n,best_map[n]); 
                } endfor
        }

        G.set_partition_count(config.k);

        PRINT(std::cout << "initial partitioning took " << t.elapsed()                << std::endl;)
        PRINT(std::cout << "log>"                       << "current initial balance " << qm.balance(G) << std::endl;)

        if(config.initial_partition_optimize || config.combine) {
                initial_refinement iniref;
                iniref.optimize(config, G, best_cut);
        }

        PRINT(std::cout << "log>" << "final current initial partitiong from " << best_cut
                        << " to " << best_cut                                 << std::endl;)

        if(!(config.graph_allready_partitioned && config.no_new_initial_partitioning)) {
                PRINT(std::cout << "finalinitialcut " << best_cut                         << std::endl;)
                PRINT(std::cout << "log>"             << "final current initial balance " << qm.balance(G) << std::endl;)
        }

        ASSERT_TRUE(graph_partition_assertions::assert_graph_has_kway_partition(config, G));

        delete[] partition_map;
        delete[] best_map;
        delete partition;
}

void initial_partitioning::perform_initial_partitioning_separator(const PartitionConfig & config, graph_access &  G) {
        initial_node_separator ipns;
        ipns.compute_node_separator(config,G);
}

