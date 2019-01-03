/******************************************************************************
 * gal_combine.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <fstream>
#include <iostream>
#include <sstream>

#include "construct_partition.h"
#include "gal_combine.h"
#include "graph_partitioner.h"
#include "random_functions.h"
#include "uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "uncoarsening/refinement/mixed_refinement.h"
#include "uncoarsening/refinement/refinement.h"
#include "uncoarsening/refinement/tabu_search/tabu_search.h"

gal_combine::gal_combine() {
                
}

gal_combine::~gal_combine() {
                
}


// implements our version of gal combine
// compute a matching between blocks (greedily) 
// extend to partition
// apply our refinements and tabu search
void gal_combine::perform_gal_combine( PartitionConfig & config, graph_access & G) {
        //first greedily compute a matching of the partitions
        std::vector< std::unordered_map<PartitionID, unsigned> > counters(config.k);
        forall_nodes(G, node) {
                //boundary_pair bp;
                if(counters[G.getPartitionIndex(node)].find(G.getSecondPartitionIndex(node)) != counters[G.getPartitionIndex(node)].end()) {
                        counters[G.getPartitionIndex(node)][G.getSecondPartitionIndex(node)] += 1;
                } else {
                        counters[G.getPartitionIndex(node)][G.getSecondPartitionIndex(node)] = 1;
                }
        } endfor

        std::vector< PartitionID > permutation(config.k);
        for( unsigned i = 0; i < permutation.size(); i++) {
                permutation[i] = i;
        }

        random_functions::permutate_vector_good_small(permutation); 
        std::vector<bool> rhs_matched(config.k, false);
        std::vector<PartitionID> bipartite_matching(config.k);
        for( unsigned i = 0; i < permutation.size(); i++) {
                PartitionID cur_partition   = permutation[i];
                PartitionID best_unassigned = config.k;
                NodeWeight  best_value      = 0;

                for( std::unordered_map<PartitionID, unsigned>::iterator it = counters[cur_partition].begin(); 
                     it != counters[cur_partition].end(); ++it) {
                        if( rhs_matched[it->first] == false && it->second > best_value ) {
                                best_unassigned = it->first; 
                                best_value = it->second;
                        }
                }

                bipartite_matching[cur_partition] = best_unassigned;
                if( best_unassigned != config.k ) {
                        rhs_matched[best_unassigned] = true;
                }
        }

        std::vector<bool> blocked_vertices(G.number_of_nodes(), false);
        forall_nodes(G, node) {
                if( bipartite_matching[G.getPartitionIndex(node)] == G.getSecondPartitionIndex(node) ){
                        blocked_vertices[node] = true;
                } else {
                        // we will reassign this vertex since the partitions do not agree on it
                        G.setPartitionIndex(node, config.k); 
                }
        } endfor

        construct_partition cp;
        cp.construct_starting_from_partition( config, G );

        refinement* refine = new mixed_refinement();

        double real_epsilon = config.imbalance/100.0;
        double epsilon = random_functions::nextDouble(real_epsilon+0.005,real_epsilon+config.kabaE_internal_bal); 
        PartitionConfig copy = config;
        copy.upper_bound_partition = (1+epsilon)*ceil(config.work_load/(double)config.k);

        complete_boundary boundary(&G);
        boundary.build();

        tabu_search ts;
        ts.perform_refinement( copy, G, boundary);
        
        //now obtain the quotient graph
        complete_boundary boundary2(&G);
        boundary2.build();

        copy = config;
        copy.upper_bound_partition = (1+epsilon)*ceil(config.work_load/(double)config.k);

        refine->perform_refinement( copy, G, boundary2);

        copy = config;
        cycle_refinement cr;
        cr.perform_refinement(config, G, boundary2);
        delete refine;

        
}
