/******************************************************************************
 * wcycle_partitioner.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "coarsening/coarsening_configurator.h"
#include "coarsening/contraction.h"
#include "coarsening/edge_rating/edge_ratings.h"
#include "coarsening/matching/gpa/gpa_matching.h"
#include "coarsening/matching/random_matching.h"
#include "coarsening/stop_rules/stop_rules.h"
#include "data_structure/graph_hierarchy.h"
#include "definitions.h"
#include "graph_partition_assertions.h"
#include "initial_partitioning/initial_partitioning.h"
#include "misc.h"
#include "random_functions.h"
#include "uncoarsening/refinement/mixed_refinement.h"
#include "uncoarsening/refinement/label_propagation_refinement/label_propagation_refinement.h"
#include "uncoarsening/refinement/refinement.h"
#include "wcycle_partitioner.h"

int wcycle_partitioner::perform_partitioning(const PartitionConfig & config, graph_access & G) {
        PartitionConfig  cfg = config; 

        if(config.stop_rule == STOP_RULE_SIMPLE) {
                m_coarsening_stop_rule = new simple_stop_rule(cfg, G.number_of_nodes());
        } else {
                m_coarsening_stop_rule = new multiple_k_stop_rule(cfg, G.number_of_nodes());
        }

        int improvement = (int) perform_partitioning_recursive(cfg, G, NULL); 
        delete m_coarsening_stop_rule;

        return improvement;
}

int wcycle_partitioner::perform_partitioning_recursive( PartitionConfig & partition_config, 
                                                        graph_access & G, 
                                                        complete_boundary ** c_boundary) {

        //if graph not small enough
        //      perform matching two times
        //      perform coarsening two times
        //      call rekursive
        //else 
        //      initial partitioning
        //
        //refinement
        NodeID no_of_coarser_vertices = G.number_of_nodes();
        NodeID no_of_finer_vertices   = G.number_of_nodes();
        int improvement = 0;

        edge_ratings rating(partition_config);
        CoarseMapping* coarse_mapping =  new CoarseMapping();

        graph_access* finer                      = &G;
        matching* edge_matcher                   = NULL;
        contraction* contracter                  = new contraction();
        PartitionConfig copy_of_partition_config = partition_config;
        graph_access* coarser                    = new graph_access();

        Matching edge_matching;
        NodePermutationMap permutation;

        coarsening_configurator coarsening_config;
        coarsening_config.configure_coarsening(partition_config, &edge_matcher, m_level);
        
        rating.rate(*finer, m_level);

        edge_matcher->match(partition_config, *finer, edge_matching, *coarse_mapping, no_of_coarser_vertices, permutation);
        delete edge_matcher; 

        if(partition_config.graph_allready_partitioned) {
                contracter->contract_partitioned(partition_config, *finer, 
                                                 *coarser, edge_matching, 
                                                 *coarse_mapping, no_of_coarser_vertices, 
                                                 permutation);
        } else {
                contracter->contract(partition_config, *finer, 
                                     *coarser, edge_matching, 
                                     *coarse_mapping, no_of_coarser_vertices, 
                                     permutation);
        }

        coarser->set_partition_count(partition_config.k);
        complete_boundary* coarser_boundary =  NULL;
        refinement* refine = NULL;

        if(!partition_config.label_propagation_refinement) {
                coarser_boundary = new complete_boundary(coarser);
                refine = new mixed_refinement();
        } else {
                refine = new label_propagation_refinement();
        }

        if(!m_coarsening_stop_rule->stop(no_of_finer_vertices, no_of_coarser_vertices)) {

                PartitionConfig cfg; cfg = partition_config;

                double factor = partition_config.balance_factor;
                cfg.upper_bound_partition = (factor +1.0)*partition_config.upper_bound_partition;

	        initial_partitioning init_part;
		init_part.perform_initial_partitioning(cfg, *coarser);

                if(!partition_config.label_propagation_refinement) coarser_boundary->build();

                //PRINT(std::cout <<  "upper bound " <<  cfg.upper_bound_partition  << std::endl;)
                improvement += refine->perform_refinement(cfg, *coarser, *coarser_boundary);
                m_deepest_level = m_level + 1;
        } else {
                m_level++;

                improvement += perform_partitioning_recursive( partition_config, *coarser, &coarser_boundary); 
                partition_config.graph_allready_partitioned = true;

                if(m_level % partition_config.level_split == 0 ) {

                        if(!partition_config.use_fullmultigrid 
                        || m_have_been_level_down.find(m_level) == m_have_been_level_down.end())  { 

                                if(!partition_config.label_propagation_refinement) {
                                        delete coarser_boundary;

                                        coarser_boundary                = new complete_boundary(coarser);
                                }
                                m_have_been_level_down[m_level] = true;

                                // configurate the algorithm to use the same amount
                                // of imbalance as was allowed on this level 
                                PartitionConfig cfg;
                                cfg = partition_config;
                                cfg.set_upperbound = false;

                                double cur_factor = partition_config.balance_factor/(m_deepest_level-m_level);
                                cfg.upper_bound_partition = ( (m_level != 0) * cur_factor+1.0)*partition_config.upper_bound_partition;

                                // do the next arm of the F-cycle
                                improvement += perform_partitioning_recursive( cfg, *coarser, &coarser_boundary); 
                        }
                }

                m_level--;

        }

        if(partition_config.use_balance_singletons && !partition_config.label_propagation_refinement) {
                coarser_boundary->balance_singletons( partition_config, *coarser );
        }
        
        //project
        graph_access& fRef = *finer;
        graph_access& cRef = *coarser;
        forall_nodes(fRef, n) {
                NodeID coarser_node              = (*coarse_mapping)[n];
                PartitionID coarser_partition_id = cRef.getPartitionIndex(coarser_node);
                fRef.setPartitionIndex(n, coarser_partition_id);
        } endfor

        finer->set_partition_count(coarser->get_partition_count());
        complete_boundary* current_boundary = NULL;
        if(!partition_config.label_propagation_refinement) {
                current_boundary = new complete_boundary(finer);
                current_boundary->build_from_coarser(coarser_boundary, no_of_coarser_vertices, coarse_mapping ); 
        }

        PartitionConfig cfg; cfg = partition_config;
        double cur_factor = partition_config.balance_factor/(m_deepest_level-m_level);

        //only set the upperbound if it is the first time 
        //we go down the F-cycle
        if( partition_config.set_upperbound ) {
                cfg.upper_bound_partition = ( (m_level != 0) * cur_factor+1.0)*partition_config.upper_bound_partition;
        } else {
                cfg.upper_bound_partition = partition_config.upper_bound_partition;
        }

        //PRINT(std::cout <<  "upper bound " <<  cfg.upper_bound_partition  << std::endl;)
        improvement += refine->perform_refinement(cfg, *finer, *current_boundary);

        if(c_boundary != NULL) {
                delete *c_boundary;
                *c_boundary = current_boundary;
        } else {
		if( current_boundary != NULL ) delete current_boundary;
	}

        //std::cout <<  "finer " <<  no_of_finer_vertices  << std::endl;
        delete contracter;
        delete coarse_mapping;
        delete coarser_boundary;
        delete coarser;
        delete refine;

        return improvement;
}
