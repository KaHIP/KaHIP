/******************************************************************************
 * cycle_refinement.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>

#include "algorithms/cycle_search.h"
#include "augmented_Qgraph_fabric.h"
#include "cycle_refinement.h"
#include "quality_metrics.h"

cycle_refinement::cycle_refinement() {

}

cycle_refinement::~cycle_refinement() {

}

EdgeWeight cycle_refinement::perform_refinement(PartitionConfig & partition_config, 
                graph_access & G, 
                complete_boundary & boundary) {
        Gain overall_gain = 0;
        PartitionConfig copy = partition_config;

        switch(partition_config.cycle_refinement_algorithm) {
                case CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL:
                        overall_gain = greedy_ultra_model(copy, G, boundary);
                        break;
                case CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL_PLUS:
                        overall_gain = greedy_ultra_model_plus(copy, G, boundary);
                        break;
                case CYCLE_REFINEMENT_ALGORITHM_PLAYFIELD:
                        //dropbox for new algorithms
                        overall_gain = playfield_algorithm(copy, G, boundary);
                        break;
        }

        return overall_gain;
}


EdgeWeight cycle_refinement::playfield_algorithm(PartitionConfig & partition_config, 
                                graph_access & G, 
                                complete_boundary & boundary) {
        greedy_ultra_model(partition_config, G, boundary);
        greedy_ultra_model_plus(partition_config, G, boundary);
        return 0;
}

EdgeWeight cycle_refinement::greedy_ultra_model(PartitionConfig & partition_config, 
                graph_access & G, 
                complete_boundary & boundary) {

        augmented_Qgraph_fabric augmented_fabric;
        unsigned s             = partition_config.kaba_internal_no_aug_steps_aug;
        bool something_changed = false;
        bool overloaded        = false;
        unsigned unsucc_count  = 0;

        do {
                augmented_Qgraph aqg;
                augmented_fabric.build_augmented_quotient_graph(partition_config, G, boundary, aqg, s, false);
                something_changed = m_advanced_modelling.compute_vertex_movements_ultra_model(partition_config, 
                                                                                              G,
                                                                                              boundary, 
                                                                                              aqg, 
                                                                                              s, 
                                                                                              false);
                if( something_changed ) {
                        unsucc_count = 0;
                } else {
                        unsucc_count++;
                }

                if(unsucc_count > 2 
                && unsucc_count <= partition_config.kaba_unsucc_iterations 
                && partition_config.kaba_enable_zero_weight_cycles) {
                        something_changed = m_advanced_modelling.compute_vertex_movements_ultra_model(partition_config, 
                                                                                                      G,
                                                                                                      boundary, 
                                                                                                      aqg, 
                                                                                                      s, 
                                                                                                      true);
                }

                if(unsucc_count >= partition_config.kaba_unsucc_iterations ) {
                        graph_access G_bar;
                        boundary.getUnderlyingQuotientGraph(G_bar); 
                        overloaded = false;
                        forall_nodes(G_bar, block) {
                                if(boundary.getBlockWeight(block) > partition_config.upper_bound_partition ) {
                                        overloaded = true;
                                        break;
                                }
                        } endfor
                  
                        if(overloaded) {
                                augmented_Qgraph aqg_rebal;
                                bool movs_allready_performed = augmented_fabric.build_augmented_quotient_graph(partition_config, 
                                                                                                               G, 
                                                                                                               boundary, 
                                                                                                               aqg_rebal, 
                                                                                                               s, true);
                                if(!movs_allready_performed) {
                                        m_advanced_modelling.compute_vertex_movements_rebalance(partition_config, 
                                                                                                G, 
                                                                                                boundary, 
                                                                                                aqg_rebal, s);
                                } // else the fall back solution has been applied
                        }
                }
        } while(unsucc_count < partition_config.kaba_unsucc_iterations || (overloaded));

        return 0; 
}

EdgeWeight cycle_refinement::greedy_ultra_model_plus(PartitionConfig & partition_config, 
                graph_access & G, 
                complete_boundary & boundary) {
        unsigned s             = partition_config.kaba_internal_no_aug_steps_aug;
        bool something_changed = false;
        bool overloaded        = false;
        
       
        augmented_Qgraph_fabric augmented_fabric;
        bool first_level = true;
        forall_nodes(G, node) {
                if(G.getNodeWeight(node) != 1) {
                        first_level = false;
                        break;
                }
        } endfor

        int unsucc_count = 0;
        do {
                augmented_Qgraph aqg;
                augmented_fabric.build_augmented_quotient_graph(partition_config, G, boundary, aqg, s, false, true);
                something_changed = m_advanced_modelling.compute_vertex_movements_ultra_model(partition_config, 
                                                                                              G, 
                                                                                              boundary, 
                                                                                              aqg, 
                                                                                              s,
                                                                                              false);

                if( something_changed ) {
                        unsucc_count = 0;
                } else {
                        unsucc_count++;
                }

                if(unsucc_count > 2 && unsucc_count < 19) {
                        something_changed = m_advanced_modelling.compute_vertex_movements_ultra_model(partition_config, 
                                                                                                      G, 
                                                                                                      boundary, 
                                                                                                      aqg, 
                                                                                                      s, true);
                }

                if(unsucc_count > 19 && first_level) {
                        graph_access G_bar;
                        boundary.getUnderlyingQuotientGraph(G_bar); 
                        overloaded = false;
                        forall_nodes(G_bar, block) {
                                if(boundary.getBlockWeight(block) > partition_config.upper_bound_partition ) {
                                        overloaded = true;
                                        break;
                                }
                        } endfor
                  
                        if(overloaded) {
                                augmented_Qgraph aqg_rebal;
                                bool moves_performed = augmented_fabric.build_augmented_quotient_graph(partition_config, 
                                                                                                       G, 
                                                                                                       boundary, 
                                                                                                       aqg_rebal, 
                                                                                                       s, true, true);
                                if(!moves_performed) {
                                        m_advanced_modelling.compute_vertex_movements_rebalance(partition_config, 
                                                                                                G, boundary, 
                                                                                                aqg_rebal, s);
                                } // else the fall back solution has been applied
                        }
                       
                }
        } while(unsucc_count < 20 || overloaded);

        return 0; 
}
