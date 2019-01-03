/******************************************************************************
 * cycle_refinement.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef CYCLE_REFINEMENT_JEPIS3F0
#define CYCLE_REFINEMENT_JEPIS3F0

#include "advanced_models.h"
#include "cycle_definitions.h"
#include "definitions.h"
#include "greedy_neg_cycle.h"
#include "random_functions.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "uncoarsening/refinement/refinement.h"

class cycle_refinement : public refinement{
        public:
                cycle_refinement();
                virtual ~cycle_refinement();

                EdgeWeight perform_refinement(PartitionConfig & config, 
                                              graph_access & G, 
                                              complete_boundary & boundary);

        private:
                EdgeWeight greedy_ultra_model(PartitionConfig & partition_config, 
                                              graph_access & G, 
                                              complete_boundary & boundary);

                EdgeWeight greedy_ultra_model_plus(PartitionConfig & partition_config, 
                                                   graph_access & G, 
                                                   complete_boundary & boundary);

                EdgeWeight playfield_algorithm(PartitionConfig & partition_config, 
                                               graph_access & G, 
                                               complete_boundary & boundary);

                advanced_models m_advanced_modelling;
};




#endif /* end of include guard: CYCLE_REFINEMENT_JEPIS3F0 */
