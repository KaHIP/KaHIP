/******************************************************************************
 * cycle_refinement.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
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
