/******************************************************************************
 * coarsening_configurator.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
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

#ifndef COARSENING_CONFIGURATOR_8UJ78WYS
#define COARSENING_CONFIGURATOR_8UJ78WYS

#include "contraction.h"
#include "data_structure/graph_hierarchy.h"
#include "definitions.h"
#include "edge_rating/edge_ratings.h"
#include "matching/gpa/gpa_matching.h"
#include "matching/random_matching.h"
#include "clustering/size_constraint_label_propagation.h"
#include "stop_rules/stop_rules.h"

class coarsening_configurator {
        public:
                coarsening_configurator( ) {};
                virtual ~coarsening_configurator() {};

                void configure_coarsening(const PartitionConfig & partition_config, 
                                          matching** edge_matcher, 
                                          unsigned level); 
};

inline void coarsening_configurator::configure_coarsening( const PartitionConfig & partition_config, 
                                                           matching** edge_matcher, 
                                                           unsigned level) {

        switch(partition_config.matching_type) {
                case MATCHING_RANDOM: 
                        *edge_matcher = new random_matching();
                        break; 
                case MATCHING_GPA:
                        *edge_matcher = new gpa_matching();
                        PRINT(std::cout <<  "gpa matching"  << std::endl;)
                        break;
                case MATCHING_RANDOM_GPA:
                        PRINT(std::cout <<  "random gpa matching"  << std::endl;)
                        *edge_matcher = new gpa_matching();
                        break;
               case CLUSTER_COARSENING:
                        PRINT(std::cout <<  "cluster_coarsening"  << std::endl;)
                        *edge_matcher = new size_constraint_label_propagation();
                        break;

        }

        if( partition_config.matching_type == MATCHING_RANDOM_GPA 
            && level < partition_config.aggressive_random_levels) {

                delete *edge_matcher;
                PRINT(std::cout <<  "random matching"  << std::endl;)
                *edge_matcher = new random_matching();
        }  
}

#endif /* end of include guard: COARSENING_CONFIGURATOR_8UJ78WYS */
