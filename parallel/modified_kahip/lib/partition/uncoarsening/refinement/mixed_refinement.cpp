/******************************************************************************
 * mixed_refinement.cpp
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

#include "cycle_improvements/cycle_refinement.h"
#include "kway_graph_refinement/kway_graph_refinement.h"
#include "kway_graph_refinement/multitry_kway_fm.h"
#include "mixed_refinement.h"
#include "quotient_graph_refinement/quotient_graph_refinement.h"

mixed_refinement::mixed_refinement() {

}

mixed_refinement::~mixed_refinement() {

}

EdgeWeight mixed_refinement::perform_refinement(PartitionConfig & config, graph_access & G, complete_boundary & boundary) {
        refinement* refine              = new quotient_graph_refinement();
        refinement* kway                = new kway_graph_refinement();
        multitry_kway_fm* multitry_kway = new multitry_kway_fm();
        cycle_refinement* cycle_refine  = new cycle_refinement();

        EdgeWeight overall_improvement = 0; 
        //call refinement
        if(config.no_change_convergence) {
                bool sth_changed = true;
                while(sth_changed) {
                        EdgeWeight improvement = 0;
                        if(config.corner_refinement_enabled) {
                                improvement += kway->perform_refinement(config, G, boundary);
                        }

                        if(!config.quotient_graph_refinement_disabled) {
                                improvement += refine->perform_refinement(config, G, boundary);
                        }

                        overall_improvement += improvement;
                        sth_changed = improvement != 0;
                }

        } else {
                if(config.corner_refinement_enabled) {
                        overall_improvement += kway->perform_refinement(config, G, boundary);
                } 

                if(!config.quotient_graph_refinement_disabled) {
                        overall_improvement += refine->perform_refinement(config, G, boundary);
                }

                if(config.kaffpa_perfectly_balanced_refinement) {
                        overall_improvement += cycle_refine->perform_refinement(config, G, boundary);
                }
        }

        delete refine;
        delete kway;
        delete multitry_kway;
        delete cycle_refine;

        return overall_improvement;
}

