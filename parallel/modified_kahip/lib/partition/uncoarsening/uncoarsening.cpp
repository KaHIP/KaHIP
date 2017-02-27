/******************************************************************************
 * uncoarsening.cpp
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

#include "graph_partition_assertions.h"
#include "misc.h"
#include "quality_metrics.h"
#include "refinement/mixed_refinement.h"
#include "refinement/label_propagation_refinement/label_propagation_refinement.h"
#include "refinement/refinement.h"
#include "separator/vertex_separator_algorithm.h"
#include "uncoarsening.h"


uncoarsening::uncoarsening() {

}

uncoarsening::~uncoarsening() {

}

int uncoarsening::perform_uncoarsening(const PartitionConfig & config, graph_hierarchy & hierarchy) {
        int improvement = 0;

        PartitionConfig cfg     = config;
        refinement* refine      = NULL;

        if(config.label_propagation_refinement) {
                refine      = new label_propagation_refinement();
        } else {
                refine      = new mixed_refinement();
        }

        graph_access * coarsest = hierarchy.get_coarsest();
        PRINT(std::cout << "log>" << "unrolling graph with " << coarsest->number_of_nodes() << std::endl;)

        complete_boundary* finer_boundary   = NULL;
        complete_boundary* coarser_boundary = NULL;
        if(!config.label_propagation_refinement) {
                coarser_boundary = new complete_boundary(coarsest);
                coarser_boundary->build();
        }
        double factor = config.balance_factor;
        cfg.upper_bound_partition = ((!hierarchy.isEmpty()) * factor +1.0)*config.upper_bound_partition;
        improvement += (int)refine->perform_refinement(cfg, *coarsest, *coarser_boundary);

        NodeID coarser_no_nodes = coarsest->number_of_nodes();
        graph_access* finest    = NULL;
        graph_access* to_delete = NULL;
        unsigned int hierarchy_deepth = hierarchy.size();

        while(!hierarchy.isEmpty()) {
                graph_access* G = hierarchy.pop_finer_and_project();

                PRINT(std::cout << "log>" << "unrolling graph with " << G->number_of_nodes()<<  std::endl;)
                
                if(!config.label_propagation_refinement) {
                        finer_boundary = new complete_boundary(G); 
                        finer_boundary->build_from_coarser(coarser_boundary, coarser_no_nodes, hierarchy.get_mapping_of_current_finer());
                }

                //call refinement
                double cur_factor = factor/(hierarchy_deepth-hierarchy.size());
                cfg.upper_bound_partition = ((!hierarchy.isEmpty()) * cur_factor+1.0)*config.upper_bound_partition;
                PRINT(std::cout <<  "cfg upperbound " <<  cfg.upper_bound_partition  << std::endl;)
                improvement += (int)refine->perform_refinement(cfg, *G, *finer_boundary);
                ASSERT_TRUE(graph_partition_assertions::assert_graph_has_kway_partition(config, *G));

                if(config.use_balance_singletons && !config.label_propagation_refinement) {
                        finer_boundary->balance_singletons( config, *G );
                }

                // update boundary pointers
                if(!config.label_propagation_refinement) delete coarser_boundary;
                coarser_boundary = finer_boundary;
                coarser_no_nodes = G->number_of_nodes();

		//clean up 
		if(to_delete != NULL) {
			delete to_delete;
		}
		if(!hierarchy.isEmpty()) {
			to_delete = G;
		}

                finest = G;
        }

        if(config.compute_vertex_separator) {
               PRINT(std::cout <<  "now computing a vertex separator from the given edge separator"  << std::endl;)
               vertex_separator_algorithm vsa;
               vsa.compute_vertex_separator(config, *finest, *finer_boundary); 
        }

        delete refine;
        if(finer_boundary != NULL) delete finer_boundary;
	delete coarsest;

        return improvement;
}


