/******************************************************************************
 * mixed_refinement.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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

