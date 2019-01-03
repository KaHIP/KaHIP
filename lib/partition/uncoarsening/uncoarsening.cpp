/******************************************************************************
 * uncoarsening.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "graph_partition_assertions.h"
#include "misc.h"
#include "quality_metrics.h"
#include "refinement/mixed_refinement.h"
#include "refinement/node_separators/greedy_ns_local_search.h"
#include "refinement/node_separators/fm_ns_local_search.h"
#include "refinement/node_separators/localized_fm_ns_local_search.h"
#include "refinement/label_propagation_refinement/label_propagation_refinement.h"
#include "refinement/refinement.h"
#include "separator/vertex_separator_algorithm.h"
#include "tools/random_functions.h"
#include "uncoarsening.h"


uncoarsening::uncoarsening() {

}

uncoarsening::~uncoarsening() {

}

int uncoarsening::perform_uncoarsening(const PartitionConfig & config, graph_hierarchy & hierarchy) {

        if(config.mode_node_separators) {
                if( config.faster_ns ) {
                        return perform_uncoarsening_nodeseparator_fast(config, hierarchy);
                } else {
                        return perform_uncoarsening_nodeseparator(config, hierarchy);
                }
        } else {
                return perform_uncoarsening_cut(config, hierarchy);
        }
}

int uncoarsening::perform_uncoarsening_cut(const PartitionConfig & config, graph_hierarchy & hierarchy) {
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

int uncoarsening::perform_uncoarsening_nodeseparator(const PartitionConfig & config, graph_hierarchy & hierarchy) {

        std::cout <<  "log> starting uncoarsening ---------------"  << std::endl;
        PartitionConfig cfg     = config;
        graph_access * coarsest = hierarchy.get_coarsest();
        quality_metrics qm;
        std::cout << "log>" << "unrolling graph with " << coarsest->number_of_nodes() << std::endl;

        if( !config.sep_fm_disabled ) {
                for( int i = 0; i < config.sep_num_fm_reps; i++) {
                        fm_ns_local_search fmnsls;
                        fmnsls.perform_refinement(config, (*coarsest));

                        int rnd_block = random_functions::nextInt(0,1);
                        fmnsls.perform_refinement(config, (*coarsest),true, rnd_block);
                        fmnsls.perform_refinement(config, (*coarsest),true, rnd_block == 0 ? 1 : 0);
                }
        }

        if( !config.sep_flows_disabled ) {
                for( int i = 0; i < config.max_flow_improv_steps; i++) {

                        vertex_separator_algorithm vsa;

                        std::vector<NodeID> separator;
                        forall_nodes((*coarsest), node) {
                                if( coarsest->getPartitionIndex(node) == 2) {
                                        separator.push_back(node);
                                }
                        } endfor

                        std::vector<NodeID> output_separator;
                        NodeWeight improvement = vsa.improve_vertex_separator(config, *coarsest, separator, output_separator);
                        if(improvement == 0) break;
                }
        }

        graph_access* to_delete = NULL;
        while(!hierarchy.isEmpty()) {
                graph_access* G = hierarchy.pop_finer_and_project();
                std::cout << "log>" << "unrolling graph with " << G->number_of_nodes() << std::endl;

                if( !config.sep_fm_disabled) {
                        for( int i = 0; i < config.sep_num_fm_reps; i++) {
                                fm_ns_local_search fmnsls;
                                fmnsls.perform_refinement(config, (*G));

                                int rnd_block = random_functions::nextInt(0,1);
                                fmnsls.perform_refinement(config, (*G), true, rnd_block);
                                fmnsls.perform_refinement(config, (*G), true, rnd_block == 0? 1 : 0);
                        }
                }

                if( !config.sep_loc_fm_disabled) {
                        for( int i = 0; i < config.sep_num_loc_fm_reps; i++) {
                                localized_fm_ns_local_search fmnsls;
                                fmnsls.perform_refinement(config, (*G));

                                int rnd_block = random_functions::nextInt(0,1);
                                fmnsls.perform_refinement(config, (*G), true, rnd_block);
                                fmnsls.perform_refinement(config, (*G), true, rnd_block == 0? 1 : 0);
                        }
                }


                if( !config.sep_flows_disabled ) {
                        for( int i = 0; i < config.max_flow_improv_steps; i++) {
                                vertex_separator_algorithm vsa;

                                std::vector<NodeID> separator;
                                forall_nodes((*G), node) {
                                        if( G->getPartitionIndex(node) == 2) {
                                                separator.push_back(node);
                                        }
                                } endfor

                                std::vector<NodeID> output_separator;
                                NodeWeight improvement = vsa.improve_vertex_separator(config, *G, separator, output_separator);
                                if(improvement == 0) break;
                        }
                }
                if(to_delete != NULL) {
			delete to_delete;
		}
		if(!hierarchy.isEmpty()) {
			to_delete = G;
		}
        }
	delete coarsest;

        return 0;
}

int uncoarsening::perform_uncoarsening_nodeseparator_fast(const PartitionConfig & config, graph_hierarchy & hierarchy) {

        std::cout <<  "log> starting uncoarsening ---------------"  << std::endl;
        PartitionConfig cfg     = config;
        graph_access * coarsest = hierarchy.get_coarsest();
        std::cout << "log>" << "unrolling graph with " << coarsest->number_of_nodes() << std::endl;

        std::vector< NodeWeight > block_weights(3,0); PartialBoundary current_separator;
        //compute coarsest block weights and separator
        forall_nodes((*coarsest), node) {
                block_weights[coarsest->getPartitionIndex(node)] += coarsest->getNodeWeight(node);
                if( coarsest->getPartitionIndex(node) == 2) {
                        current_separator.insert(node);
                }
        } endfor
        

        std::vector< bool > moved_out_of_S(coarsest->number_of_nodes(), false);
        if( !config.sep_fm_disabled ) {
                for( int i = 0; i < config.sep_num_fm_reps; i++) {
                        fm_ns_local_search fmnsls;
                        fmnsls.perform_refinement(config, (*coarsest), block_weights, moved_out_of_S, current_separator);

                        int rnd_block = random_functions::nextInt(0,1);
                        fmnsls.perform_refinement(config, (*coarsest), block_weights, moved_out_of_S, current_separator, true, rnd_block);
                        fmnsls.perform_refinement(config, (*coarsest), block_weights, moved_out_of_S, current_separator, true, rnd_block == 0? 1 : 0);
                }
        }

        if( !config.sep_flows_disabled ) {
                for( int i = 0; i < config.max_flow_improv_steps; i++) {

                        vertex_separator_algorithm vsa;
                        std::vector<NodeID> output_separator;
                        NodeWeight improvement = vsa.improve_vertex_separator(config, *coarsest, block_weights, current_separator);
                        if(improvement == 0) break;
                }
        }

        graph_access* to_delete = NULL;
        while(!hierarchy.isEmpty()) {
                graph_access* G = hierarchy.pop_finer_and_project_ns(current_separator);
                std::cout << "log>" << "unrolling graph with " << G->number_of_nodes() << std::endl;

                std::vector< bool > moved_out_of_S(G->number_of_nodes(), false);
                if( !config.sep_fm_disabled) {
                        for( int i = 0; i < config.sep_num_fm_reps; i++) {
                                
                                fm_ns_local_search fmnsls;
                                NodeWeight improvement = 0;
                                improvement += fmnsls.perform_refinement(config, (*G), block_weights, moved_out_of_S, current_separator);

                                int rnd_block = random_functions::nextInt(0,1);
                                improvement += fmnsls.perform_refinement(config, (*G), block_weights, moved_out_of_S, current_separator, true, rnd_block);
                                improvement += fmnsls.perform_refinement(config, (*G), block_weights, moved_out_of_S, current_separator, true, rnd_block == 0 ? 1 : 0);
                                if( improvement == 0 ) break;
                        }
                }

                if( !config.sep_loc_fm_disabled) {
                        for( int i = 0; i < config.sep_num_loc_fm_reps; i++) {
                                localized_fm_ns_local_search fmnsls;
                                NodeWeight improvement = 0;
                                improvement += fmnsls.perform_refinement(config, (*G), block_weights, moved_out_of_S, current_separator);

                                int rnd_block = random_functions::nextInt(0,1);
                                improvement += fmnsls.perform_refinement(config, (*G), block_weights, moved_out_of_S, current_separator, true, rnd_block);
                                improvement += fmnsls.perform_refinement(config, (*G), block_weights, moved_out_of_S, current_separator, true, rnd_block == 0 ? 1 : 0);
                        }
                }

                if( !config.sep_flows_disabled ) {
                        for( int i = 0; i < config.max_flow_improv_steps; i++) {
                                vertex_separator_algorithm vsa;
                                std::vector<NodeID> output_separator;
                                NodeWeight improvement = vsa.improve_vertex_separator(config, *G, block_weights, current_separator);
                                if(improvement == 0) break;
                        }
                }
                if(to_delete != NULL) {
			delete to_delete;
		}
		if(!hierarchy.isEmpty()) {
			to_delete = G;
		}
        }
	delete coarsest;

        return 0;
}

