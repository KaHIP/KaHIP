/******************************************************************************
 * uncoarsening.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <queue>
#include <limits>

#include "graph_partition_assertions.h"
#include "misc.h"
#include "quality_metrics.h"
#include "refinement/connectivity_check.h"
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

// Eliminate disconnected components by moving them to the lightest adjacent block,
// then greedy-rebalance by moving non-articulation boundary nodes from overweight blocks.
void eliminate_and_rebalance(const PartitionConfig & config, graph_access & G) {
        // Step 1: Eliminate disconnected components
        for(PartitionID block = 0; block < config.k; block++) {
                std::vector<int> comp_id(G.number_of_nodes(), -1);
                int num_comps = 0;
                std::vector<NodeWeight> comp_weight;
                std::vector<std::vector<NodeID>> comp_nodes;

                forall_nodes(G, node) {
                        if(G.getPartitionIndex(node) == block && comp_id[node] == -1) {
                                NodeWeight cw = 0;
                                std::vector<NodeID> nodes;
                                std::queue<NodeID> q;
                                q.push(node); comp_id[node] = num_comps;
                                while(!q.empty()) {
                                        NodeID v = q.front(); q.pop();
                                        cw += G.getNodeWeight(v);
                                        nodes.push_back(v);
                                        forall_out_edges(G, e, v) {
                                                NodeID u = G.getEdgeTarget(e);
                                                if(G.getPartitionIndex(u) == block && comp_id[u] == -1) {
                                                        comp_id[u] = num_comps; q.push(u);
                                                }
                                        } endfor
                                }
                                comp_weight.push_back(cw);
                                comp_nodes.push_back(nodes);
                                num_comps++;
                        }
                } endfor

                if(num_comps <= 1) continue;

                int largest = 0;
                for(int c = 1; c < num_comps; c++) {
                        if(comp_weight[c] > comp_weight[largest]) largest = c;
                }

                for(int c = 0; c < num_comps; c++) {
                        if(c == largest) continue;
                        std::vector<NodeWeight> bw(config.k, 0);
                        forall_nodes(G, n) { bw[G.getPartitionIndex(n)] += G.getNodeWeight(n); } endfor

                        PartitionID target = block;
                        NodeWeight min_w = std::numeric_limits<NodeWeight>::max();
                        for(NodeID n : comp_nodes[c]) {
                                forall_out_edges(G, e, n) {
                                        PartitionID ab = G.getPartitionIndex(G.getEdgeTarget(e));
                                        if(ab != block && bw[ab] < min_w) { min_w = bw[ab]; target = ab; }
                                } endfor
                        }
                        if(target != block) {
                                for(NodeID n : comp_nodes[c]) G.setPartitionIndex(n, target);
                        }
                }
        }

        // Step 2: Greedy rebalance
        bool progress = true;
        while(progress) {
                progress = false;
                std::vector<NodeWeight> bw(config.k, 0);
                forall_nodes(G, n) { bw[G.getPartitionIndex(n)] += G.getNodeWeight(n); } endfor

                forall_nodes(G, node) {
                        PartitionID from = G.getPartitionIndex(node);
                        if(bw[from] <= config.upper_bound_partition) continue;

                        PartitionID target = from;
                        NodeWeight min_w = std::numeric_limits<NodeWeight>::max();
                        forall_out_edges(G, e, node) {
                                PartitionID ab = G.getPartitionIndex(G.getEdgeTarget(e));
                                if(ab != from && bw[ab] + G.getNodeWeight(node) <= config.upper_bound_partition && bw[ab] < min_w) {
                                        min_w = bw[ab]; target = ab;
                                }
                        } endfor

                        if(target == from) continue;
                        if(would_disconnect_block(G, node, from)) continue;

                        bw[from] -= G.getNodeWeight(node);
                        bw[target] += G.getNodeWeight(node);
                        G.setPartitionIndex(node, target);
                        progress = true;
                } endfor
        }
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

        // Disable connectivity guard during regular refinement (METIS approach)
        bool connected_blocks_saved = cfg.connected_blocks;
        cfg.connected_blocks = false;

        improvement += (int)refine->perform_refinement(cfg, *coarsest, *coarser_boundary);

        // Checkpoint 1: coarsest level - repair + guarded refinement
        if(connected_blocks_saved) {
                eliminate_and_rebalance(cfg, *coarsest);
                if(!config.label_propagation_refinement) {
                        delete coarser_boundary;
                        coarser_boundary = new complete_boundary(coarsest);
                        coarser_boundary->build();
                }
                cfg.connected_blocks = true;
                improvement += (int)refine->perform_refinement(cfg, *coarsest, *coarser_boundary);
                cfg.connected_blocks = false;
        }

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

                cfg.connected_blocks = false; // unconstrained FM
                improvement += (int)refine->perform_refinement(cfg, *G, *finer_boundary);
                ASSERT_TRUE(graph_partition_assertions::assert_graph_has_kway_partition(config, *G));

                // Checkpoint at midpoint and finest level
                unsigned int current_level = hierarchy_deepth - hierarchy.size();
                if(connected_blocks_saved && (current_level == hierarchy_deepth/2 || hierarchy.isEmpty())) {
                        eliminate_and_rebalance(cfg, *G);
                        if(!config.label_propagation_refinement) {
                                delete finer_boundary;
                                finer_boundary = new complete_boundary(G);
                                finer_boundary->build();
                        }
                        cfg.connected_blocks = true; // guarded refinement after repair
                        improvement += (int)refine->perform_refinement(cfg, *G, *finer_boundary);
                        cfg.connected_blocks = false;
                }

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

        PRINT(std::cout <<  "log> starting uncoarsening ---------------"  << std::endl;)
        PartitionConfig cfg     = config;
        graph_access * coarsest = hierarchy.get_coarsest();
        quality_metrics qm;
        PRINT(std::cout << "log>" << "unrolling graph with " << coarsest->number_of_nodes() << std::endl;)

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
                PRINT(std::cout << "log>" << "unrolling graph with " << G->number_of_nodes() << std::endl;)

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

        PRINT(std::cout <<  "log> starting uncoarsening ---------------"  << std::endl;)
        PartitionConfig cfg     = config;
        graph_access * coarsest = hierarchy.get_coarsest();
        PRINT(std::cout << "log>" << "unrolling graph with " << coarsest->number_of_nodes() << std::endl;)

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
                PRINT(std::cout << "log>" << "unrolling graph with " << G->number_of_nodes() << std::endl;)

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

