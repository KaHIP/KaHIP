/******************************************************************************
 * vertex_separator_algorithm.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <sstream>

#include "area_bfs.h"
#include "algorithms/push_relabel.h"
#include "graph_io.h"
#include "most_balanced_minimum_cuts/most_balanced_minimum_cuts.h"
#include "tools/random_functions.h"
#include "tools/graph_extractor.h"
#include "tools/quality_metrics.h"
#include "tools/graph_extractor.h"
#include "data_structure/union_find.h"
#include "uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/simple_quotient_graph_scheduler.h"
#include "vertex_separator_algorithm.h"
#include "vertex_separator_flow_solver.h"

vertex_separator_algorithm::vertex_separator_algorithm() {

}

vertex_separator_algorithm::~vertex_separator_algorithm() {

}

void vertex_separator_algorithm::build_flow_problem(const PartitionConfig & config, 
                                                          graph_access & G, 
                                                          std::vector< NodeID > & lhs_nodes,
                                                          std::vector< NodeID > & rhs_nodes,
                                                          std::vector< NodeID > & separator_nodes,
                                                          flow_graph & rG, 
                                                          std::vector< NodeID > & forward_mapping, 
                                                          NodeID & source, NodeID & sink) {
        // *************** determine outer boundary of boundary stripes ***************
        std::vector<NodeID> outer_lhs_boundary_nodes;
        std::vector<NodeID> outer_rhs_boundary_nodes;

        for( NodeID node : lhs_nodes ) {
                G.setPartitionIndex(node, 3);
        }

        for( NodeID node : rhs_nodes ) {
                G.setPartitionIndex(node, 3);
        }

        for( unsigned i = 0; i < separator_nodes.size(); i++) {
                G.setPartitionIndex(separator_nodes[i], 3);
        }

        for( NodeID node : lhs_nodes ) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(G.getPartitionIndex(target) == 0 || G.getPartitionIndex(target) == 1) {
                                // outer boundary node
                                outer_lhs_boundary_nodes.push_back(node);
                                break;
                        }
                } endfor
        }
        
        if(lhs_nodes.size() > 0) {
                for( NodeID node : separator_nodes) {
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if(G.getPartitionIndex(target) == 0) {
                                        // outer boundary node
                                        outer_lhs_boundary_nodes.push_back(node);
                                        break;
                                }
                        } endfor
                }
        } else {
                outer_lhs_boundary_nodes = separator_nodes;
        }

        for( NodeID node : rhs_nodes ) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(G.getPartitionIndex(target) == 0 || G.getPartitionIndex(target) == 1) {
                                // outer boundary node
                                outer_rhs_boundary_nodes.push_back(node);
                                break;
                        }
                } endfor
        }

        if(rhs_nodes.size() > 0) {
                for( NodeID node : separator_nodes) {
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if(G.getPartitionIndex(target) == 1) {
                                        // outer boundary node
                                        outer_rhs_boundary_nodes.push_back(node);
                                        break;
                                }
                        } endfor
                }
        } else {
                outer_rhs_boundary_nodes = separator_nodes;
        }

        NodeID n = 2*(lhs_nodes.size() + rhs_nodes.size() + separator_nodes.size()) + 2; // source and sink

        // find forward and backward mapping
        std::unordered_map< NodeID, NodeID > backward_mapping;
        forward_mapping.clear();
        forward_mapping.resize(n-2);

        NodeID node_count = 0;
        for( NodeID v : lhs_nodes ) {
                backward_mapping[v] = node_count;
                forward_mapping[node_count++] = v;
                forward_mapping[node_count++] = v;
        }

        for( NodeID v : rhs_nodes ) {
                backward_mapping[v] = node_count;
                forward_mapping[node_count++] = v;
                forward_mapping[node_count++] = v;
        }

        for( NodeID v : separator_nodes) {
                backward_mapping[v] = node_count;
                forward_mapping[node_count++] = v;
                forward_mapping[node_count++] = v;
        }


        source = n-2;
        sink   = n-1;
        FlowType infinite = std::numeric_limits<FlowType>::max()/2;
        rG.start_construction(n, 0);

        for( NodeID v : outer_lhs_boundary_nodes) {
                rG.new_edge(source, backward_mapping[v], infinite);
        }

        for( NodeID v : outer_rhs_boundary_nodes) {
                rG.new_edge(backward_mapping[v]+1, sink, infinite);
        }

        for( NodeID v : lhs_nodes) {
                rG.new_edge(backward_mapping[v], backward_mapping[v]+1, G.getNodeWeight(v));
                forall_out_edges(G, e, v) {
                        NodeID target = G.getEdgeTarget(e);
                        if(G.getPartitionIndex(target) != 3) continue; // not part of the flow problem
                        rG.new_edge( backward_mapping[target]+1, backward_mapping[v], infinite);
                } endfor
        }

        for( NodeID v : rhs_nodes) {
                rG.new_edge(backward_mapping[v], backward_mapping[v]+1, G.getNodeWeight(v));
                forall_out_edges(G, e, v) {
                        NodeID target = G.getEdgeTarget(e);
                        if(G.getPartitionIndex(target) != 3) continue; // not part of the flow problem
                        rG.new_edge( backward_mapping[target]+1, backward_mapping[v], infinite);
                } endfor
        }

        for( NodeID v : separator_nodes) {
                rG.new_edge(backward_mapping[v], backward_mapping[v]+1, G.getNodeWeight(v));
                forall_out_edges(G, e, v) {
                        NodeID target = G.getEdgeTarget(e);
                        if(G.getPartitionIndex(target) != 3) continue; // not part of the flow problem
                        rG.new_edge( backward_mapping[target]+1, backward_mapping[v], infinite);
                } endfor
        }

        rG.finish_construction();
}

NodeWeight vertex_separator_algorithm::improve_vertex_separator(const PartitionConfig & config, 
							        graph_access & G, 
								std::vector< NodeWeight > & block_weights,
						 	        PartialBoundary & separator) {

        G.set_partition_count(3);
        
        double current_region_factor    = config.region_factor_node_separators;
        double iteration                = 0;
        NodeWeight prev_solution_value  = block_weights[2];

        NodeWeight solution_value;
        bool solution_imbalanced;
        std::vector< NodeWeight > old_block_weights = block_weights;

        do {
                PartitionConfig cfg               = config;
                cfg.region_factor_node_separators = 1+current_region_factor;
                solution_imbalanced               = false;
	
                std::vector< NodeID > old_lhs;
                std::vector< NodeID > old_rhs;
                std::vector< NodeID > old_sep;

                solution_value = improve_vertex_separator_internal( cfg , G, block_weights, separator, old_lhs, old_rhs, old_sep);

                if( solution_value == prev_solution_value ) {
                        NodeWeight cur_diff = abs( (int)(block_weights[1]) - (int)(block_weights[0]) );
                        NodeWeight old_diff = abs( (int)(old_block_weights[1]) - (int)(old_block_weights[0]) );
                        if( cur_diff > old_diff ) {
                                //reconstruct old solution and return since solution value will not be better
                                apply_vectors( G, old_lhs, old_rhs, old_sep);
                                separator.clear();
                                for( NodeID node : old_sep ) {
                                        separator.insert(node);
                                }
                                block_weights = old_block_weights;
                        }
                        return 0;
                }
	
                if( block_weights[0] > config.upper_bound_partition || block_weights[1] > config.upper_bound_partition ) {
                        solution_imbalanced = true;
                        current_region_factor /= 2;

                        //reconstruct old solution
                        apply_vectors( G, old_lhs, old_rhs, old_sep);
			separator.clear();
                        for( NodeID node : old_sep ) {
                                separator.insert(node);
                        }
			block_weights = old_block_weights;
                }
                iteration++;
        } while ( solution_imbalanced && iteration < 10);

        if( solution_imbalanced ) {
                PartitionConfig cfg = config;
                cfg.region_factor_node_separators = 1;

                std::vector< NodeID > old_lhs;
                std::vector< NodeID > old_rhs;
                std::vector< NodeID > old_sep;
                solution_value = improve_vertex_separator_internal(cfg , G, block_weights, separator, old_lhs, old_rhs, old_sep);
        }

        return prev_solution_value-solution_value;
}

NodeWeight vertex_separator_algorithm::improve_vertex_separator_internal(const PartitionConfig & config, 
	                                                                 graph_access & G, 
									 std::vector< NodeWeight > & block_weights,
			 						 PartialBoundary & separator,
                                                                         std::vector< NodeID > & lhs_nodes, 
                                                                         std::vector< NodeID > & rhs_nodes,
                                                                         std::vector< NodeID > & start_nodes) {
        NodeWeight lhs_part_weight  = block_weights[0];
        NodeWeight rhs_part_weight  = block_weights[1];
        NodeWeight separator_weight = block_weights[2];

        area_bfs abfs;
        // perform BFS into one side
	forall_boundary_nodes( separator, node ) {
		start_nodes.push_back(node);
	} endfor
        abfs.perform_bfs(config, G, start_nodes, 0, block_weights, lhs_nodes);

        // perform BFS into other side
        abfs.perform_bfs(config, G, start_nodes, 1, block_weights, rhs_nodes);

        // now build the flow problem
        flow_graph rG; NodeID source, sink;
        std::vector< NodeID > forward_mapping; // maps a node from rG to original G
        build_flow_problem(config, G, lhs_nodes, rhs_nodes, start_nodes, rG, forward_mapping, source, sink);

	push_relabel mfmc_solver; std::vector<NodeID> source_set;
        bool compute_source_set = !config.most_balanced_minimum_cuts_node_sep;
	FlowType value =  mfmc_solver.solve_max_flow_min_cut(rG, source, sink, compute_source_set, source_set);

        std::vector< bool > is_in_source_set( rG.number_of_nodes());
        bool start_value = config.most_balanced_minimum_cuts_node_sep;
        forall_nodes(rG, node) {
                is_in_source_set[node] = start_value;
        } endfor

        // most balanced minimum cuts
        if(!config.most_balanced_minimum_cuts_node_sep) {
                for( NodeID v : source_set ) {
                        is_in_source_set[v] = true;
                }
        } else {
                graph_access residualGraph; 
                convert_residualGraph( G, forward_mapping, source, sink, rG, residualGraph );

                NodeWeight rhs_stripe_weight = 0;
                for( NodeID v : rhs_nodes ) {
                        rhs_stripe_weight += G.getNodeWeight(v);
                }

                NodeWeight overall_weight = lhs_part_weight + separator_weight + rhs_part_weight;
                NodeWeight ideal_new_block_weight = (overall_weight - value)/2;
                NodeWeight amount_to_be_added = ideal_new_block_weight - rhs_part_weight;
                NodeWeight perfect_rhs_stripe_weight = std::max((NodeWeight)(2*amount_to_be_added+value),(NodeWeight)0);
                 
                PartitionConfig tmpconfig = config;
                tmpconfig.mode_node_separators = true;

                most_balanced_minimum_cuts mbmc;
                std::vector<NodeID> rhs_set;
                mbmc.compute_good_balanced_min_cut(residualGraph, tmpconfig, perfect_rhs_stripe_weight, rhs_set);

                for( NodeID v : rhs_set) {
                        is_in_source_set[v] = false ;
                }
        }
        

        //reconstruct old pidx for weight computations
        apply_vectors( G, lhs_nodes, rhs_nodes, start_nodes );

        //reconstruct the partition IDs
        forall_nodes(rG, node) {
                if(node == sink || node == source) continue;
		NodeID v              = forward_mapping[node];
		PartitionID to_set    = is_in_source_set[node] ? 0 : 1;
		bool is_frontier_node = is_in_source_set[node] && !is_in_source_set[node+1];
		if(!is_frontier_node) {
			block_weights[G.getPartitionIndex(v)] -= G.getNodeWeight(v);
			G.setPartitionIndex(v, to_set);
			block_weights[G.getPartitionIndex(v)] += G.getNodeWeight(v);
		} else {
			block_weights[G.getPartitionIndex(v)] -= G.getNodeWeight(v);
		}
		node++;
        } endfor

        separator.clear();
	block_weights[2] = value;
        forall_nodes(rG, node) {
                if(node == sink || node == source) continue;
                if( is_in_source_set[node] && !is_in_source_set[node+1]) {
			NodeID v = forward_mapping[node];
                        separator.insert(v);
                        G.setPartitionIndex(v, 2);
                }
                node++;
        } endfor

        return value; // return flow value 
}


NodeWeight vertex_separator_algorithm::improve_vertex_separator(const PartitionConfig & config, 
                                              graph_access & G, 
                                              std::vector<NodeID> & input_separator,
                                              std::vector<NodeID> & output_separator) {

        std::vector< NodeID > current_solution(G.number_of_nodes(), 0);
        forall_nodes(G, node) {
                current_solution[node] = G.getPartitionIndex(node);
        } endfor
        
        bool solution_imbalanced;
        NodeWeight solution_value;
        double current_region_factor = config.region_factor_node_separators;
        double iteration = 0;
        quality_metrics qm;
        do {
                PartitionConfig cfg = config;
                cfg.region_factor_node_separators = 1+current_region_factor;
                solution_imbalanced = false;
                solution_value = improve_vertex_separator_internal( cfg , G, input_separator, output_separator);
                G.set_partition_count(3);
                double balance = qm.balance_separator(G);
                if( balance > (1+config.epsilon/(double)100) ) {
                        solution_imbalanced = true;
                        current_region_factor /= 2;

                        forall_nodes(G, node) {
                                G.setPartitionIndex(node, current_solution[node]);
                        } endfor
                }
                iteration++;
        } while ( solution_imbalanced && iteration < 10);

        if( solution_imbalanced ) {
                PartitionConfig cfg = config;
                cfg.region_factor_node_separators = 1;
                solution_value = improve_vertex_separator_internal( cfg , G, input_separator, output_separator);
        }

        return solution_value;
}

NodeWeight vertex_separator_algorithm::improve_vertex_separator_internal(const PartitionConfig & config, 
                                              graph_access & G, 
                                              std::vector<NodeID> & input_separator,
                                              std::vector<NodeID> & output_separator) {

        NodeWeight lhs_part_weight = 0;
        NodeWeight rhs_part_weight = 0;
        NodeWeight separator_weight = 0;
        forall_nodes(G, node) {
                if( G.getPartitionIndex(node) == 1 ) {
                       rhs_part_weight += G.getNodeWeight(node);
                } else if ( G.getPartitionIndex(node) == 0) {
                       lhs_part_weight += G.getNodeWeight(node); 
                } else if (G.getPartitionIndex(node) == 2) {
                       separator_weight += G.getNodeWeight(node); 
                } 
        } endfor

        NodeWeight old_separator_weight = separator_weight;

        area_bfs abfs;
        std::vector< NodeWeight > block_weights(3, 0);
        block_weights[0] = lhs_part_weight;
        block_weights[1] = rhs_part_weight;
        block_weights[2] = old_separator_weight;
        // perform BFS into one side
        std::vector<NodeID> lhs_nodes;
        abfs.perform_bfs(config, G, input_separator, 0, block_weights, lhs_nodes);

        // perform BFS into other side
        std::vector<NodeID> rhs_nodes;
        abfs.perform_bfs(config, G, input_separator, 1, block_weights, rhs_nodes);

        // now build the flow problem
        flow_graph rG; NodeID source, sink;
        std::vector< NodeID > forward_mapping; // maps a node from rG to original G
        build_flow_problem(config, G, lhs_nodes, rhs_nodes, input_separator, rG, forward_mapping, source, sink);

	push_relabel mfmc_solver; std::vector<NodeID> source_set;
        bool compute_source_set = !config.most_balanced_minimum_cuts_node_sep;
	FlowType value =  mfmc_solver.solve_max_flow_min_cut(rG, source, sink, compute_source_set, source_set);

        std::vector< bool > is_in_source_set( rG.number_of_nodes());
        bool start_value = config.most_balanced_minimum_cuts_node_sep;
        forall_nodes(rG, node) {
                is_in_source_set[node] = start_value;
        } endfor

        // most balanced minimum cuts
        if(!config.most_balanced_minimum_cuts_node_sep) {
                for( NodeID v : source_set ) {
                        is_in_source_set[v] = true;
                }
        } else {
                graph_access residualGraph; 
                convert_residualGraph( G, forward_mapping, source, sink, rG, residualGraph );

                NodeWeight rhs_stripe_weight = 0;
                for( NodeID v : rhs_nodes ) {
                        rhs_stripe_weight += G.getNodeWeight(v);
                }

                NodeWeight overall_weight = lhs_part_weight + separator_weight + rhs_part_weight;
                NodeWeight ideal_new_block_weight = (overall_weight - value)/2;
                NodeWeight amount_to_be_added = ideal_new_block_weight - rhs_part_weight;
                NodeWeight perfect_rhs_stripe_weight = std::max((NodeWeight)(2*amount_to_be_added+value),(NodeWeight)0);
                 
                PartitionConfig tmpconfig = config;
                tmpconfig.mode_node_separators = true;

                most_balanced_minimum_cuts mbmc;
                std::vector<NodeID> rhs_set;
                mbmc.compute_good_balanced_min_cut(residualGraph, tmpconfig, perfect_rhs_stripe_weight, rhs_set);

                for( NodeID v : rhs_set) {
                        is_in_source_set[v] = false ;
                }
        }
        
        //reconstruct the partition IDs
        forall_nodes(rG, node) {
                if(node == sink || node == source) continue;
                if( is_in_source_set[node]) {
                        G.setPartitionIndex(forward_mapping[node],0);
                } else {
                        G.setPartitionIndex(forward_mapping[node],1);
                }
        } endfor

        output_separator.clear();
        forall_nodes(rG, node) {
                if(node == sink || node == source) continue;
                if( is_in_source_set[node] && !is_in_source_set[node+1]) {
                        output_separator.push_back(forward_mapping[node]);
                        G.setPartitionIndex(forward_mapping[node],2);
                }
                node++;
        } endfor


        std::unordered_map<NodeID, bool> allready_separator;
        for( unsigned int i = 0; i < output_separator.size(); i++) {
                allready_separator[output_separator[i]] = true;
        }

        //TODO remove them later on
        //is_vertex_separator(G, allready_separator);

        return old_separator_weight - value; // return improvement
}

void vertex_separator_algorithm::compute_vertex_separator(const PartitionConfig & config, 
                                                          graph_access & G, 
                                                          complete_boundary & boundary, 
                                                          std::vector<NodeID> & overall_separator) {

        PartitionConfig cfg     = config;
        cfg.bank_account_factor = 1;

        std::unordered_map<NodeID, bool> allready_separator;

        QuotientGraphEdges qgraph_edges;
        boundary.getQuotientGraphEdges(qgraph_edges);

        if(qgraph_edges.size() == 0) {
                //is_vertex_separator(G, allready_separator);         
                return;
        }

        quotient_graph_scheduling* scheduler = new simple_quotient_graph_scheduler(cfg, qgraph_edges,qgraph_edges.size()); 

        do {
                boundary_pair & bp = scheduler->getNext();
                PartitionID lhs = bp.lhs;
                PartitionID rhs = bp.rhs;

                boundary_starting_nodes start_nodes_lhs;
                boundary_starting_nodes start_nodes_rhs;

                PartialBoundary & lhs_b = boundary.getDirectedBoundary(lhs, lhs, rhs);
                PartialBoundary & rhs_b = boundary.getDirectedBoundary(rhs, lhs, rhs);

                forall_boundary_nodes(lhs_b, cur_bnd_node) {
                        if(allready_separator.find(cur_bnd_node) == allready_separator.end()) {
                                start_nodes_lhs.push_back(cur_bnd_node);
                        }
                } endfor

                forall_boundary_nodes(rhs_b, cur_bnd_node) {
                        if(allready_separator.find(cur_bnd_node) == allready_separator.end()) {
                                start_nodes_rhs.push_back(cur_bnd_node);
                        }
                } endfor

                vertex_separator_flow_solver vsfs;
                std::vector<NodeID> separator;
                vsfs.find_separator(config, G, lhs, rhs, start_nodes_lhs, start_nodes_rhs, separator);
                for( unsigned i = 0; i < separator.size(); i++) {
                        allready_separator[separator[i]] = true;
                }
                //*************************** end **************************************** 
        } while(!scheduler->hasFinished());


        // now print the computed vertex separator to disk
        std::unordered_map<NodeID, bool>::iterator it;
        for( it = allready_separator.begin(); it != allready_separator.end(); ++it) {
                overall_separator.push_back(it->first);
                G.setPartitionIndex(it->first, G.getSeparatorBlock());
        }
        delete scheduler;
        //std::cout <<  "performing check "  << std::endl;
        //is_vertex_separator(G, allready_separator);         
        //std::cout <<  "performing check done "  << std::endl;
}

void vertex_separator_algorithm::compute_vertex_separator(const PartitionConfig & config, 
                                                          graph_access & G, 
                                                          complete_boundary & boundary) {

        std::vector<NodeID> overall_separator;
        compute_vertex_separator(config, G, boundary, overall_separator); 

        // write the partition to the disc 
        std::stringstream filename;
        filename << "tmpseparator" << config.k;
        graph_io::writeVector(overall_separator, filename.str());
}

void vertex_separator_algorithm::compute_vertex_separator_simpler(const PartitionConfig & config, 
                                                          graph_access & G, 
                                                          complete_boundary & boundary, 
                                                          std::vector<NodeID> & overall_separator) {

        PartitionConfig cfg     = config;
        cfg.bank_account_factor = 1;

        QuotientGraphEdges qgraph_edges;
        boundary.getQuotientGraphEdges(qgraph_edges);

        if(qgraph_edges.size() == 0) {
                //is_vertex_separator(G, allready_separator);         
                return;
        }

        quotient_graph_scheduling* scheduler = new simple_quotient_graph_scheduler(cfg, qgraph_edges,qgraph_edges.size()); 

        std::unordered_map<NodeID, bool> allready_separator;
        do {
                boundary_pair & bp = scheduler->getNext();
                PartitionID lhs = bp.lhs;
                PartitionID rhs = bp.rhs;

                boundary_starting_nodes start_nodes_lhs;
                boundary_starting_nodes start_nodes_rhs;

                PartialBoundary & lhs_b = boundary.getDirectedBoundary(lhs, lhs, rhs);
                PartialBoundary & rhs_b = boundary.getDirectedBoundary(rhs, lhs, rhs);

                        forall_boundary_nodes(lhs_b, cur_bnd_node) {
                                if(allready_separator.find(cur_bnd_node) == allready_separator.end()) {
                                        allready_separator[cur_bnd_node] = true;
                                }
                        } endfor
                        forall_boundary_nodes(rhs_b, cur_bnd_node) {
                                if(allready_separator.find(cur_bnd_node) == allready_separator.end()) {
                                        allready_separator[cur_bnd_node] = true;
                                }
                        } endfor

                //*************************** end **************************************** 
        } while(!scheduler->hasFinished());

        // now print the computed vertex separator to disk
        std::unordered_map<NodeID, bool>::iterator it;
        for( it = allready_separator.begin(); it != allready_separator.end(); ++it) {
                overall_separator.push_back(it->first);
                G.setPartitionIndex(it->first, G.getSeparatorBlock());
        }
        //is_vertex_separator(G, allready_separator);         
        delete scheduler;
}

void vertex_separator_algorithm::compute_vertex_separator_simple(const PartitionConfig & config, 
                                                          graph_access & G, 
                                                          complete_boundary & boundary, 
                                                          std::vector<NodeID> & overall_separator) {

        PartitionConfig cfg     = config;
        cfg.bank_account_factor = 1;

        QuotientGraphEdges qgraph_edges;
        boundary.getQuotientGraphEdges(qgraph_edges);

        if(qgraph_edges.size() == 0) {
                //is_vertex_separator(G, allready_separator);         
                return;
        }

        quotient_graph_scheduling* scheduler = new simple_quotient_graph_scheduler(cfg, qgraph_edges,qgraph_edges.size()); 

        std::unordered_map<NodeID, bool> allready_separator;
        do {
                boundary_pair & bp = scheduler->getNext();
                PartitionID lhs = bp.lhs;
                PartitionID rhs = bp.rhs;

                boundary_starting_nodes start_nodes_lhs;
                boundary_starting_nodes start_nodes_rhs;

                PartialBoundary & lhs_b = boundary.getDirectedBoundary(lhs, lhs, rhs);
                PartialBoundary & rhs_b = boundary.getDirectedBoundary(rhs, lhs, rhs);

                if(lhs_b.size() < rhs_b.size()) {
                        forall_boundary_nodes(lhs_b, cur_bnd_node) {
                                if(allready_separator.find(cur_bnd_node) == allready_separator.end()) {
                                        allready_separator[cur_bnd_node] = true;
                                }
                        } endfor
                } else {
                        forall_boundary_nodes(rhs_b, cur_bnd_node) {
                                if(allready_separator.find(cur_bnd_node) == allready_separator.end()) {
                                        allready_separator[cur_bnd_node] = true;
                                }
                        } endfor
                }

                //*************************** end **************************************** 
        } while(!scheduler->hasFinished());

        // now print the computed vertex separator to disk
        std::unordered_map<NodeID, bool>::iterator it;
        for( it = allready_separator.begin(); it != allready_separator.end(); ++it) {
                overall_separator.push_back(it->first);
                G.setPartitionIndex(it->first, G.getSeparatorBlock());
        }
        //is_vertex_separator(G, allready_separator);         
        delete scheduler;
}



bool vertex_separator_algorithm::is_vertex_separator(graph_access & G, std::unordered_map<NodeID, bool> & separator) {
         forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( G.getPartitionIndex(node) != G.getPartitionIndex(target)) {

                                // in this case one of them has to be a separator
                                if( separator.find(node)   == separator.end() && 
                                    separator.find(target) == separator.end()) {
                                        std::cout <<  "not a separator! " <<  node <<  " " <<  target << std::endl;
                                        std::cout <<  "PartitionIndex node " << G.getPartitionIndex(node)  << std::endl;
                                        std::cout <<  "PartitionIndex target " << G.getPartitionIndex(target)  << std::endl;
                                        ASSERT_TRUE(false);
                                        exit(0);
                                } 
                        } 
                } endfor
         } endfor
         return true;
}

