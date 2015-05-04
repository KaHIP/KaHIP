/******************************************************************************
 * vertex_separator_algorithm.cpp 
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

#include <sstream>

#include "algorithms/push_relabel.h"
#include "graph_io.h"
#include "tools/random_functions.h"
#include "uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/simple_quotient_graph_scheduler.h"
#include "vertex_separator_algorithm.h"
#include "vertex_separator_flow_solver.h"

vertex_separator_algorithm::vertex_separator_algorithm() {

}

vertex_separator_algorithm::~vertex_separator_algorithm() {

}

void vertex_separator_algorithm::region_bfs(const PartitionConfig & config, 
                graph_access & G, 
                std::vector< NodeID > & input_separator, 
                PartitionID block, 
                std::vector< NodeID > & reached_nodes) {
        std::queue<NodeID> node_queue;
	std::vector<int> deepth(G.number_of_nodes(), -1);
	int cur_deepth = 0;
       
        random_functions::permutate_vector_good(input_separator, false);
	/***************************
	 * Initialize the Queue
	 * *************************/
	for(unsigned int i = 0; i < input_separator.size(); i++) {
		node_queue.push(input_separator[i]);
		deepth[input_separator[i]] = cur_deepth;
	}
	++cur_deepth;

        NodeWeight size_lhs = 0;
        NodeWeight size_rhs = 0;
        forall_nodes(G, node) {
                if(G.getPartitionIndex(node) == 0) {
                        size_lhs += G.getNodeWeight(node);
                } else if( G.getPartitionIndex(node) == 1) {
                        size_rhs += G.getNodeWeight(node);
                }
        } endfor
        
        NodeWeight accumulated_weight = 0;
        NodeID upper_bound_no_nodes;

        if( block == 0 ) {
                upper_bound_no_nodes = config.upper_bound_partition - size_rhs;
        } else {
                upper_bound_no_nodes = config.upper_bound_partition - size_lhs;
        }

        //upper_bound_no_nodes *= 2;

        /***************************
	 * Do the BFS
	 ***************************/
	while (!node_queue.empty()) {
		if(accumulated_weight >= upper_bound_no_nodes) break;
		NodeID n = node_queue.front();
		node_queue.pop();

		if (deepth[n] == cur_deepth) {
			cur_deepth++;
		}
		forall_out_edges(G,e,n) {
			NodeID target = G.getEdgeTarget(e);
			if(deepth[target] == -1 && G.getPartitionIndex(target) == block 
                           && accumulated_weight + G.getNodeWeight(target) <= upper_bound_no_nodes) {
				deepth[target] = cur_deepth;
				node_queue.push(target);
				reached_nodes.push_back(target);
                                accumulated_weight += G.getNodeWeight(target);
			}
		} endfor
	}
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
                G.setPartitionIndex(separator_nodes[i],3);
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

void vertex_separator_algorithm::improve_vertex_separator(const PartitionConfig & config, 
                                              graph_access & G, 
                                              std::vector<NodeID> & input_separator,
                                              std::vector<NodeID> & output_separator) {
        // perform BFS into one side
        std::vector<NodeID> lhs_nodes;
        region_bfs(config, G, input_separator, 0, lhs_nodes);

        // perform BFS into other side
        std::vector<NodeID> rhs_nodes;
        region_bfs(config, G, input_separator, 1, rhs_nodes);

        std::cout <<  "lhs nodes size " <<  lhs_nodes.size()  << std::endl;
        std::cout <<  "rhs nodes size " <<  rhs_nodes.size()  << std::endl;

        flow_graph rG;
        NodeID source, sink;
        // now build the flow problem
        std::vector< NodeID > forward_mapping;
        build_flow_problem(config, G, lhs_nodes, rhs_nodes, input_separator, rG, forward_mapping, source, sink);


        std::vector<NodeID> source_set;
	push_relabel mfmc_solver;
	FlowType value =  mfmc_solver.solve_max_flow_min_cut(rG, source, sink, true, source_set);

        std::cout <<  "source set size " <<  source_set.size()  << std::endl;
        std::cout <<  "flow value " <<  value  << std::endl;
        
        std::vector< bool > is_in_source_set( rG.number_of_nodes(), false);
        std::vector< bool > is_in_separator( G.number_of_nodes(), false);
        for( NodeID v : source_set) {
                is_in_source_set[v] = true;
        }

        output_separator.clear();
        forall_nodes(rG, node) {
                if(node == sink || node == source) continue;
                if( is_in_source_set[node] && !is_in_source_set[node+1]) {
                        output_separator.push_back(forward_mapping[node]);
                }
                node++;
        } endfor

        std::cout <<  "improvement achieved " <<  (input_separator.size()-value)  << std::endl;
        std::cout <<  "relative improvement achieved " <<  (input_separator.size()/(double)value)  << std::endl;
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
                is_vertex_separator(G, allready_separator);         
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
        }
        is_vertex_separator(G, allready_separator);         
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

void vertex_separator_algorithm::compute_vertex_separator_simple(const PartitionConfig & config, 
                                                          graph_access & G, 
                                                          complete_boundary & boundary, 
                                                          std::vector<NodeID> & overall_separator) {

        PartitionConfig cfg     = config;
        cfg.bank_account_factor = 1;

        QuotientGraphEdges qgraph_edges;
        boundary.getQuotientGraphEdges(qgraph_edges);

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
        }
        is_vertex_separator(G, allready_separator);         
}

bool vertex_separator_algorithm::is_vertex_separator(graph_access & G, std::unordered_map<NodeID, bool> & separator) {
         forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( G.getPartitionIndex(node) != G.getPartitionIndex(target)) {
                                // in this case one of them has to be a separator
                                if( separator.find(node)   == separator.end() && 
                                    separator.find(target) == separator.end()) {
                                        std::cout <<  "not a separator!"  << std::endl;
                                        ASSERT_TRUE(false);
                                } 
                        } 
                } endfor
         } endfor
         return true;
}
