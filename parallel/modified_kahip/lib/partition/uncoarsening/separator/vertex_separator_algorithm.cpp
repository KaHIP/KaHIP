/******************************************************************************
 * vertex_separator_algorithm.cpp
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

#include <sstream>

#include "graph_io.h"
#include "uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/simple_quotient_graph_scheduler.h"
#include "vertex_separator_algorithm.h"
#include "vertex_separator_flow_solver.h"

vertex_separator_algorithm::vertex_separator_algorithm() {

}

vertex_separator_algorithm::~vertex_separator_algorithm() {

}

void vertex_separator_algorithm::compute_vertex_separator(const PartitionConfig & config, 
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
                                        //overall_separator.push_back(cur_bnd_node);
                                        allready_separator[cur_bnd_node] = true;
                                }
                        } endfor
                } else {
                        forall_boundary_nodes(rhs_b, cur_bnd_node) {
                                if(allready_separator.find(cur_bnd_node) == allready_separator.end()) {
                                        //overall_separator.push_back(cur_bnd_node);
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
