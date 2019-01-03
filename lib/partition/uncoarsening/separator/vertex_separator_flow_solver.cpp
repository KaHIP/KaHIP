/******************************************************************************
 * vertex_separator_flow_solver.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>
#include <map>
#include <math.h>
#include <unordered_map>

#include "algorithms/push_relabel.h"
#include "data_structure/flow_graph.h"
#include "vertex_separator_flow_solver.h"

vertex_separator_flow_solver::vertex_separator_flow_solver() {

}

vertex_separator_flow_solver::~vertex_separator_flow_solver() {

}

void vertex_separator_flow_solver::find_separator(const PartitionConfig & config, 
                                                  graph_access & G, 
                                                  PartitionID lhs, 
                                                  PartitionID rhs,  
                                                  boundary_starting_nodes lhs_nodes,
                                                  boundary_starting_nodes rhs_nodes,
                                                  std::vector<NodeID> & separator) {
        if(lhs_nodes.size() == 0 || rhs_nodes.size() == 0) return;

        std::vector<NodeID> new_to_old_ids; flow_graph fG;
        build_flow_pb(config, G, lhs, rhs, lhs_nodes, rhs_nodes, new_to_old_ids, fG);

        push_relabel pr;
        NodeID source = fG.number_of_nodes() - 2;
        NodeID sink   = fG.number_of_nodes() - 1;

        std::vector<NodeID> S_tmp;
        pr.solve_max_flow_min_cut( fG, source, sink, true, S_tmp);

        std::vector<NodeID> S;
        for( unsigned i = 0; i < S_tmp.size(); i++) {
                if( S_tmp[i] != source && S_tmp[i] != sink) 
                        S.push_back(new_to_old_ids[S_tmp[i]]);
        }
        std::sort(lhs_nodes.begin(), lhs_nodes.end());
        std::sort(rhs_nodes.begin(), rhs_nodes.end());
        std::sort(S.begin(), S.end());
 
        std::vector<int> separator_tmp(lhs_nodes.size() + rhs_nodes.size(), -1);
        std::vector<int>::iterator it;
        it = std::set_intersection(rhs_nodes.begin(), rhs_nodes.end(), S.begin(), S.end(), separator_tmp.begin()); 

        for( unsigned i = 0; i < separator_tmp.size(); i++) {
                if(separator_tmp[i] != -1) {
                        separator.push_back(separator_tmp[i]);                
                }
        }

        std::vector<int>::iterator it2;
        std::vector<int> separator_tmp2(lhs_nodes.size() + rhs_nodes.size(), -1);
        it2 = std::set_difference(lhs_nodes.begin(), lhs_nodes.end(), S.begin(), S.end(), separator_tmp2.begin()); 
        for( unsigned i = 0; i < separator_tmp2.size(); i++) {
                if(separator_tmp2[i] != -1) {
                        separator.push_back(separator_tmp2[i]);                
                }
        }

}

bool vertex_separator_flow_solver::build_flow_pb( const PartitionConfig & config, 
                               graph_access & G, 
                               PartitionID & lhs, 
                               PartitionID & rhs, 
                               std::vector<NodeID> & lhs_nodes,
                               std::vector<NodeID> & rhs_nodes,
                               std::vector<NodeID> & new_to_old_ids,              
                               flow_graph & fG) {

        unsigned no_edges = 0;
        for( unsigned i = 0; i < lhs_nodes.size(); i++) {
                NodeID node = lhs_nodes[i];
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(rhs == G.getPartitionIndex(target)) {
                                ++no_edges;
                        }
                } endfor
        }

        //build mappings from old to new node ids and reverse
        NodeID idx = 0;
        new_to_old_ids.resize(lhs_nodes.size() + rhs_nodes.size());
        std::unordered_map<NodeID, NodeID> old_to_new;
        for( unsigned i = 0; i < lhs_nodes.size(); i++) {
                new_to_old_ids[idx] = lhs_nodes[i];
                old_to_new[lhs_nodes[i]] = idx++ ;
        }
        for( unsigned i = 0; i < rhs_nodes.size(); i++) {
                new_to_old_ids[idx] = rhs_nodes[i];
                old_to_new[rhs_nodes[i]] = idx++;
        }

        NodeID n = lhs_nodes.size() + rhs_nodes.size() + 2; //+source and target
        if(n == 2) return false;

        NodeID source = n - 2;
        NodeID sink   = n - 1;

        idx = 0;
        FlowType max_capacity = std::numeric_limits<FlowType>::max();

        fG.start_construction(n);
        //insert directed edges from L to R
        for( unsigned i = 0; i < lhs_nodes.size(); i++, idx++) {
                NodeID node = lhs_nodes[i];
                NodeID sourceID = idx;
                forall_out_edges(G, e, node) {
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == rhs)  {
                                NodeID targetID = old_to_new[G.getEdgeTarget(e)];
                                fG.new_edge(sourceID, targetID, max_capacity);
                        }
                } endfor
        }

        //connect source and target with outer boundary nodes 
        for(unsigned i = 0; i < lhs_nodes.size(); i++) {
                NodeID targetID = old_to_new[lhs_nodes[i]];
                fG.new_edge(source, targetID, G.getNodeWeight(lhs_nodes[i]));
        }

        for(unsigned i = 0; i < rhs_nodes.size(); i++) {
                NodeID sourceID = old_to_new[rhs_nodes[i]];
                fG.new_edge(sourceID, sink, G.getNodeWeight(rhs_nodes[i]));
        }

        return true;
}

