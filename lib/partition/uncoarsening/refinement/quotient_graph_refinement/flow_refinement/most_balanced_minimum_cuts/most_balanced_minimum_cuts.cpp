/******************************************************************************
 * most_balanced_minimum_cuts.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <unordered_map>

#include "algorithms/strongly_connected_components.h"
#include "algorithms/topological_sort.h"
#include "most_balanced_minimum_cuts.h"

most_balanced_minimum_cuts::most_balanced_minimum_cuts() {

}

most_balanced_minimum_cuts::~most_balanced_minimum_cuts() {

}

void most_balanced_minimum_cuts::compute_good_balanced_min_cut( graph_access & residualGraph, 
                                                                const PartitionConfig & config,
                                                                NodeWeight & perfect_rhs_weight, 
                                                                std::vector< NodeID > & new_rhs_nodes ) {

        strongly_connected_components scc;
        std::vector<int> components(residualGraph.number_of_nodes());
        int comp_count = scc.strong_components(residualGraph,components);
        
        std::vector< std::vector<NodeID> > comp_nodes(comp_count);
        std::vector< NodeWeight > comp_weights(comp_count, 0);
        
        forall_nodes(residualGraph, node) {
               comp_nodes[components[node]].push_back(node);
               comp_weights[components[node]] += residualGraph.getNodeWeight(node); 
        } endfor
       
        NodeID s = residualGraph.number_of_nodes()-2;
        NodeID t = residualGraph.number_of_nodes()-1;
        int comp_of_s = components[s];
        int comp_of_t = components[t];
        
        graph_access scc_graph;
        build_internal_scc_graph( residualGraph, components, comp_count, scc_graph);
         
        std::vector<int> comp_for_rhs;
        compute_new_rhs(scc_graph, config, comp_weights, comp_of_s, comp_of_t, perfect_rhs_weight, comp_for_rhs);
   
        //add comp_for_rhs nodes to new rhs   
        if(!config.mode_node_separators) {
                for( unsigned i = 0; i < comp_for_rhs.size(); i++) {    
                        int cur_component = comp_for_rhs[i];
                        if(cur_component != comp_of_s && cur_component != comp_of_t) {
                                for( unsigned j = 0; j < comp_nodes[cur_component].size(); j++) {
                                        new_rhs_nodes.push_back(comp_nodes[cur_component][j]);
                                }
                        }
                }
        } else {
                for( unsigned i = 0; i < comp_for_rhs.size(); i++) {    
                        int cur_component = comp_for_rhs[i];
                        if(cur_component != comp_of_s ) {
                                for( unsigned j = 0; j < comp_nodes[cur_component].size(); j++) {
                                        if( comp_nodes[cur_component][j] != t )
                                                new_rhs_nodes.push_back(comp_nodes[cur_component][j]);
                                }
                        }
                }

        }
} 

void most_balanced_minimum_cuts::compute_new_rhs( graph_access & scc_graph, 
                                                  const PartitionConfig & config,
                                                  std::vector< NodeWeight > & comp_weights,
                                                  int comp_of_s,
                                                  int comp_of_t,
                                                  NodeWeight optimal_rhs_stripe_weight,
                                                  std::vector<int> & comp_for_rhs) {
        //all successors of s cant be in any closure so they are marked invalid / this is basically a bfs in a DAG
        std::vector<bool> valid_to_add(scc_graph.number_of_nodes(), true);
        std::queue<NodeID> node_queue;
        node_queue.push(comp_of_s);
        valid_to_add[comp_of_s] = false;

        while (!node_queue.empty()) {
                NodeID node = node_queue.front();
                node_queue.pop();
                forall_out_edges(scc_graph, e, node) {
                        NodeID target = scc_graph.getEdgeTarget(e);
                        if(valid_to_add[target] == true) {
                                valid_to_add[target] = false;
                                node_queue.push(target);
                        }
                } endfor
        }

        std::vector<int> tmp_comp_for_rhs;
        int best_diff = std::numeric_limits<int>::max();
        for(unsigned i = 0; i < config.toposort_iterations; i++) {
                topological_sort ts;
                std::vector<NodeID> sorted_sequence;
                ts.sort(scc_graph, sorted_sequence); 

                tmp_comp_for_rhs.clear();

                bool t_contained   = false;
                int cur_rhs_weight = 0;
                int diff           = std::numeric_limits<int>::max();
                for( unsigned idx = 0; idx < sorted_sequence.size(); idx++) {
                        int cur_component = sorted_sequence[idx];

                        if( cur_component == comp_of_t ) {
                                t_contained = true;
                                tmp_comp_for_rhs.push_back(cur_component);
                                continue;
                        }

                        if(valid_to_add[cur_component]) {
                                int tmpdiff = optimal_rhs_stripe_weight - cur_rhs_weight - comp_weights[cur_component];
                                bool would_break = tmpdiff <= 0 && t_contained;
                                if(!would_break) {
                                        tmp_comp_for_rhs.push_back(cur_component);
                                        cur_rhs_weight += comp_weights[cur_component];
                                } else {
                                        //decide wether we should add this component now
                                        if(abs(tmpdiff) < abs(diff)) {
                                                //add it 
                                                tmp_comp_for_rhs.push_back(cur_component);
                                                cur_rhs_weight += comp_weights[cur_component];
                                                diff = optimal_rhs_stripe_weight - cur_rhs_weight;
                                        }       
                                        break;                                
                                }

                        }       

                        diff = optimal_rhs_stripe_weight - cur_rhs_weight;
                        if( diff <= 0 && t_contained) {
                                break;
                        }
                       
                }

                if(abs(diff) < best_diff) {
                        best_diff    = abs(diff);
                        comp_for_rhs = tmp_comp_for_rhs;
                }
        }
}

void most_balanced_minimum_cuts::build_internal_scc_graph( graph_access & residualGraph,  
                                                           std::vector<int> & components, 
                                                           int comp_count, 
                                                           graph_access & scc_graph) {
        std::vector< std::vector<int> > edges(comp_count);
        unsigned edge_count = 0;
        forall_nodes(residualGraph, node) {
                forall_out_edges(residualGraph, e, node) {
                        NodeID target = residualGraph.getEdgeTarget(e);
                        if(components[node] != components[target]) {
                                edges[components[node]].push_back(components[target]);
                                edge_count++;
                        }
                } endfor
        } endfor

        //build_scc_graph 
        scc_graph.start_construction(comp_count, edge_count);
        for( unsigned i = 0; i < (unsigned) comp_count; i++) {
                NodeID node = scc_graph.new_node();
                std::unordered_map<NodeID, bool> allready_contained;
                for(unsigned j = 0; j < edges[i].size(); j++) {
                        if(allready_contained.find(edges[i][j]) == allready_contained.end()) {
                                scc_graph.new_edge(node, edges[i][j]);
                                allready_contained[edges[i][j]] = true;
                        }
                }
        }

        scc_graph.finish_construction();
}
 
