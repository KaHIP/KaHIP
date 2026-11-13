/******************************************************************************
 * most_balanced_minimum_cuts.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <unordered_map>

#include <algorithm>
#include <queue>
#include "algorithms/topological_sort.h"
#include "random_functions.h"
#include "most_balanced_minimum_cuts.h"

most_balanced_minimum_cuts::most_balanced_minimum_cuts() {

}

most_balanced_minimum_cuts::~most_balanced_minimum_cuts() {

}

// Strongly connected components of the graph with the edges targets[first[v] .. first[v+1]) at v:
// the iterative Tarjan DFS of strongly_connected_components (same visiting order, so the same
// component numbers), with vector stacks kept between calls.
static int csr_strong_components( const std::vector<EdgeID> & first, const std::vector<NodeID> & targets,
                                  NodeID n, std::vector<int> & comp_num) {
        static std::vector<int> dfsnum;
        static std::vector<NodeID> unfinished;
        static std::vector<NodeID> roots;
        static std::vector< std::pair<NodeID,EdgeID> > iteration_stack;
        dfsnum.assign(n, -1);
        comp_num.assign(n, -1);
        int dfscount   = 0;
        int comp_count = 0;

        for( NodeID start = 0; start < n; start++) {
                if(dfsnum[start] != -1) continue;
                iteration_stack.push_back( std::pair<NodeID,EdgeID>( start, first[start] ) );
                dfsnum[start] = dfscount++;
                unfinished.push_back(start);
                roots.push_back(start);

                while( !iteration_stack.empty() ) {
                        NodeID current_node = iteration_stack.back().first;
                        EdgeID current_edge = iteration_stack.back().second;
                        iteration_stack.pop_back();

                        for( EdgeID e = current_edge, end = first[current_node+1]; e < end; ++e) {
                                NodeID target = targets[e];
                                if(dfsnum[target] == -1) {
                                        iteration_stack.push_back( std::pair<NodeID,EdgeID>( current_node, e ) );
                                        iteration_stack.push_back( std::pair<NodeID,EdgeID>( target, first[target] ) );
                                        dfsnum[target] = dfscount++;
                                        unfinished.push_back(target);
                                        roots.push_back(target);
                                        break;
                                } else if( comp_num[target] == -1) {
                                        while( dfsnum[roots.back()] > dfsnum[target] ) roots.pop_back();
                                }
                        }

                        if(current_node == roots.back()) {
                                NodeID w = 0;
                                do {
                                        w = unfinished.back();
                                        unfinished.pop_back();
                                        comp_num[w] = comp_count;
                                } while( w != current_node );
                                comp_count++;
                                roots.pop_back();
                        }
                }
        }
        return comp_count;
}

void most_balanced_minimum_cuts::compute_good_balanced_min_cut( graph_access & residualGraph, 
                                                                const PartitionConfig & config,
                                                                NodeWeight & perfect_rhs_weight, 
                                                                std::vector< NodeID > & new_rhs_nodes ) {
        std::vector<EdgeID> first(residualGraph.number_of_nodes()+1);
        std::vector<NodeID> targets(residualGraph.number_of_edges());
        std::vector<NodeWeight> weights(residualGraph.number_of_nodes());
        forall_nodes(residualGraph, node) {
                first[node]   = residualGraph.get_first_edge(node);
                weights[node] = residualGraph.getNodeWeight(node);
                forall_out_edges(residualGraph, e, node) {
                        targets[e] = residualGraph.getEdgeTarget(e);
                } endfor
        } endfor
        first[residualGraph.number_of_nodes()] = residualGraph.number_of_edges();
        compute_good_balanced_min_cut(first, targets, weights, config, perfect_rhs_weight, new_rhs_nodes);
}

// topological_sort::sort on CSR arrays: the same random start order (same draws), the same
// depth-first search with an explicit stack (same finishing order), reversed
static void csr_topological_sort( const std::vector<NodeID> & dag_first, const std::vector<NodeID> & dag_targets,
                                  std::vector<NodeID> & sorted_sequence) {
        const NodeID n = dag_first.size() - 1;
        static std::vector<char> visited;
        static std::vector<std::pair<NodeID, NodeID> > stack;
        visited.assign(n, 0);
        std::vector<NodeID> nodes(n);
        random_functions::permutate_vector_good(nodes, true);
        sorted_sequence.reserve(n);
        for( NodeID i = 0; i < n; i++) {
                NodeID start = nodes[i];
                if(visited[start]) continue;
                visited[start] = 1;
                stack.clear();
                stack.push_back(std::make_pair(start, dag_first[start]));
                while(!stack.empty()) {
                        std::pair<NodeID, NodeID> & top = stack.back();
                        if(top.second < dag_first[top.first+1]) {
                                NodeID target = dag_targets[top.second++];
                                if(!visited[target]) {
                                        visited[target] = 1;
                                        stack.push_back(std::make_pair(target, dag_first[target]));
                                }
                        } else {
                                sorted_sequence.push_back(top.first);
                                stack.pop_back();
                        }
                }
        }
        std::reverse(sorted_sequence.begin(), sorted_sequence.end());
}

static void compute_new_rhs_csr( const std::vector<NodeID> & dag_first, const std::vector<NodeID> & dag_targets, 
                                                  const PartitionConfig & config,
                                                  std::vector< NodeWeight > & comp_weights,
                                                  int comp_of_s,
                                                  int comp_of_t,
                                                  NodeWeight optimal_rhs_stripe_weight,
                                                  std::vector<int> & comp_for_rhs) {
        //all successors of s cant be in any closure so they are marked invalid / this is basically a bfs in a DAG
        const NodeID dag_n = dag_first.size() - 1;
        std::vector<bool> valid_to_add(dag_n, true);
        std::queue<NodeID> node_queue;
        node_queue.push(comp_of_s);
        valid_to_add[comp_of_s] = false;

        while (!node_queue.empty()) {
                NodeID node = node_queue.front();
                node_queue.pop();
                for( NodeID e = dag_first[node]; e < dag_first[node+1]; e++) {
                        NodeID target = dag_targets[e];
                        if(valid_to_add[target] == true) {
                                valid_to_add[target] = false;
                                node_queue.push(target);
                        }
                }
        }

        std::vector<int> tmp_comp_for_rhs;
        int best_diff = std::numeric_limits<int>::max();
        for(unsigned i = 0; i < config.toposort_iterations; i++) {
                std::vector<NodeID> sorted_sequence;
                csr_topological_sort(dag_first, dag_targets, sorted_sequence);

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
                                int tmpdiff = (int)optimal_rhs_stripe_weight - (int)cur_rhs_weight - (int)comp_weights[cur_component];
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

void most_balanced_minimum_cuts::compute_good_balanced_min_cut( const std::vector<EdgeID> & first,
                                                                const std::vector<NodeID> & targets,
                                                                const std::vector<NodeWeight> & weights,
                                                                const PartitionConfig & config,
                                                                NodeWeight & perfect_rhs_weight, 
                                                                std::vector< NodeID > & new_rhs_nodes ) {

        NodeID n = weights.size();
        static std::vector<int> components;
        int comp_count = csr_strong_components(first, targets, n, components);

        // nodes of each component in increasing id order: comp_nodes[comp_begin[c] .. comp_begin[c+1])
        static std::vector<NodeID> comp_begin;
        static std::vector<NodeID> comp_nodes;
        static std::vector<NodeID> pos;
        std::vector< NodeWeight > comp_weights(comp_count, 0);
        comp_begin.assign(comp_count+1, 0);
        for( NodeID node = 0; node < n; node++) {
               comp_begin[components[node]+1]++;
               comp_weights[components[node]] += weights[node]; 
        }
        for( int c = 0; c < comp_count; c++) comp_begin[c+1] += comp_begin[c];
        pos.assign(comp_begin.begin(), comp_begin.end()-1);
        comp_nodes.resize(n);
        for( NodeID node = 0; node < n; node++) comp_nodes[pos[components[node]]++] = node;
       
        NodeID s = n-2;
        NodeID t = n-1;
        int comp_of_s = components[s];
        int comp_of_t = components[t];
        
        // component graph: the edges between different components, grouped by source component
        // in node and edge order, each target component once per source component (first occurrence)
        static std::vector<NodeID> cedge_begin;
        static std::vector<int> cedges;
        static std::vector<int> contained;
        cedge_begin.assign(comp_count+1, 0);
        unsigned edge_count = 0;
        for( NodeID node = 0; node < n; node++) {
                for( EdgeID e = first[node]; e < first[node+1]; e++) {
                        if(components[node] != components[targets[e]]) {
                                cedge_begin[components[node]+1]++;
                                edge_count++;
                        }
                }
        }
        for( int c = 0; c < comp_count; c++) cedge_begin[c+1] += cedge_begin[c];
        pos.assign(cedge_begin.begin(), cedge_begin.end()-1);
        cedges.resize(edge_count);
        for( NodeID node = 0; node < n; node++) {
                for( EdgeID e = first[node]; e < first[node+1]; e++) {
                        if(components[node] != components[targets[e]]) {
                                cedges[pos[components[node]]++] = components[targets[e]];
                        }
                }
        }

        // the component DAG as CSR arrays (the same edges in the same order as the graph built before)
        static std::vector<NodeID> dag_first;
        static std::vector<NodeID> dag_targets;
        dag_first.resize(comp_count + 1);
        dag_targets.clear();
        contained.assign(comp_count, -1);
        for( int i = 0; i < comp_count; i++) {
                dag_first[i] = dag_targets.size();
                for( NodeID j = cedge_begin[i]; j < cedge_begin[i+1]; j++) {
                        if(contained[cedges[j]] != i) {
                                dag_targets.push_back(cedges[j]);
                                contained[cedges[j]] = i;
                        }
                }
        }
        dag_first[comp_count] = dag_targets.size();
         
        std::vector<int> comp_for_rhs;
        compute_new_rhs_csr(dag_first, dag_targets, config, comp_weights, comp_of_s, comp_of_t, perfect_rhs_weight, comp_for_rhs);
   
        //add comp_for_rhs nodes to new rhs   
        if(!config.mode_node_separators) {
                for( unsigned i = 0; i < comp_for_rhs.size(); i++) {    
                        int cur_component = comp_for_rhs[i];
                        if(cur_component != comp_of_s && cur_component != comp_of_t) {
                                for( NodeID j = comp_begin[cur_component]; j < comp_begin[cur_component+1]; j++) {
                                        new_rhs_nodes.push_back(comp_nodes[j]);
                                }
                        }
                }
        } else {
                for( unsigned i = 0; i < comp_for_rhs.size(); i++) {    
                        int cur_component = comp_for_rhs[i];
                        if(cur_component != comp_of_s ) {
                                for( NodeID j = comp_begin[cur_component]; j < comp_begin[cur_component+1]; j++) {
                                        if( comp_nodes[j] != t )
                                                new_rhs_nodes.push_back(comp_nodes[j]);
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
                                int tmpdiff = (int)optimal_rhs_stripe_weight - (int)cur_rhs_weight - (int)comp_weights[cur_component];
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
