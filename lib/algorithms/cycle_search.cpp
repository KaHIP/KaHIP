/******************************************************************************
 * cycle_search.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>

#include "algorithms/strongly_connected_components.h"
#include "cycle_search.h"
#include "random_functions.h"
#include "timer.h"

double cycle_search::total_time = 0;

cycle_search::cycle_search() {

}

cycle_search::~cycle_search() {

}

void cycle_search::find_random_cycle(graph_access & G, std::vector<NodeID> & cycle) {
	//first perform a bfs starting from a random node and build the parent array
        std::deque<NodeID>* bfsqueue = new std::deque<NodeID>;
	NodeID v = random_functions::nextInt(0, G.number_of_nodes()-1);
	bfsqueue->push_back(v); 

	std::vector<bool>   touched(G.number_of_nodes(),false);
	std::vector<bool>   is_leaf(G.number_of_nodes(),false);
	std::vector<NodeID> parent(G.number_of_nodes(),0);
	std::vector<NodeID> leafes;
	touched[v] = true;
	parent[v]  = v;

	while(!bfsqueue->empty()) {
		NodeID source = bfsqueue->front();
		bfsqueue->pop_front();

		bool is_leaf = true;
		forall_out_edges(G, e, source) {
			NodeID target = G.getEdgeTarget(e);
			if(!touched[target]) {
				is_leaf         = false;
				touched[target] = true; 
				parent[target]  = source;
				bfsqueue->push_back(target);
			}
		} endfor

		if(is_leaf) 
			leafes.push_back(source);

	}

	std::vector<NodeID> sources(G.number_of_edges(), 0);
	std::vector<NodeID> targets(G.number_of_edges(), 0);
	forall_nodes(G, node) {
		forall_out_edges(G, e, node) {
			NodeID target = G.getEdgeTarget(e);
			sources[e] = node;
			targets[e] = target;
		} endfor
	} endfor


	//now find two random leafes 
	NodeID v_1, v_2;
	unsigned r_idx = random_functions::nextInt(0, G.number_of_edges()-1);
	while(true) {
		NodeID source = sources[r_idx];
		NodeID target = targets[r_idx];
		if( parent[source] != target && parent[target] != source) {
			//found a non-tree edge
			v_1 = source;
			v_2 = target;
			break;
		} 

		r_idx = random_functions::nextInt(0, G.number_of_edges()-1);
	}

	// NodeID now climb up the parent array step wise left and right
	std::vector<NodeID> lhs_path, rhs_path;
	lhs_path.push_back(v_1);
	rhs_path.push_back(v_2);

	std::vector<bool> touched_nodes(G.number_of_nodes(),false);
	std::vector<unsigned> index(G.number_of_nodes(),0);
	index[v_1] = 0;
	index[v_2] = 0;
	touched_nodes[v_1] = true;
	touched_nodes[v_2] = true;

	NodeID cur_lhs = v_1, cur_rhs = v_2; 
	NodeID counter = 0;
	bool break_lhs = false;
	while(true) {
		counter++;
		if(cur_lhs != parent[cur_lhs]) {
			if(touched_nodes[parent[cur_lhs]] == true) {
				break_lhs = true;
				lhs_path.push_back(parent[cur_lhs]);
				break;
			} else {
				cur_lhs                = parent[cur_lhs];
				touched_nodes[cur_lhs] = true;
				lhs_path.push_back(cur_lhs);
				index[cur_lhs] = counter;
			}
		}
		if(cur_rhs != parent[cur_rhs]) {
			if(touched_nodes[parent[cur_rhs]] == true) {
				rhs_path.push_back(parent[cur_rhs]);
				break;
			} else {
				cur_rhs                = parent[cur_rhs];
				touched_nodes[cur_rhs] = true;
				rhs_path.push_back(cur_rhs);
				index[cur_rhs] = counter;
			}
		}

	}

	if(break_lhs) {
		for( unsigned i = 0; i < lhs_path.size(); i++) {
			cycle.push_back(lhs_path[i]);
		}

		NodeID connecting_vertice = cycle[cycle.size()-1];
		for( int i = index[connecting_vertice]-1; i >= 0; i--) {
			cycle.push_back(rhs_path[i]);
		}
	} else {
		for( unsigned i = 0; i < rhs_path.size(); i++) {
			cycle.push_back(rhs_path[i]);
		}

                NodeID connecting_vertice = cycle[cycle.size()-1];
                for( int i = index[connecting_vertice]-1; i >= 0; i--) {
                        cycle.push_back(lhs_path[i]);
                }
        }

        cycle.push_back(cycle[0]);
}

bool cycle_search::find_shortest_path(graph_access & G, 
                                      NodeID & start, 
                                      NodeID & dest, 
                                      std::vector<NodeID> & cycle) {

        std::vector<EdgeWeight> distance(G.number_of_nodes(), std::numeric_limits<EdgeWeight>::max()/2);
        std::vector<NodeID> parent(G.number_of_nodes(), std::numeric_limits<NodeID>::max());

        bool negative_cycle_detected = negative_cycle_detection(G, start, distance, parent, cycle);

        if( !negative_cycle_detected) {
                //if there is no negative cycle then we should return a shortest path from s to t
                cycle.clear();
                cycle.push_back(dest);
                NodeID cur = dest;
                while(cur != start) {
                        cur = parent[cur];
                        cycle.push_back(cur);
                }
                std::reverse(cycle.begin(), cycle.end());
        } 
        return negative_cycle_detected;
}

bool cycle_search::find_negative_cycle(graph_access & G, NodeID & start, std::vector<NodeID> & cycle) {
	//simplest bellman ford algorithm

	std::vector<EdgeWeight> distance(G.number_of_nodes(), std::numeric_limits<EdgeWeight>::max()/2);
        std::vector<NodeID> parent(G.number_of_nodes(), std::numeric_limits<NodeID>::max());

        return negative_cycle_detection(G, start, distance, parent, cycle);
}

int cycle_search::bellman_ford_with_subtree_disassembly_and_updates(graph_access & G, 
                                                                    NodeID & start, 
                                                                    std::vector<EdgeWeight> & distance, 
                                                                    std::vector<NodeID> & parent, 
                                                                    std::vector<NodeID> & cycle) {
        // Goldberg spc-1.2 similar implementation using our data structures 
        int NULL_NODE = -1;
        distance[start] = 0;
        std::queue<NodeID> L;

        //doubly linked list of shortest path tree in preorder
        short OUT_OF_QUEUE = 0;
        short INACTIVE     = 1;
        short ACTIVE       = 2;
        short IN_QUEUE     = 2;

        std::vector<int> before(G.number_of_nodes(), NULL_NODE);
        std::vector<int> after(G.number_of_nodes());
        std::vector<int> degree(G.number_of_nodes());
        std::vector<short> status(G.number_of_nodes(), OUT_OF_QUEUE);

        L.push(start);

        after[start]  = start;
        before[start] = start;
        degree[start] = -1;
        status[start] = IN_QUEUE;

        while( ! L.empty() ) {
                NodeID v = L.front();
                L.pop();

                short current_status = status[v];
                status[v]            = OUT_OF_QUEUE;

                if( current_status == INACTIVE ) continue;

                forall_out_edges(G, e, v) {
                        NodeID w = G.getEdgeTarget(e);
                        int delta = distance[w] - distance[v] - G.getEdgeWeight(e);
                        if(delta > 0) {
                                // dissassemble subtree
                                // in this case we disassemble the subtree looking for v
                                int new_distance = distance[w] - delta; 

                                int x = before[w];
                                int y = w;
                                if( x != NULL_NODE) {
                                        // in this case w is allready in the tree and we remove it / disassemble the subtree
                                        for( int total_degree = 0; total_degree >= 0; y = after[y]) {
                                                //disassemble the subtree 
                                                //w <-> .... <-> y <-> ...
                                                if( y == (int) v ) { 
                                                        parent[w]   = v;
                                                        return w; // since parent[w] = v  

                                                }  else{
                                                        distance[y]   = distance[y] - delta;
                                                        before[y]     = NULL_NODE;
                                                        total_degree += degree[y];

                                                        if( status[y] == ACTIVE ) {
                                                                status[y] = INACTIVE;
                                                        }
                                                }
                                        }

                                        // since we removed w from the shortest path tree
                                        degree[parent[w]]--;

                                        // the old subtree
                                        // x <-> [w <-> subtree]_is cut <-> y
                                        // y is the vertex after the subtree of w
                                        // afterwards: x <-> y
                                        after[x]         = y;
                                        before[y]        = x;
                                }
                                distance[w] = new_distance;
                                parent[w]   = v;

                        }

                        // negative cycle is not found
                        // ===================================================
                        // take care of the rest of the shortest path Tree T
                        // ===================================================
                        if( before[w] == NULL_NODE && parent[w] == v) {
                                // in this case w was not in the shortest path tree
                                // so we integrate it 
                                degree[v] ++;
                                degree[w] = -1;

                                NodeID after_v = after[v];

                                // integrate w into the tree
                                after[v]         = w;
                                before[w]        = v;
                                after[w]         = after_v;
                                before[after_v]  = w;
                                // we now have v <-> w <-> after_v in the preorder tree representation

                                // handle the queue
                                if( status[w] == OUT_OF_QUEUE ) {
                                        L.push(w);
                                        status[w] = IN_QUEUE;
                                } else {
                                        status[w] = ACTIVE;
                                }


                        } 
                } endfor
        }
        return NULL_NODE;
}



bool cycle_search::negative_cycle_detection(graph_access & G, 
                                            NodeID & start, 
                                            std::vector<EdgeWeight> & distance, 
                                            std::vector<NodeID> & parent, 
                                            std::vector<NodeID> & cycle) {
        timer timeR;
        
        int w = bellman_ford_with_subtree_disassembly_and_updates(G, start, distance, parent, cycle);
        
        if(w >= 0) { // found a cycle 
                // the edge yielding the cycle was (t,w)
                NodeID t = parent[w];
                NodeID u = t;

                std::vector<bool> seen(G.number_of_nodes(), false); //use hashing?
                seen[u] = true;
                NodeID predecessor = parent[u];
                NodeID start_vertex;

                while( true  ) {
                        if( seen[predecessor] ) {
                                start_vertex = predecessor;
                                break;
                        }
                        seen[predecessor] = true;
                        predecessor = parent[predecessor];
                }

                cycle.push_back(start_vertex);
                predecessor = parent[start_vertex];
                while( predecessor != start_vertex) {
                        cycle.push_back(predecessor);
                        predecessor = parent[predecessor];
                }
                cycle.push_back(start_vertex);
                std::reverse(cycle.begin(), cycle.end());

                total_time += timeR.elapsed();
                return true;

        } 

        total_time += timeR.elapsed();
	return false;

}

//preconditition: no negative cycles 
bool cycle_search::find_zero_weight_cycle(graph_access & G, NodeID & start, std::vector<NodeID> & cycle) {

        std::vector<EdgeWeight> distance(G.number_of_nodes(), std::numeric_limits<EdgeWeight>::max()/2);
	std::vector<NodeID> parent(G.number_of_nodes(), std::numeric_limits<NodeID>::max());
        bool negative_weight_cycle = negative_cycle_detection(G, start, distance, parent, cycle);
        if(!negative_weight_cycle) {
               //now we try to return a random directed zero weight gain cycle
               //therefore we use W(e) = d(u) + w(e) - d(v)
               //and keep edges with weight 0
               graph_access W;
               W.start_construction(G.number_of_nodes(), G.number_of_edges());

               forall_nodes(G, node) {
                       NodeID shadow_node              = W.new_node();
                       W.setNodeWeight(shadow_node, G.getNodeWeight(node));
                       forall_out_edges(G, e, node) {
                               NodeID target                   = G.getEdgeTarget(e);
                               EdgeWeight modified_edge_weight = G.getEdgeWeight(e) + distance[node] - distance[target]; 
                               ASSERT_GEQ(modified_edge_weight, 0);
                               if(modified_edge_weight == 0) {
                                        EdgeID shadow_edge              = W.new_edge(shadow_node, target);
                                        W.setEdgeWeight(shadow_edge, 0);
                               }
                       } endfor
               } endfor
               W.finish_construction();

               strongly_connected_components scc;
               std::vector<int> comp_num(W.number_of_nodes());
               scc.strong_components(W, comp_num);

               //first check wether there are components with more then one vertex
               std::vector<unsigned> comp_count(W.number_of_nodes(), 0);
               forall_nodes(W, node) {
                       comp_count[comp_num[node]]++;
               } endfor

               std::vector<NodeID> candidates;
               forall_nodes(W, node) {
                       if(comp_count[comp_num[node]] > 1) {
                               candidates.push_back(node);
                       }
               } endfor 

               if(candidates.size() == 0) {return false;}

               //now pick a random start vertex
               NodeID start_vertex_idx = random_functions::nextInt(0, candidates.size()-1);
               NodeID start_vertex = candidates[start_vertex_idx];
               std::vector<bool> seen(W.number_of_nodes(), false);
               std::vector<NodeID> list;

               NodeID successor = start_vertex;
               NodeID comp_of_sv = comp_num[start_vertex];
               do {

                       seen[successor] = true;
                       list.push_back(successor);

                       std::vector<NodeID> same_comp_neighbors;
                       forall_out_edges(W, e, successor) {
                               NodeID target = W.getEdgeTarget(e);
                               if(comp_num[target] == (int)comp_of_sv) {
                                        same_comp_neighbors.push_back(target);
                               }
                       } endfor
                       NodeID succ_id = random_functions::nextInt(0, same_comp_neighbors.size()-1);
                       
                       successor = same_comp_neighbors[succ_id];
               } while(seen[successor] == false);

               NodeID start_idx = 0;
               for( unsigned i = 0; i < list.size(); i++) {
                       if(list[i] == successor) {
                               start_idx = i;
                               break;
                       }
               }

               for( unsigned i = start_idx; i < list.size(); i++) {
                       cycle.push_back(list[i]);
               }
               cycle.push_back(successor);

               return true;
        }
        return false;

}
