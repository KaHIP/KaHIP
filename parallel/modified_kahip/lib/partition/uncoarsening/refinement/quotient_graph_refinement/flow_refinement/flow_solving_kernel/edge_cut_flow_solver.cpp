/******************************************************************************
 * edge_cut_flow_solver.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/


#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <unordered_map>

#include "edge_cut_flow_solver.h"
#include "flow_macros.h"
#include "most_balanced_minimum_cuts/most_balanced_minimum_cuts.h"


edge_cut_flow_solver::edge_cut_flow_solver() {
}

edge_cut_flow_solver::~edge_cut_flow_solver() {
}

EdgeID edge_cut_flow_solver::regions_no_edges( graph_access & G,
                                               std::vector<NodeID> & lhs_boundary_stripe,
                                               std::vector<NodeID> & rhs_boundary_stripe,
                                               PartitionID & lhs, 
                                               PartitionID & rhs,
                                               std::vector<NodeID> & outer_lhs_boundary_nodes,
                                               std::vector<NodeID> & outer_rhs_boundary_nodes ) {

        EdgeID no_of_edges = 0;
        unsigned idx = 0;
        for( unsigned i = 0; i < lhs_boundary_stripe.size(); i++, idx++) {
                NodeID node = lhs_boundary_stripe[i];
                bool is_outer_boundary = false;
                forall_out_edges(G, e, node) {
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == BOUNDARY_STRIPE_NODE) no_of_edges++;
                        else is_outer_boundary = true;
                } endfor
                if(is_outer_boundary) {
                        outer_lhs_boundary_nodes.push_back(idx);
                }
        }

        for( unsigned i = 0; i < rhs_boundary_stripe.size(); i++, idx++) {
                NodeID node = rhs_boundary_stripe[i];
                bool is_outer_boundary = false;
                forall_out_edges(G, e, node) {
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == BOUNDARY_STRIPE_NODE) no_of_edges++;
                        else is_outer_boundary = true;
                } endfor
                if(is_outer_boundary) {
                        outer_rhs_boundary_nodes.push_back(idx);
                }
        }

        return no_of_edges; 
}

EdgeWeight edge_cut_flow_solver::convert_ds( const PartitionConfig & config, 
                                             graph_access & G, 
                                             PartitionID & lhs, 
                                             PartitionID & rhs, 
                                             std::vector<NodeID> & lhs_boundary_stripe,
                                             std::vector<NodeID> & rhs_boundary_stripe,
                                             std::vector<NodeID> & new_to_old_ids,              
                                             long *n_ad, 
                                             long* m_ad, 
                                             node** nodes_ad, 
                                             arc** arcs_ad, 
                                             long ** cap_ad,
                                             node** source_ad, 
                                             node** sink_ad, 
                                             long* node_min_ad,
                                             EdgeID & no_edges_in_flow_graph) {

        //should soon be refactored
        #include "convert_ds_variables.h"

        //building up the graph as in parse.h of hi_pr code
        NodeID idx = 0;
        new_to_old_ids.resize(lhs_boundary_stripe.size() + rhs_boundary_stripe.size());
        std::unordered_map<NodeID, NodeID> old_to_new;
        for( unsigned i = 0; i < lhs_boundary_stripe.size(); i++) {
                G.setPartitionIndex(lhs_boundary_stripe[i], BOUNDARY_STRIPE_NODE);
                new_to_old_ids[idx]                = lhs_boundary_stripe[i];
                old_to_new[lhs_boundary_stripe[i]] = idx++ ;
        }
        for( unsigned i = 0; i < rhs_boundary_stripe.size(); i++) {
                G.setPartitionIndex(rhs_boundary_stripe[i], BOUNDARY_STRIPE_NODE); 
                new_to_old_ids[idx]                = rhs_boundary_stripe[i];
                old_to_new[rhs_boundary_stripe[i]] = idx++;
        }

        std::vector<NodeID>  outer_lhs_boundary;
        std::vector<NodeID>  outer_rhs_boundary;
        EdgeID no_edges = regions_no_edges(G, lhs_boundary_stripe, rhs_boundary_stripe, 
                                              lhs, rhs, 
                                              outer_lhs_boundary, outer_rhs_boundary);
        no_edges_in_flow_graph = no_edges;
        
        if(outer_lhs_boundary.size() == 0 || outer_rhs_boundary.size() == 0) return false;
        n = lhs_boundary_stripe.size() + rhs_boundary_stripe.size() + 2; //+source and target
        m = no_edges + outer_lhs_boundary.size() + outer_rhs_boundary.size(); 

        nodes    = (node*) calloc ( n+2, sizeof(node) );
        arcs     = (arc*)  calloc ( 2*m+1, sizeof(arc) );
        arc_tail = (long*) calloc ( 2*m,   sizeof(long) ); 
        arc_first= (long*) calloc ( n+2, sizeof(long) );
        acap     = (long*) calloc ( 2*m, sizeof(long) );
        arc_current = arcs;

        node_max = 0;
        node_min = n;

        unsigned nodeoffset = 1;
        source              = n - 2 + nodeoffset;
        sink                = source+1;
        idx                 = 0;
        for( unsigned i = 0; i < lhs_boundary_stripe.size(); i++, idx++) {
                NodeID node = lhs_boundary_stripe[i];
                NodeID sourceID = idx + nodeoffset;
                forall_out_edges(G, e, node) {
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == BOUNDARY_STRIPE_NODE)  {
                                NodeID targetID     = old_to_new[G.getEdgeTarget(e)] + nodeoffset;
                                EdgeWeight capacity = G.getEdgeWeight(e);
                                tail                = sourceID;
                                head                = targetID;
                                cap                 = capacity;
                        
                                createEdge()                        
                        }
                } endfor
        }

        for( unsigned i = 0; i < rhs_boundary_stripe.size(); i++, idx++) {
                NodeID node = rhs_boundary_stripe[i];
                NodeID sourceID = idx + nodeoffset;
                forall_out_edges(G, e, node) {
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == BOUNDARY_STRIPE_NODE)  {
                                NodeID targetID     = old_to_new[G.getEdgeTarget(e)] + nodeoffset;
                                EdgeWeight capacity = G.getEdgeWeight(e);
                                tail                = sourceID;
                                head                = targetID;
                                cap                 = capacity;

                                createEdge()
                        }              
                } endfor
        }

        //connect source and target with outer boundary nodes 
        long max_capacity = std::numeric_limits<long>::max();
        for(unsigned i = 0; i < outer_lhs_boundary.size(); i++) {
                NodeID targetID = outer_lhs_boundary[i]+ nodeoffset;
                tail            = source;
                head            = targetID;
                cap             = max_capacity;

                createEdge()
        }

        for(unsigned i = 0; i < outer_rhs_boundary.size(); i++) {
                NodeID sourceID = outer_rhs_boundary[i]+ nodeoffset;
                tail            = sourceID;
                head            = sink;
                cap             = max_capacity;

                createEdge()
        }

        //this is so dirty ;)
        #include "linear_ordering_n_assign.h"
        
        /* Thanks God! all is done */
        return true;
}

EdgeWeight edge_cut_flow_solver::get_min_flow_max_cut(const PartitionConfig & config, 
                                                      graph_access & G, 
                                                      PartitionID & lhs, 
                                                      PartitionID & rhs, 
                                                      std::vector<NodeID> & lhs_boundary_stripe,
                                                      std::vector<NodeID> & rhs_boundary_stripe,
                                                      std::vector<NodeID> & new_to_old_ids,
                                                      EdgeWeight & initial_cut,
                                                      NodeWeight & rhs_part_weight,
                                                      NodeWeight & rhs_stripe_weight,
                                                      std::vector<NodeID> & new_rhs_nodes) {

        node *j = NULL;
        int  cc;
        bucket *l;


        globUpdtFreq = GLOB_UPDT_FREQ;

        EdgeID no_edges_in_flow_graph = 0;
        bool do_something = convert_ds(config, G, lhs, rhs, lhs_boundary_stripe, rhs_boundary_stripe, new_to_old_ids, &n, 
                        &m,
                        &nodes,
                        &arcs, 
                        &cap,
                        &source, 
                        &sink,
                        &nMin, 
                        no_edges_in_flow_graph );

        if(!do_something) return initial_cut; 

        cc = internal_allocDS();
        if ( cc ) { fprintf ( stderr, "Allocation error\n"); exit ( 1 ); }

        internal_init();
        internal_stage_one( );
        
        if(config.most_balanced_minimum_cuts) {
                internal_stage_two();
        }

        /* check if mincut is saturated */
        aMax = dMax = 0;
        for (l = buckets; l < buckets + n; l++) {
                l->firstActive = sentinelNode;
                l->firstInactive = sentinelNode;
        }
        internal_global_update();

        if(!config.most_balanced_minimum_cuts) {
                forAllNodes(j) {
                        if (j->d < n) {
                                new_rhs_nodes.push_back(nNode(j)-1);
                        }      
                }
        } else {
                node *i;
                node *t;
                long ni, na, flow_value, back_flow_value;
                arc *a;
                arc *innerStopA;

                NodeID no_nodes_flow_graph = lhs_boundary_stripe.size() + rhs_boundary_stripe.size()+2;
                graph_access residualGraph;
                residualGraph.start_construction(no_nodes_flow_graph,2*no_edges_in_flow_graph);

                //(u.v) \in E iff (u,v) in E and f_uv < c(u,v)
                //             or (v,u) in E and f_vu > 0
                forAllNodes(i) {
                        ni = nNode(i);
                        NodeID node = residualGraph.new_node(); // for each node here create a new node 
                        
                        if( (unsigned)(ni - 1) < new_to_old_ids.size()) { //note: unsigned has been introduced without testing
                                residualGraph.setNodeWeight( node, G.getNodeWeight(new_to_old_ids[ni-1]));
                        }

                        forAllArcs(i,a) {
                                na = nArc(a);
                                if ( cap[na] > 0 ) {
                                        flow_value = cap[na] - a->resCap;

                                        if( flow_value < cap[na] ) {
                                                //create that edge
                                                residualGraph.new_edge(node, nNode(a->head)-1);
                                        } else {
                                                //check wether backwards edge has positive flow
                                                t = a->head;
                                                arc* outarc; 
                                                bool prev_found = false;
                                                //we cannot use the makro here because it would overwrite stopA!
                                                for (outarc = t->first, innerStopA = (t+1)->first; outarc != innerStopA; outarc++) {
                                                        if(nNode(outarc->head) == ni) {

                                                                back_flow_value = cap[nArc(outarc)] - outarc->resCap;
                                                                if(back_flow_value > 0) {
                                                                        residualGraph.new_edge(node, nNode(a->head)-1); 
                                                                        break;
                                                                }
                                                                if(prev_found) {
                                                                        break;
                                                                }
                                                                prev_found = true;
                                                        }
                                                }

                                        }
                                }
                        }
                }

                residualGraph.finish_construction();
                NodeWeight average_partition_weight = ceil(config.largest_graph_weight / config.k);
                NodeWeight perfect_rhs_stripe_weight = abs((int)average_partition_weight - (int)rhs_part_weight+(int) rhs_stripe_weight);
                
                most_balanced_minimum_cuts mbmc;
                mbmc.compute_good_balanced_min_cut(residualGraph, config, perfect_rhs_stripe_weight, new_rhs_nodes);
        }
        return flow;
}

