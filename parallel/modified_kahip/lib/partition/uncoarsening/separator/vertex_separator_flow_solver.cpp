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

#include "vertex_separator_flow_solver.h"
#include "flow_solving_kernel/flow_macros.h"

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

        node *j = NULL;
        int  cc;
        bucket *l;


        globUpdtFreq = GLOB_UPDT_FREQ;
        std::vector<NodeID> new_to_old_ids;

        EdgeID no_edges_in_flow_graph = 0;
        bool success = construct_flow_pb(config, G, lhs, rhs, lhs_nodes, rhs_nodes, new_to_old_ids, 
                        &n, 
                        &m,
                        &nodes,
                        &arcs, 
                        &cap,
                        &source, 
                        &sink,
                        &nMin, 
                        no_edges_in_flow_graph );


        cc = internal_allocDS();
        if(!success) return; 
        if ( cc ) { fprintf ( stderr, "Allocation error\n"); exit ( 1 ); }

        internal_init();
        internal_stage_one( );

        internal_stage_two();

        /* check if mincut is saturated */
        aMax = dMax = 0;
        for (l = buckets; l < buckets + n; l++) {
                l->firstActive = sentinelNode;
                l->firstInactive = sentinelNode;
        }
        internal_global_update();

        std::vector<NodeID> S;
        forAllNodes(j) {
                if (!(j->d < n)) {
                        if((unsigned) (nNode(j) - 1) < (int)lhs_nodes.size() + rhs_nodes.size()) { //Note: unsigned has been introduced without testing
                                S.push_back(new_to_old_ids[nNode(j) -1]);
                        } 
                }      
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

bool vertex_separator_flow_solver::construct_flow_pb( const PartitionConfig & config, 
                                                      graph_access & G, 
                                                      PartitionID & lhs, 
                                                      PartitionID & rhs, 
                                                      std::vector<NodeID> & lhs_nodes,
                                                      std::vector<NodeID> & rhs_nodes,
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
        
        //very dirty for loading variables :). some time this should all be refactored. for now we can focus on the important stuff.
        #include "../refinement/quotient_graph_refinement/flow_refinement/flow_solving_kernel/convert_ds_variables.h"

        //building up the graph as in parse.h of hi_pr code
        //first we have to count the number of edges 
        // s to lhs + rhs to t + lhs to rhs 
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

        n = lhs_nodes.size() + rhs_nodes.size() + 2; //+source and target
        m = no_edges + lhs_nodes.size() + rhs_nodes.size(); 

        nodes    = (node*) calloc ( n+2, sizeof(node) );
        arcs     = (arc*)  calloc ( 2*m+1, sizeof(arc) );
        arc_tail = (long*) calloc ( 2*m,   sizeof(long) ); 
        arc_first= (long*) calloc ( n+2, sizeof(long) );
        acap     = (long*) calloc ( 2*m, sizeof(long) );
        arc_current = arcs;

        node_max = 0;
        node_min = n;

        if(n == 2) return false;

        unsigned nodeoffset = 1; 
        source = n - 2 + nodeoffset;
        sink = source+1;

        idx = 0;
        long max_capacity = std::numeric_limits<long>::max();
        //insert directed edges from L to R
        for( unsigned i = 0; i < lhs_nodes.size(); i++, idx++) {
                NodeID node = lhs_nodes[i];
                NodeID sourceID = idx + nodeoffset;
                forall_out_edges(G, e, node) {
                        if(G.getPartitionIndex(G.getEdgeTarget(e)) == rhs)  {
                                NodeID targetID = old_to_new[G.getEdgeTarget(e)] + nodeoffset;
                                tail = sourceID;
                                head = targetID;
                                cap = max_capacity; 

                                createEdge()                                
                        }
                } endfor
        }

        //connect source and target with outer boundary nodes 
        for(unsigned i = 0; i < lhs_nodes.size(); i++) {
                NodeID targetID = old_to_new[lhs_nodes[i]]+nodeoffset;
                tail = source;
                head = targetID;
                cap = G.getNodeWeight(lhs_nodes[i]); 

                createEdge()                                
        }

        for(unsigned i = 0; i < rhs_nodes.size(); i++) {
                NodeID sourceID = old_to_new[rhs_nodes[i]]+ nodeoffset;
                tail = sourceID;
                head = sink;
                cap = G.getNodeWeight(rhs_nodes[i]); 

                createEdge()        
        }

        //very dirty
        #include "../refinement/quotient_graph_refinement/flow_refinement/flow_solving_kernel/linear_ordering_n_assign.h"
        /* Thanks God! all is done */

        return true;



}


