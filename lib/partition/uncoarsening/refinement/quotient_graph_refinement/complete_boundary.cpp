/******************************************************************************
 * complete_boundary.cpp 
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

#include "complete_boundary.h"
#include "quality_metrics.h"

complete_boundary::complete_boundary(graph_access * G) {
        m_graph_ref   = G;
        m_pb_lhs_lazy = 0;
        m_pb_rhs_lazy = 0;
        m_last_pair   = 0;
        m_last_key    = -1;
        m_block_infos.resize(G->get_partition_count());
        delete Q.graphref;
        Q.graphref    = NULL;
}

complete_boundary::~complete_boundary() {
}

void complete_boundary::updatePeripheralCutEdges(NodeID node, boundary_pair * pair) {
  graph_access & G = *m_graph_ref;
        PartitionID to   = m_graph_ref->getPartitionIndex(node);
        PartitionID from = to == pair->lhs ? pair->rhs : pair->lhs;
        ASSERT_NEQ(from, to);
        forall_out_edges(G, e, node) {
                NodeID target = G.getEdgeTarget(e);
                PartitionID targetPartition = G.getPartitionIndex(target);
                EdgeWeight edge_weight = G.getEdgeWeight(e);

                if (targetPartition != from && targetPartition != to){
                //   bool target_is_still_incident_to_from = false;
                //   //this should only be delete if there is other incident partition
                //   forall_out_edges(G, t_e, target) {
                //     NodeID targets_target = G.getEdgeTarget(t_e);
                //     NodeID targets_target_partition = G.getPartitionIndex(targets_target);
                //     if(targets_target_partition == from) {
                //       //since partition index of node is to it cant be node, and this edge is
                //       //a widness that target can remain in this boundary
                //       target_is_still_incident_to_from = true;
                //       break;
                //     }
                //   } endfor

                        boundary_pair delete_bp;
                  delete_bp.k   = m_graph_ref->get_partition_count();
                  delete_bp.lhs = from;
                  delete_bp.rhs = targetPartition;
                  // deleteNode(node, from, &delete_bp);

                  // if (!target_is_still_incident_to_from) {
                  //    deleteNode(target, targetPartition, &delete_bp);
                  // }

                  boundary_pair insert_bp;
                  insert_bp.k   = m_graph_ref->get_partition_count();
                  insert_bp.lhs = to;
                  insert_bp.rhs = targetPartition;

                  // insert(node, to, &insert_bp);
                  // insert(target, targetPartition, &insert_bp);

                  m_pairs[delete_bp].edge_cut -= edge_weight;
                  m_pairs[insert_bp].edge_cut += edge_weight;


                  }
        } endfor

}




void complete_boundary::postMovedBoundaryNodeUpdates(NodeID node, boundary_pair * pair, 
                                                     bool update_edge_cuts, bool update_all_boundaries) {

        graph_access & G = *m_graph_ref;
        PartitionID to   = m_graph_ref->getPartitionIndex(node);
        PartitionID from = to == pair->lhs ? pair->rhs : pair->lhs;
        ASSERT_NEQ(from, to);
        //First delete this node from all incidient partition boudnary and decreas the edgecut (from, target_partition != to)
        //then insert it in the right target
        forall_out_edges(G, e, node) {

                NodeID target = G.getEdgeTarget(e);
                PartitionID targetPartition = G.getPartitionIndex(target);


                if(update_all_boundaries || targetPartition != to ) {
                        //delete 
                        boundary_pair delete_bp;
                        delete_bp.k   = m_graph_ref->get_partition_count();
                        delete_bp.lhs = from;
                        delete_bp.rhs = targetPartition;

                        EdgeWeight edge_weight = G.getEdgeWeight(e);
                        if(targetPartition != from) {
                                deleteNode(node, from, &delete_bp);

                                bool target_is_still_incident = false;
                                //this should only be delete if there is other incident partition
                                forall_out_edges(G, t_e, target) {
                                        NodeID targets_target = G.getEdgeTarget(t_e); 
                                        NodeID targets_target_partition = G.getPartitionIndex(targets_target); 
                                        if(targets_target_partition == from) {
                                                //since partition index of node is to it cant be node, and this edge is 
                                                //a widness that target can remain in this boundary
                                                target_is_still_incident = true;
                                                break;
                                        }
                                } endfor

                                      // std::cout << "target is still incident to lhs" <<  target_is_still_incident << std::endl;

                                if(!target_is_still_incident)
                                        deleteNode(target, targetPartition, &delete_bp);

                                if(update_edge_cuts) {
                                        m_pairs[delete_bp].edge_cut -= edge_weight;    
                                }
                        }

                        if(targetPartition != to) {
                                //insert
                                boundary_pair insert_bp;
                                insert_bp.k   = m_graph_ref->get_partition_count();
                                insert_bp.lhs = to;
                                insert_bp.rhs = targetPartition;

                                insert(node, to, &insert_bp);
                                insert(target, targetPartition, &insert_bp); 

                                if(update_edge_cuts) {
                                        m_pairs[insert_bp].edge_cut += edge_weight;    
                                }
                        }
                } 
        } endfor
}

void complete_boundary::balance_singletons(const PartitionConfig & config, graph_access & G) {
        for( unsigned i = 0; i < m_singletons.size(); i++) {
                NodeWeight min = m_block_infos[0].block_weight;
                PartitionID p  = 0;
                for( unsigned j = 0; j < m_block_infos.size(); j++) {
                        if( m_block_infos[j].block_weight < min ) {
                                min = m_block_infos[j].block_weight;
                                p = j;
                        }
                }

                NodeID node = m_singletons[i];
                if( m_block_infos[p].block_weight + G.getNodeWeight(node) <= config.upper_bound_partition) {
                        m_block_infos[G.getPartitionIndex(node)].block_weight -= G.getNodeWeight(node);
                        m_block_infos[p].block_weight += G.getNodeWeight(node);
                        G.setPartitionIndex(node, p);
                }
        }
}
