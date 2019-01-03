/******************************************************************************
 * construct_mapping.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <cstdio>

#include "construct_mapping.h"
#include "tools/random_functions.h"
#include "fast_construct_mapping.h"
#include "data_structure/priority_queues/maxNodeHeap.h"


construct_mapping::construct_mapping() {

}

construct_mapping::~construct_mapping() {

}

void construct_mapping::construct_initial_mapping( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank) {
        switch( config.construction_algorithm ) {
                case MAP_CONST_IDENTITY:
                        PRINT(std::cout <<  "running identity mapping"  << std::endl;)
                        construct_identity( config, C, D, perm_rank);
                        break;
                case MAP_CONST_RANDOM:
                        PRINT(std::cout <<  "running random initial mapping"  << std::endl;)
                        construct_random( config, C, D, perm_rank);
                        break;
                case MAP_CONST_OLDGROWING:
                        PRINT(std::cout <<  "running old growing"  << std::endl;)
                        construct_old_growing( config, C, D, perm_rank);
                        break;
                case MAP_CONST_OLDGROWING_FASTER:
                        PRINT(std::cout <<  "running faster growing"  << std::endl;)
                        construct_old_growing_faster( config, C, D, perm_rank);
                        break;
                case MAP_CONST_FASTHIERARCHY_BOTTOMUP:
                        PRINT(std::cout <<  "running fast hierarchy bottom up"  << std::endl;)
                        construct_fast_hierarchy_bottomup( config, C, D, perm_rank);
                        break;
                case MAP_CONST_FASTHIERARCHY_TOPDOWN:
                        PRINT(std::cout <<  "running fast hierarchy top down"  << std::endl;)
                        construct_fast_hierarchy_topdown( config, C, D, perm_rank);
                        break;
                default: 
                        PRINT(std::cout <<  "running identity mapping"  << std::endl;)
                        construct_identity( config, C, D, perm_rank);
        }
}

void construct_mapping::construct_old_growing_matrix( PartitionConfig & config, matrix & C, matrix & D, std::vector< NodeID > & perm_rank) {
        std::cout <<  "constructing initial mapping matrix version of growing"  << std::endl;

        //initialze perm rank
        //interpretation task 'node' is assinged to perm_rank[node] 
        for( unsigned int i = 0; i < perm_rank.size(); i++) {
                perm_rank[i] = UNASSIGNED;
        }

        NodeWeight max_vol      = 0;
        NodeWeight max_vol_elem = 0;
        for( unsigned int i = 0; i < C.get_x_dim(); i++) {
                NodeWeight cur_vol = 0;
                for( unsigned int j = 0; j < C.get_x_dim(); j++) {
                        cur_vol += C.get_xy(i,j);
                }

                if( cur_vol > max_vol ) {
                        max_vol = cur_vol;
                        max_vol_elem = i;
                }
        }

        NodeWeight min_dist      = std::numeric_limits< NodeWeight >::max();
        NodeWeight min_dist_elem = 0;
        for( unsigned int cpu = 0; cpu < D.get_x_dim(); cpu++) {
                NodeWeight cur_dist = 0;
                for( unsigned int cpu_bar = 0; cpu_bar < D.get_y_dim(); cpu_bar++) {
                        cur_dist += D.get_xy( cpu, cpu_bar );
                }

                if( cur_dist < min_dist ) {
                        min_dist = cur_dist;
                        min_dist_elem = cpu;
                }
        }

        std::vector< NodeWeight > core_assigned( C.get_x_dim(), UNASSIGNED);
        std::vector< NodeWeight > total_vol( C.get_x_dim(), 0); // store volume to assigned nodes
        std::vector< NodeWeight > total_dist( C.get_x_dim(), 0);

        std::vector< NodeWeight > unassigned_PEs;
        std::vector< NodeWeight > unassigned_tasks;
        for( unsigned int i = 0; i < C.get_x_dim(); i++) {
                unassigned_PEs.push_back(i);
                unassigned_tasks.push_back(i);
        }

        perm_rank[max_vol_elem]      = min_dist_elem;
        core_assigned[min_dist_elem] = ASSIGNED;

        std::swap(unassigned_PEs[min_dist_elem], unassigned_PEs[unassigned_PEs.size()-1 ]);
        unassigned_PEs.pop_back();

        std::swap(unassigned_tasks[max_vol_elem], unassigned_tasks[unassigned_tasks.size()-1 ]);
        unassigned_tasks.pop_back();

        for( unsigned int i = 0; i < unassigned_PEs.size(); i++) {
                NodeID cpu      = unassigned_PEs[i];
                total_dist[cpu] += D.get_xy( min_dist_elem, cpu );
        }

        for( unsigned int i = 0; i < C.get_x_dim(); i++) {
                total_vol[i] += C.get_xy(max_vol_elem, i);
        }

        while( unassigned_tasks.size() > 0 ) {
                max_vol      = 0;
                max_vol_elem = 0;
                int idx_task = 0;
                for( unsigned int i = 0; i < unassigned_tasks.size(); i++) {
                        NodeID task = unassigned_tasks[i];
                        if( total_vol[task] > max_vol ) {
                                max_vol = total_vol[task];
                                max_vol_elem = task;
                                idx_task = i;
                        }
                }

                NodeID cur_task = max_vol_elem;

                min_dist      = std::numeric_limits< NodeWeight >::max();
                min_dist_elem = std::numeric_limits< NodeWeight >::max();
                unsigned int idx = 0;
                for( unsigned int i = 0; i < unassigned_PEs.size(); i++) {
                        NodeID cpu      = unassigned_PEs[i];
                        if( core_assigned[cpu] == ASSIGNED ) { 
                                continue;
                        }
                        if( total_dist[cpu] <= min_dist ) {
                                min_dist = total_dist[cpu];
                                min_dist_elem = cpu;
                                idx = i;
                        }
                }
                NodeID cur_PE = min_dist_elem;

                perm_rank[cur_task]   = cur_PE;
                core_assigned[cur_PE] = ASSIGNED;

                std::swap(unassigned_tasks[idx_task], unassigned_tasks[unassigned_tasks.size()-1 ]);
                unassigned_tasks.pop_back();

                std::swap(unassigned_PEs[idx], unassigned_PEs[unassigned_PEs.size()-1 ]);
                unassigned_PEs.pop_back();

                // update the priority of the other PEs
                for( unsigned int i = 0; i < unassigned_PEs.size(); i++) {
                        NodeID cpu      = unassigned_PEs[i];
                        total_dist[cpu] += D.get_xy( cur_PE, cpu );
                }

                ////update priorities
                for( unsigned int i = 0; i < unassigned_tasks.size(); i++) {
                        NodeID target_task = unassigned_tasks[i];
                        total_vol[target_task] += C.get_xy(cur_task, target_task);
                }


        } 
}

void construct_mapping::construct_old_growing( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank) {
        std::cout <<  "constructing initial mapping"  << std::endl;

        std::vector< NodeWeight > total_dist( C.number_of_nodes(), 0);
        std::vector< NodeWeight > total_vol( C.number_of_nodes(), 0);
        std::vector< NodeWeight > core_assigned( C.number_of_nodes(), UNASSIGNED);

        NodeWeight max_vol      = 0;
        NodeWeight max_vol_elem = 0;
        forall_nodes(C, node) {
                forall_out_edges(C, e, node) {
                        total_vol[node] += C.getEdgeWeight(e);
                } endfor
                if( total_vol[node] > max_vol ) {
                        max_vol = total_vol[node];
                        max_vol_elem = node;
                }
        } endfor

        NodeWeight min_dist      = std::numeric_limits< NodeWeight >::max();
        NodeWeight min_dist_elem = 0;
        for( unsigned int cpu = 0; cpu < C.number_of_nodes(); cpu++) {
                total_dist[cpu] = 0;
                for( unsigned int cpu_bar = 0; cpu_bar < C.number_of_nodes(); cpu_bar++) {
                        total_dist[cpu] += D.get_xy( cpu, cpu_bar );
                }

                if( total_dist[cpu] < min_dist ) {
                        min_dist = total_dist[cpu];
                        min_dist_elem = cpu;
                }
        }

        //initialze perm rank
        //interpretation task 'node' is assinged to perm_rank[node] 
        for( unsigned int i = 0; i < perm_rank.size(); i++) {
                perm_rank[i] = UNASSIGNED;
        }
        perm_rank[max_vol_elem]      = min_dist_elem;
        core_assigned[min_dist_elem] = ASSIGNED;
        //initialization, now assign the rest of the ranks to the tasks
        for( unsigned i = 0; i < C.number_of_nodes()-1; i++) {
                max_vol      = 0;
                max_vol_elem = std::numeric_limits< NodeWeight >::max();
                forall_nodes(C, node) {
                        if(perm_rank[node] != UNASSIGNED) continue;

                        total_vol[node] = 0;
                        forall_out_edges(C, e, node) {
                                NodeID target = C.getEdgeTarget(e);
                                if( perm_rank[target] != UNASSIGNED ) {
                                        total_vol[node] += C.getEdgeWeight(e);
                                }
                        } endfor
                        if( total_vol[node] >= max_vol ) {
                                max_vol      = total_vol[node];
                                max_vol_elem = node;
                        }
                } endfor

                min_dist      = std::numeric_limits< NodeWeight >::max();
                min_dist_elem = std::numeric_limits< NodeWeight >::max();
                for( unsigned int cpu = 0; cpu < C.number_of_nodes(); cpu++) {
                        total_dist[cpu] = 0;
                        if( core_assigned[cpu] == ASSIGNED ) continue;
                        for( unsigned int cpu_bar = 0; cpu_bar < C.number_of_nodes(); cpu_bar++) {
                                if( core_assigned[cpu_bar] == ASSIGNED ) {
                                        total_dist[cpu] += D.get_xy( cpu, cpu_bar );
                                }
                        }

                        if( total_dist[cpu] <= min_dist ) {
                                min_dist = total_dist[cpu];
                                min_dist_elem = cpu;
                        }
                }

                perm_rank[max_vol_elem]      = min_dist_elem;
                core_assigned[min_dist_elem] = ASSIGNED;
        }
}

void construct_mapping::construct_old_growing_faster( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank) {
        std::cout <<  "constructing initial mapping with faster growing"  << std::endl;

        //initialze perm rank
        //interpretation task 'node' is assinged to perm_rank[node] 
        for( unsigned int i = 0; i < perm_rank.size(); i++) {
                perm_rank[i] = UNASSIGNED;
        }

        maxNodeHeap unassigned_tasks; // contains unassigned tasks and their priority
        NodeWeight max_vol      = 0;
        NodeWeight max_vol_elem = 0;
        forall_nodes(C, node) {
                NodeWeight cur_vol = 0;
                forall_out_edges(C, e, node) {
                        cur_vol += C.getEdgeWeight(e);
                } endfor

                unassigned_tasks.insert( node, 0);

                if( cur_vol > max_vol ) {
                        max_vol = cur_vol;
                        max_vol_elem = node;
                }
        } endfor

        NodeWeight min_dist      = std::numeric_limits< NodeWeight >::max();
        NodeWeight min_dist_elem = 0;
        for( unsigned int cpu = 0; cpu < C.number_of_nodes(); cpu++) {
                NodeWeight cur_dist = 0;
                for( unsigned int cpu_bar = 0; cpu_bar < C.number_of_nodes(); cpu_bar++) {
                        cur_dist += D.get_xy( cpu, cpu_bar );
                }

                if( cur_dist < min_dist ) {
                        min_dist = cur_dist;
                        min_dist_elem = cpu;
                }
        }

        std::vector< NodeWeight > core_assigned( C.number_of_nodes(), UNASSIGNED);
        std::vector< NodeWeight > total_vol( C.number_of_nodes(), 0); // store volume to assigned nodes
        std::vector< NodeWeight > total_dist( C.number_of_nodes(), 0);

        std::vector< NodeWeight > unassigned_PEs;
        for( unsigned int cpu = 0; cpu < C.number_of_nodes(); cpu++) {
                unassigned_PEs.push_back(cpu);
        }

        perm_rank[max_vol_elem]      = min_dist_elem;
        core_assigned[min_dist_elem] = ASSIGNED;

        std::swap(unassigned_PEs[min_dist_elem], unassigned_PEs[unassigned_PEs.size()-1 ]);
        unassigned_PEs.pop_back();
        unassigned_tasks.deleteNode(max_vol_elem);

        for( unsigned int i = 0; i < unassigned_PEs.size(); i++) {
                NodeID cpu      = unassigned_PEs[i];
                total_dist[cpu] += D.get_xy( min_dist_elem, cpu );
        }

        forall_out_edges(C, e, max_vol_elem) {
                NodeID target_task = C.getEdgeTarget(e);
                if( unassigned_tasks.contains(target_task) ) {
                        total_vol[target_task] += C.getEdgeWeight(e);
                        unassigned_tasks.changeKey(target_task, total_vol[target_task]);
                }
        } endfor

        while( unassigned_tasks.size() > 0 ) {
                NodeID cur_task = unassigned_tasks.deleteMax();

                min_dist      = std::numeric_limits< NodeWeight >::max();
                min_dist_elem = std::numeric_limits< NodeWeight >::max();
                unsigned int idx = 0;
                for( unsigned int i = 0; i < unassigned_PEs.size(); i++) {
                        NodeID cpu      = unassigned_PEs[i];
                        if( core_assigned[cpu] == ASSIGNED ) { 
                                continue;
                        }
                        if( total_dist[cpu] <= min_dist ) {
                                min_dist = total_dist[cpu];
                                min_dist_elem = cpu;
                                idx = i;
                        }
                }
                NodeID cur_PE = min_dist_elem;

                perm_rank[cur_task]   = cur_PE;
                core_assigned[cur_PE] = ASSIGNED;

                //update priorities
                forall_out_edges(C, e, cur_task) {
                        NodeID target_task = C.getEdgeTarget(e);
                        if( unassigned_tasks.contains(target_task) ) {
                                total_vol[target_task] += C.getEdgeWeight(e);
                                unassigned_tasks.changeKey(target_task, total_vol[target_task]);
                        }
                } endfor

                std::swap(unassigned_PEs[idx], unassigned_PEs[unassigned_PEs.size()-1 ]);
                unassigned_PEs.pop_back();

                // update the priority of the other PEs
                for( unsigned int i = 0; i < unassigned_PEs.size(); i++) {
                        NodeID cpu      = unassigned_PEs[i];
                        total_dist[cpu] += D.get_xy( cur_PE, cpu );
                }
        } 
}

void construct_mapping::construct_identity( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank) {
        for( unsigned int i = 0; i < perm_rank.size(); i++) {
                perm_rank[i] = i;
        }
}

void construct_mapping::construct_random( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank) {
        construct_identity( config, C, D, perm_rank);
        random_functions::permutate_vector_good( perm_rank, false);
}

void construct_mapping::construct_fast_hierarchy_bottomup( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank) {
        fast_construct_mapping fcm;
        fcm.construct_initial_mapping_bottomup( config, C, D, perm_rank);
}

void construct_mapping::construct_fast_hierarchy_topdown( PartitionConfig & config, graph_access & C, matrix & D, std::vector< NodeID > & perm_rank) {
        fast_construct_mapping fcm;
        fcm.construct_initial_mapping_topdown( config, C, D, perm_rank);
}


