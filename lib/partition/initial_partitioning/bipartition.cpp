/******************************************************************************
 * bipartition.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "bipartition.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "uncoarsening/refinement/quotient_graph_refinement/quotient_graph_refinement.h"
#include "uncoarsening/refinement/refinement.h"

bipartition::bipartition() {

}

bipartition::~bipartition() {

}

void bipartition::initial_partition( const PartitionConfig & config, 
                                     const unsigned int seed, 
                                     graph_access & G, 
                                     int* partition_map) {

        timer t;
        t.restart();
        unsigned iterations = config.bipartition_tries;
        EdgeWeight best_cut = std::numeric_limits<EdgeWeight>::max();
        int best_load       = std::numeric_limits<int>::max();

        for( unsigned i = 0; i < iterations; i++) {
                if(config.bipartition_algorithm == BIPARTITION_BFS)  {
                        grow_regions_bfs(config, G);
                } else if( config.bipartition_algorithm == BIPARTITION_FM) {
                        grow_regions_fm(config, G);
                } 

                G.set_partition_count(2);

                post_fm(config, G);

                quality_metrics qm;
                EdgeWeight curcut = qm.edge_cut(G); 

                int lhs_block_weight = 0;
                int rhs_block_weight = 0;

                forall_nodes(G, node) {
                        if(G.getPartitionIndex(node) == 0) {
                                lhs_block_weight += G.getNodeWeight(node);
                        } else {
                                rhs_block_weight += G.getNodeWeight(node);
                        }
                } endfor

                int lhs_overload = std::max(lhs_block_weight - config.target_weights[0],0);
                int rhs_overload = std::max(rhs_block_weight - config.target_weights[1],0);

                if(curcut < best_cut || (curcut == best_cut && lhs_overload + rhs_block_weight < best_load) ) {
                        //store it
                        best_cut  = curcut;
                        best_load = lhs_overload + rhs_overload;

                        forall_nodes(G, n) {
                                partition_map[n] =  G.getPartitionIndex(n);
                        } endfor
                } 

        }
        PRINT(std::cout <<  "bipartition took " <<  t.elapsed()  << std::endl;)
}

void bipartition::initial_partition( const PartitionConfig & config, 
                                     const unsigned int seed,  
                                     graph_access & G, 
                                     int* xadj,
                                     int* adjncy, 
                                     int* vwgt, 
                                     int* adjwgt,
                                     int* partition_map) {

        std::cout <<  "not implemented yet"  << std::endl;

}

void bipartition::post_fm(const PartitionConfig & config, graph_access & G) {
                refinement* refine          = new quotient_graph_refinement();
                complete_boundary* boundary = new complete_boundary(&G);
                boundary->build();

                PartitionConfig initial_cfg                 = config;
                initial_cfg.fm_search_limit                 = config.bipartition_post_fm_limits;
                initial_cfg.refinement_type                 = REFINEMENT_TYPE_FM;
                initial_cfg.refinement_scheduling_algorithm = REFINEMENT_SCHEDULING_ACTIVE_BLOCKS;
                initial_cfg.bank_account_factor             = 5;
                initial_cfg.rebalance                       = true;
                initial_cfg.softrebalance                   = true;
                initial_cfg.upper_bound_partition           = 100000000;
                initial_cfg.initial_bipartitioning          = true;
                refine->perform_refinement(initial_cfg, G, *boundary);

                delete refine;
                delete boundary;

}

NodeID bipartition::find_start_node( const PartitionConfig & config, graph_access & G) {
        NodeID startNode = random_functions::nextInt(0, G.number_of_nodes()-1);
        NodeID lastNode = startNode;

        int counter = G.number_of_nodes();
        while( G.getNodeDegree(startNode) == 0 && --counter > 0) {
                startNode = random_functions::nextInt(0, G.number_of_nodes()-1);
        }

        //now perform a bfs to get a partition
        for( unsigned i = 0; i < 3; i++) {
                std::vector<bool> touched(G.number_of_nodes(), false);
                startNode = lastNode;
                touched[startNode]  = true;

                std::queue<NodeID>* bfsqueue = new std::queue<NodeID>;
                bfsqueue->push(startNode);
                while(!bfsqueue->empty()) {
                        NodeID source = bfsqueue->front();
                        lastNode = source;
                        bfsqueue->pop();

                        forall_out_edges(G, e, source) {
                                NodeID target = G.getEdgeTarget(e);
                                if(!touched[target]) {
                                        touched[target] = true;
                                        bfsqueue->push(target);
                                }
                        } endfor
                }
                delete bfsqueue;

        }
        return lastNode;
}

void bipartition::grow_regions_bfs(const PartitionConfig & config, graph_access & G) {
        if(G.number_of_nodes() == 0) return;

        NodeID startNode = random_functions::nextInt(0, G.number_of_nodes()-1);
        if(config.buffoon) { startNode = find_start_node(config, G);  } // more likely to produce connected partitions

        std::vector<bool> touched(G.number_of_nodes(), false);
        touched[startNode]              = true;
        NodeWeight cur_partition_weight = 0;

        forall_nodes(G, node) {
                G.setPartitionIndex(node, 1);
        } endfor

        NodeID nodes_left = G.number_of_nodes()-1;

        //now perform a bfs to get a partition
        std::queue<NodeID>* bfsqueue = new std::queue<NodeID>;
        bfsqueue->push(startNode);
        for(;;) {
                if( nodes_left == 1 ) {
                        //only one node left --> we have to break
                        break;
                }

                if(bfsqueue->empty() && nodes_left > 0) {
                        //disconnected graph -> find a new start node among those that havent been touched
                        NodeID k = random_functions::nextInt(0, nodes_left-1);
                        NodeID start_node = 0;
                        forall_nodes(G, node) {
                                if(!touched[node]) {
                                        if(k == 0) {
                                                if( G.getNodeDegree(node) != 0) {
                                                        start_node = node;
                                                        break;
                                                } else {
                                                        G.setPartitionIndex(node, 0);
                                                        nodes_left--;
                                                        cur_partition_weight += G.getNodeWeight(node);
                                                        touched[node] = true;

                                                        if(cur_partition_weight >= (NodeWeight) config.grow_target) break;
                                                }
                                        } else {
                                                k--;
                                        }      
                                }
                        } endfor
                        
                        if(cur_partition_weight >= (NodeWeight) config.grow_target) break;
                        
                        bfsqueue->push(start_node);
                        touched[start_node] = true;
                } else if (bfsqueue->empty() && nodes_left == 0) {
                        break;
                }

                NodeID source = bfsqueue->front();
                bfsqueue->pop();
                G.setPartitionIndex(source, 0);

                nodes_left--;
                cur_partition_weight += G.getNodeWeight(source);

                if(cur_partition_weight >= (NodeWeight) config.grow_target) break;

                forall_out_edges(G, e, source) {
                        NodeID target = G.getEdgeTarget(e);
                        if(!touched[target]) {
                                touched[target] = true;
                                bfsqueue->push(target);
                        }
                } endfor
        }
        delete bfsqueue;
}


void bipartition::grow_regions_fm(const PartitionConfig & config, graph_access & G) {
        if(G.number_of_nodes() == 0) return;

        //NodeID startNode        = random_functions::nextInt(0, G.number_of_nodes()-1);
        NodeID startNode        = find_start_node(config, G);

        std::vector<bool> touched(G.number_of_nodes(), false);
        touched[startNode]              = true;
        NodeWeight cur_partition_weight = 0;

        forall_nodes(G, node) {
                G.setPartitionIndex(node, 1);
        } endfor

        NodeID nodes_left = G.number_of_nodes()-1;

        //now perform a pseudo dijkstra to get a partition
        maxNodeHeap* queue = new maxNodeHeap();
        queue->insert(startNode, 0); // in this case the gain doesn't really matter

        for(;;) {
                if( nodes_left == 1 ) {
                        //only one node left --> we have to break
                        break;
                }

                if(queue->empty() && nodes_left > 0) {
                        //disconnected graph -> find a new start node among those that havent been touched
                        NodeID k = random_functions::nextInt(0, nodes_left-1);
                        NodeID start_node = 0;
                        forall_nodes(G, node) {
                                if(!touched[node]) {
                                        if(k == 0) {
                                                start_node = node;
                                                break;
                                        } else {
                                                k--;
                                        }      
                                }
                        } endfor
                        
                        queue->insert(start_node, 0);
                        touched[start_node] = true;
                } else if (queue->empty() && nodes_left == 0) {
                        break;
                }

                NodeID source = queue->deleteMax();
                G.setPartitionIndex(source, 0);

                nodes_left--;
                cur_partition_weight += G.getNodeWeight(source);

                if(cur_partition_weight >= (NodeWeight)config.grow_target) break;

                forall_out_edges(G, e, source) {
                        NodeID target = G.getEdgeTarget(e);
                        if(G.getPartitionIndex(target) == 1) { //then we might need to update the gain!
                                Gain gain = compute_gain(G, target, 0);
                                touched[target] = true;

                                if(queue->contains(target)) {
                                        //change the gain
                                        queue->changeKey(target, gain);
                                } else {
                                        //insert
                                        queue->insert(target, gain);
                                }

                        }
                } endfor
        }
        delete queue;
}
