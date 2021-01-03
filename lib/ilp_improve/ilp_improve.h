/******************************************************************************
 * ilp_improve.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#pragma once

#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <vector>
#include <algorithm>
#include <stdlib.h>

#include "gurobi_c++.h"

#include "balance_configuration.h"
#include "ilp_improve/ilp_helpers.h"
#include "partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "partition/coarsening/clustering/size_constraint_label_propagation.h"
#include "partition/graph_partitioner.h"

class ilp_improve {
public:
    size_t computeBFS(graph_access &G, std::unordered_set<NodeID> &nodesAvailable,
                    PartitionConfig partition_config, size_t limit_nonzeroes) {
        // BFS to get available nodes
        std::queue<std::vector<NodeID> > queue;
        ilp_helpers help;

        std::vector<std::pair<NodeID, Gain>> gains;

        switch(partition_config.ilp_mode) {
            case OptimizationMode::TREES :
                gains = help.bestStartNodes(G, nodesAvailable, queue);
                break;
            case OptimizationMode::GAIN :
                help.gainBFSStartNodes(G, nodesAvailable, queue,
                                  partition_config.ilp_min_gain);
                break;
            case OptimizationMode::BOUNDARY :
                help.cutBFSStartNodes(G, nodesAvailable, queue);
                break;
            case OptimizationMode::OVERLAP :
                //DO NOTHING!
                break;
        }

        std::vector<size_t> degrees(G.number_of_nodes(),0);

        size_t n = nodesAvailable.size() + G.get_partition_count();
        //we ignore edges between blocks but those can only be up to k^2
        //we also count edges twice between start vertices
        size_t m = 0;
        size_t k = G.get_partition_count();

        m = help.numEdgesBetweenNonStarters(G, nodesAvailable);

        for(auto q : nodesAvailable) {
            m += help.edgesInCoarse(G,nodesAvailable,q,degrees);
        }


        size_t trees = queue.size();



        if (partition_config.ilp_mode == OptimizationMode::TREES) {

            nodesAvailable.clear();

            n = G.get_partition_count();
            m = help.numEdgesBetweenNonStarters(G, nodesAvailable);


            Gain maxgain = gains[0].second;
            Gain curgain = maxgain;

            size_t tree_index = 0;
            size_t i = 0;

            //std::vector<size_t> degrees(G.number_of_nodes(),0);
            while(tree_index < gains.size()) {
                std::vector<NodeID> equalGain;

                while(gains[tree_index].second == curgain) {
                    equalGain.push_back(gains[tree_index].first);
                    tree_index++;
                }


                curgain++;

                std::random_device rd;
                std::mt19937 g(rd());
                std::shuffle(equalGain.begin(), equalGain.end(), g);

                for (NodeID el : equalGain) {
                    std::queue<std::vector<size_t> > bfs_queue;
                    bfs_queue.push(std::vector<size_t>(1,el));
                    std::unordered_set<NodeID> inMyTree;
                    nodesAvailable.emplace(el);

                    n++;
                    m+=help.edgesInCoarse(G,nodesAvailable,el,degrees);
                    i++;
                    if (help.numNonzeroes(n, m, k) > limit_nonzeroes) {
                        std::cout << "building trees - nonzeroes: " << help.numNonzeroes(n,m,k) << std::endl;
                        std::cout << "n: " << n << " m: " << m << " i: " << i << std::endl;
                        std::cout << " treeindex: " << tree_index << std::endl;
                        std::cout << "nodes available: " << nodesAvailable.size() << std::endl;
                        return i;
                    }

                    std::vector<bool> new_vertex(G.number_of_nodes(), false);
                    while(!bfs_queue.empty()) {
                        std::vector<size_t> curVec = bfs_queue.front();
                        bfs_queue.pop();
                        auto curNode = (NodeID) curVec.back();
                        if (curVec.size() < (size_t) partition_config.ilp_bfs_depth) {
                            if (new_vertex[curNode])
                                m -= degrees[curNode];

                            forall_out_edges(G, e, curNode)
                                    {
                                        NodeID newNode = G.getEdgeTarget(e);
                                        if (inMyTree.count(newNode) == 0) {
                                            inMyTree.emplace(newNode);
                                            bool inSet = (nodesAvailable.count(newNode) != 0);
                                            std::vector<bool> neighbors(G.get_partition_count(), false);
                                            std::vector<size_t> newVec = curVec;
                                            newVec.push_back(newNode);
                                            nodesAvailable.emplace(newNode);
                                            bfs_queue.emplace(newVec);


                                            if (!inSet) {
                                                new_vertex[newNode] = true;
                                                //guess for nonzeroes:
                                                n++;

                                                bool test = true;

                                                forall_out_edges(G, e2, newNode)
                                                        {

                                                            PartitionID p = G.getPartitionIndex(G.getEdgeTarget(e2));
                                                            if (nodesAvailable.count(G.getEdgeTarget(e2)) == 0 &&
                                                                neighbors[p] == false) {
                                                                test = false;

                                                                neighbors[p] = true;
                                                                m++;
                                                                degrees[newNode]++;
                                                            }

                                                            if (nodesAvailable.count(G.getEdgeTarget(e2)) != 0) {
                                                                m++;
                                                            }
                                                        }
                                                endfor

                                                if (test)
                                                    m--;
                                            }
                                        }
                                    }endfor

                            if (help.numNonzeroes(n, m, k) > limit_nonzeroes) {

                                std::cout << "building trees - nonzeroes: " << help.numNonzeroes(n,m,k) << std::endl;
                                std::cout << "n: " << n << " m: " << m << " i: " << i << std::endl;
                                return i;
                            }
                        }
                    }
                }

            }

        } else {

            std::vector<NodeID> shuffleq;

            while(!queue.empty()) {
                auto el = queue.front();
                queue.pop();
                shuffleq.push_back(el[0]);
            }


            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(shuffleq.begin(), shuffleq.end(), g);

            for (auto el : shuffleq) {
                queue.push(std::vector<NodeID>(1, el));
            }

            if (help.numNonzeroes(n,m,k) > limit_nonzeroes) {
                // just the boundaries have too much nonzeroes already. remove some vertices
                n=G.get_partition_count();

                nodesAvailable.clear();

                m=help.numEdgesBetweenNonStarters(G,nodesAvailable);

                while (!queue.empty()) {
                    NodeID curNode = queue.front()[0];
                    queue.pop();
                    nodesAvailable.emplace(curNode);
                    n++;

                    m += help.edgesInCoarse(G,nodesAvailable,curNode,degrees);

                    if (help.numNonzeroes(n,m,k) > limit_nonzeroes) {
                        std::cout << "early break - nonzeroes: " << help.numNonzeroes(n,m,k) << std::endl;
                        std::cout << "n: " << n << " m: " << m << std::endl;
                        return n;
                    }
                }
            }

            while (!queue.empty()) {
                std::vector<NodeID> curVec = queue.front();
                queue.pop();
                NodeID curNode = curVec.back();
                if (curVec.size() < (size_t) partition_config.ilp_bfs_depth) {
                    m -= degrees[curNode];
                    forall_out_edges(G, e, curNode){
                                if (nodesAvailable.count(G.getEdgeTarget(e)) == 0) {

                                    std::vector<bool> neighbors(G.get_partition_count(), false);

                                    std::vector<NodeID> newVec(curVec);
                                    NodeID newNode = G.getEdgeTarget(e);
                                    newVec.push_back(newNode);
                                    nodesAvailable.emplace(newNode);
                                    queue.push(newVec);

                                    //guess for nonzeroes:
                                    n++;
                                    bool test = true;

                                    forall_out_edges(G, e2, newNode){

                                                PartitionID p = G.getPartitionIndex(G.getEdgeTarget(e2));
                                                if (nodesAvailable.count(G.getEdgeTarget(e2)) == 0 &&
                                                    neighbors[p] == false) {
                                                    test = false;

                                                    neighbors[p] = true;
                                                    m++;
                                                    degrees[newNode]++;
                                                }

                                                if (nodesAvailable.count(G.getEdgeTarget(e2)) != 0) {
                                                    m++;
                                                }
                                            }
                                    endfor

                                    if (test)
                                        m--;
                                }
                            }endfor

                    if (help.numNonzeroes(n, m, k) > limit_nonzeroes)
                        break;
                }

            }
        }
        return trees;
    }

    PartitionID computeBlocks(graph_access &G,
                              std::unordered_set<NodeID> &nodesAvailable, size_t k) {
        // create blocks
        std::unordered_map<PartitionID, std::unordered_set<NodeID> > blocks;

        std::vector<size_t> not_found(k, 0);
        forall_nodes(G, node) {
            if (nodesAvailable.count(node) == 0) {
                blocks[G.getPartitionIndex(node)].emplace(node);
                not_found[G.getPartitionIndex(node)]++;
            }
        } endfor

        // compress block indexes if there are less than k blocks
        PartitionID idx = 0;
        if (blocks.size() < k) {
            for (auto pair: blocks) {

                for (NodeID node: pair.second) {
                    G.setPartitionIndex(node, idx);
                }

                idx++;
            }
        }

        return blocks.size();
    }

    void setPartitionIDs(graph_access &G,
                         std::unordered_set<NodeID> &nodesAvailable,
                         PartitionID numBlocks) {
        PartitionID count = numBlocks;
        for (NodeID n: nodesAvailable) {

            G.setPartitionIndex(n, count++);
        }
    }

    void transferIlp(graph_access &G, graph_access &coarser) {

        forall_nodes(G, n) {
            PartitionID new_id = G.getPartitionIndex(n);
            G.setPartitionIndex(n, coarser.getPartitionIndex(new_id));
        } endfor
    }

    void graphInfo (quality_metrics &qm, graph_access &G) {
        std::cout << "cut \t\t\t" << qm.edge_cut(G) << std::endl;
        std::cout << "bnd \t\t\t" << qm.boundary_nodes(G) << std::endl;
        std::cout << "balance \t\t" << qm.balance(G) << std::endl;
        std::cout << "max_comm_vol \t"
                  << qm.max_communication_volume(G) << std::endl;
    }

    void computeIlp(graph_access &coarser, PartitionConfig pc,
                    PartitionID numBlocks, double timelimit, std::vector<PartitionID> & pid,
                    std::vector<bool> & coarse_presets, bool finalrun = false) {

        double imbalance = pc.imbalance / 100 + 1;

        try {
        GRBEnv* env = new GRBEnv();
        GRBModel model = GRBModel(*env);


        NodeID numNodes = coarser.number_of_nodes();
        NodeID numEdges = coarser.number_of_edges() / 2;

        // create index for start edges of nodes
        std::unordered_map<EdgeID, NodeID> edge_source;
        forall_nodes(coarser, node) {
                    forall_out_edges(coarser, e, node) {
                                edge_source[e] = node;
                            } endfor
                } endfor

        double ceilNodes = imbalance * std::ceil((double) pc.largest_graph_weight / (double) pc.k);
        auto upper = (long) ceilNodes;

        GRBVar* nodes = new GRBVar [numNodes*pc.k];
        GRBVar* edges = new GRBVar [numEdges];

        model.set(GRB_StringAttr_ModelName, "Partition");
        model.set(GRB_DoubleParam_TimeLimit, timelimit);
        model.set(GRB_DoubleParam_MIPGap, 0);
#if (defined(MODE_KAFFPAE) || defined(MODE_DEVEL))
        if(!finalrun) {
            //in KaffpaE, we have multiple programs running on the same machine, thus we only use a single thread per program
            model.set(GRB_IntParam_Threads, 1);
            //output is also unfavourable if we run 50+ instances in parallel :)
            model.set(GRB_IntParam_LogToConsole, 0);
            model.set(GRB_IntParam_PoolSearchMode, 0);
            //model.set(GRB_IntParam_PoolSolutions, 3);
        }
#endif

        // Set decision variables for nodes
        for (size_t q = 0; q < pc.k; q++) {
            GRBLinExpr nodeTot = 0;
            for (NodeID i = 0; i < numNodes; i++) {
                if (i < numBlocks  //mode IS NOT Overlap
                    || coarse_presets[i] //mode IS overlap
                        ) {// constrain indexes of contracted nodes
                    if (coarser.getPartitionIndex(i) == q) {
                        nodes[i + q*numNodes] = model.addVar(1.0, 1.0, 0, GRB_BINARY);
                        nodes[i + q*numNodes].set(GRB_DoubleAttr_Start, 1.0);
                    } else {
                        nodes[i + q*numNodes] = model.addVar(0.0, 0.0, 0, GRB_BINARY);
                        nodes[i + q*numNodes].set(GRB_DoubleAttr_Start, 0.0);
                    }
                } else {
                    nodes[i + q*numNodes] = model.addVar(0.0, 1.0, 0, GRB_BINARY);
                    nodes[i + q*numNodes].set(GRB_DoubleAttr_Start, (coarser.getPartitionIndex(i) == q));
                }
                nodeTot += coarser.getNodeWeight(i) * nodes[i + q*numNodes];
            }
            // Add constraint: Balanced partition
            model.addConstr(nodeTot <= upper, "upper bound " + pc.k);
        }


        size_t j = 0;
        // Decision variables for edges
        GRBLinExpr edgeTot = 0;
        for (NodeID i = 0; i < (2 * numEdges); i++) {
            if (edge_source[i] > coarser.getEdgeTarget(i)) {
                edges[j] = model.addVar(0.0, 1.0, 0, GRB_BINARY);
                double start = (coarser.getPartitionIndex(edge_source[i]) != coarser.getPartitionIndex(coarser.getEdgeTarget(i))) ? 1 : 0;
                edges[j].set(GRB_DoubleAttr_Start, start);
                edgeTot += edges[j] * coarser.getEdgeWeight(i);
                for (size_t q = 0; q < pc.k; q++) {
                    GRBLinExpr cons = nodes[edge_source[i] + q * numNodes] - nodes[coarser.getEdgeTarget(i) + q * numNodes];
                    // Add constraint: valid partiton
                    model.addConstr(edges[j] >= cons, "valid partition");
                    model.addConstr(edges[j] >= -cons, "valid partition");
                }
                j++;
            }
        }

        // set model objective: minimize sum of weights of cut edges
        model.setObjective(edgeTot, GRB_MINIMIZE);

        // Add constraint: sum of all decision variables for 1 node is 1
        for (size_t i = 0; i < numNodes; i++) {
            GRBLinExpr sumCons = 0;
            for (size_t q = 0; q < pc.k; q++) {
                sumCons += nodes[i + q*numNodes];
            }
            model.addConstr(sumCons == 1);
        }

        // Optimize model
        model.optimize();

        // if solution is found
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL ||
            model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
            // set partition
            for (PartitionID q = 0; q < pc.k; q++) {
                for (size_t i = 0; i < numNodes; i++) {
                    auto val = (long) nodes[i + q*numNodes].get(GRB_DoubleAttr_X);
                    if (val == 1) coarser.setPartitionIndex((NodeID) i, q);
                }
            }
        } else {
            std::cout << "No solution" << std::endl;
        }

        delete[] nodes;
        delete[] edges;
        } catch(GRBException e) {
            std::cout << "Gurobi exception occurred:" << " ERROR " << e.getErrorCode() << ": " << e.getMessage() << std::endl;
	    exit (EXIT_FAILURE);
        }
    }

    std::vector<std::vector<NodeID>> createPartitions(graph_access & G, PartitionConfig pc) {
        std::vector<std::vector<NodeID>> partitions;
        G.set_partition_count(pc.k);

        balance_configuration bc;
        bc.configurate_balance(pc, G);

        graph_partitioner partitioner;


        for (int i = 0; i < pc.ilp_overlap_runs; ++i) {
            pc.graph_allready_partitioned = false;
            random_functions::setSeed(i);
            partitioner.perform_partitioning(pc, G);
            if( pc.kaffpa_perfectly_balance ) {
                double epsilon                         = pc.imbalance/100.0;
                pc.upper_bound_partition = (NodeWeight) ((1 + epsilon) * ceil(pc.largest_graph_weight / (double)pc.k));

                complete_boundary boundary(&G);
                boundary.build();

                cycle_refinement cr;
                cr.perform_refinement(pc, G, boundary);
            }

            std::vector<NodeID> partition;
            forall_nodes(G, n){
                        partition.push_back(G.getPartitionIndex(n));
                    }endfor

            partitions.push_back(partition);
        }
        return partitions;
    }

    std::tuple<size_t, size_t, size_t> buildOverlapGraph(graph_access & G, PartitionConfig pc, std::vector<NodeID> & best,
                                                         std::vector<std::vector<NodeID>> & partitions) {
        return buildOverlapGraph(G, pc, best, partitions, pc.ilp_overlap_runs);
    };

    std::tuple<size_t, size_t, size_t> buildOverlapGraph(graph_access & G, PartitionConfig pc, std::vector<NodeID> & best,
                                                         std::vector<std::vector<NodeID>> & partitions, int overlap_runs) {

        quality_metrics qm;

        size_t bestval = std::numeric_limits<size_t>::max();

        PartitionID  no_vertices = pc.k;

        std::vector<NodeID> part;

        size_t bestindex = 0;

        for (int i = 0; i < overlap_runs; ++i) {
            std::vector<NodeID> part_new = partitions[i];

            forall_nodes(G, n){
                        G.setPartitionIndex(n, part_new[n]);
                    } endfor

            auto val = (size_t) qm.edge_cut(G);

            if (val < bestval) {
                best = part_new;
                bestval = val;

                bestindex = i;
            }
            size_constraint_label_propagation sclp;
            
            if (i) {
                sclp.ensemble_two_clusterings(G, part, part_new, part, no_vertices);
            } else {
                part = part_new;
            }
        }

        forall_nodes(G, n) {
                    G.setPartitionIndex(n, part[n]);
                } endfor

        return std::make_tuple(no_vertices, bestval, bestindex);
    }

    std::vector<bool> findPresets(graph_access & G, PartitionConfig pc, std::vector<std::vector<NodeID>> & partitions,
                                  size_t best_index, graph_access & coarser) {

        std::vector<bool> presets(G.number_of_nodes(), false);

        if (pc.ilp_overlap_presets == OverlapPresets::NOEQUAL) {
            std::vector<bool> lockable_vertices(G.number_of_nodes(), true);
            size_t num_lockable = G.number_of_nodes();

            std::vector<bool> used_up(G.get_partition_count(), false);
            for (size_t i = 0; i < pc.k; ++i) {
                if (num_lockable > 0) {
                    size_t index;
                    do {
                        index = random_functions::nextInt(0, G.number_of_nodes() - 1);
                    } while (!lockable_vertices[index] || used_up[partitions[best_index][index]]);

                    used_up[partitions[best_index][index]] = true;

                    for (auto &partition : partitions) {
                        PartitionID part_index = partition[index];
                        forall_nodes(G, n)
                                {
                                    if (lockable_vertices[n] && partition[n] == part_index) {
                                        lockable_vertices[n] = false;
                                        --num_lockable;
                                    }
                                }
                        endfor
                    }
                    presets[index] = true;
                }
            }
        } else if (pc.ilp_overlap_presets == OverlapPresets::CENTER) {
            std::vector<std::queue<NodeID>> queue(pc.k);
            std::vector<bool> vtx_discovered(G.number_of_nodes(), false);
            std::vector<size_t> discovered(pc.k, 0);
            std::vector<size_t> block_size(pc.k, 0);

            forall_nodes(G, n)
                    {
                        PartitionID p = partitions[best_index][n];
                        ++block_size[p];
                        forall_out_edges(G, e, n)
                                {
                                    NodeID tgt = G.getEdgeTarget(e);


                                    if (partitions[best_index][tgt] != p) {
                                        vtx_discovered[n] = true;
                                        ++discovered[p];
                                        queue[p].push(n);
                                        break;
                                    }
                                } endfor

                    } endfor

            for(size_t p = 0; p < pc.k; ++p) {
                while(!queue[p].empty()) {
                    NodeID n = queue[p].front();
                    queue[p].pop();
                    forall_out_edges(G, e, n){
                                NodeID tgt = G.getEdgeTarget(e);
                                if (!vtx_discovered[tgt] && (discovered[p] < block_size[p])) {
                                    ++discovered[p];
                                    queue[p].push(tgt);
                                    vtx_discovered[tgt] = true;

                                    if (discovered[p] == block_size[p]) {
                                        presets[tgt] = true;
                                    }
                                }
                            }endfor
                }
            }

        } else if (pc.ilp_overlap_presets == OverlapPresets::RANDOM) {
            size_t set = 0;
            std::vector<bool> already_set(pc.k, false);

            while(set < pc.k) {
                int n = random_functions::nextInt(0, G.number_of_nodes() - 1);
                if (already_set[partitions[best_index][n]] == false) {
                    already_set[partitions[best_index][n]] = true;
                    set++;
                    presets[n] = true;
                }
            }
        } else if (pc.ilp_overlap_presets == OverlapPresets::HEAVIEST) {
            std::vector<std::pair<NodeID, NodeWeight>> nodeWeights;
            size_t set = 0;
            std::vector<bool> already_set(pc.k, false);

            for(size_t i = 0; i < coarser.number_of_nodes(); ++i) {
                nodeWeights.emplace_back(i, 0);
            }

            std::vector<int> reverse(coarser.number_of_nodes(), -1);

            forall_nodes(G, n) {
                        NodeID p = G.getPartitionIndex(n);
                        reverse[p] = n;
                        nodeWeights[p].second++;
                    } endfor


            std::sort(nodeWeights.begin(), nodeWeights.end(),
                      [](const std::pair<NodeID, NodeID>& n1,
                         const std::pair<NodeID, NodeID>& n2) {
                          return n1.second > n2.second;
                      });

            size_t i = 0;
            while (set < pc.k) {

                NodeID n = reverse[nodeWeights[i].first];

                i++;

                if (already_set[partitions[best_index][n]] == false) {
                    already_set[partitions[best_index][n]] = true;
                    set++;
                    presets[n] = true;
                }
            }
        }

        std::vector<bool> coarse_presets(coarser.number_of_nodes(), false);
        forall_nodes(G, node) {
                    NodeID c_node = G.getPartitionIndex(node);
                    if(presets[node]) {
                        coarse_presets[c_node] = true;
                    }
                } endfor

        return coarse_presets;
    }
};
