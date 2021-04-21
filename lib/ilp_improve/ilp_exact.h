/******************************************************************************
 * ilp_exact.h 
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

class ilp_exact {
public:
    void computeIlp(graph_access &G, PartitionConfig & config) {

        double imbalance = config.imbalance / 100 + 1;

        try {
        GRBEnv* env = new GRBEnv();
        GRBModel model = GRBModel(*env);


        //NodeID numNodes = coarser.number_of_nodes();
        //NodeID numEdges = coarser.number_of_edges() / 2;

        //// create index for start edges of nodes
        //std::unordered_map<EdgeID, NodeID> edge_source;
        //forall_nodes(G, node) {
                    //forall_out_edges(G, e, node) {
                                //edge_source[e] = node;
                    //} endfor
        //} endfor

        double ceilNodes = imbalance * std::ceil((double) config.largest_graph_weight / (double) config.k);
        auto upper = (long) ceilNodes;

        GRBVar* nodes = new GRBVar [G.number_of_nodes()*config.k];
        GRBVar* edges = new GRBVar [G.number_of_edges()];

        model.set(GRB_StringAttr_ModelName, "Partition");
        model.set(GRB_DoubleParam_TimeLimit, config.ilp_timeout);
        model.set(GRB_DoubleParam_MIPGap, 0);


        //// Set decision variables for nodes
        for (size_t q = 0; q < config.k; q++) {
            GRBLinExpr nodeTot = 0;
            for (NodeID i = 0; i < G.number_of_nodes(); i++) {
                    nodes[i + q*G.number_of_nodes()] = model.addVar(0.0, 1.0, 0, GRB_BINARY);
                    nodes[i + q*G.number_of_nodes()].set(GRB_DoubleAttr_Start, 0.0);
                    nodeTot += G.getNodeWeight(i) * nodes[i + q*G.number_of_nodes()];
            }
            //// Add constraint: Balanced partition
            model.addConstr(nodeTot <= upper, "upper bound " + config.k);
            model.addConstr(nodeTot >= 1, "lower bound " + config.k);
        }

        //// Decision variables for edges
        GRBLinExpr edgeTot = 0;
        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        edges[e] = model.addVar(0.0, 1.0, 0, GRB_BINARY);
                        edges[e].set(GRB_DoubleAttr_Start, 0);

                        edgeTot += edges[e] * G.getEdgeWeight(e);
                        for (size_t q = 0; q < config.k; q++) {
                                GRBLinExpr cons = nodes[node + q * G.number_of_nodes()] - nodes[target + q * G.number_of_nodes()];
                                // Add constraint: valid partiton
                                model.addConstr(edges[e] >= cons, "valid partition");
                                model.addConstr(edges[e] >= -cons, "valid partition");
                        }
                } endfor
        } endfor
        
        //// set model objective: minimize sum of weights of cut edges
        model.setObjective(edgeTot/2, GRB_MINIMIZE);

        //// Add constraint: sum of all decision variables for 1 node is 1
        for (size_t i = 0; i < G.number_of_nodes(); i++) {
            GRBLinExpr sumCons = 0;
            for (size_t q = 0; q < config.k; q++) {
                sumCons += nodes[i + q*G.number_of_nodes()];
            }
            model.addConstr(sumCons == 1);
        }

        //// Optimize model
        model.optimize();

        //// if solution is found
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL ||
            model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
            // set partition
            for (PartitionID q = 0; q < config.k; q++) {
                for (size_t i = 0; i < G.number_of_nodes(); i++) {
                    auto val = (long) nodes[i + q*G.number_of_nodes()].get(GRB_DoubleAttr_X);
                    if (val == 1) G.setPartitionIndex((NodeID) i, q);
                }
            }
        } else {
            std::cout << "No solution" << std::endl;

            for (PartitionID q = 0; q < config.k; q++) {
                for (size_t i = 0; i < G.number_of_nodes(); i++) {
                    auto val = (long) nodes[i + q*G.number_of_nodes()].get(GRB_DoubleAttr_X);
                    if (val == 1) G.setPartitionIndex((NodeID) i, q);
                }
            }

        }

        delete[] nodes;
        delete[] edges;
        } catch(GRBException e) {
            std::cout << "Gurobi exception occurred:" << " ERROR " << e.getErrorCode() << ": " << e.getMessage() << std::endl;
	    exit (EXIT_FAILURE);
        }
    }

};
