/******************************************************************************
 * kaHIP_interface_mpi.cpp
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *****************************************************************************/

#include <fstream>
#include <iostream>

#include "kaHIP_interface_mpi.h"
#include "../lib/data_structure/graph_access.h"
#include "../lib/parallel_mh/parallel_mh_async.h"
#include "../lib/partition/partition_config.h"
#include "../lib/tools/quality_metrics.h"
#include "../lib/tools/random_functions.h"
#include "../app/balance_configuration.h"
#include "../app/configuration.h"

static_assert(sizeof(kahip_idx) == sizeof(EdgeID),
              "kahip_idx and EdgeID size mismatch -- check KAHIP_64BIT vs MODE64BITEDGES");

void kaffpaE_comm(int* n, int* vwgt, kahip_idx* xadj,
                  kahip_idx* adjcwgt, kahip_idx* adjncy, int* nparts,
                  double* imbalance, bool suppress_output,
                  bool graph_partitioned, int time_limit, int seed,
                  int mode, bool single_individual,
                  MPI_Comm communicator,
                  kahip_idx* edgecut, double* balance, int* part) {

        std::streambuf* backup = std::cout.rdbuf();
        std::ofstream ofs;
        if(suppress_output) {
                ofs.open("/dev/null");
                std::cout.rdbuf(ofs.rdbuf());
        }

        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = *nparts;
        switch( mode ) {
                case FAST:         cfg.fast(partition_config);         break;
                case ECO:          cfg.eco(partition_config);          break;
                case STRONG:       cfg.strong(partition_config);       break;
                case ULTRA:        cfg.ultra(partition_config);        break;
                case FASTSOCIAL:   cfg.fastsocial(partition_config);   break;
                case ECOSOCIAL:    cfg.ecosocial(partition_config);    break;
                case STRONGSOCIAL: cfg.strongsocial(partition_config); break;
                default:           cfg.eco(partition_config);          break;
        }

        // as in app/kaffpaE.cpp
        partition_config.seed        = seed;
        partition_config.imbalance   = 100*(*imbalance);
        partition_config.time_limit  = time_limit;
        partition_config.kaffpaE     = true;
        partition_config.repetitions = 1;
        partition_config.kabapE      = partition_config.imbalance < 1;

        graph_access G;
        G.build_from_metis(*n, xadj, adjncy);
        G.set_partition_count(partition_config.k);
        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);
        if(vwgt != NULL) {
                forall_nodes(G, node) {
                        G.setNodeWeight(node, vwgt[node]);
                } endfor
        }
        if(adjcwgt != NULL) {
                forall_edges(G, e) {
                        G.setEdgeWeight(e, adjcwgt[e]);
                } endfor
        }

        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        if(graph_partitioned) {
                forall_nodes(G, node) {
                        G.setPartitionIndex(node, part[node]);
                } endfor
        }
        partition_config.graph_allready_partitioned  = graph_partitioned;
        partition_config.no_new_initial_partitioning = graph_partitioned;

        // own communicator: messages of the evolutionary algorithm never match messages of the caller
        MPI_Comm comm;
        MPI_Comm_dup(communicator, &comm);
        {
                parallel_mh_async mh(comm);
                mh.set_single_individual(single_individual);
                mh.perform_partitioning(partition_config, G);
        }
        MPI_Comm_free(&comm);

        forall_nodes(G, node) {
                part[node] = G.getPartitionIndex(node);
        } endfor

        quality_metrics qm;
        *edgecut = qm.edge_cut(G);
        *balance = qm.balance(G);

        if(suppress_output) {
                ofs.close();
                std::cout.rdbuf(backup);
        }
}
