/******************************************************************************
 * dspac.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Author: Daniel Seemaier <daniel.seemaier@student.kit.edu>
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <mpi.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>

#include "communication/mpi_tools.h"
#include "communication/dummy_operations.h"
#include "data_structure/parallel_graph_access.h"
#include "distributed_partitioning/distributed_partitioner.h"
#include "io/parallel_graph_io.h"
#include "io/parallel_vector_io.h"
#include "macros_assertions.h"
#include "parse_dspac_parameters.h"
#include "partition_config.h"
#include "random_functions.h"
#include "timer.h"
#include "tools/distributed_quality_metrics.h"
#include "dspac/dspac.h"
#include "dspac/edge_balanced_graph_io.h"

static void executeParhip(parallel_graph_access &G, PPartitionConfig &partitionConfig);

int main(int argn, char **argv) {
    MPI_Init(&argn, &argv);

    int rank, size;
    MPI_Comm communicator = MPI_COMM_WORLD;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    PPartitionConfig partition_config;
    DspacConfig dspac_config;
    std::string graph_filename;
    std::string partition_filename;

    int ret_code = parse_dspac_parameters(argn, argv, partition_config, dspac_config, graph_filename, partition_filename);

    if (ret_code) {
        MPI_Finalize();
        return 0;
    }
    if (rank == ROOT) {
        std::cout << "graph: " << graph_filename << "\n"
                  << "infinity edge weight: " << dspac_config.infinity << "\n"
                  << "seed: " << partition_config.seed << "\n"
                  << "k: " << partition_config.k << "\n"
                  << "ncores: " << size << std::endl;
    }

    timer t;
    MPI_Barrier(MPI_COMM_WORLD);
    {
        t.restart();
        if (rank == ROOT) std::cout << "running collective dummy operations ";
        dummy_operations dop;
        dop.run_collective_dummy_operations();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == ROOT) {
        std::cout << "took " << t.elapsed() << std::endl;
    }

    // load input graph
    t.restart();
    parallel_graph_access input_graph(communicator);
    edge_balanced_graph_io::read_binary_graph_edge_balanced(input_graph, graph_filename, partition_config, rank, size);
    if (rank == ROOT) {
        std::cout << "input IO took " << t.elapsed() << "\n"
                  << "n(input): " << input_graph.number_of_global_nodes() << "\n"
                  << "m(input): " << input_graph.number_of_global_edges() << std::endl;
    }
    MPI_Barrier(communicator);

    // construct split graph
    t.restart();
    parallel_graph_access split_graph(communicator);
    dspac splitter(input_graph, communicator, dspac_config.infinity);
    splitter.construct(split_graph);
    if (rank == ROOT) {
        std::cout << "split graph construction took " << t.elapsed() << "\n"
                  << "n(split): " << split_graph.number_of_global_nodes() << "\n"
                  << "m(split): " << split_graph.number_of_global_edges() << std::endl;
    }

    // partition split graph
    t.restart();
    executeParhip(split_graph, partition_config);
    if (rank == ROOT) {
        std::cout << "parhip took " << t.elapsed() << std::endl;
    }

    // evaluate edge partition
    t.restart();
    splitter.fix_cut_dominant_edges(split_graph);
    std::vector<PartitionID> edge_partition = splitter.project_partition(split_graph);
    EdgeWeight vertex_cut = splitter.calculate_vertex_cut(partition_config.k, edge_partition);
    if (rank == ROOT) {
        std::cout << "evaluation took " << t.elapsed() << "\n"
                  << "vertex cut: " << vertex_cut << std::endl;
    }

    if( partition_config.save_partition ) {
            parallel_vector_io pvio;
            std::string filename = partition_filename.empty() ? "tmpedgepartition.txtp" : partition_filename;
            pvio.writePartitionSimpleParallel(split_graph, filename);
    }

    if( partition_config.save_partition_binary ) {
            parallel_vector_io pvio;
            std::string filename = partition_filename.empty() ? "tmpedgepartition.binp" : partition_filename;
            pvio.writePartitionBinaryParallelPosix(partition_config, split_graph, filename);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

static void executeParhip(parallel_graph_access &G, PPartitionConfig &partitionConfig) {
    timer t;
    int rank, size;
    MPI_Comm communicator = MPI_COMM_WORLD;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &size);

    MPI_Barrier(MPI_COMM_WORLD);

    if (communicator != MPI_COMM_NULL) {
        MPI_Comm_rank(communicator, &rank);
        MPI_Comm_size(communicator, &size);

        if (rank == ROOT) {
            PRINT(std::cout << "log> cluster coarsening factor is set to "
                            << partitionConfig.cluster_coarsening_factor << std::endl;)
        }

        partitionConfig.stop_factor /= partitionConfig.k;
        if (rank != 0) partitionConfig.seed = partitionConfig.seed * size + rank;
        srand(static_cast<unsigned int>(partitionConfig.seed));

        random_functions::setSeed(partitionConfig.seed);
        parallel_graph_access::set_comm_rounds(partitionConfig.comm_rounds / size);
        parallel_graph_access::set_comm_rounds_up(partitionConfig.comm_rounds / size);
        distributed_partitioner::generate_random_choices(partitionConfig);

        G.printMemoryUsage(std::cout);

        //compute some stats
        EdgeWeight interPEedges = 0;
        EdgeWeight localEdges = 0;
        forall_local_nodes(G, node)
        {
            forall_out_edges(G, e, node)
            {
                NodeID target = G.getEdgeTarget(e);
                if (!G.is_local_node(target)) {
                    interPEedges++;
                } else {
                    localEdges++;
                }
            }
            endfor
        }
        endfor

        EdgeWeight globalInterEdges = 0;
        EdgeWeight globalIntraEdges = 0;
        MPI_Reduce(&interPEedges, &globalInterEdges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, communicator);
        MPI_Reduce(&localEdges, &globalIntraEdges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, communicator);

        if (rank == ROOT) {
            std::cout << "log> ghost edges " << globalInterEdges / (double) G.number_of_global_edges() << std::endl;
            std::cout << "log> local edges " << globalIntraEdges / (double) G.number_of_global_edges() << std::endl;
        }

        t.restart();
        double epsilon = (partitionConfig.inbalance) / 100.0;
        if (partitionConfig.vertex_degree_weights) {
            throw std::logic_error("not allowed to overwrite vertex degrees");
        } else {
            partitionConfig.number_of_overall_nodes = G.number_of_global_nodes();
            partitionConfig.upper_bound_partition =
                    (1 + epsilon) * ceil(G.number_of_global_nodes() / (double) partitionConfig.k);
        }


        distributed_partitioner dpart;
        dpart.perform_partitioning(communicator, partitionConfig, G);

        MPI_Barrier(communicator);

        double running_time = t.elapsed();
        distributed_quality_metrics qm;
        EdgeWeight edge_cut = qm.edge_cut(G, communicator);
        double balance = qm.balance(partitionConfig, G, communicator);
        PRINT(double
        balance_load = qm.balance_load(partitionConfig, G, communicator);)
        PRINT(double
        balance_load_dist = qm.balance_load_dist(partitionConfig, G, communicator);)

        if (rank == ROOT) {
            std::cout << "log>" << "=====================================" << std::endl;
            std::cout << "log>" << "============AND WE R DONE============" << std::endl;
            std::cout << "log>" << "=====================================" << std::endl;
            std::cout << "log>total partitioning time elapsed " << running_time << std::endl;
            std::cout << "log>final edge cut " << edge_cut << std::endl;
            std::cout << "log>final balance " << balance << std::endl;
            PRINT(std::cout << "log>final balance load " << balance_load << std::endl;)
            PRINT(std::cout << "log>final balance load dist " << balance_load_dist << std::endl;)
        }
        PRINT(qm.comm_vol(partitionConfig, G, communicator);)
        PRINT(qm.comm_vol_dist(G, communicator);)
    }

    MPI_Status st;
    int flag;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &st);
    while (flag) {
        std::cout << "attention: still incoming messages! rank " << rank << " from " << st.MPI_SOURCE << std::endl;
        int message_length;
        MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
        MPI_Status rst;
        std::vector <NodeID> message;
        message.resize(message_length);
        MPI_Recv(&message[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, st.MPI_TAG, MPI_COMM_WORLD, &rst);
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &st);
    };
    MPI_Barrier(MPI_COMM_WORLD);
}
