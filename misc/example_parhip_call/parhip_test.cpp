/******************************************************************************
 * parhip_test.cpp
 *
 * Example of how to use the ParHIP parallel partitioning interface.
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#include <iostream>
#include <vector>
#include <mpi.h>

#include "parhip_interface.h"

// Small example graph (5 nodes, 12 directed edges = 6 undirected edges):
//
//   0 --- 1 --- 2
//   |     |     |
//   +---- 4 --- 3
//
// Adjacency (CSR, global):
//   node 0: neighbors {1, 4}
//   node 1: neighbors {0, 2, 4}
//   node 2: neighbors {1, 3}
//   node 3: neighbors {2, 4}
//   node 4: neighbors {0, 1, 3}

int main(int argc, char **argv) {
        MPI_Init(&argc, &argv);

        int rank, size;
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);

        // ---------------------------------------------------------------
        // vtxdist: distributed CSR vertex distribution array.
        //   vtxdist[p] is the global index of the first node on process p.
        //   vtxdist[size] is the total number of nodes.
        //   Process p owns nodes vtxdist[p] .. vtxdist[p+1]-1.
        //
        // Each process provides xadj/adjncy for its LOCAL nodes only,
        // but adjncy values are GLOBAL node indices.
        // ---------------------------------------------------------------

        int n_global = 5;

        // Distribute nodes across processes as evenly as possible.
        std::vector<idxtype> vtxdist(size + 1);
        for (int p = 0; p <= size; p++) {
                vtxdist[p] = (idxtype)p * n_global / size;
        }

        idxtype local_n = vtxdist[rank + 1] - vtxdist[rank];

        // Full graph adjacency (global CSR)
        idxtype full_xadj[]  = {0, 2, 5, 7, 9, 12};
        idxtype full_adjncy[] = {1,4,  0,2,4,  1,3,  2,4,  0,1,3};

        // Extract local portion of CSR
        idxtype global_from = vtxdist[rank];
        std::vector<idxtype> xadj(local_n + 1);
        idxtype edge_offset = full_xadj[global_from];
        for (idxtype i = 0; i <= local_n; i++) {
                xadj[i] = full_xadj[global_from + i] - edge_offset;
        }

        idxtype local_edges = xadj[local_n];
        std::vector<idxtype> adjncy(local_edges);
        for (idxtype i = 0; i < local_edges; i++) {
                adjncy[i] = full_adjncy[edge_offset + i];
        }

        // Partition into 2 blocks with 3% imbalance
        int nparts = 2;
        double imbalance = 0.03;
        int edgecut = 0;
        std::vector<idxtype> part(local_n);

        ParHIPPartitionKWay(vtxdist.data(), xadj.data(), adjncy.data(),
                            NULL,  // vwgt: NULL means unit weights
                            NULL,  // adjwgt: NULL means unit weights
                            &nparts, &imbalance,
                            true,  // suppress_output
                            42,    // seed
                            FASTSOCIAL,
                            &edgecut, part.data(), &comm);

        if (rank == 0) {
                std::cout << "ParHIP edge cut: " << edgecut << std::endl;
        }

        // Gather and print the full partition vector on rank 0
        std::vector<int> recvcounts(size), displs(size);
        for (int p = 0; p < size; p++) {
                recvcounts[p] = (int)(vtxdist[p + 1] - vtxdist[p]);
                displs[p] = (int)vtxdist[p];
        }

        std::vector<idxtype> full_part;
        if (rank == 0) full_part.resize(n_global);

        // idxtype is unsigned long long, use MPI_UNSIGNED_LONG_LONG
        MPI_Gatherv(part.data(), (int)local_n, MPI_UNSIGNED_LONG_LONG,
                     full_part.data(), recvcounts.data(), displs.data(),
                     MPI_UNSIGNED_LONG_LONG, 0, comm);

        if (rank == 0) {
                std::cout << "Partition: ";
                for (int i = 0; i < n_global; i++) {
                        std::cout << full_part[i] << " ";
                }
                std::cout << std::endl;
        }

        MPI_Finalize();
        return 0;
}
