/******************************************************************************
 * parhip_test.cpp
 *
 * Example of how to use the ParHIP parallel partitioning interface.
 * Each MPI process reads the full METIS graph file but only keeps
 * its own portion of the distributed CSR structure.
 *
 * Usage: mpirun -np <P> ./parhip_test <graph.metis> [<k>]
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <mpi.h>

#include "parhip_interface.h"

int main(int argc, char **argv) {
        MPI_Init(&argc, &argv);

        int rank, size;
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);

        if (argc < 2) {
                if (rank == 0) {
                        std::cerr << "Usage: mpirun -np <P> " << argv[0]
                                  << " <graph.metis> [<k>]" << std::endl;
                }
                MPI_Finalize();
                return 1;
        }

        std::string filename = argv[1];
        int nparts = (argc >= 3) ? atoi(argv[2]) : 2;

        // -----------------------------------------------------------------
        // Step 1: Every PE reads the full METIS graph file.
        //         METIS format: first non-comment line is "n m [fmt]"
        //         followed by n lines, one per node, listing neighbors
        //         (1-indexed). fmt=1 means edge weights, fmt=11 means
        //         node+edge weights, fmt=10 means node weights.
        // -----------------------------------------------------------------
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << filename << std::endl;
                MPI_Finalize();
                return 1;
        }

        std::string line;

        // skip comments
        while (std::getline(in, line)) {
                if (line[0] != '%') break;
        }

        idxtype n_global, m_global;
        int fmt = 0;
        {
                std::stringstream ss(line);
                ss >> n_global >> m_global >> fmt;
        }

        bool read_ew = (fmt == 1 || fmt == 11);
        bool read_nw = (fmt == 10 || fmt == 11);

        // Read full graph into CSR arrays (0-indexed)
        std::vector<idxtype> full_xadj(n_global + 1, 0);
        std::vector<idxtype> full_adjncy;
        std::vector<idxtype> full_adjwgt;
        std::vector<idxtype> full_vwgt(n_global, 1);

        idxtype node_counter = 0;
        while (std::getline(in, line)) {
                if (line[0] == '%') continue;
                if (node_counter >= n_global) break;

                std::stringstream ss(line);

                if (read_nw) {
                        idxtype w;
                        ss >> w;
                        full_vwgt[node_counter] = w;
                }

                idxtype target;
                while (ss >> target) {
                        full_adjncy.push_back(target - 1);  // convert to 0-indexed

                        if (read_ew) {
                                idxtype ew;
                                ss >> ew;
                                full_adjwgt.push_back(ew);
                        }
                }

                node_counter++;
                full_xadj[node_counter] = full_adjncy.size();
        }
        in.close();

        // -----------------------------------------------------------------
        // Step 2: Build the distributed CSR structure.
        //         vtxdist[p] = first global node index owned by process p.
        //         Each PE extracts its local portion of xadj/adjncy.
        // -----------------------------------------------------------------
        std::vector<idxtype> vtxdist(size + 1);
        for (int p = 0; p <= size; p++) {
                vtxdist[p] = (idxtype)p * n_global / size;
        }

        idxtype local_from = vtxdist[rank];
        idxtype local_to   = vtxdist[rank + 1];
        idxtype local_n    = local_to - local_from;

        // Extract local xadj (re-based to start at 0)
        std::vector<idxtype> xadj(local_n + 1);
        idxtype edge_offset = full_xadj[local_from];
        for (idxtype i = 0; i <= local_n; i++) {
                xadj[i] = full_xadj[local_from + i] - edge_offset;
        }

        // Extract local adjncy and adjwgt
        idxtype local_edges = xadj[local_n];
        std::vector<idxtype> adjncy(full_adjncy.begin() + edge_offset,
                                     full_adjncy.begin() + edge_offset + local_edges);

        idxtype* adjwgt_ptr = NULL;
        std::vector<idxtype> adjwgt;
        if (read_ew) {
                adjwgt.assign(full_adjwgt.begin() + edge_offset,
                              full_adjwgt.begin() + edge_offset + local_edges);
                adjwgt_ptr = adjwgt.data();
        }

        // Extract local vwgt
        std::vector<idxtype> vwgt(full_vwgt.begin() + local_from,
                                   full_vwgt.begin() + local_to);
        idxtype* vwgt_ptr = read_nw ? vwgt.data() : NULL;

        // -----------------------------------------------------------------
        // Step 3: Call ParHIP to partition the graph.
        // -----------------------------------------------------------------
        double imbalance = 0.03;
        int edgecut = 0;
        std::vector<idxtype> part(local_n);

        ParHIPPartitionKWay(vtxdist.data(), xadj.data(), adjncy.data(),
                            vwgt_ptr, adjwgt_ptr,
                            &nparts, &imbalance,
                            false,  // suppress_output
                            42,     // seed
                            FASTSOCIAL,
                            &edgecut, part.data(), &comm);

        if (rank == 0) {
                std::cout << "edge cut: " << edgecut << std::endl;
        }

        MPI_Finalize();
        return 0;
}
