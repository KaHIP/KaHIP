/******************************************************************************
 * parhip_test.cpp
 *
 * Example of how to use the ParHIP parallel partitioning interface.
 * Each MPI process reads the METIS graph file but only stores its
 * own portion of the distributed CSR structure in memory.
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
        // Step 1: Read header to get n and fmt, then compute vtxdist.
        // -----------------------------------------------------------------
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << filename << std::endl;
                MPI_Finalize();
                return 1;
        }

        std::string line;
        while (std::getline(in, line)) {
                if (line[0] != '%') break;
        }

        idxtype n_global;
        int fmt = 0;
        {
                idxtype m_global;
                std::stringstream ss(line);
                ss >> n_global >> m_global >> fmt;
        }

        bool read_ew = (fmt == 1 || fmt == 11);
        bool read_nw = (fmt == 10 || fmt == 11);

        // Compute which nodes this PE owns
        std::vector<idxtype> vtxdist(size + 1);
        for (int p = 0; p <= size; p++) {
                vtxdist[p] = (idxtype)p * n_global / size;
        }

        idxtype local_from = vtxdist[rank];
        idxtype local_to   = vtxdist[rank + 1];
        idxtype local_n    = local_to - local_from;

        // -----------------------------------------------------------------
        // Step 2: Scan through the file, only storing data for local nodes.
        // -----------------------------------------------------------------
        std::vector<idxtype> xadj(local_n + 1);
        std::vector<idxtype> adjncy;
        std::vector<idxtype> adjwgt;
        std::vector<idxtype> vwgt(local_n, 1);

        idxtype node_counter = 0;
        idxtype local_idx = 0;
        xadj[0] = 0;

        while (std::getline(in, line)) {
                if (line[0] == '%') continue;

                bool is_local = (node_counter >= local_from && node_counter < local_to);

                if (is_local) {
                        std::stringstream ss(line);

                        if (read_nw) {
                                idxtype w;
                                ss >> w;
                                vwgt[local_idx] = w;
                        }

                        idxtype target;
                        while (ss >> target) {
                                adjncy.push_back(target - 1);  // 0-indexed
                                if (read_ew) {
                                        idxtype ew;
                                        ss >> ew;
                                        adjwgt.push_back(ew);
                                }
                        }

                        local_idx++;
                        xadj[local_idx] = adjncy.size();
                }

                node_counter++;
                if (node_counter >= local_to) break;
        }
        in.close();

        // -----------------------------------------------------------------
        // Step 3: Call ParHIP to partition the graph.
        // -----------------------------------------------------------------
        double imbalance = 0.03;
        int edgecut = 0;
        std::vector<idxtype> part(local_n);

        ParHIPPartitionKWay(vtxdist.data(), xadj.data(), adjncy.data(),
                            read_nw ? vwgt.data() : NULL,
                            read_ew ? adjwgt.data() : NULL,
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
