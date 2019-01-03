/******************************************************************************
 * dummy_operations.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <mpi.h>
#include <vector>
#include "dummy_operations.h"

dummy_operations::dummy_operations() {
                
}

dummy_operations::~dummy_operations() {
                        
}

void dummy_operations::run_collective_dummy_operations() {
        int rank, size;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);
        
        // Run Broadcast
        {
                int x;
                MPI_Comm_rank( MPI_COMM_WORLD, &x);
                MPI_Bcast(&x, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
        // Run Allgather.
        {
                int x, size;
                MPI_Comm_rank( MPI_COMM_WORLD, &x);
                MPI_Comm_size( MPI_COMM_WORLD, &size);
                
                std::vector<int> rcv(size);
                MPI_Allgather(&x, 1, MPI_INT, &rcv[0], 1, MPI_INT, MPI_COMM_WORLD);
        }

        // Run Allreduce.
        {
                int x;
                MPI_Comm_rank( MPI_COMM_WORLD, &x);
                
                int y = 0;
                MPI_Allreduce(&x, &y, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }

        // Dummy Prefix Sum
        {
                int x  = 1;
                int y  = 0;

                MPI_Scan(&x, &y, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
        }

        // Run Alltoallv.
        {
                std::vector<int> snd(size);
                std::vector<int> rcv(size);
                std::vector<int> scounts(size, 1);
                std::vector<int> rcounts(size, 1);
                std::vector<int> sdispls(size);
                std::vector<int> rdispls(size);
                for (int i = 0, iend = sdispls.size(); i < iend; ++i) {
                        sdispls[i] = rdispls[i] = i;
                }
                MPI_Alltoallv(&snd[0], &scounts[0], &sdispls[0], MPI_INT,
                              &rcv[0], &rcounts[0], &rdispls[0], MPI_INT, MPI_COMM_WORLD);
        }
        

}
