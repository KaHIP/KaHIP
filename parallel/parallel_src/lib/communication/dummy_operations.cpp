/******************************************************************************
 * dummy_operations.cpp
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
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
                MPI_Allreduce(&x, &y, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
        }

        // Dummy Prefix Sum
        {
                int x  = 1;
                int y  = 0;

                MPI_Scan(&x, &y, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD); 
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
                MPI_Alltoallv(&snd[0], &scounts[0], &sdispls[0], MPI_INTEGER,
                              &rcv[0], &rcounts[0], &rdispls[0], MPI_INTEGER, MPI_COMM_WORLD);
        }
        

}
