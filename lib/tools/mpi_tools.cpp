/******************************************************************************
 * mpi_tools.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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
#include <unistd.h>

#include "mpi_tools.h"

mpi_tools::mpi_tools() {
                
}

mpi_tools::~mpi_tools() {
                

}

void mpi_tools::non_active_wait_for_root() {
        int rank   = MPI::COMM_WORLD.Get_rank();
        int size   = MPI::COMM_WORLD.Get_size();
        int MASTER = 0;

        if(rank == MASTER) {
                //wake up call
                bool wakeup = true;
                for( int to = 1; to < size; to++) {
                        MPI::COMM_WORLD.Send(&wakeup, 1, MPI::BOOL, to, 0);
                }
        } else {
                //non-busy waiting:
                bool stop = false;
                do {
                        usleep(5000);
                        stop = MPI::COMM_WORLD.Iprobe(MASTER,0);
                } while(!stop);

                bool wakeup = true;
                MPI::COMM_WORLD.Recv(&wakeup, 1, MPI::BOOL, MASTER, 0);
        }
}

