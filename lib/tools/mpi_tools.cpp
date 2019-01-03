/******************************************************************************
 * mpi_tools.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <mpi.h>
#include <unistd.h>

#include "mpi_tools.h"

mpi_tools::mpi_tools() {
                
}

mpi_tools::~mpi_tools() {
                

}

//void mpi_tools::non_active_wait_for_root() {
        //int rank, size;
        //MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        //MPI_Comm_size( MPI_COMM_WORLD, &size);
        
        //int MASTER = 0;

        //if(rank == MASTER) {
                ////wake up call
                //bool wakeup = true;
                //for( int to = 1; to < size; to++) {
                        //MPI_Send(&wakeup, 1, MPI_BOOL, to, 0, MPI_COMM_WORLD);
                //}
        //} else {
                ////non-busy waiting:
                //bool stop = false;
                //do {
                        //usleep(5000);
                        //stop = MPI::COMM_WORLD.Iprobe(MASTER,0);
                //} while(!stop);

                //bool wakeup = true;
                //MPI::COMM_WORLD.Recv(&wakeup, 1, MPI::BOOL, MASTER, 0);
        //}
//}

