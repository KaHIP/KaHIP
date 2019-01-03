/******************************************************************************
 * mpi_tools.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/



#ifndef MPI_TOOLS_HMESDXF2
#define MPI_TOOLS_HMESDXF2

#include "data_structure/parallel_graph_access.h"
#include "partition_config.h"

class mpi_tools {
public:
        mpi_tools();
        virtual ~mpi_tools();

        void collect_and_write_labels( MPI_Comm communicator, PPartitionConfig & config, 
                                       parallel_graph_access & G);

        void collect_parallel_graph_to_local_graph( MPI_Comm communicator, 
                                                    PPartitionConfig & config, 
                                                    parallel_graph_access & G,
                                                    complete_graph_access & Q);

        // G is input (only on ROOT)
        // G is output (on every other PE)
        void distribute_local_graph( MPI_Comm communicator, PPartitionConfig & config, complete_graph_access & G);

        // alltoallv that can send more than int-count elements
        void alltoallv( void * sendbuf, 
                        ULONG sendcounts[], ULONG displs[], 
                        const MPI_Datatype & sendtype, void * recvbuf,
                        ULONG recvcounts[], ULONG rdispls[],
                        const MPI_Datatype & recvtype  ) {
                alltoallv( sendbuf, sendcounts, displs, sendtype, recvbuf, recvcounts, rdispls, recvtype, MPI_COMM_WORLD);
        };

        void alltoallv( void * sendbuf, 
                       ULONG sendcounts[], ULONG displs[], 
                       const MPI_Datatype & sendtype, void * recvbuf,
                       ULONG recvcounts[], ULONG rdispls[],
                       const MPI_Datatype & recvtype, MPI_Comm communicator  );


};


#endif /* end of include guard: MPI_TOOLS_HMESDXF2 */
