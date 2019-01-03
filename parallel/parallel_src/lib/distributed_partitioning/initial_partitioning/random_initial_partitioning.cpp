/******************************************************************************
 * random_initial_partitioning.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "random_initial_partitioning.h"
#include "data_structure/parallel_graph_access.h"
#include "tools/random_functions.h"
#include "tools/distributed_quality_metrics.h"

random_initial_partitioning::random_initial_partitioning() {
                
}

random_initial_partitioning::~random_initial_partitioning() {
                
}


void random_initial_partitioning::perform_partitioning( MPI_Comm communicator, PPartitionConfig & config, parallel_graph_access & G) {

        forall_local_nodes(G, node) {
                G.setNodeLabel(node, random_functions::nextInt(0ULL, config.k-1));
        } endfor
        
        G.update_ghost_node_data_global(); // exchange the labels of ghost nodes

        distributed_quality_metrics qm;
        EdgeWeight edgecut = qm.edge_cut(G, communicator );
        double balance = qm.balance(config, G, communicator );

        PEID rank;
        MPI_Comm_rank( communicator, &rank);
        
        if( rank == (int)ROOT) {
                std::cout <<  "log>initial edge edge cut " <<  edgecut  << std::endl;
                std::cout <<  "log>initial imbalance "     <<  balance << std::endl;
        }

}
