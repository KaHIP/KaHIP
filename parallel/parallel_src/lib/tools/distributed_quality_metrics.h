/******************************************************************************
 * distributed_quality_metrics.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef DISTRIBUTED_QUALITY_METRICS_UAVSEXBT
#define DISTRIBUTED_QUALITY_METRICS_UAVSEXBT

#include "pdefinitions.h"
#include "data_structure/parallel_graph_access.h"
#include "data_structure/processor_tree.h"
#include "ppartition_config.h"

class distributed_quality_metrics {
public:
        distributed_quality_metrics();
        virtual ~distributed_quality_metrics();

        EdgeWeight local_edge_cut( parallel_graph_access & G, int * partition_map, MPI_Comm communicator );
        EdgeWeight edge_cut( parallel_graph_access & G, MPI_Comm communicator );
        EdgeWeight edge_cut_second( parallel_graph_access & G, MPI_Comm communicator  );
        NodeWeight local_max_block_weight( PPartitionConfig & config, parallel_graph_access & G, int * partition_map, MPI_Comm communicator  );
        double balance( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator  );
        double balance_load( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator  );
        double balance_load_dist( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator  );
        double balance_second( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator  );
	EdgeWeight comm_vol( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator  );
	EdgeWeight comm_bnd( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator );
	EdgeWeight comm_vol_dist( parallel_graph_access & G, MPI_Comm communicator );

	EdgeWeight total_qap( parallel_graph_access & G, const processor_tree &PEtree, MPI_Comm communicator);


	
	//private: 

	EdgeWeight initialobjective;
	double blc;
        int m_num_levels;
        timer t_c;
	timer t_i;
	timer t_p;
	timer t_r;



};



#endif /* end of include guard: DISTRIBUTED_QUALITY_METRICS_UAVSEXBT */
