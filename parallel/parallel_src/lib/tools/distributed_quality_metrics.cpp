/******************************************************************************
 * distributed_quality_metrics.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>
#include <numeric>
#include "distributed_quality_metrics.h"

distributed_quality_metrics::distributed_quality_metrics() : initial_qap(0), initial_cut(0), ml_time(3,0.0){  
}

distributed_quality_metrics::distributed_quality_metrics(EdgeWeight qap, EdgeWeight cut) : initial_qap(qap), initial_cut(cut), ml_time(3,0.0){
}


distributed_quality_metrics::~distributed_quality_metrics() {
                        
}

EdgeWeight distributed_quality_metrics::edge_cut_second( parallel_graph_access & G, MPI_Comm communicator ) {
        EdgeWeight local_cut = 0;
        forall_local_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( G.getSecondPartitionIndex( node ) != G.getSecondPartitionIndex(target)) {
                                local_cut += G.getEdgeWeight(e);
                        }
                } endfor
        } endfor

        EdgeWeight global_cut = 0;
        MPI_Allreduce(&local_cut, &global_cut, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        return global_cut/2;
}

EdgeWeight distributed_quality_metrics::local_edge_cut( parallel_graph_access & G, int* partition_map, MPI_Comm communicator ) {
        EdgeWeight local_cut = 0;
        forall_local_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( partition_map[ node ] != partition_map[ target ]) {
                                local_cut += G.getEdgeWeight(e);
                        }
                } endfor
        } endfor

        return local_cut/2;
}

EdgeWeight distributed_quality_metrics::edge_cut( parallel_graph_access & G, MPI_Comm communicator ) {
        EdgeWeight local_cut = 0;
        forall_local_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( G.getNodeLabel( node ) != G.getNodeLabel(target)) {
                                local_cut += G.getEdgeWeight(e);
                        }
                } endfor
        } endfor

        EdgeWeight global_cut = 0;
        MPI_Allreduce(&local_cut, &global_cut, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        return global_cut/2;
}


EdgeWeight distributed_quality_metrics::total_qap( parallel_graph_access & G, const processor_tree &PEtree, MPI_Comm communicator ) {
        EdgeWeight local_qap = 0;
	// int clz = __builtin_clzll(G.number_of_global_nodes());
	// const int label_size = 8*sizeof(unsigned long long int) - clz;

    //probably the tree is not initialized
    if( PEtree.get_numPUs()<=1){
        return 0;
    }

	/* JUST FOR TESTING */
	int rank;
        MPI_Comm_rank( communicator, &rank);
	// if( rank == ROOT ) 
	//   std::cout <<  " PRINT : label_size " << label_size << std::endl;

	
        forall_local_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( G.getNodeLabel( node ) != G.getNodeLabel(target)) {
			  // if( rank == ROOT ) {
			  //   std::cout << "(" << node << ", " << target << ") -> ( "
			  // 	      << G.getNodeLabel( node ) << ", " << G.getNodeLabel( target )
			  // 	      << ")" << std::endl;
			  // }

			  local_qap += G.getEdgeWeight(e) * PEtree.getDistance_PxPy(G.getNodeLabel( node ),G.getNodeLabel( target));
                        }
                } endfor
        } endfor

        EdgeWeight global_qap = 0;
        MPI_Allreduce(&local_qap, &global_qap, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

	// if( rank == ROOT ) 
	//   std::cout <<  " PRINT : total_global_qap " << global_qap << std::endl;
	
        return global_qap;
}



void distributed_quality_metrics::add_timing(std::vector<double> vec) {
  assert(vec.size() == ml_time.size());
  for( int i = 0; i < vec.size(); i++) {
    ml_time[i] += vec[i];
  }
}

void distributed_quality_metrics::print() {
  std::cout << "log>qm: initial_qap " <<  initial_qap << std::endl;
  std::cout << "log>qm: initial_cut " <<  initial_cut << std::endl;
  std::cout << "log>qm: coarse_time " <<  ml_time[0] << std::endl;
  std::cout << "log>qm: inpart_time " <<  ml_time[1] << std::endl;
  std::cout << "log>qm: refine_time " <<  ml_time[2] << std::endl;
}


NodeWeight distributed_quality_metrics::local_max_block_weight( PPartitionConfig & config, parallel_graph_access & G, int * partition_map, MPI_Comm communicator ) {
        std::vector<PartitionID> block_weights(config.k, 0);

        NodeWeight graph_vertex_weight = 0;

        forall_local_nodes(G, n) {
                PartitionID curPartition     = partition_map[n];
                block_weights[curPartition] += G.getNodeWeight(n);
                graph_vertex_weight   += G.getNodeWeight(n);
        } endfor

        NodeWeight cur_max = 0;

        for( PartitionID block = 0; block < config.k; block++) {
                NodeWeight cur_weight = block_weights[block];
                if (cur_weight > cur_max) {
                        cur_max = cur_weight;
                }
        }

        return cur_max;
}

double distributed_quality_metrics::balance( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator ) {
        std::vector<PartitionID> block_weights(config.k, 0);

        NodeWeight local_graph_vertex_weight = 0;

        forall_local_nodes(G, n) {
                PartitionID curPartition     = G.getNodeLabel(n);
                block_weights[curPartition] += G.getNodeWeight(n);
                local_graph_vertex_weight   += G.getNodeWeight(n);
        } endfor

        std::vector<PartitionID> overall_weights(config.k, 0);
        MPI_Allreduce(&block_weights[0], &overall_weights[0], config.k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        NodeWeight graph_vertex_weight = 0;
        MPI_Allreduce(&local_graph_vertex_weight, &graph_vertex_weight, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        double balance_part_weight = ceil(graph_vertex_weight / (double)config.k);
        double cur_max             = -1;

        for( PartitionID block = 0; block < config.k; block++) {
                double cur = overall_weights[block];
                if (cur > cur_max) {
                        cur_max = cur;
                }
        }

        double percentage = cur_max/balance_part_weight;
        return percentage;
}

double distributed_quality_metrics::balance_second( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator ) {
        std::vector<PartitionID> block_weights(config.k, 0);

        NodeWeight local_graph_vertex_weight = 0;

        forall_local_nodes(G, n) {
                PartitionID curPartition     = G.getSecondPartitionIndex(n);
                block_weights[curPartition] += G.getNodeWeight(n);
                local_graph_vertex_weight   += G.getNodeWeight(n);
        } endfor

        std::vector<PartitionID> overall_weights(config.k, 0);
        MPI_Allreduce(&block_weights[0], &overall_weights[0], config.k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        NodeWeight graph_vertex_weight = 0;
        MPI_Allreduce(&local_graph_vertex_weight, &graph_vertex_weight, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        double balance_part_weight = ceil(graph_vertex_weight / (double)config.k);
        double cur_max             = -1;

        for( PartitionID block = 0; block < config.k; block++) {
                double cur = overall_weights[block];
                if (cur > cur_max) {
                        cur_max = cur;
                }
        }

        double percentage = cur_max/balance_part_weight;
        return percentage;
}


double distributed_quality_metrics::balance_load( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator ) {
        std::vector<NodeWeight> block_weights(config.k, 0);

        NodeWeight local_weight = 0;
        forall_local_nodes(G, n) {
                PartitionID curPartition     = G.getNodeLabel(n);
                block_weights[curPartition] += G.getNodeWeight(n)+G.getNodeDegree(n);
                local_weight   += G.getNodeWeight(n)+G.getNodeDegree(n);
        } endfor

        std::vector<PartitionID> overall_weights(config.k, 0);
        MPI_Allreduce(&block_weights[0], &overall_weights[0], config.k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        NodeWeight total_weight = 0;
        MPI_Allreduce(&local_weight, &total_weight, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        double balance_part_weight = ceil(total_weight / (double)config.k);
        double cur_max             = -1;

        for( PartitionID block = 0; block < config.k; block++) {
                double cur = overall_weights[block];
                if (cur > cur_max) {
                        cur_max = cur;
                }
        }

        double percentage = cur_max/balance_part_weight;
        return percentage;

}

double distributed_quality_metrics::balance_load_dist( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator ) {

        int rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        std::vector<NodeWeight> block_weights(size, 0);
        
        NodeWeight local_weight = 0;
        forall_local_nodes(G, n) {
                block_weights[rank] += G.getNodeWeight(n)+G.getNodeDegree(n);
                local_weight   += G.getNodeWeight(n)+G.getNodeDegree(n);
        } endfor

        std::vector<PartitionID> overall_weights(size, 0);
        MPI_Allreduce(&block_weights[0], &overall_weights[0], size, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        NodeWeight total_weight = 0;
        MPI_Allreduce(&local_weight, &total_weight, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        double balance_part_weight = ceil(total_weight / (double)size);
        double cur_max             = -1;

        for( PartitionID block = 0; block < (PartitionID)size; block++) {
                double cur = overall_weights[block];
                if (cur > cur_max) {
                        cur_max = cur;
                }
        }

        double percentage = cur_max/balance_part_weight;
        return percentage;

}

// measure the communication volume of the current graph distribution
EdgeWeight distributed_quality_metrics::comm_vol( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator ) {
        EdgeWeight local_comm_vol = 0; int rank;
        MPI_Comm_rank( communicator, &rank);

        std::vector<PartitionID> block_volume(config.k, 0);
        forall_local_nodes(G, node) {
                std::vector<bool> block_incident(config.k, false);
                PartitionID block = G.getNodeLabel( node );
                block_incident[block]    = true;
                int num_incident_blocks = 0;
                
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        PartitionID target_block = G.getNodeLabel( target );
                        if(!block_incident[target_block]) {
                                block_incident[target_block] = true;
                                num_incident_blocks++;
                        }
                } endfor
                block_volume[block] += num_incident_blocks;
        } endfor

        std::vector<PartitionID> overall_weights(config.k, 0);
        MPI_Allreduce(&block_volume[0], &overall_weights[0], config.k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        if( rank == ROOT ) {
                EdgeWeight total_comm_vol = 0;
                for( PEID i = 0; i < (PEID)overall_weights.size(); i++) {
                        total_comm_vol += overall_weights[i];
                }
                EdgeWeight max_comm_vol = *(std::max_element(overall_weights.begin(), overall_weights.end()));
                EdgeWeight min_comm_vol = *(std::min_element(overall_weights.begin(), overall_weights.end()));

                std::cout <<  "log> total vol part " <<  total_comm_vol << std::endl;
                std::cout <<  "log> max vol part " <<  max_comm_vol << std::endl;
                std::cout <<  "log> min vol part " <<  min_comm_vol << std::endl;
                std::cout <<  "log> vol part ratio " <<  max_comm_vol/(double)min_comm_vol  << std::endl << std::endl;
        }

        return local_comm_vol;
}


EdgeWeight distributed_quality_metrics::comm_bnd( PPartitionConfig & config, parallel_graph_access & G, MPI_Comm communicator ) {
        EdgeWeight local_comm_vol = 0; int rank;
        MPI_Comm_rank( communicator, &rank);

        std::vector<EdgeWeight> num_inner_nodes(config.k, 0);
        std::vector<EdgeWeight> num_bnd_nodes(config.k, 0);
        
        forall_local_nodes(G, node) {                
                bool is_bnd = false;
                PartitionID my_block = G.getNodeLabel( node );
                
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        PartitionID target_block = G.getNodeLabel( target );
                        if( my_block!=target_block ) {     // found neighbor in different block
                            is_bnd = true;
                            break;
                        }
                } endfor

                if( is_bnd ){
                    num_bnd_nodes[my_block]++;
                }
                else{
                    num_inner_nodes[my_block]++;
                }
                
        } endfor

        std::vector<EdgeWeight> global_num_inner_nodes(config.k, 0);
        std::vector<EdgeWeight> global_num_bnd_nodes(config.k, 0);

        MPI_Allreduce(&num_inner_nodes[0], &global_num_inner_nodes[0], config.k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);
        MPI_Allreduce(&num_bnd_nodes[0], &global_num_bnd_nodes[0], config.k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);

        
        EdgeWeight total_bnd_nodes = std::accumulate( global_num_bnd_nodes.begin(), global_num_bnd_nodes.end(), 0);
        EdgeWeight total_inner_nodes = std::accumulate( global_num_inner_nodes.begin(), global_num_inner_nodes.end(), 0);
        
        EdgeWeight max_bnd_nodes = *std::max_element( global_num_bnd_nodes.begin(), global_num_bnd_nodes.end() );
        
        EdgeWeight total_nodes = total_bnd_nodes + total_inner_nodes;
        //TODO assertion to check if total sum of nodes is correct
                
        if( rank == ROOT ) {
                std::vector<double> percent_bnd_nodes(config.k, 0.0);    
                for( PEID i = 0; i < PEID( config.k); i++) {
                        percent_bnd_nodes[i] = (double( global_num_bnd_nodes[i]))/(global_num_bnd_nodes[i]+global_num_inner_nodes[i]);
                }
                
                assert( total_inner_nodes+total_bnd_nodes==total_nodes );
                std::cout <<  "log> total bnd nodes " <<  total_bnd_nodes << std::endl;
                std::cout <<  "log> max bnd nodes " <<  max_bnd_nodes << std::endl;
                std::cout <<  "log> max percentage of bnd nodes " << *std::max_element(percent_bnd_nodes.begin(), percent_bnd_nodes.end())  << std::endl;
                std::cout <<  "log> average percentage of bnd nodes " << std::accumulate(percent_bnd_nodes.begin(), percent_bnd_nodes.end(), 0.0)/(double(config.k))  << std::endl<< std::endl;
        }

        return local_comm_vol;
}

// measure the communication volume of the current graph distribution
EdgeWeight distributed_quality_metrics::comm_vol_dist( parallel_graph_access & G, MPI_Comm communicator ) {
        EdgeWeight local_comm_vol = 0;
        int rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);

        forall_local_nodes(G, node) {
                std::vector<bool> block_incident(size, false);
                block_incident[rank]    = true;
                int num_incident_blocks = 0;
                
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if( !G.is_local_node( target ) ) {
                                PartitionID target_block = G.getTargetPE( target );
                                if(!block_incident[target_block]) {
                                        block_incident[target_block] = true;
                                        num_incident_blocks++;
                                }
                        }
                } endfor
                local_comm_vol += num_incident_blocks;
        } endfor

        EdgeWeight total_comm_vol = 0;
        EdgeWeight max_comm_vol   = 0;
        EdgeWeight min_comm_vol   = 0;

        MPI_Reduce(&local_comm_vol, &total_comm_vol, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, communicator);
        MPI_Reduce(&local_comm_vol, &max_comm_vol, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, ROOT, communicator);
        MPI_Reduce(&local_comm_vol, &min_comm_vol, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, ROOT, communicator);

        if( rank == ROOT ) {
                std::cout <<  "log> total vol currentdist " <<  total_comm_vol  << std::endl;
                std::cout <<  "log> max vol currentdist " <<  max_comm_vol  << std::endl;
                std::cout <<  "log> min vol currentdist " <<  min_comm_vol  << std::endl;
                std::cout <<  "log> vol dist currentratio " <<  max_comm_vol/(double)min_comm_vol  << std::endl;
        }

        return local_comm_vol;
}

