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
	distributed_quality_metrics(EdgeWeight qap, EdgeWeight cut);
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

	void set_initial_qap(EdgeWeight qap) {initial_qap = qap;};
	void set_initial_cut(EdgeWeight cut) {initial_cut = cut;};
	void add_timing(std::vector<double> vec);
	EdgeWeight get_initial_qap() { return initial_qap; };
	EdgeWeight get_initial_cut() { return initial_cut; };
        std::vector< double > get_cycle_time() { return ml_time; };
        double get_coarse_time() { return ml_time[0]; };
	double get_inpart_time() { return ml_time[1]; };
	double get_refine_time() { return ml_time[2]; };
	void print();


/******************************************************/
/*                  evaluateMapping                   */
/******************************************************/
	void evaluateMapping(parallel_graph_access & C, const processor_tree & PEtree, MPI_Comm communicator) {
	  

		    

  	  int rank, comm_size;
  	  MPI_Comm_rank( communicator, &rank);
	  MPI_Comm_size( communicator, &comm_size);

	  
	  unsigned k = PEtree.get_numPUs();//number of nodes in C
	  
	  // pe p reads the lines p*ceil(n/size) to (p+1)floor(n/size) lines of that file
	  ULONG from  = rank     * ceil(k / (double)comm_size);
	  ULONG to    = (rank+1) * ceil(k / (double)comm_size) - 1;
	  to = std::min<unsigned long>(to, k-1);

	  ULONG local_no_nodes = to - from + 1;



	  unsigned ksq = k * k;
	  parallel_graph_access cg(communicator);
	  
	  //find number of edges in C and edge weights
	  EdgeWeight edgeWeights[ksq];
	  for(unsigned i = 0; i < ksq ; i ++) {
	    edgeWeights[i] = 0;
	  }
	  unsigned nmbEdges = 0;
	  forall_local_nodes(C,u) {
	    PartitionID uBlock = C.getNodeLabel(u);
	    forall_out_edges(C, e, u) {
	      NodeID v = C.getEdgeTarget(e);
	      PartitionID vBlock = C.getNodeLabel(v);
	      //cout << "uBlock and vBlock are " << uBlock << " and " << vBlock << endl;
	      if(uBlock != vBlock) {
		unsigned indexC = (unsigned) ((uBlock * k) + vBlock);
		if(edgeWeights[indexC] == 0) {
		  nmbEdges++;
		}
		(edgeWeights[indexC])++;
	      }
	    } endfor    
		} endfor

   for(PartitionID i = 0; i < ksq; i++) {
     std::cout <<  " PRINT " << rank << " : edgeWeights[" << i << "] = "
	       << edgeWeights[i] << std::endl;
   }

   MPI_Barrier(communicator);


  cg.start_construction((NodeID) local_no_nodes, nmbEdges,(NodeID) k, nmbEdges);
  cg.set_range(from, to);
  
  std::vector< NodeID > vertex_dist( comm_size+1, 0 );
  for( PEID peID = 0; peID <= comm_size; peID++) {
    vertex_dist[peID] = peID * ceil(k / (double)comm_size); // from positions
  }
  cg.set_range_array(vertex_dist);
  
  for (NodeID i = 0; i < local_no_nodes; ++i) {
    NodeID node = cg.new_node();
    cg.setNodeWeight(node, 1);
    cg.setNodeLabel(node, from+node);
    cg.setSecondPartitionIndex(node, 0);
    for(unsigned j = 0; j < k; j ++) {
      unsigned indexC = ((i * k) + j);
      //std::cout <<  " indexC " << indexC << std::endl;
      if(edgeWeights[indexC] != 0) {
	//cout << i << ", " << j<< "  \n";
	EdgeID e = cg.new_edge(i, j);
	//cg.setEdgeWeight(e, 1);
	cg.setEdgeWeight(e, edgeWeights[indexC]);
      }
    }
  }

  cg.finish_construction(); 
  MPI_Barrier(communicator);

  cout << "cg.number_of_global_nodes() = " << cg.number_of_global_nodes()
     << ", cg.number_of_local_nodes() = " << cg.number_of_global_edges() << "\n";


   /* forall_local_nodes(cg,i) { */
   /*   forall_out_edges(cg, edgeG, i) { */
   /*     unsigned int start = i; */
   /*     unsigned int target = cg.getEdgeTarget(edgeG); */
   /*     cout << "edge[" << start << "][" << target << "]: "  << */
   /* 	 cg.getEdgeWeight(edgeG) << "\n"; */

   /*   } endfor */
   /* 	 } endfor */

 
	  
  PartitionID n = C.number_of_global_nodes();
  cout << "C.number_of_global_nodes() = " << n << ", C.number_of_local_nodes() = " << C.number_of_local_nodes() << "\n";
  
  	  //compute maximum congestion, maximum dilation, average dilation
  	  int local_sumDilation = 0;
  	  int local_maxDilation = 0;
	  
  	  int np = PEtree.get_numPUs();
	  cout << "np = " << np  << "\n";
  	  vector<int> local_qap(np, 0);

	  vector<int> congestion(ksq, 0);
	  
	  int maxCongestion = 0;
	  int local_maxCongestion = 0;
  	  double avgDilation = 0.0;

	  // create and build processor parallel_graph_access object
	  parallel_graph_access P(communicator);
	  vector< vector<int>> predecessorMatrix;
	  PEtree.build_parallelPGraph(P, communicator);
	  predecessorMatrix = PEtree.build_predecessorMatrix(P);

	  forall_local_nodes(cg, i) {
  	    forall_out_edges(cg, edgeC, i) {
	      if( i < cg.getEdgeTarget(edgeC) ) {
  	  /* 	//only one edge direction considered */
  	  /* 	//find dilation of edgeC, update sumDilation, maxDilation */
  	  /* 	// the mapping is already included within C, don't need to map vertices */
  	  	unsigned int start = i;
  	  	unsigned int target = cg.getEdgeTarget(edgeC);
		if( cg.getNodeLabel( start ) != cg.getNodeLabel(target)) {
  	  	  int distance = PEtree.getDistance_PxPy(cg.getNodeLabel( start ),cg.getNodeLabel( target));
  	  	  int currDilation = distance * (cg.getEdgeWeight(edgeC));
  	  	  cout << "D[" << start << "][" << target << "]: "  << distance  << " -> currDilation"  << currDilation << "\n";
		  local_sumDilation += currDilation;
		  if(currDilation > local_maxDilation)
		    local_maxDilation = currDilation;
		  
		  //update congestion[]
		  unsigned int current = target;
		  unsigned int next = target;
		  while (current != start) {
		    current = predecessorMatrix[start][current];
		    if (next >= current) {
		      
		      /* for(PartitionID j = 0; j < np; j++) { */
		      /* 	if( j == next) { */
		      /* 	  PartitionID uBlock = C.getNodeLabel(current); */
		      /* 	  unsigned indexP = (unsigned) ((uBlock * k) + j); */
		      /* 	  congestion[indexP] += cg.getEdgeWeight(edgeC); */
		      /* 	} */
		      /* } */
		      forall_out_edges(P, edgeP, current) {
		      	if(P.getEdgeTarget(edgeP) == next) {
		      	  congestion[edgeP] += cg.getEdgeWeight(edgeC);
		      	}
		      } endfor
			  
		    } else {
		      /* for(PartitionID j = 0; j < np; j++) { */
		      /* 	if( j == current) { */
		      /* 	  PartitionID uBlock = C.getNodeLabel(next); */
		      /* 	  unsigned indexP = (unsigned) ((uBlock * k) + j); */
		      /* 	  congestion[indexP] += cg.getEdgeWeight(edgeC); */
		      /* 	} */
		      /* } */
		      
		      forall_out_edges(P, edgeP, next) {
		      	if (P.getEdgeTarget(edgeP) == current) {
		      	  congestion[edgeP] += cg.getEdgeWeight(edgeC);
		      	}
		      } endfor
		    }
		    next = current;
		  } // while
		  //update congestion[]		  
		}
	      }
  	    } endfor
        } endfor
		    

	 /* //update congestion[] */
	 /*  for(PartitionID j = 0; j < np; j++) { */
	 /*    for(PartitionID l = 0; l < np; l++) { */
	 /*      // TODO: cureful here division with 0!! */
	 /*      unsigned indexP = (unsigned) ((j * k) + l); */
	 /*      (congestion[indexP]) /= PEtree.getDistance_PxPy(j,l); */
	 /*    } */
	 /*  } */
	 forall_local_edges(P, edgeP) {
	    (congestion[edgeP]) /= P.getEdgeWeight(edgeP);//recall that edge weight indicates bandwidth
	  } endfor

	  /* for(PartitionID j = 0; j < ksq; j++) { */
	  /*   if (congestion[j] > local_maxCongestion) { */
	  /*     local_maxCongestion = congestion[j]; */
	  /*   } */
	  /* } */
	 forall_local_edges(P, edgeP) {
	    if (congestion[edgeP] > local_maxCongestion) {
	      local_maxCongestion = congestion[edgeP];
	    }
	  } endfor

	  MPI_Allreduce(&local_maxCongestion, &maxCongestion, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, communicator);	    
	  //update congestion[]
	  
  	  /* forall_local_nodes(C, i) { */
  	  /*   forall_out_edges(C, edgeC, i) { */
	  /*     //  	      if( i < C.getEdgeTarget(edgeC) ) { */
  	  /* /\* 	//only one edge direction considered *\/ */
  	  /* /\* 	//find dilation of edgeC, update sumDilation, maxDilation *\/ */
  	  /* /\* 	// the mapping is already included within C, don't need to map vertices *\/ */
  	  /* 	unsigned int start = i; */
  	  /* 	unsigned int target = C.getEdgeTarget(edgeC); */
	  /* 	if( C.getNodeLabel( start ) != C.getNodeLabel(target)) { */
  	  /* 	  //int distance = PEtree.getDistance_PxPy(C.getNodeLabel( start ),C.getNodeLabel( target)); */
  	  /* 	  //int currDilation = distance * (C.getEdgeWeight(edgeC)); */
  	  /* 	  //cout << "D[" << start << "][" << target << "]: "  << distance  << " -> currDilation"  << currDilation << "\n"; */
  	  /* 	  local_qap[C.getNodeLabel( target)] += C.getEdgeWeight(edgeC); */
	  /* 	} */
	  /* 	//} */
  	  /*   } endfor */
	  /* 	} endfor */
		    

	  /* MPI_Barrier(communicator); */
	  
  	  /* for(PartitionID i = 0; i < np; i++) { */
	  /*   std::cout <<  " PRINT " << rank << " : local_qap[" << i << "] = " << local_qap[i] << std::endl; */
	  /* } */
	  
  	  /* for(PartitionID i = 0; i < np; i++) { */
  	  /*   //if( i != rank) { */
	  /*     if( i < rank ) { */
  	  /*     int distance = PEtree.getDistance_PxPy(i,rank); */
  	  /*     int currDilation = distance * local_qap[i]; */
  	  /*     cout << "D[" << i << "][" << rank << "]: "  << distance  << ", C.getEdgeWeight = " << local_qap[i] << " -> currDilation = "  << currDilation << "\n"; */
  	  /*     local_sumDilation += currDilation; */
  	  /*     if(currDilation > local_maxDilation) { */
  	  /* 	local_maxDilation = currDilation; */
  	  /*     } */
  	  /*   } */
  	  /* } */
  //cout << rank << " local_maxDilation::" << local_maxDilation << "\n";
	  cout << rank << " local_sumDilation::" << local_sumDilation
	       << " num edges " << (cg.number_of_global_edges() / 2)<< "\n";

  	  int global_maxDilation = 0;
  	  MPI_Allreduce(&local_maxDilation, &global_maxDilation, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, communicator);
  	  int global_sumDilation = 0;
  	  MPI_Allreduce(&local_sumDilation, &global_sumDilation, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);


  if( rank == ROOT ) {
    cout << "max_congestion::" << maxCongestion << "\n";
    cout << "max_dilation::" << global_maxDilation << "\n";
    //cout << "Sum dilation: " << 2*sumDilation << "\n";
    cout << "sum_dilation::" <<  global_sumDilation << "\n";
    avgDilation = ((double) global_sumDilation) / ((double) (cg.number_of_global_edges() / 2));
    cout << "avg_dilation::" << avgDilation << "\n";
  }

}
	
	
private: 

	EdgeWeight initial_qap, initial_cut;
        std::vector< double > ml_time;


};



#endif /* end of include guard: DISTRIBUTED_QUALITY_METRICS_UAVSEXBT */
