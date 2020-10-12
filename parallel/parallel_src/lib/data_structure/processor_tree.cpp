#include "processor_tree.h"


/** @brief Default constructor. **/
processor_tree::processor_tree(){

};

/**	@brief Constructor to create a processor tree based on leaves.
	@param[in] distances[vector]: each vector element corresponds to the communication costs on each level from the leaves
	@param[in] descendants[vector]: each vector element corresponds to the number of descendants on each level
	@brief convention that ascending element positions (0,1,2, ...) corresponds to levels in the tree from higher to lower. 
	*/
processor_tree::processor_tree(const vector<int> &distances, const vector<int> &descendants ) {
  assert( distances.size() == descendants.size());
  traversalDistances = distances.size()==0 ? std::vector<int>{0} : distances;
  traversalDescendants = descendants.size()==0 ? std::vector<int>{1}  : descendants; //TODO: should this be 0 or 1
  numOfLevels = distances.size()==0 ? 1 : distances.size();
  for( unsigned int i = 0; i < get_numOfLevels(); i++)
    numPUs *= traversalDescendants[i];
  setCompactRep();
}

/* processor_tree::processor_tree(const vector<int> &distances, const vector<int> &descendants , MPI_Comm communicator ) { */
/*   assert( distances.size() == descendants.size()); */
/*   traversalDistances = distances.size()==0 ? std::vector<int>{0} : distances; */
/*   traversalDescendants = descendants.size()==0 ? std::vector<int>{1}  : descendants; //TODO: should this be 0 or 1 */
/*   numOfLevels = distances.size()==0 ? 1 : distances.size(); */
/*   for( unsigned int i = 0; i < get_numOfLevels(); i++) */
/*     numPUs *= traversalDescendants[i]; */
/*   buildProcGraph(communicator); */
/*   build_predecessorMatrix(); */
/* } */


void processor_tree::setCompactRep() {

        int groups_size = traversalDescendants.size();
	assert ( groups_size == numOfLevels);
	compact_bin_id = new std::vector<unsigned int> (numPUs,0);
	
	bit_sec_len = 1;
	for( unsigned k = 0; k < groups_size; k++) {
	  int tmp = ceil(log2(traversalDescendants[k]));
	  if (tmp > bit_sec_len) {
	    bit_sec_len = tmp;
	  }
	}

	for (unsigned i = 0; i < numPUs; i++) {
	  unsigned int lay_id = i;
	  for(int k=0; k < groups_size; k++) {
	    int remainder = lay_id % traversalDescendants[k];
	    lay_id = lay_id / traversalDescendants[k];
	    (*compact_bin_id)[i] += remainder << (k*bit_sec_len);
	  }
	}
}

void processor_tree::print() const {
        assert( traversalDistances.size()==get_numOfLevels() );
	assert( traversalDescendants.size()==get_numOfLevels() );
	std::cout << " ===== Printing Tree Information ===== " << std::endl;
	for( unsigned int i = 0; i < get_numOfLevels(); i++) {
	  std::cout << "Level ==" << i << "== Distance : "
		    <<   traversalDistances[i] << " Descedants : "
		    <<   traversalDescendants[i]  << std::endl;
	}
	std::cout << "Total number of processors = " << get_numPUs() << std::endl;
	std::cout << " ===================================== " << std::endl;
}

	
void processor_tree::print_allPairDistances() const {
        std::cout << " ========== Distance Matrix ==========" << std::endl;
	std::cout << " ===================================== " << std::endl;
	for( unsigned int i = 0; i < get_numPUs(); i++) {
	  for( unsigned int j = 0; j < get_numPUs(); j++) {
	    std::cout << getDistance_PxPy(i, j) << "  "; // << std::endl; //
	  }
	  std::cout  << std::endl;
	}
	std::cout << " ===================================== " << std::endl;
	//printDistance_PxPy(1,1);
}



 int processor_tree::printDistance_PxPy(int x, int y) const {
        assert((x <= numPUs) and (y <= numPUs) );
	int groups_size = traversalDescendants.size();
	assert ( groups_size == numOfLevels);
	std::vector<unsigned int>  * compact_bin_id = new std::vector<unsigned int> (numPUs,0);
	
	int bit_sec_len = 1;
	for( unsigned k = 0; k < groups_size; k++) {
	  int tmp = ceil(log2(traversalDescendants[k]));
	  if (tmp > bit_sec_len) {
	    bit_sec_len = tmp;
	  }
	}
	
	for (unsigned i = 0; i < numPUs; i++) {
	  unsigned int lay_id = i;
	  for(int k=0; k < groups_size; k++) {
	    int remainder = lay_id % traversalDescendants[k];
	    lay_id = lay_id / traversalDescendants[k];
	    (*compact_bin_id)[i] += remainder << (k*bit_sec_len);
	  }
	}
		
	for (unsigned i = 0; i < numPUs; i++) {
	  std::cout << (*compact_bin_id)[i] << ", "
		    << std::bitset<8>( (*compact_bin_id)[i])<< std::endl;
	}
	
	int k = 0;
	unsigned long long int xor_x_y = (*compact_bin_id)[x] ^ (*compact_bin_id)[y];
	if (!xor_x_y) return 0;
	int count_leading_zeros = __builtin_clzll(xor_x_y);
	int total_n_bits = 8*sizeof(unsigned long long int);
	int clz = total_n_bits - count_leading_zeros -1;
	if (clz >= 0) {
	  k = (int)floor(clz / bit_sec_len);
	  return traversalDistances[k];
	} else  {
	  return 0;
	}       		
 }	




// create predecessorMatrix (replicated in all procs).
// TODO: do not store whole matrix for all procs only the corresponding lines
void processor_tree::create_predecessorMatrix(parallel_graph_access & P, vector< vector<int>> & predecessorMatrix) const {
        int nodesNo = P.number_of_global_nodes();
	forall_local_nodes(P, i) {
	  for (int j = 0; j < nodesNo; j++) {
	    //forall_out_edges(P, e, j) {
	    bool done = false;
	    int start = i;
	    int target = j;
	    int tempsave;
	    while (start != target) {
	      if (start > target) {
		tempsave = start;
		start = P.getEdgeTarget(P.get_first_edge(start));
		if ((start == j) && !done) {
		  predecessorMatrix[i][j] = tempsave;
		}
	      } else {
		target = P.getEdgeTarget(P.get_first_edge(target));
		if (!done) {
		  predecessorMatrix[i][j] = P.getEdgeTarget(P.get_first_edge(j));
		  done = true;
		}
	      }
	    } // while
	  } // for
	} endfor  // for
}


/* void processor_tree::build_predecessorMatrix() const { */
/*   int nodesNo = proc.number_of_nodes();	   */
/*   for (int i = 0; i < nodesNo; i++) { */
/*     vector<int> row(nodesNo); */
/*     predecessorMatrix.push_back(row); */
/*   } */
/*   for (int i = 0; i < nodesNo; i++) { */
/*     for (int j = 0; j < nodesNo; j++) { */
/*       bool done = false; */
/*       int start = i; */
/*       int target = j; */
/*       int tempsave; */
/*       while (start != target) { */
/* 	if (start > target) { */
/* 	  tempsave = start;	   */
/* 	  start = proc.getEdgeTarget(proc.get_first_edge(start)); */
/* 	  if ((start == j) && !done) { */
/* 	    proc.predecessorMatrix[i][j] = tempsave; */
/* 	  } */
/* 	} else { */
/* 	  target = proc.getEdgeTarget(proc.get_first_edge(target)); */
/* 	  if (!done) { */
/* 	    proc.predecessorMatrix[i][j] = proc.getEdgeTarget(proc.get_first_edge(j)); */
/* 	    done = true; */
/* 	  } */
/* 	} */
/*       } // while */
/*     } */
/*   } */
/* } */



void processor_tree::print_predecessorMatrix(parallel_graph_access & P, vector< vector<int>> & predecessorMatrix) const {
        int nodesNo = P.number_of_global_nodes();
	std::cout << " PredecessorMatrix: "<< std::endl;
	for (int i = 0; i < nodesNo; i++) {
	  for (int j = 0; j < nodesNo; j++) {
	    std::cout << predecessorMatrix[i][j] << " ";
	  }
	  std::cout <<  std::endl;
	}
}
	  
// function to build the processor graph.
// Processor graph is not really constructed in parallel.
void processor_tree::create_procGraph(parallel_graph_access & cg, MPI_Comm communicator) const {


        vector<int> bandwidths(traversalDistances.size(),0);
	vector<int> extensions(traversalDescendants.size(),0);
  
	for (unsigned int i = 0; i < traversalDescendants.size(); i++) {    
	  bandwidths[traversalDistances.size() -1 -i] = traversalDistances[i];
	  extensions[traversalDescendants.size() -1 -i] = traversalDescendants[i];
	}
	
	  vector<int> vec(numOfLevels + 1, 1);
	  NodeID k = 0;
	  int j;
	  for (j = 0; j < numOfLevels; j++) {
	    vec[j+1] = extensions[j]*vec[j];
	    k += vec[j];
	  }
	  k += vec[j];
	  
	  EdgeID nmbEdges = (EdgeID)(k-1) * 2; // todo maybe *2

  	  int rank, comm_size;
  	  MPI_Comm_rank( communicator, &rank);
	  MPI_Comm_size( communicator, &comm_size);
	  
	  ULONG from  = 0;
	  ULONG to    = k - 1;
	  ULONG local_no_nodes = to - from + 1;


	  vector<int> nodes;
	  vector<int> nodeWeights;
	  vector<int> edges;
	  vector<int> edgeWeights;
	  
	  int sizeOfThisLevel = 1;
	  int startOfLastLevel = 0;
	  int startOfThisLevel = 0;
	  int startOfNextLevel = 1;
	  int numEdges = 0;
	  int numNodes = 1;
  
	  nodes.push_back(0);
	  nodeWeights.push_back(0);
  

	  for (unsigned int i = 0; i < extensions.size(); i++) {    
	    //add edges for all nodes in level i
	    for (int j = 0; j < sizeOfThisLevel; j++) {
	      //edge to parent
	      if (i > 0) {
		edges.push_back(startOfLastLevel + (j / extensions[i-1]));
		edgeWeights.push_back(bandwidths[i-1]);
	      }
	      
	      for (int l = 0; l < extensions[i]; l++) {
		edges.push_back(startOfNextLevel + l + j * extensions[i]);
		edgeWeights.push_back(bandwidths[i]);
		numEdges++;
	      }
	    }
	    //add new nodes for level i+1
	    for (int j = 0; j < (sizeOfThisLevel * extensions[i]) ; j++) {
	      if (i + 1 == extensions.size()) {
		nodes.push_back(numEdges);
		nodeWeights.push_back(1);
	      } else {
		nodes.push_back(numEdges + j * extensions[i + 1]);
		nodeWeights.push_back(0);
	      }
	      numNodes++;
	      numEdges++;
	      
	    }
	    startOfLastLevel = startOfThisLevel;
	    startOfThisLevel += sizeOfThisLevel;
	    sizeOfThisLevel *= extensions[i];
	    startOfNextLevel += sizeOfThisLevel;
	  }
	  
	  for (int j = 0; j < sizeOfThisLevel; j++) {
	    edges.push_back(startOfLastLevel + (j / extensions[extensions.size() - 1]));
	    edgeWeights.push_back(bandwidths[bandwidths.size() - 1]);
	  }
	  nodes.push_back(numEdges);

	  //std::cout <<  "numNodes: " << numNodes << " numEdges: " << numEdges << std::endl;
	  //std::cout <<  "k: " << k << " nmbEdges: " << nmbEdges << std::endl;
	  //std::cout <<  "proc R " << rank << " from: " << from << " to: " << to << std::endl;
	  
	  /* for (int i = 0; i < local_edge_lists.size(); i++) { */
	  /*   std::cout << "R:" << rank << " node: " << i << " "; */
	  /*   for (int j = 0; j < local_edge_lists[i].size(); j++) */
	  /*     std::cout <<  local_edge_lists[i][j] -1 << "  "; */
	  /*   std::cout << std::endl; */
	  /* } */
	  
	  cg.start_construction((NodeID) local_no_nodes, nmbEdges,(NodeID) k, nmbEdges);
	  cg.set_range(from, to);


	  std::vector< NodeID > vertex_dist( comm_size+1, 0 );
	  for( PEID peID = 0; peID <= comm_size; peID++) {
	    vertex_dist[peID] = peID * ceil(k / (double)comm_size); // from positions
	    //std::cout <<  " R " << rank << " proc vertex_dist[" << peID << "] = " << vertex_dist[peID] << std::endl;
	  }
	  cg.set_range_array(vertex_dist);

	  
	  for (NodeID i = 0; i < local_no_nodes; ++i) {
	    NodeID node = cg.new_node();
	    cg.setNodeWeight(node, 1);
	    cg.setNodeLabel(node, from+node);
	    cg.setSecondPartitionIndex(node, 0);
	    for(unsigned j = nodes[i]; j < nodes[i+1]; j++) {
	    	EdgeID e = cg.new_edge(node, edges[j]);
	    	cg.setEdgeWeight(e, edgeWeights[e]);
	    } // for
	  }
	  cg.finish_construction(); 
	  // TODO: DO I REALLY NEED IT HERE?
	  MPI_Barrier(communicator);

	  /* if (rank == ROOT) { */
	  /*   std::cout << " proc_graph nodes " << cg.number_of_global_nodes() */
	  /* 	      << " proc_graph edges " << cg.number_of_global_edges() << std::endl; */
	  /* } */

	    /* forall_local_nodes(cg,u) { */
	    /*   forall_out_edges(cg, edgeP, u) { */
	    /* 	unsigned int start = u; */
	    /* 	unsigned int target = cg.getEdgeTarget(edgeP); */
	    /* 	cout << "edgeP[" << start << "][" << target << "]: "  << */
	    /* 	  cg.getEdgeWeight(edgeP) << "\n"; */

	    /*   } endfor */
	    /* 	  } endfor */


	}



void processor_tree::create_commGraph(parallel_graph_access & C, parallel_graph_access & cg, MPI_Comm communicator) const {

	  int rank, comm_size;
  	  MPI_Comm_rank( communicator, &rank);
	  MPI_Comm_size( communicator, &comm_size);
	  unsigned k = numPUs; // number of nodes in comm graph
	  
	  // pe p reads the lines p*ceil(n/size) to (p+1)floor(n/size) lines of that file
	  ULONG from  = rank     * ceil(k / (double)comm_size);
	  ULONG to    = (rank+1) * ceil(k / (double)comm_size) - 1;
	  to = std::min<unsigned long>(to, k-1);
	  ULONG local_no_nodes = 0;
	  if (from <= to)
	    local_no_nodes = to - from + 1;
	  
	  unsigned ksq = k * k;
	  
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
	      if(uBlock != vBlock) {
		unsigned indexC = (unsigned) ((uBlock * k) + vBlock);
		if(edgeWeights[indexC] == 0) {
		  nmbEdges++;
		}
		(edgeWeights[indexC])++;
	      }
	    } endfor    
	  } endfor

	  EdgeWeight global_edgeWeights[ksq];
	  for(unsigned i = 0; i < ksq ; i ++)
	    global_edgeWeights[i] = 0;
	  for(unsigned i = 0; i < ksq ; i ++)
	    MPI_Allreduce(&edgeWeights[i], &global_edgeWeights[i], 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);
	  
	  if (rank == ROOT) {
	  for(PartitionID i = 0; i < ksq; i++) {
	    std::cout <<  " R " << rank << " : g_edgeWeights[" << i << "] = "
	  	      << global_edgeWeights[i] << std::endl;
	  }
	  }
	  
	  // TODO: Do I really need it? Since I have an allreduce command just above.. 
	  MPI_Barrier(communicator);
	  std::cout <<  " R " << rank << " local_no_nodes = "  << local_no_nodes  << "( from:" << from << ", to:" << to << ") nmbEdges = " << nmbEdges << " k = " << k << std::endl;
	  cg.start_construction((NodeID) local_no_nodes, nmbEdges,(NodeID) k, nmbEdges);
	  cg.set_range(from, to);
  
	  std::vector< NodeID > vertex_dist( comm_size+1, 0 );
	  for( PEID peID = 0; peID <= comm_size; peID++) {
	    vertex_dist[peID] = peID * ceil(k / (double)comm_size); // from positions
	    std::cout <<  " R " << rank << " comm vertex_dist[" << peID << "] = " << vertex_dist[peID] << std::endl;
	  }
	  cg.set_range_array(vertex_dist);

    
	  for (NodeID i = 0; i < local_no_nodes; ++i) {
	    NodeID node = cg.new_node();
	    cg.setNodeWeight(node, 1);
	    cg.setNodeLabel(node, from+node);
	    cg.setSecondPartitionIndex(node, 0);
	    for(unsigned j = 0; j < k; j ++) {
	      unsigned indexC = (((i + vertex_dist[rank]) * k) + j);
	      if(global_edgeWeights[indexC] != 0) {
		EdgeID e = cg.new_edge(i, j);
		cg.setEdgeWeight(e, global_edgeWeights[indexC]);
	      }
	    }
	  }

	  cg.finish_construction(); 
	  MPI_Barrier(communicator);
}




void processor_tree::create_parallelprocGraph(parallel_graph_access & cg, MPI_Comm communicator)
  const {

	  vector<int> vec(numOfLevels + 1, 1);
	  NodeID k = 0;
	  int j;
	  for (j = 0; j < numOfLevels; j++) {
	    vec[j+1] = traversalDescendants[j]*vec[j];
	    k += vec[j];
	  }
	  k += vec[j];
	  
	  EdgeID nmbEdges = (EdgeID)(k-1) * 2; // todo maybe *2

  	  int rank, comm_size;
  	  MPI_Comm_rank( communicator, &rank);
	  MPI_Comm_size( communicator, &comm_size);
	  
	  // pe p reads the lines p*ceil(n/size) to (p+1)floor(n/size) lines of that file
	  ULONG from  = rank     * ceil(k / (double)comm_size);
	  ULONG to    = (rank+1) * ceil(k / (double)comm_size) - 1;
	  to = std::min<unsigned long>(to, k-1);

	  ULONG local_no_nodes = to - from + 1;
	  
	  vector<int> nodes;
	  vector<int> nodeWeights;
	  vector<int> edges;
	  vector<int> edgeWeights;
	  
	  int sizeOfThisLevel = 1;
	  int startOfLastLevel = 0;
	  int startOfThisLevel = 0;
	  int startOfNextLevel = 1;
	  int numEdges = 0;
	  int numNodes = 1;
  
	  nodes.push_back(0);
	  nodeWeights.push_back(0);
  

	  for (unsigned int i = 0; i < traversalDescendants.size(); i++) {    
	    //add edges for all nodes in level i
	    for (int j = 0; j < sizeOfThisLevel; j++) {
	      //edge to parent
	      if (i > 0) {
		edges.push_back(startOfLastLevel + (j / traversalDescendants[i-1]));
		edgeWeights.push_back(traversalDistances[i-1]);
	      }
	      
	      for (int l = 0; l < traversalDescendants[i]; l++) {
		edges.push_back(startOfNextLevel + l + j * traversalDescendants[i]);
		edgeWeights.push_back(traversalDistances[i]);
		numEdges++;
	      }
	    }
	    //add new nodes for level i+1
	    for (int j = 0; j < (sizeOfThisLevel * traversalDescendants[i]) ; j++) {
	      if (i + 1 == traversalDescendants.size()) {
		nodes.push_back(numEdges);
		nodeWeights.push_back(1);
	      } else {
		nodes.push_back(numEdges + j * traversalDescendants[i + 1]);
		nodeWeights.push_back(0);
	      }
	      numNodes++;
	      numEdges++;
	      
	    }
	    startOfLastLevel = startOfThisLevel;
	    startOfThisLevel += sizeOfThisLevel;
	    sizeOfThisLevel *= traversalDescendants[i];
	    startOfNextLevel += sizeOfThisLevel;
	  }
	  
	  for (int j = 0; j < sizeOfThisLevel; j++) {
	    edges.push_back(startOfLastLevel + (j / traversalDescendants[traversalDescendants.size() - 1]));
	    edgeWeights.push_back(traversalDistances[traversalDistances.size() - 1]);
	  }
	  nodes.push_back(numEdges);

	  std::cout <<  "numNodes: " << numNodes << " numEdges: " << numEdges << std::endl;
	  std::cout <<  "k: " << k << " nmbEdges: " << nmbEdges << std::endl;

	  std::cout <<  " R " << rank << " from: " << from << " to: " << to << std::endl;
	  
	  MPI_Barrier(communicator);

	  /* for (NodeID i = 0; i < local_no_nodes; ++i) { */
	  /*   std::cout << "R:" << rank << " node: " << i << " "; */
	  /*   for(unsigned j = nodes[i]; j < nodes[i+1]; j++) { */
	  /*     std::cout <<  edges[j] << "  "; */
	  /*   } */
	  /*     std::cout << std::endl; */
	  /* } */
	  
	  cg.start_construction((NodeID) local_no_nodes, nmbEdges,(NodeID) k, nmbEdges);
	  cg.set_range(from, to);

	  std::vector< NodeID > vertex_dist( comm_size+1, 0 );
	  for( PEID peID = 0; peID <= comm_size; peID++) {
	    vertex_dist[peID] = peID * ceil(k / (double)comm_size); // from positions
	    //std::cout <<  " R " << rank << " proc_vertex_dist[" << peID << "] = " << vertex_dist[peID] << std::endl;
	  }
	  cg.set_range_array(vertex_dist);

	  
	  for (NodeID i = 0; i < local_no_nodes; ++i) {
	    NodeID node = cg.new_node();
	    cg.setNodeWeight(node, 1);
	    cg.setNodeLabel(node, from+node);
	    cg.setSecondPartitionIndex(node, 0);
	    for(unsigned j = nodes[i]; j < nodes[i+1]; j++) {
	    	EdgeID e = cg.new_edge(node, edges[j]);
	    	cg.setEdgeWeight(e, edgeWeights[e]);
	    } // for
	  }
	  cg.finish_construction(); 
	  MPI_Barrier(communicator);
	  /* if (rank == ROOT) { */
	  /*   std::cout << " proc_graph nodes " << cg.number_of_global_nodes() */
	  /* 	      << " proc_graph edges " << cg.number_of_global_edges() << std::endl; */
	  /* } */

	    /* forall_local_nodes(cg,u) { */
	    /*   forall_out_edges(cg, edgeP, u) { */
	    /* 	unsigned int start = u; */
	    /* 	unsigned int target = cg.getEdgeTarget(edgeP); */
	    /* 	cout << "edgeP[" << start << "][" << target << "]: "  << */
	    /* 	  cg.getEdgeWeight(edgeP) << "\n"; */

	    /*   } endfor */
	    /* 	  } endfor */


}
