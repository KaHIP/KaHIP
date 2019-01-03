/******************************************************************************
 * definitions.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef DEFINITIONS_H_CHRA
#define DEFINITIONS_H_CHRA

#include <limits>
#include <queue>
#include <vector>

#include "limits.h"
#include "macros_assertions.h"
#include "stdio.h"

// allows us to disable most of the output during partitioning
#ifndef NOOUTPUT 
        #define PRINT(x) x
#else
        #define PRINT(x) do {} while (false);
#endif

/**********************************************
 * Constants
 * ********************************************/
//Types needed for the parallel graph ds
//we use long since we want to partition huge graphs
typedef unsigned long long ULONG;
typedef unsigned int UINT;
typedef unsigned long long NodeID;
typedef unsigned long long EdgeID;
typedef unsigned long long PartitionID;
typedef unsigned long long NodeWeight;
typedef unsigned long long EdgeWeight;
typedef int PEID; 

const PEID ROOT = 0;

typedef enum {
        PERMUTATION_QUALITY_NONE, 
	PERMUTATION_QUALITY_FAST,  
	PERMUTATION_QUALITY_GOOD
} PermutationQuality;

typedef enum {
        KAFFPAESTRONG,
        KAFFPAEECO,
        KAFFPAEFAST,
        KAFFPAEULTRAFASTSNW,
        KAFFPAEFASTSNW,
        KAFFPAEECOSNW,
        KAFFPAESTRONGSNW,
        RANDOMIP
} InitialPartitioningAlgorithm;

struct source_target_pair {
        NodeID source;
        NodeID target;
};

typedef enum {
        RANDOM_NODEORDERING, 
        DEGREE_NODEORDERING,
	LEASTGHOSTNODESFIRST_DEGREE_NODEODERING,
	DEGREE_LEASTGHOSTNODESFIRST_NODEODERING
} NodeOrderingType;


#endif

  //Tag Listing of Isend Operations(they should be unique per level) 
  //**************************************************************************
  //rank +   size                projection algorithm
  //rank + 2*size                projection algorithm
  //rank + 3*size                update labels global
  //rank + 4*size                contraction algorithm / label mapping
  //rank + 5*size                 --  ""  --
  //rank + 6*size                contracion algorithm / get nodes to cnodes
  //rank + 7*size                redist hashed graph
  //rank + 8*size                redist hashed graph
  //rank + 9*size                communicate node weights
  //rank + 10*size               down propagation
  //rank + 11*size               down propagation
  //rank + 12*size               MPI Tools
  //rank + 13*size               MPI Tools
  //rank + 100*size + x          Label Isends  
 


