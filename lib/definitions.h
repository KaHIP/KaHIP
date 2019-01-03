/******************************************************************************
 * definitions.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef DEFINITIONS_H_CHR
#define DEFINITIONS_H_CHR

#include <limits>
#include <queue>
#include <vector>

#include "limits.h"
#include "macros_assertions.h"
#include "stdio.h"


// allows us to disable most of the output during partitioning
#ifdef KAFFPAOUTPUT
        #define PRINT(x) x
#else
        #define PRINT(x) do {} while (false);
#endif

/**********************************************
 * Constants
 * ********************************************/
//Types needed for the graph ds
typedef unsigned int 	NodeID;
typedef double 		EdgeRatingType;
typedef unsigned int 	PathID;
typedef unsigned int 	PartitionID;
typedef unsigned int 	NodeWeight;
typedef int 		EdgeWeight;
typedef EdgeWeight 	Gain;
#ifdef MODE64BITEDGES
typedef uint64_t 	EdgeID;
#else
typedef unsigned int 	EdgeID;
#endif
typedef int 		Color;
typedef unsigned int 	Count;
typedef std::vector<NodeID> boundary_starting_nodes;
typedef long FlowType;

const EdgeID UNDEFINED_EDGE            = std::numeric_limits<EdgeID>::max();
const NodeID UNDEFINED_NODE            = std::numeric_limits<NodeID>::max();
const NodeID UNASSIGNED                = std::numeric_limits<NodeID>::max();
const NodeID ASSIGNED                  = std::numeric_limits<NodeID>::max()-1;
const PartitionID INVALID_PARTITION    = std::numeric_limits<PartitionID>::max();
const PartitionID BOUNDARY_STRIPE_NODE = std::numeric_limits<PartitionID>::max();
const int NOTINQUEUE 		       = std::numeric_limits<int>::max();
const int ROOT 			       = 0;

//for the gpa algorithm
struct edge_source_pair {
        EdgeID e;
        NodeID source;       
};

struct source_target_pair {
        NodeID source;       
        NodeID target;       
};

//matching array has size (no_of_nodes), so for entry in this table we get the matched neighbor
typedef std::vector<NodeID> CoarseMapping;
typedef std::vector<NodeID> Matching;
typedef std::vector<NodeID> NodePermutationMap;

typedef double ImbalanceType;
//Coarsening
typedef enum {
        EXPANSIONSTAR, 
        EXPANSIONSTAR2, 
 	WEIGHT, 
 	REALWEIGHT, 
	PSEUDOGEOM, 
	EXPANSIONSTAR2ALGDIST, 
        SEPARATOR_MULTX,
        SEPARATOR_ADDX,
        SEPARATOR_MAX,
        SEPARATOR_LOG,
        SEPARATOR_R1,
        SEPARATOR_R2,
        SEPARATOR_R3,
        SEPARATOR_R4,
        SEPARATOR_R5,
        SEPARATOR_R6,
        SEPARATOR_R7,
        SEPARATOR_R8
} EdgeRating;

typedef enum {
        PERMUTATION_QUALITY_NONE, 
	PERMUTATION_QUALITY_FAST,  
	PERMUTATION_QUALITY_GOOD
} PermutationQuality;

typedef enum {
        MATCHING_RANDOM, 
	MATCHING_GPA, 
	MATCHING_RANDOM_GPA,
        CLUSTER_COARSENING
} MatchingType;

typedef enum {
	INITIAL_PARTITIONING_RECPARTITION, 
	INITIAL_PARTITIONING_BIPARTITION
} InitialPartitioningType;

typedef enum {
        REFINEMENT_SCHEDULING_FAST, 
	REFINEMENT_SCHEDULING_ACTIVE_BLOCKS, 
	REFINEMENT_SCHEDULING_ACTIVE_BLOCKS_REF_KWAY
} RefinementSchedulingAlgorithm;

typedef enum {
        REFINEMENT_TYPE_FM, 
	REFINEMENT_TYPE_FM_FLOW, 
	REFINEMENT_TYPE_FLOW
} RefinementType;

typedef enum {
        STOP_RULE_SIMPLE, 
	STOP_RULE_MULTIPLE_K, 
	STOP_RULE_STRONG 
} StopRule;

typedef enum {
        BIPARTITION_BFS, 
	BIPARTITION_FM
} BipartitionAlgorithm ;

typedef enum {
        KWAY_SIMPLE_STOP_RULE, 
	KWAY_ADAPTIVE_STOP_RULE
} KWayStopRule;

typedef enum {
        COIN_RNDTIE, 
	COIN_DIFFTIE, 
	NOCOIN_RNDTIE, 
	NOCOIN_DIFFTIE 
} MLSRule;

typedef enum {
        CYCLE_REFINEMENT_ALGORITHM_PLAYFIELD, 
        CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL, 
	CYCLE_REFINEMENT_ALGORITHM_ULTRA_MODEL_PLUS
} CycleRefinementAlgorithm;

typedef enum {
        RANDOM_NODEORDERING, 
        DEGREE_NODEORDERING
} NodeOrderingType;

typedef enum {
        NSQUARE, 
        NSQUAREPRUNED, 
        COMMUNICATIONGRAPH
} LsNeighborhoodType;

typedef enum {
        MAP_CONST_RANDOM, 
        MAP_CONST_IDENTITY,
        MAP_CONST_OLDGROWING,
        MAP_CONST_OLDGROWING_FASTER,
        MAP_CONST_OLDGROWING_MATRIX,
        MAP_CONST_FASTHIERARCHY_BOTTOMUP,
        MAP_CONST_FASTHIERARCHY_TOPDOWN
} ConstructionAlgorithm;

typedef enum {
        DIST_CONST_RANDOM, 
        DIST_CONST_IDENTITY,
        DIST_CONST_HIERARCHY,
        DIST_CONST_HIERARCHY_ONLINE
} DistanceConstructionAlgorithm;

typedef enum {
        PRE_CONFIG_MAPPING_FAST, 
        PRE_CONFIG_MAPPING_ECO,
        PRE_CONFIG_MAPPING_STRONG
} PreConfigMapping;


#endif

