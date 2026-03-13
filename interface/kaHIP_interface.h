/******************************************************************************
 * kaffpa_interface.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *****************************************************************************/


#ifndef KAFFPA_INTERFACE_RYEEZ6WJ
#define KAFFPA_INTERFACE_RYEEZ6WJ

#include <stdint.h>

#ifdef KAHIP_64BIT
typedef int64_t kahip_idx;
#else
typedef int32_t kahip_idx;
#endif

#ifdef __cplusplus

extern "C"
{
#endif

// returns the size of kahip_idx in bytes (4 for 32-bit, 8 for 64-bit)
int kahip_sizeof_idx();

const int FAST           = 0;
const int ECO            = 1;
const int STRONG         = 2;
const int FASTSOCIAL     = 3;
const int ECOSOCIAL      = 4;
const int STRONGSOCIAL   = 5;

const int MAPMODE_MULTISECTION = 0;
const int MAPMODE_BISECTION = 1;

// same data structures as in metis
// edgecut and part are output parameters
// part has to be an array of n ints
void kaffpa(int* n, int* vwgt, kahip_idx* xadj,
                   kahip_idx* adjcwgt, kahip_idx* adjncy, int* nparts,
                   double* imbalance, bool suppress_output, int seed, int mode,
                   kahip_idx* edgecut, int* part);

// same as kaffpa, provides an additional parameter for perfect balance
void kaffpa_balance(int* n, int* vwgt, kahip_idx* xadj,
                   kahip_idx* adjcwgt, kahip_idx* adjncy, int* nparts,
                   double* imbalance,
                   bool perfectly_balance,
                   bool suppress_output, int seed, int mode,
                   kahip_idx* edgecut, int* part);

// balance constraint on nodes and edges
void kaffpa_balance_NE(int* n, int* vwgt, kahip_idx* xadj,
                kahip_idx* adjcwgt, kahip_idx* adjncy, int* nparts,
                double* imbalance,  bool suppress_output, int seed, int mode,
                kahip_idx* edgecut, int* part);

// same data structures as in metis
// edgecut and part and qap are output parameters
// part has to be an array of n ints
void process_mapping(int* n, int* vwgt, kahip_idx* xadj,
                   kahip_idx* adjcwgt, kahip_idx* adjncy,
                   int* hierarchy_parameter,  int* distance_parameter, int hierarchy_depth,
                   int mode_partitioning, int mode_mapping,
                   double* imbalance,
                   bool suppress_output, int seed,
                   kahip_idx* edgecut, int* qap, int* part);


void node_separator(int* n, int* vwgt, kahip_idx* xadj,
                    kahip_idx* adjcwgt, kahip_idx* adjncy, int* nparts,
                    double* imbalance,  bool suppress_output, int seed, int mode,
                    int* num_separator_vertices, int** separator);

// takes an unweighted graph and performs reduced nested dissection
// ordering is the output parameter, an array of n ints
void reduced_nd(int* n, kahip_idx* xadj, kahip_idx* adjncy,
                bool suppress_output, int seed, int mode,
                int* ordering);

void edge_partitioning(int* n, int* vwgt, kahip_idx* xadj,
                   kahip_idx* adjcwgt, kahip_idx* adjncy, int* nparts,
                   double* imbalance, bool suppress_output, int seed, int mode,
                   int* vertexcut, int* part, kahip_idx infinity_edge_weight = 1000);

#ifdef USEMETIS
// reduced nested dissection with metis
void reduced_nd_fast(int* n, kahip_idx* xadj, kahip_idx* adjncy,
                      bool suppress_output, int seed, int* ordering);
#endif

#ifdef __cplusplus
}
#endif

#endif /* end of include guard: KAFFPA_INTERFACE_RYEEZ6WJ */
