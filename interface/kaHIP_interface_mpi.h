/******************************************************************************
 * kaHIP_interface_mpi.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *****************************************************************************/

#ifndef KAHIP_INTERFACE_MPI_H
#define KAHIP_INTERFACE_MPI_H

#include <mpi.h>
#include "kaHIP_interface.h"

#ifdef __cplusplus
extern "C"
{
#endif

// distributed evolutionary partitioning (kaffpaE) on the processes of communicator
// every process passes the same graph; part is an input (graph_partitioned) and output array of n ints
// mode: FAST, ECO, STRONG, FASTSOCIAL, ECOSOCIAL, STRONGSOCIAL, ULTRA
// single_individual: every process computes one partition with mode, the best one is returned (no evolution)
void kaffpaE_comm(int* n, int* vwgt, kahip_idx* xadj,
                  kahip_idx* adjcwgt, kahip_idx* adjncy, int* nparts,
                  double* imbalance, bool suppress_output,
                  bool graph_partitioned, int time_limit, int seed,
                  int mode, bool single_individual,
                  MPI_Comm communicator,
                  kahip_idx* edgecut, double* balance, int* part);

#ifdef __cplusplus
}
#endif

#endif
