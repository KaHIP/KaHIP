/*
 *  * ParHIP.h
 *  *
 *  Author: Christian Schulz <christian.schulz.phone@gmail.com>
 *  */

#ifndef PARHIP_INTERFACE
#define PARHIP_INTERFACE
#include <mpi.h>
#ifdef __cplusplus

extern "C"
{
#endif
typedef unsigned long long idxtype;

const int ULTRAFASTMESH   = 0;
const int FASTMESH        = 1;
const int ECOMESH         = 2;
const int ULTRAFASTSOCIAL = 3;
const int FASTSOCIAL      = 4;
const int ECOSOCIAL       = 5;

void ParHIPPartitionKWay(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt,
                         int *nparts, double* imbalance, bool suppress_output, int seed, int mode, int *edgecut, idxtype *part, 
                         MPI_Comm *comm);
#ifdef __cplusplus
}
#endif


#endif /* end of include guard: PARHIP_INTERFAVE */
