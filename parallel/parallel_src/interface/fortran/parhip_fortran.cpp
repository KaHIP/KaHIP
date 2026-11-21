/******************************************************************************
 * parhip_fortran.cpp
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 * C helper for the Fortran module parhip (parhip.F90): Fortran passes MPI
 * communicators as integer handles, which are converted to MPI_Comm here.
 *****************************************************************************/

#include <mpi.h>

#include "parhip_interface.h"

extern "C" void parhip_fortran_partition_kway(idxtype* vtxdist, idxtype* xadj, idxtype* adjncy,
                                              idxtype* vwgt, idxtype* adjwgt, int* nparts,
                                              double* imbalance, bool suppress_output, int seed, int mode,
                                              int* edgecut, idxtype* part, MPI_Fint* fcomm) {
        MPI_Comm comm = MPI_Comm_f2c(*fcomm);
        ParHIPPartitionKWay(vtxdist, xadj, adjncy, vwgt, adjwgt, nparts, imbalance,
                            suppress_output, seed, mode, edgecut, part, &comm);
}
