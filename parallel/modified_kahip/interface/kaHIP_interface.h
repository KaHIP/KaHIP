/******************************************************************************
 * kaHIP_interface.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef KAFFPA_INTERFACE_RYEEZ6WJ
#define KAFFPA_INTERFACE_RYEEZ6WJ

#ifdef __cplusplus

#include <mpi.h>

extern "C"
{
#endif

const int FAST            = 0;
const int ECO             = 1;
const int STRONG          = 2;
const int FASTSOCIAL      = 3;
const int ECOSOCIAL       = 4;
const int STRONGSOCIAL    = 5;
const int ULTRAFASTSOCIAL = 6;

// same data structures as in metis 
// edgecut and part are output parameters
// part has to be an array of n ints
void kaffpa(int* n, int* vwgt, int* xadj, 
                   int* adjcwgt, int* adjncy, int* nparts, 
                   double* imbalance,  bool suppress_output, int seed, int mode,
                   int* edgecut, int* part);


void node_separator(int* n, int* vwgt, int* xadj, 
                    int* adjcwgt, int* adjncy, int* nparts, 
                    double* imbalance,  bool suppress_output, int seed, int mode,
                    int* num_separator_vertices, int** separator); 

void kaffpaE(int* n, int* vwgt, int* xadj, int* adjcwgt, 
                    int* adjncy, int* nparts, 
                    double* inbalance,  bool suppress_output, 
                    bool graph_partitioned, int time_limit, int seed, 
                    int mode, 
                    MPI_Comm communicator, 
                    int* edgecut, double* balance, int* part); 


#ifdef __cplusplus
}
#endif

#endif /* end of include guard: KAFFPA_INTERFACE_RYEEZ6WJ */
