/******************************************************************************
 * kaffpa_interface.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
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

extern "C"
{
#endif

// same data structures as in metis 
// edgecut and part are output parameters
// part has to be an array of n ints
void kaffpa_strong(int* n, int* vwgt, int* xadj, 
                   int* adjcwgt, int* adjncy, int* nparts, 
                   double* inbalance,  bool suppress_output, int seed,
                   int* edgecut, int* part);

void kaffpa_eco(int* n, int* vwgt, int* xadj, 
                int* adjcwgt, int* adjncy,  
                int* nparts, double* inbalance,  
                bool suppress_output, int seed,
                int* edgecut, int* part);

void kaffpa_fast(int* n, int* vwgt, int* xadj, int* adjcwgt, 
                 int* adjncy, int* nparts, 
                 double* inbalance,  bool suppress_output, int seed,
                 int* edgecut, int* part); 

#ifdef __cplusplus
}
#endif

#endif /* end of include guard: KAFFPA_INTERFACE_RYEEZ6WJ */
