/******************************************************************************
 * kahip_fortran.cpp
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 * C helpers for the Fortran module kahip (kahip.F90). Fortran cannot free
 * memory that the C++ library allocates with new[], so node_separator is
 * wrapped here: the separator is copied into an array provided by the caller.
 *****************************************************************************/

#include <algorithm>

#include "kaHIP_interface.h"

extern "C" void kahip_fortran_node_separator(int* n, int* vwgt, kahip_idx* xadj,
                                             kahip_idx* adjcwgt, kahip_idx* adjncy, int* nparts,
                                             double* imbalance, bool suppress_output, int seed, int mode,
                                             int* num_separator_vertices, int* separator) {
        int* internal_separator = nullptr;
        node_separator(n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, suppress_output, seed, mode,
                       num_separator_vertices, &internal_separator);
        std::copy(internal_separator, internal_separator + *num_separator_vertices, separator);
        delete[] internal_separator;
}
