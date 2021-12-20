/******************************************************************************
 * kway_graph_refinement_commons.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include "kway_graph_refinement_commons.h"

kway_graph_refinement_commons::kway_graph_refinement_commons( PartitionConfig & config ) {
        init(config);
}

kway_graph_refinement_commons::~kway_graph_refinement_commons() {
}

