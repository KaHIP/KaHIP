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

std::vector<kway_graph_refinement_commons*>* kway_graph_refinement_commons::m_instances = NULL;

kway_graph_refinement_commons::kway_graph_refinement_commons() {

}

kway_graph_refinement_commons::~kway_graph_refinement_commons() {
}


kway_graph_refinement_commons* kway_graph_refinement_commons::getInstance( PartitionConfig & config ) {
        bool created = false;
        #ifdef USE_OPENMP
        int max_threads = omp_get_max_threads();
        #pragma omp critical
        #else
        int  max_threads = 1;
        #endif
        {
                if( m_instances == NULL ) {
                        m_instances = new std::vector< kway_graph_refinement_commons*>(max_threads, NULL);
                }
        }
        #ifdef USE_OPENMP
        int id = omp_get_thread_num();
        #else
        int id = 0;
        #endif
        if((*m_instances)[id] == NULL) {
                (*m_instances)[id] = new kway_graph_refinement_commons();
                (*m_instances)[id]->init(config);
                created = true;
        }

        if(created == false) {
                if(config.k != (*m_instances)[id]->getUnderlyingK()) {
                        //should be a very rare case
                        (*m_instances)[id]->init(config);
                }
        }

        return  (*m_instances)[id];
}
