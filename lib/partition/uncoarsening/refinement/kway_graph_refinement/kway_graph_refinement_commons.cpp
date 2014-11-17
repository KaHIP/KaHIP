/******************************************************************************
 * kway_graph_refinement_commons.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#include <omp.h>

#include "kway_graph_refinement_commons.h"

std::vector<kway_graph_refinement_commons*>* kway_graph_refinement_commons::m_instances = NULL;

kway_graph_refinement_commons::kway_graph_refinement_commons() {

}

kway_graph_refinement_commons::~kway_graph_refinement_commons() {
}


kway_graph_refinement_commons* kway_graph_refinement_commons::getInstance( PartitionConfig & config ) {

        bool created = false;
        #pragma omp critical 
        {
                if( m_instances == NULL ) {
                        m_instances = new std::vector< kway_graph_refinement_commons*>(omp_get_max_threads(), reinterpret_cast<kway_graph_refinement_commons*>(NULL));
                }
        } 

        int id = omp_get_thread_num();
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
