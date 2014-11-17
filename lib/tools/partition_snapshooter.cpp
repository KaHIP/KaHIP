/******************************************************************************
 * partition_snapshooter.cpp 
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

#include <sstream>

#include "definitions.h"
#include "graph_io.h"
#include "partition_snapshooter.h"

partition_snapshooter* partition_snapshooter::m_instance = NULL;

partition_snapshooter::partition_snapshooter() {
        m_buffer_size = 500;
        m_idx         = 0;
}

partition_snapshooter::~partition_snapshooter() {
        flush_buffer();
}

partition_snapshooter * partition_snapshooter::getInstance() {
        if( m_instance == NULL ) {
                m_instance = new partition_snapshooter();
        } 
        return m_instance;
}

void partition_snapshooter::addSnapshot(graph_access & G) {
        std::cout <<  "idx " << m_partition_map_buffer.size() << std::endl;
        std::vector<PartitionID>* partition_map = new std::vector<PartitionID>();
        m_partition_map_buffer.push_back(partition_map);

        forall_nodes(G, node) {
               partition_map->push_back(G.getPartitionIndex(node)); 
        } endfor

        if( m_partition_map_buffer.size() > m_buffer_size) {
                flush_buffer();
        }
}

void partition_snapshooter::addSnapshot(graph_access & G, std::vector<PartitionID> & ext_partition_map) {
        std::vector<PartitionID>* partition_map = new std::vector<PartitionID>();
        m_partition_map_buffer.push_back(partition_map);

        forall_nodes(G, node) {
               partition_map->push_back(ext_partition_map[node]); 
        } endfor

        if( m_partition_map_buffer.size() > m_buffer_size) {
                flush_buffer();
        }

}

void partition_snapshooter::flush_buffer() {

        for( unsigned i = 0; i < m_partition_map_buffer.size(); i++) {
                std::stringstream snapshot_name;
                snapshot_name << "snapshot_" << m_idx;

                graph_io::writeVector(*(m_partition_map_buffer[i]), snapshot_name.str());

                m_idx++;
        }

        //flush buffer
        for( int i = m_partition_map_buffer.size()-1; i >= 0; i--) {
                delete m_partition_map_buffer[i];
                m_partition_map_buffer.pop_back();
        }

}

void partition_snapshooter::set_buffer_size( unsigned int new_buffer_size ) {
        m_buffer_size = new_buffer_size;
}

