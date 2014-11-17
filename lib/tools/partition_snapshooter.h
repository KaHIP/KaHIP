/******************************************************************************
 * partition_snapshooter.h 
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

#ifndef PARTITION_SNAPSHOOTER_LGCUMS2I
#define PARTITION_SNAPSHOOTER_LGCUMS2I

#include "data_structure/graph_access.h"

//buffered partition snapshooter (singleton)
class partition_snapshooter {
        public: 
                static partition_snapshooter * getInstance();

                void addSnapshot(graph_access & G); 
                void addSnapshot(graph_access & G, std::vector<PartitionID> & partition_map); 

                //flushes buffer to disk
                void flush_buffer();
                void set_buffer_size( unsigned int new_buffer_size );
        private: 
                partition_snapshooter(); 
                partition_snapshooter(const partition_snapshooter&) {}          

                virtual ~partition_snapshooter();
                static partition_snapshooter* m_instance;

                unsigned int m_buffer_size;
                unsigned int m_idx;

                std::vector< std::vector< PartitionID >* > m_partition_map_buffer;
};


#endif /* end of include guard: PARTITION_SNAPSHOOTER_LGCUMS2I */
