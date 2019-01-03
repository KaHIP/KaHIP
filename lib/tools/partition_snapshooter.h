/******************************************************************************
 * partition_snapshooter.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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
