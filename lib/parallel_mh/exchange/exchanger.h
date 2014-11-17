/******************************************************************************
 * exchanger.h 
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

#ifndef EXCHANGER_YPB6QKNL
#define EXCHANGER_YPB6QKNL

#include "data_structure/graph_access.h"
#include "parallel_mh/population.h"
#include "partition_config.h"
#include "tools/quality_metrics.h"

class exchanger {
public:
        exchanger();
        virtual ~exchanger();

        void diversify_population( PartitionConfig & config, graph_access & G, population & island, bool replace );
        void quick_start( PartitionConfig & config,  graph_access & G, population & island );
        void push_best( PartitionConfig & config,  graph_access & G, population & island );
        void recv_incoming( PartitionConfig & config,  graph_access & G, population & island );

private:
        void exchange_individum(const PartitionConfig & config, 
                                graph_access & G, 
                                int & from, 
                                int & rank, 
                                int & to, 
                                Individuum & in, Individuum & out);

        std::vector< int* >          m_partition_map_buffers;
        std::vector< MPI::Request* > m_request_pointers;
        std::vector<bool>            m_allready_send_to;

        int m_prev_best_objective;
        int m_max_num_pushes;
        int m_cur_num_pushes;

        quality_metrics m_qm;
};



#endif /* end of include guard: EXCHANGER_YPB6QKNL */
