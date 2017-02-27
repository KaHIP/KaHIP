/******************************************************************************
 * parallel_vector_io.h
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

#ifndef PARALLEL_VECTOR_IO_BZVNZ570A
#define PARALLEL_VECTOR_IO_BZVNZ570A

#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "parallel_graph_io.h"
#include "partition_config.h"

const ULONG fileTypeVersionNumberPartition = 1;
const ULONG header_count_partition         = 2;

class parallel_vector_io {
public:
        parallel_vector_io();
        virtual ~parallel_vector_io();

        template<typename vectortype> 
        void readVectorSequentially(std::vector<vectortype> & vec, std::string filename);

        template<typename vectortype> 
        void writeVectorSequentially(std::vector<vectortype> & vec, std::string filename);

        void writePartitionSimpleParallel(parallel_graph_access & G, std::string filename);

        void writePartitionBinaryParallel(PPartitionConfig & config, parallel_graph_access & G, std::string filename);
        void writePartitionBinaryParallelPosix(PPartitionConfig & config, parallel_graph_access & G, std::string filename);

        void readPartitionSimpleParallel(parallel_graph_access & G, std::string filename);

        void readPartitionBinaryParallel(PPartitionConfig & config, parallel_graph_access & G, std::string filename);

        void readPartition(PPartitionConfig & config, parallel_graph_access & G, std::string filename);
};

template<typename vectortype> 
void parallel_vector_io::writeVectorSequentially(std::vector<vectortype> & vec, std::string filename) {
        std::ofstream f(filename.c_str());
        for( ULONG i = 0; i < vec.size(); ++i) {
                f << vec[i] <<  std::endl;
        }

        f.close();
}

template<typename vectortype> 
void parallel_vector_io::readVectorSequentially(std::vector<vectortype> & vec, std::string filename) {

        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening vectorfile" << filename << std::endl;
                return;
        }

        ULONG pos = 0;
        std::getline(in, line);
        while( !in.eof() ) {
                if (line[0] == '%') { //Comment
                        continue;
                }

                vectortype value = (vectortype) atof(line.c_str());
                vec[pos++] = value;
                std::getline(in, line);
        }

        in.close();
}

#endif /* end of include guard: PARALLEL_VECTOR_IO_BZVNZ570 */
