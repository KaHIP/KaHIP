/******************************************************************************
 * parallel_vector_io.cpp
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

#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <unistd.h>
#include "parallel_vector_io.h"
#include "tools/helpers.h"

parallel_vector_io::parallel_vector_io() {
                
}

parallel_vector_io::~parallel_vector_io() {
                
}

void parallel_vector_io::writePartitionBinaryParallelPosix(PPartitionConfig & config,
                                                      parallel_graph_access & G, 
                                                      std::string filename) {

        PEID rank, size;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);
        
        if( rank == ROOT ) {
                ULONG n = G.number_of_global_nodes();
                int output_fd = open(filename.c_str(), O_WRONLY | O_CREAT, 0644);
                write(output_fd, (char*)(&fileTypeVersionNumberPartition), sizeof( ULONG ));
                write(output_fd, (char*)(&n), sizeof( ULONG ));
                close(output_fd);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        PEID window_size = std::min(config.binary_io_window_size, size);
        PEID lowPE = 0;
        PEID highPE = window_size;
        while ( lowPE < size ) {
                if( rank >= lowPE && rank < highPE ) {
                        int output_fd = open(filename.c_str(), O_WRONLY, 0644);

                        ULONG from = G.get_from_range();
                        ULONG start_pos = (header_count_partition + from)*(sizeof(ULONG));
                        lseek( output_fd, start_pos, SEEK_SET);

                        std::vector< ULONG > partition_ids(G.number_of_local_nodes(),0);
                        forall_local_nodes(G, node) {
                               ULONG block = G.getNodeLabel(node);
                               partition_ids[node]= block;
                        } endfor

                        write(output_fd, (char*)(&partition_ids[0]), G.number_of_local_nodes()*sizeof( ULONG ));
                        close(output_fd);
                }
                lowPE  += window_size;
                highPE += window_size;
                MPI_Barrier(MPI_COMM_WORLD);
        }
         
        MPI_Barrier(MPI_COMM_WORLD);
}

void parallel_vector_io::writePartitionBinaryParallel(PPartitionConfig & config,
                                                      parallel_graph_access & G, 
                                                      std::string filename) {

        PEID rank, size;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);
        
        if( rank == ROOT ) {
                // ROOT writes head
                ULONG n = G.number_of_global_nodes();
                std::ofstream outfile;
                outfile.open(filename.c_str(), std::ios::binary | std::ios::out);
                outfile.write((char*)(&fileTypeVersionNumberPartition), sizeof( ULONG ));
                outfile.write((char*)(&n), sizeof( ULONG ));
                outfile.close();
                
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        PEID window_size = 1;
        PEID lowPE = 0;
        PEID highPE = window_size;
        while ( lowPE < size ) {
                if( rank >= lowPE && rank < highPE ) {
                        std::ofstream file;
                        file.open(filename.c_str(), std::ios::binary | std::ios::out | std::ios::app);

                        ULONG from = G.get_from_range();
                        ULONG start_pos = (header_count_partition + from)*(sizeof(ULONG));
                        file.seekp(start_pos);

                        std::vector< ULONG > partition_ids(G.number_of_local_nodes(),0);
                        forall_local_nodes(G, node) {
                               ULONG block = G.getNodeLabel(node);
                               partition_ids[node]= block;
                        } endfor
                        file.write((char*)(&partition_ids[0]), G.number_of_local_nodes()*sizeof( ULONG ));
                        
                        file.close();
                }
                lowPE  += window_size;
                highPE += window_size;
                MPI_Barrier(MPI_COMM_WORLD);
        }
         
        MPI_Barrier(MPI_COMM_WORLD);
}

void parallel_vector_io::readPartitionBinaryParallel(PPartitionConfig & config,
                                                     parallel_graph_access & G, 
                                                     std::string filename) {

        PEID rank, size;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);
        
        if( rank == ROOT ) {
                // ROOT reades head
                std::cout <<  "reading binary partition"  << std::endl;
                std::ifstream file;
                file.open(filename.c_str(), std::ios::binary | std::ios::in);
                std::vector< ULONG > buffer(2, 0);
                if(file) {
                        file.read((char*)(&buffer[0]), 2*sizeof(ULONG));
                        if( buffer[0] != fileTypeVersionNumberPartition ) {
                                std::cout <<  "filetype version mismatch " <<  buffer[0] << "!=" <<  fileTypeVersionNumberPartition  << std::endl;
                                exit(0);
                        }
                        if( buffer[1] != G.number_of_global_nodes()) {
                                std::cout <<  "wrong number of nodes in partition file"  << std::endl;
                                exit(0);
                        }
                }
                file.close();
        }
        
        PEID window_size = std::min(config.binary_io_window_size, size);
        PEID lowPE = 0;
        PEID highPE = window_size;
        while ( lowPE < size ) {
                if( rank >= lowPE && rank < highPE ) {
                        std::ifstream file;
                        file.open(filename.c_str(), std::ios::binary | std::ios::in);

                        ULONG from = G.get_from_range();
                        ULONG ids_to_read = G.number_of_local_nodes();
                        ULONG start_pos = (header_count_partition + from)*(sizeof(ULONG));
                        file.seekg(start_pos);

                        std::vector< ULONG > partition_ids(ids_to_read, 0);
                        file.read((char*)(&partition_ids[0]), ids_to_read*sizeof( ULONG ));
                        file.close();
                        forall_local_nodes(G, node) {
                                G.setNodeLabel(node, partition_ids[node]);
                        } endfor
                        
                }
                lowPE  += window_size;
                highPE += window_size;
                MPI_Barrier(MPI_COMM_WORLD);
        }
         
        MPI_Barrier(MPI_COMM_WORLD);
        G.update_ghost_node_data_global();
        MPI_Barrier(MPI_COMM_WORLD);
}


void parallel_vector_io::writePartitionSimpleParallel(parallel_graph_access & G, 
                                                      std::string filename) {
        PEID rank, size;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);
        
        if( rank == ROOT ) {
                std::ofstream f(filename.c_str());

                forall_local_nodes(G, node) {
                        f <<  G.getNodeLabel(node) ;
                        f <<  std::endl;
                } endfor

                f.close();
        } 

        for( int i = 1; i < size; i++) {
                MPI_Barrier(MPI_COMM_WORLD);
                
                if( rank == i ) {
                        std::ofstream f;
                        f.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
                        forall_local_nodes(G, node) {
                                f <<  G.getNodeLabel(node) ;
                                f <<  std::endl;
                        } endfor
                        f.close();
                }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
}
void parallel_vector_io::readPartition(PPartitionConfig & config, parallel_graph_access & G, 
                                       std::string filename) {
        std::string text_ending(".txtp");
        std::string bin_ending(".binp");

        if( hasEnding(filename, text_ending) ) {
                return readPartitionSimpleParallel(G, filename);
        }

        if( hasEnding(filename, bin_ending) ) {
                return readPartitionBinaryParallel(config, G, filename);
        }
}

void parallel_vector_io::readPartitionSimpleParallel(parallel_graph_access & G, 
                                                      std::string filename) {
        PEID rank, size;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size);
        
        MPI_Barrier(MPI_COMM_WORLD);
        if( rank == ROOT ) { std::cout <<  "reading text partition"  << std::endl; }

        std::string line;
        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening file" << filename << std::endl;
                return;
        }

        NodeID counter  = 0;
        NodeID from = G.get_from_range();
        NodeID to   = G.get_to_range();

        std::getline(in, line);
        while( !in.eof() ) {
                if( counter > to ) {
                        break;
                }

                if( counter >= from ) {
                        PartitionID block_id = (PartitionID) atof(line.c_str());
                        G.setNodeLabel(counter-from, block_id);
                }
                counter++;
                std::getline(in, line);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        G.update_ghost_node_data_global();
        MPI_Barrier(MPI_COMM_WORLD);
        
}

