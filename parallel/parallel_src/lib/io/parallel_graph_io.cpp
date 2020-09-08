/******************************************************************************
 * parallel_graph_io.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#define _FILE_OFFSET_BITS 64

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <mpi.h>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <vector>

#include "parallel_graph_io.h"
#include "tools/helpers.h"

const ULONG fileTypeVersionNumber = 3;
const ULONG header_count          = 3;


parallel_graph_io::parallel_graph_io() {
                
}

parallel_graph_io::~parallel_graph_io() {
                
}

int parallel_graph_io::readGraphWeighted(PPartitionConfig & config, 
                                         parallel_graph_access & G, 
                                         std::string filename, 
                                         PEID peID, PEID comm_size, MPI_Comm communicator) {

        std::string metis_ending(".graph");
        std::string bin_ending(".bgf");

        if( hasEnding(filename, metis_ending) ) {
                std::stringstream ss; 
                ss << filename << bin_ending;
                if(file_exists(ss.str())) {
                        return readGraphBinary(config, G, ss.str(), peID, comm_size, communicator);
                } else {
                        return readGraphWeightedFlexible(G, filename, peID, comm_size, communicator);
                }
        }

        if( hasEnding(filename, bin_ending) ) {
                return readGraphBinary(config, G, filename, peID, comm_size, communicator);
        }

        //non of both is true -- try metis format
        return readGraphWeightedFlexible(G, filename, peID, comm_size, communicator);
}

int parallel_graph_io::readGraphWeightedFlexible(parallel_graph_access & G, 
                                                 std::string filename, 
                                                 PEID peID, PEID comm_size, MPI_Comm communicator) {
        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << filename << std::endl;
                return 1;
        }

        NodeID nmbNodes;
        EdgeID nmbEdges;

        std::getline(in,line);
        //skip comments
        while( line[0] == '%' ) {
                std::getline(in, line);
        }

        int ew = 0;
        std::stringstream ss(line);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;

        if(ew != 0) {
                in.close();
                if(peID == 0) std::cout <<  "graph is weighted --> using a different IO routine"  << std::endl;
                return readGraphWeightedMETIS_fixed(G, filename, peID, comm_size, communicator);
        }
        
        // pe p reads the lines p*ceil(n/size) to (p+1)floor(n/size) lines of that file
        ULONG from  = peID     * ceil(nmbNodes / (double)comm_size);
        ULONG to    = (peID+1) * ceil(nmbNodes / (double)comm_size) - 1;
        to = std::min(to, nmbNodes-1);

        ULONG local_no_nodes = to - from + 1;
        PRINT(std::cout <<  "peID " <<  peID <<  " from " <<  from <<  " to " <<  to  <<  " amount " <<  local_no_nodes << std::endl;);

        std::vector< std::vector< NodeID > > local_edge_lists;
        local_edge_lists.resize(local_no_nodes);
        
        ULONG counter  = 0;
        NodeID node_counter = 0;
        EdgeID edge_counter = 0;

        char *oldstr, *newstr;
        while( std::getline(in, line) ) {
                if( counter > to ) {
                        break;
                }
                if (line[0] == '%') { // a comment in the file
                        continue;
                }

                if( counter >= from ) {
                        oldstr = &line[0];
                        newstr = 0;

                        for (;;) {
                                NodeID target;
                                target = (NodeID) strtol(oldstr, &newstr, 10);

                                if (target == 0) {
                                        break;
                                }

                                oldstr = newstr; 

                                local_edge_lists[node_counter].push_back(target);
                                edge_counter++;

                        }

                        node_counter++;
                }

                counter++;

                if( in.eof() ) {
                        break;
                }
        }

        MPI_Barrier(communicator);
        
        G.start_construction(local_no_nodes, 2*edge_counter, nmbNodes, 2*nmbEdges);
        G.set_range(from, to);

        std::vector< NodeID > vertex_dist( comm_size+1, 0 );
        for( PEID peID = 0; peID <= comm_size; peID++) {
                vertex_dist[peID] = peID * ceil(nmbNodes / (double)comm_size); // from positions
        }
        G.set_range_array(vertex_dist);

        for (NodeID i = 0; i < local_no_nodes; ++i) {
                NodeID node = G.new_node();
                G.setNodeWeight(node, 1);
                G.setNodeLabel(node, from+node);
                G.setSecondPartitionIndex(node, 0);
        
                for( ULONG j = 0; j < local_edge_lists[i].size(); j++) {
                        NodeID target = local_edge_lists[i][j]-1; // -1 since there are no nodes with id 0 in the file
                        EdgeID e = G.new_edge(node, target);
                        G.setEdgeWeight(e, 1);
                }
        }

        G.finish_construction();
        MPI_Barrier(communicator);

        return 0;
}

//int parallel_graph_io::readGraphWeightedMETISFast(parallel_graph_access & G, 
                                         //std::string filename, 
                                         //PEID peID, PEID comm_size, MPI_Comm communicator) {
        //std::string line;

        //// open file for reading
        //std::ifstream in(filename.c_str());
        //if (!in) {
                //std::cerr << "Error opening " << filename << std::endl;
                //return 1;
        //}

        //NodeID nmbNodes;
        //EdgeID nmbEdges;

        //std::getline(in,line);
        ////skip comments
        //while( line[0] == '%' ) {
                //std::getline(in, line);
        //}

        //int ew = 0;
        //std::stringstream ss(line);
        //ss >> nmbNodes;
        //ss >> nmbEdges;
        //ss >> ew;

        //// pe p reads the lines p*ceil(n/size) to (p+1)floor(n/size) lines of that file
        //ULONG from  = peID     * ceil(nmbNodes / (double)comm_size);
        //ULONG to    = (peID+1) * ceil(nmbNodes / (double)comm_size) - 1;
        //to = std::min(to, nmbNodes-1);

        //ULONG local_no_nodes = to - from + 1;
        //std::cout <<  "peID " <<  peID <<  " from " <<  from <<  " to " <<  to  <<  " amount " <<  local_no_nodes << std::endl;

        //std::vector< std::vector< NodeID > > local_edge_lists;
        //local_edge_lists.resize(local_no_nodes);
        
        //ULONG counter  = 0;
        //NodeID node_counter = 0;
        //EdgeID edge_counter = 0;

        //char *oldstr, *newstr;
        //while( std::getline(in, line) ) {
                //if( counter > to ) {
                        //break;
                //}
                //if (line[0] == '%') { // a comment in the file
                        //continue;
                //}

                //if( counter >= from ) {
                        //oldstr = &line[0];
                        //newstr = 0;

                        //for (;;) {
                                //NodeID target;
                                //target = (NodeID) strtol(oldstr, &newstr, 10);

                                //if (target == 0) {
                                        //break;
                                //}

                                //oldstr = newstr; 

                                //local_edge_lists[node_counter].push_back(target);
                                //edge_counter++;

                        //}

                        //node_counter++;
                //}

                //counter++;

                //if( in.eof() ) {
                        //break;
                //}
        //}

        //MPI_Barrier(communicator);
        
        //G.start_construction(local_no_nodes, 2*edge_counter, nmbNodes, 2*nmbEdges);
        //G.set_range(from, to);

        //for (NodeID i = 0; i < local_no_nodes; ++i) {
                //NodeID node = G.new_node();
                //G.setNodeWeight(node, 1);
                //G.setNodeLabel(node, from+node);
                //G.setSecondPartitionIndex(node, 0);
        
                //for( ULONG j = 0; j < local_edge_lists[i].size(); j++) {
                        //NodeID target = local_edge_lists[i][j]-1; // -1 since there are no nodes with id 0 in the file
                        //EdgeID e = G.new_edge(node, target);
                        //G.setEdgeWeight(e, 1);
                //}
        //}

        //G.finish_construction();
        //MPI_Barrier(communicator);
        //return 0;
//}
// we start with the simplest version of IO 
// where each process reads the graph sequentially
// 
//int parallel_graph_io::readGraphWeightedMETIS(parallel_graph_access & G, 
                                         //std::string filename, 
                                         //PEID peID, PEID comm_size, MPI_Comm communicator) {
        //std::string line;

        //// open file for reading
        //std::ifstream in(filename.c_str());
        //if (!in) {
                //std::cerr << "Error opening " << filename << std::endl;
                //return 1;
        //}

        //NodeID nmbNodes;
        //EdgeID nmbEdges;

        //std::getline(in,line);
        ////skip comments
        //while( line[0] == '%' ) {
                //std::getline(in, line);
        //}

        //int ew = 0;
        //std::stringstream ss(line);
        //ss >> nmbNodes;
        //ss >> nmbEdges;
        //ss >> ew;

        //if(ew == 1) {
                //std::cout <<  "io of weighted graphs not supported yet"  << std::endl;
                //exit(0);
        //} else if (ew == 11) {
                //std::cout <<  "io of weighted graphs not supported yet"  << std::endl;
                //exit(0);
        //} else if (ew == 10) {
                //std::cout <<  "io of weighted graphs not supported yet"  << std::endl;
                //exit(0);
        //}

        //// pe p reads the lines p*ceil(n/size) to (p+1)floor(n/size) lines of that file
        //ULONG from           = peID     * ceil(nmbNodes / (double)comm_size);
        //ULONG to             = (peID+1) * ceil(nmbNodes / (double)comm_size) - 1;
        //to = std::min(to, nmbNodes-1);

        //ULONG local_no_nodes = to - from + 1;
        //std::cout <<  "peID " <<  peID <<  " from " <<  from <<  " to " <<  to  <<  " amount " <<  local_no_nodes << std::endl;

        //std::vector< std::vector< NodeID > > local_edge_lists;
        //local_edge_lists.resize(local_no_nodes);
        

        ////std::getline(in, line);
        //ULONG counter      = 0;
        //NodeID node_counter = 0;
        //EdgeID edge_counter = 0;

        //while( std::getline(in, line) ) {
                //if( counter > to ) {
                        //break;
                //}
                //if (line[0] == '%') { // a comment in the file
                        //continue;
                //}

                //if( counter >= from ) {
                        //std::stringstream ss(line);

                        //NodeID target;
                        //while( ss >> target ) {
                                //local_edge_lists[node_counter].push_back(target);
                                //edge_counter++;
                        //}
                        //node_counter++;
                //}

                //counter++;

                //if( in.eof() ) {
                        //break;
                //}
        //}

        //MPI_Barrier(communicator);
        
        //G.start_construction(local_no_nodes, 2*edge_counter, nmbNodes, 2*nmbEdges);
        //G.set_range(from, to);

        //for (NodeID i = 0; i < local_no_nodes; ++i) {
                //NodeID node = G.new_node();
                //G.setNodeWeight(node, 1);
                //G.setNodeLabel(node, from+node);
                //G.setSecondPartitionIndex(node, 0);
        
                //for( ULONG j = 0; j < local_edge_lists[i].size(); j++) {
                        //NodeID target = local_edge_lists[i][j]-1; // -1 since there are no nodes with id 0 in the file
                        //EdgeID e = G.new_edge(node, target);
                        //G.setEdgeWeight(e, 1);
                //}
        //}

        //G.finish_construction();
        //return 0;
//}

int parallel_graph_io::writeGraphExternallyBinary(std::string input_filename, std::string output_filename) {

        std::string line;

        // open file for reading
        std::ifstream in(input_filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << input_filename << std::endl;
                return 1;
        }

        NodeID n;
        EdgeID m;

        std::getline(in,line);
        //skip comments
        while( line[0] == '%' ) {
                std::getline(in, line);
        }

        int ew = 0;
        std::stringstream ss(line);
        ss >> n;
        ss >> m;
        ss >> ew;

        m *= 2;

        std::ofstream outfile;
        outfile.open(output_filename.c_str(), std::ios::binary | std::ios::out);
        outfile.write((char*)(&fileTypeVersionNumber), sizeof( ULONG ));
        outfile.write((char*)(&n), sizeof( ULONG ));
        outfile.write((char*)(&m), sizeof( ULONG ));
        
        NodeID offset = (header_count + n + 1) * (sizeof(ULONG));

        while( std::getline(in, line) ) {
                if (line[0] == '%') { // a comment in the file
                        continue;
                }

                std::stringstream ss(line);

                EdgeID edge_counter = 0;
                NodeID target;
                while( ss >> target ) {
                        edge_counter++;
                }
                outfile.write((char*)(&offset), sizeof( ULONG ));
                offset += edge_counter*sizeof( ULONG );

                if( in.eof() ) {
                        break;
                }
        }
        outfile.write((char*)(&offset), sizeof( ULONG ));
        in.close();

        // second stream to actually write the edges
        std::ifstream second_in(input_filename.c_str());
        std::getline(second_in,line);
        //skip comments
        while( line[0] == '%' ) {
                std::getline(second_in, line);
        }

        while( std::getline(second_in, line) ) {
                if (line[0] == '%') { // a comment in the file
                        continue;
                }

                std::stringstream ss(line);

                NodeID target;
                while( ss >> target ) {
                        target -= 1;
                        outfile.write((char*)(&target), sizeof( ULONG ));
                }

                if( second_in.eof() ) {
                        break;
                }
        }
        second_in.close();
        
        return 0;

}

int parallel_graph_io::writeGraphSequentiallyBinary(complete_graph_access & G, std::string filename) {

        std::ofstream outfile;
        outfile.open(filename.c_str(), std::ios::binary | std::ios::out);
        PEID size; MPI_Comm_size( MPI_COMM_WORLD, &size);
        
        if( size > 1 ) {
                std::cout <<  "currently only one process supported."  << std::endl;
                return 0;
        }

        std::cout <<  "Writing graph " << filename  << std::endl;
        printf("Writing graph with n = %lld, m = %lld\n", G.number_of_global_nodes(), G.number_of_global_edges());

        //write version number
        outfile.write((char*)(&fileTypeVersionNumber), sizeof( ULONG ));

        //write number of nodes etc
        NodeID n = G.number_of_global_nodes();
        NodeID m = G.number_of_global_edges();

        outfile.write((char*)(&n), sizeof( ULONG ));
        outfile.write((char*)(&m), sizeof( ULONG ));

        NodeID * offset_array = new NodeID[n+1];
        ULONG pos              = 0;
        NodeID offset         = (header_count + G.number_of_global_nodes() + 1) * (sizeof(ULONG));

        forall_local_nodes(G, node) {
                offset_array[pos++] = offset;
                offset += G.getNodeDegree(node) * (sizeof(ULONG)); 
        } endfor

        offset_array[pos] = offset;
        outfile.write((char*)(offset_array), (n+1)*sizeof(ULONG));
        delete[] offset_array;

        NodeID * edge_array = new NodeID[m];
        pos = 0;

        // now write the edges 
        forall_local_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        edge_array[pos++] = target;
                } endfor
        } endfor
        outfile.write((char*)(edge_array), (m)*sizeof(ULONG));
        outfile.close();
        delete[] edge_array;


        return 0;
}

int parallel_graph_io::readGraphBinary(PPartitionConfig & config, parallel_graph_access & G, 
                                            std::string filename, 
                                            PEID peID, PEID size, MPI_Comm communicator) { 

        // read header
        std::vector< ULONG > buffer(3, 0);
        int success = 0;
        if( peID == ROOT) {
                std::cout <<  "Reading binary graph ..."  << std::endl;
                std::ifstream file;
                file.open(filename.c_str(), std::ios::binary | std::ios::in);
                if(file) {
                        success = 1;
                        file.read((char*)(&buffer[0]), 3*sizeof(ULONG));
                }
                file.close();
        }

        MPI_Bcast(&success, 1, MPI_INT, ROOT, communicator);

        if( !success ) {
                if( peID == ROOT ) std::cout <<  "problem to open the file"  << std::endl;
                MPI_Finalize();
                exit(0);
        }

        MPI_Bcast(&buffer[0], 3, MPI_LONG, ROOT, communicator);
        ULONG version = buffer[0];
        NodeID n     = buffer[1];
        NodeID m     = buffer[2];

        if(peID == ROOT) std::cout <<  "version: " <<  version <<  " n: "<<  n <<  " m: " <<  m  << std::endl;
        if( version != fileTypeVersionNumber ) {
                if(peID == ROOT) std::cout <<  "filetype version missmatch"  << std::endl;
                MPI_Finalize(); exit(0);
        }

        PEID window_size = std::min(config.binary_io_window_size, size);
        PEID lowPE = 0;
        PEID highPE = window_size;

                
        while ( lowPE < size ) {
                if( peID >= lowPE && peID < highPE ) {
                        std::ifstream file;
                        file.open(filename.c_str(), std::ios::binary | std::ios::in);

                        ULONG from = peID * ceil(n / (double)size);
                        ULONG to   = (peID +1) * ceil(n / (double)size) - 1;
                        to = std::min(to, n-1);

                        ULONG local_no_nodes = to - from + 1;
                        std::cout <<  "peID " <<  peID <<  " from " <<  from <<  " to " <<  to  <<  " amount " <<  local_no_nodes << std::endl;

                        // read the offsets
                        ULONG start_pos = (header_count + from)*(sizeof(ULONG));
                        NodeID* vertex_offsets = new NodeID[local_no_nodes+1]; // we also need the next vertex offset
                        file.seekg(start_pos);
                        file.read((char*)(vertex_offsets), (local_no_nodes+1)*sizeof(ULONG));

                        ULONG  edge_start_pos = vertex_offsets[0];
                        EdgeID num_reads             = vertex_offsets[local_no_nodes]-vertex_offsets[0];
                        EdgeID num_edges_to_read     = num_reads/sizeof(ULONG);
                        EdgeID* edges = new EdgeID[num_edges_to_read]; // we also need the next vertex offset
                        file.seekg(edge_start_pos);
                        file.read((char*)(edges), (num_edges_to_read)*sizeof(ULONG));

                        G.start_construction(local_no_nodes, num_edges_to_read, n, m);
                        G.set_range(from, to);

                        std::vector< NodeID > vertex_dist( size+1, 0 );
                        for( PEID peID = 0; peID <= size; peID++) {
                                vertex_dist[peID] = peID * ceil(n / (double)size); // from positions
                        }
                        G.set_range_array(vertex_dist);

                        ULONG pos = 0;
                        for (NodeID i = 0; i < local_no_nodes; ++i) {
                                NodeID node = G.new_node();
                                G.setNodeWeight(node, 1);
                                G.setNodeLabel(node, from+node);
                                G.setSecondPartitionIndex(node, 0);

                                NodeID degree =  (vertex_offsets[i+1] - vertex_offsets[i]) / sizeof(ULONG);
                                for( ULONG j = 0; j < degree; j++, pos++) {
                                        NodeID target = edges[pos]; 
                                        EdgeID e = G.new_edge(node, target);
                                        G.setEdgeWeight(e, 1);
                                }
                        }

                        G.finish_construction();

                        delete[] vertex_offsets;
                        delete[] edges;
                        file.close();
                }
                lowPE  += window_size;
                highPE += window_size;
                MPI_Barrier(communicator);
        }
        
        return 0;
}

int parallel_graph_io::writeGraphParallelSimple(parallel_graph_access & G, 
                                                std::string filename, MPI_Comm communicator) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        if( rank == ROOT ) {
                std::ofstream f(filename.c_str());
                f << G.number_of_global_nodes() <<  " " <<  G.number_of_global_edges()/2 <<   std::endl;

                forall_local_nodes(G, node) {
                        forall_out_edges(G, e, node) {
                                f << (G.getGlobalID(G.getEdgeTarget(e))+1) << " " ;
                        } endfor
                        f <<  "\n";
                } endfor

                f.close();
        } 

        for( int i = 1; i < size; i++) {
                MPI_Barrier(communicator);
                
                if( rank == i ) {
                        std::ofstream f;
                        f.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
                        forall_local_nodes(G, node) {
                                forall_out_edges(G, e, node) {
                                        f <<  (G.getGlobalID(G.getEdgeTarget(e))+1) << " " ;
                                } endfor 
                                f <<  "\n";
                        } endfor
                        f.close();
                }
        }

        MPI_Barrier(communicator);
        
        return 0;
}

int parallel_graph_io::writeGraphWeightedParallelSimple(parallel_graph_access & G, 
                                                               std::string filename, MPI_Comm communicator) {
        PEID rank, size;
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);
        
        if( rank == ROOT ) {
                std::ofstream f(filename.c_str());
                f << G.number_of_global_nodes() <<  " " <<  G.number_of_global_edges()/2 << " 11" <<  std::endl;

                forall_local_nodes(G, node) {
                        f <<  G.getNodeWeight(node) ;
                        forall_out_edges(G, e, node) {
                                f << " " <<   (G.getGlobalID(G.getEdgeTarget(e))+1) <<  " " <<  G.getEdgeWeight(e) ;
                        } endfor 
                        f <<  "\n";
                } endfor

                f.close();
        } 

        for( PEID i = 1; i < size; i++) {
                MPI_Barrier(communicator);
                
                if( rank == i ) {
                        std::ofstream f;
                        f.open(filename.c_str(), std::ofstream::out | std::ofstream::app);
                        forall_local_nodes(G, node) {
                                f <<  G.getNodeWeight(node) ;
                                forall_out_edges(G, e, node) {
                                        f << " " <<   (G.getGlobalID(G.getEdgeTarget(e))+1) <<  " " <<  G.getEdgeWeight(e) ;
                                } endfor 
                                f <<  "\n";
                        } endfor
                        f.close();
                }
        }

        MPI_Barrier(communicator);
        
        return 0;
}

int parallel_graph_io::writeGraphWeightedSequentially(complete_graph_access & G, std::string filename) {
        std::ofstream f(filename.c_str());
        f << G.number_of_global_nodes() <<  " " <<  G.number_of_global_edges()/2 <<  " 11" <<  std::endl;

        forall_local_nodes(G, node) {
                f <<  G.getNodeWeight(node) ;
                forall_out_edges(G, e, node) {
                        f << " " <<   (G.getEdgeTarget(e)+1) <<  " " <<  G.getEdgeWeight(e) ;
                } endfor 
                f <<  "\n";
        } endfor

        f.close();
        return 0;
}

int parallel_graph_io::writeGraphSequentially(complete_graph_access & G, std::ofstream & f) {
        f << G.number_of_global_nodes() <<  " " <<  G.number_of_global_edges()/2 <<   std::endl;

        forall_local_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        f << " " <<   (G.getEdgeTarget(e)+1)  ;
                } endfor 
                f <<  "\n";
        } endfor
        return 0;
}

int parallel_graph_io::writeGraphSequentially(complete_graph_access & G, std::string filename) {
        std::ofstream f(filename.c_str());
        writeGraphSequentially(G, f);
        f.close();
        return 0;

}

// we start with the simplest version of IO 
// where each process reads the graph sequentially
// TODO write weighted code and fully parallel io code
int parallel_graph_io::readGraphWeightedMETIS_fixed(parallel_graph_access & G, 
                                         std::string filename, 
                                         PEID peID, PEID comm_size, MPI_Comm communicator) {
        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << filename << std::endl;
                return 1;
        }

        NodeID nmbNodes;
        EdgeID nmbEdges;

        std::getline(in,line);
        //skip comments
        while( line[0] == '%' ) {
                std::getline(in, line);
        }

        int ew = 0;
        std::stringstream ss(line);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;

        // pe p reads the lines p*ceil(n/size) to (p+1)floor(n/size) lines of that file
        ULONG from           = peID     * ceil(nmbNodes / (double)comm_size);
        ULONG to             = (peID+1) * ceil(nmbNodes / (double)comm_size) - 1;
        to = std::min(to, nmbNodes-1);

        unsigned long local_no_nodes = to - from + 1;
        std::cout <<  "peID " <<  peID <<  " from " <<  from <<  " to " <<  to  <<  " amount " <<  local_no_nodes << std::endl;

        std::vector< std::vector< NodeID > > local_edge_lists;
        local_edge_lists.resize(local_no_nodes);
 
        std::vector< std::vector< NodeID > > local_edge_weights;
        local_edge_weights.resize(local_no_nodes);
       
        std::vector< NodeID > local_node_weights;

        bool read_ew = false;
        bool read_nw = false;

        if(ew == 1) {
                read_ew = true;
        } else if (ew == 11) {
                read_ew = true;
                read_nw = true;
        } else if (ew == 10) {
                read_nw = true;
        }

        //std::getline(in, line);
        unsigned long counter      = 0;
        NodeID node_counter = 0;
        EdgeID edge_counter = 0;

        while( std::getline(in, line) ) {
                if( counter > to ) {
                        break;
                }
                if (line[0] == '%') { // a comment in the file
                        continue;
                }

                if( counter >= from ) {
                        std::stringstream ss(line);

                        NodeWeight weight = 1;
                        if( read_nw ) {
                                ss >> weight;
                        }
                        local_node_weights.push_back(weight);

                        NodeID target;
                        while( ss >> target ) {
                                EdgeWeight edge_weight = 1;
                                if( read_ew ) {
                                        ss >> edge_weight;
                                }

                                local_edge_weights[node_counter].push_back(edge_weight);
                                local_edge_lists[node_counter].push_back(target);
                                edge_counter++;
                        }
                        node_counter++;
                }

                counter++;

                if( in.eof() ) {
                        break;
                }
        }

        MPI_Barrier(communicator);
        G.start_construction(local_no_nodes, 2*edge_counter, nmbNodes, 2*nmbEdges);
        G.set_range(from, to);

        std::vector< NodeID > vertex_dist( comm_size+1, 0 );
        for( PEID peID = 0; peID <= comm_size; peID++) {
                vertex_dist[peID] = peID * ceil(nmbNodes / (double)comm_size); // from positions
        }
        G.set_range_array(vertex_dist);

        for (NodeID i = 0; i < local_no_nodes; ++i) {
                NodeID node = G.new_node();
                G.setNodeWeight(node, local_node_weights[i]);
                G.setNodeLabel(node, from+node);
                G.setSecondPartitionIndex(node, 0);
        

                for( unsigned j = 0; j < local_edge_lists[i].size(); j++) {
                        NodeID target = local_edge_lists[i][j]-1; // -1 since there are no nodes with id 0 in the file
                        EdgeID e = G.new_edge(node, target);
                        G.setEdgeWeight(e, local_edge_weights[i][j]);
                }
        }

        G.finish_construction();
        MPI_Barrier(communicator);
        return 0;


}

