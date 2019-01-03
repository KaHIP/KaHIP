/******************************************************************************
 * graph_io.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <sstream>
#include "graph_io.h"

graph_io::graph_io() {
                
}

graph_io::~graph_io() {
                
}

int graph_io::writeGraphWeighted(graph_access & G, std::string filename) {
        std::ofstream f(filename.c_str());
        f << G.number_of_nodes() <<  " " <<  G.number_of_edges()/2 <<  " 11" <<  std::endl;

        forall_nodes(G, node) {
                f <<  G.getNodeWeight(node) ;
                forall_out_edges(G, e, node) {
                        f << " " <<   (G.getEdgeTarget(e)+1) <<  " " <<  G.getEdgeWeight(e) ;
                } endfor 
                f <<  std::endl;
        } endfor

        f.close();
        return 0;
}

int graph_io::writeGraph(graph_access & G, std::string filename) {
        std::ofstream f(filename.c_str());
        f << G.number_of_nodes() <<  " " <<  G.number_of_edges()/2 << std::endl;

        forall_nodes(G, node) {
                f <<  node <<  " ";
                forall_out_edges(G, e, node) {
                        f <<   G.getEdgeTarget(e) << " " ;
                } endfor 
                f <<  std::endl;
        } endfor

        f.close();
        return 0;
}

int graph_io::readPartition(graph_access & G, std::string filename) {
        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening file" << filename << std::endl;
                return 1;
        }

        PartitionID max = 0;
        forall_nodes(G, node) {
                // fetch current line
                std::getline(in, line);
                if (line[0] == '%') { //Comment
                        node--;
                        continue;
                }

                // in this line we find the block of Node node 
                G.setPartitionIndex(node, (PartitionID) atol(line.c_str()));

                if(G.getPartitionIndex(node) > max)
                        max = G.getPartitionIndex(node);
        } endfor

        G.set_partition_count(max+1);
        in.close();

        return 0;
}

int graph_io::readGraphWeighted(graph_access & G, std::string filename) {
        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << filename << std::endl;
                return 1;
        }

        long nmbNodes;
        long nmbEdges;

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

        if( 2*nmbEdges > std::numeric_limits<int>::max() || nmbNodes > std::numeric_limits<int>::max()) {
                std::cout <<  "The graph is too large. Currently only 32bit supported!"  << std::endl;
                exit(0);
        }

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
        nmbEdges *= 2; //since we have forward and backward edges
        
        NodeID node_counter   = 0;
        EdgeID edge_counter   = 0;

        G.start_construction(nmbNodes, nmbEdges);

        while(  std::getline(in, line)) {
       
                if (line[0] == '%') { // a comment in the file
                        continue;
                }

                NodeID node = G.new_node(); node_counter++;
                G.setPartitionIndex(node, 0);

                std::stringstream ss(line);

                NodeWeight weight = 1;
                if( read_nw ) {
                        ss >> weight;
                }
                G.setNodeWeight(node, weight);

                NodeID target;
                while( ss >> target ) {
                        EdgeWeight edge_weight = 1;
                        if( read_ew ) {
                                ss >> edge_weight;
                        }
                        edge_counter++;
                        EdgeID e = G.new_edge(node, target-1);
                        G.setEdgeWeight(e, edge_weight);
                }

                if(in.eof()) {
                        break;
                }
        }

        if( edge_counter != (EdgeID)nmbEdges ) {
                std::cout <<  "number of specified edges mismatch"  << std::endl;
                std::cout <<  edge_counter <<  " " <<  nmbEdges  << std::endl;
                exit(0);
        }

        if( node_counter != (NodeID)nmbNodes) {
                std::cout <<  "number of specified nodes mismatch"  << std::endl;
                std::cout <<  node_counter <<  " " <<  nmbNodes  << std::endl;
                exit(0);
        }


        G.finish_construction();
        return 0;
}


void graph_io::writePartition(graph_access & G, std::string filename) {
        std::ofstream f(filename.c_str());
        std::cout << "writing partition to " << filename << " ... " << std::endl;

        forall_nodes(G, node) {
                f << G.getPartitionIndex(node) <<  std::endl;
        } endfor

        f.close();
}


