/******************************************************************************
 * graph_io.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <memory>
#include <sstream>
#include <cstring>
#include "graph_io.h"

graph_io::graph_io() {

}

graph_io::~graph_io() {

}

int graph_io::writeGraphWeighted(graph_access & G, const std::string & filename) {
        std::ofstream f(filename.c_str());
        f << G.number_of_nodes() <<  " " <<  G.number_of_edges()/2 <<  " 11" <<  std::endl;

        forall_nodes(G, node) {
                f <<  G.getNodeWeight(node) ;
                forall_out_edges(G, e, node) {
                        f << " " <<   (G.getEdgeTarget(e)+1) <<  " " <<  G.getEdgeWeight(e) ;
                } endfor
                f <<  "\n";
        } endfor

        f.close();
        return 0;
}

int graph_io::writeGraph(graph_access & G, const std::string & filename) {
        std::ofstream f(filename.c_str());
        f << G.number_of_nodes() <<  " " <<  G.number_of_edges()/2 << std::endl;

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        f <<   (G.getEdgeTarget(e)+1) << " " ;
                } endfor
                f <<  "\n";
        } endfor

        f.close();
        return 0;
}

int graph_io::readPartition(graph_access & G, const std::string & filename) {
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

// Reads the whole file into memory and parses it with a hand-written scanner. Semantics are those of the
// former std::getline/std::stringstream reader: lines starting with '%' are comments, every other line
// (including an empty one) is a node, numbers are separated by whitespace.
static inline bool is_space_char(char c) {
        return c == ' ' || c == '\t' || c == '\r' || c == '\v' || c == '\f';
}

int graph_io::readGraphWeighted(graph_access & G, const std::string & filename) {
        FILE * fp = fopen(filename.c_str(), "rb");
        if (!fp) {
                std::cerr << "Error opening " << filename << std::endl;
                return 1;
        }
        // file contents plus a sentinel, in a buffer that is not zeroed before fread overwrites it
        std::unique_ptr<char[]> buf;
        size_t got = 0;
        {
#ifdef _WIN32
                // long is 32 bits on Windows: 64-bit positions for files over 2 GB
                _fseeki64(fp, 0, SEEK_END);
                long long size = _ftelli64(fp);
                _fseeki64(fp, 0, SEEK_SET);
#else
                fseek(fp, 0, SEEK_END);
                long size = ftell(fp);
                fseek(fp, 0, SEEK_SET);
#endif
                buf.reset(new char[size > 0 ? size + 1 : 1]);
                got = size > 0 ? fread(buf.get(), 1, size, fp) : 0;
                buf[got] = '\n'; // sentinel
        }
        fclose(fp);

        const char * p   = buf.get();
        const char * eof = p + got; // position of the sentinel

        // header: first line that does not start with '%'
        const char * line_end = (const char*)memchr(p, '\n', eof - p + 1);
        while( p < eof && *p == '%' ) {
                p = line_end + 1;
                line_end = (const char*)memchr(p, '\n', eof - p + 1);
        }

        long nmbNodes = 0;
        long nmbEdges = 0;
        int ew = 0;
        {
                std::string line(p, line_end);
                std::stringstream ss(line);
                ss >> nmbNodes;
                ss >> nmbEdges;
                ss >> ew;
        }
        p = line_end < eof ? line_end + 1 : eof;

        if( 2*nmbEdges > std::numeric_limits<int>::max() || nmbNodes > std::numeric_limits<int>::max()) {
                std::cerr <<  "The graph is too large. Currently only 32bit supported!"  << std::endl;
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
        long long total_nodeweight = 0;

        G.start_construction(nmbNodes, nmbEdges);

        while( p < eof ) {
                if (*p == '%') { // a comment in the file
                        p = (const char*)memchr(p, '\n', eof - p + 1) + 1;
                        continue;
                }

                if( node_counter >= (NodeID) nmbNodes ) {
                        // more node lines than announced: count them and fail below like the old reader
                        node_counter++;
                        p = (const char*)memchr(p, '\n', eof - p + 1) + 1;
                        continue;
                }

                NodeID node = G.new_node(); node_counter++;
                G.setPartitionIndex(node, 0);

                NodeWeight weight = 1;
                if( read_nw ) {
                        while(is_space_char(*p)) p++;
                        if(*p >= '0' && *p <= '9') {
                                unsigned long long w = 0;
                                while(*p >= '0' && *p <= '9') { w = w*10 + (*p - '0'); p++; }
                                weight = (NodeWeight) w;
                        }
                        total_nodeweight += weight;
                        if( total_nodeweight > (long long) std::numeric_limits<NodeWeight>::max()) {
                                std::cerr <<  "The sum of the node weights is too large (it exceeds the node weight type)."  << std::endl;
                                std::cerr <<  "Currently not supported. Please scale your node weights."  << std::endl;
                                exit(0);
                        }
                }
                G.setNodeWeight(node, weight);

                for(;;) {
                        while(is_space_char(*p)) p++;
                        if(*p < '0' || *p > '9') break;
                        NodeID target = 0;
                        while(*p >= '0' && *p <= '9') { target = target*10 + (*p - '0'); p++; }

                        //check for self-loops
                        if(target-1 == node) {
                                std::cerr <<  "The graph file contains self-loops. This is not supported. Please remove them from the file."  << std::endl;
                        }

                        EdgeWeight edge_weight = 1;
                        if( read_ew ) {
                                while(is_space_char(*p)) p++;
                                if(*p >= '0' && *p <= '9') {
                                        EdgeWeight w = 0;
                                        while(*p >= '0' && *p <= '9') { w = w*10 + (*p - '0'); p++; }
                                        edge_weight = w;
                                }
                        }
                        edge_counter++;
                        if( edge_counter > (EdgeID) nmbEdges ) {
                                continue; // counted, rejected below
                        }
                        EdgeID e = G.new_edge(node, target-1);

                        G.setEdgeWeight(e, edge_weight);
                }
                // skip the rest of the line
                if(*p != '\n') p = (const char*)memchr(p, '\n', eof - p + 1);
                p++;
        }

        if( edge_counter != (EdgeID) nmbEdges ) {
                std::cerr <<  "number of specified edges mismatch"  << std::endl;
                std::cerr <<  edge_counter <<  " " <<  nmbEdges  << std::endl;
                exit(0);
        }

        if( node_counter != (NodeID) nmbNodes) {
                std::cerr <<  "number of specified nodes mismatch"  << std::endl;
                std::cerr <<  node_counter <<  " " <<  nmbNodes  << std::endl;
                exit(0);
        }


        G.finish_construction();
        return 0;
}


void graph_io::writePartition(graph_access & G, const std::string & filename) {
        FILE * fp = fopen(filename.c_str(), "wb");
        std::cout << "writing partition to " << filename << " ... " << std::endl;
        if(!fp) return;

        std::unique_ptr<char[]> buf(new char[12 * (size_t)G.number_of_nodes() + 16]); // not zeroed
        char * out = buf.get();
        forall_nodes(G, node) {
                PartitionID v = G.getPartitionIndex(node);
                char tmp[12];
                int len = 0;
                do { tmp[len++] = '0' + (v % 10); v /= 10; } while(v);
                while(len) *out++ = tmp[--len];
                *out++ = '\n';
        } endfor
        fwrite(buf.get(), 1, out - buf.get(), fp);
        fclose(fp);
}


