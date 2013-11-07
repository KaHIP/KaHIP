/******************************************************************************
 * graph_io.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
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
        char* line = new char[MAXLINE+1];

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << filename << std::endl;
                return 1;
        }

        NodeID nmbNodes;
        EdgeID nmbEdges;

        in.getline(line, MAXLINE);
        //skip comments
        while( line[0] == '%' ) {
                in.getline(line, MAXLINE);
        }

        int ew = 0;
        sscanf(line, "%d %d %d", &nmbNodes, &nmbEdges, &ew );

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

        
        G.start_construction(nmbNodes, nmbEdges);
        char *oldstr, *newstr;

        int edge_counter = 0;
        for (NodeID i = 0; i < nmbNodes; ++i) {
                NodeID node = G.new_node();

                G.setPartitionIndex(node, 0);
                in.getline(line, MAXLINE);

                //******************************************
                // insert the edges given in this line
                // Parse the Line
                // *****************************************
                oldstr = line;
                newstr = NULL;

                // First Number in this Line is the Nodeweight
                if(read_nw) {
                        int weight = (NodeID) strtol(oldstr, &newstr, 10);
                        oldstr = newstr;
                        G.setNodeWeight(node, weight);
                } else {
                        G.setNodeWeight(node, 1);
                }

                NodeID source = i;
                NodeID target = 0;
                int edgeWeight = 0;
                for (;;) {
                        target = (NodeID) strtol(oldstr, &newstr, 10);
                        oldstr = newstr;

                        if (target == 0)
                                break;

                        if(read_ew) {
                                edgeWeight = (NodeID) strtol(oldstr, &newstr, 10);
                                oldstr = newstr;
                        }

                        target -= 1; // -1 since there are no nodes with id 0 in the file

                        ASSERT_NEQ(source, target);
                        ASSERT_LEQ(target, nmbNodes);

                        EdgeID e = G.new_edge(source, target);
                        if(read_ew) {
                                G.setEdgeWeight(e, edgeWeight);
                        } else {
                                G.setEdgeWeight(e, 1);
                        }
                        edge_counter++;
                }
        }
        std::cout <<  "edge conter " <<  edge_counter  << std::endl;

        G.finish_construction();
        delete[] line;
        return 0;
}

int graph_io::readGraphUnweighted(graph_access & G, std::string filename) {
        char* line = new char[MAXLINE+1];

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening unweighted graph " << filename << std::endl;
                return 1;
        }

        NodeID nmbNodes;
        EdgeID nmbEdges;

        in.getline(line, MAXLINE);
        //skip comments
        while( line[0] == '%' ) {
                in.getline(line, MAXLINE);
        }

        sscanf(line, "%d %d", &nmbNodes, &nmbEdges);

        nmbEdges *= 2; //since we have forward and backward edges

        G.start_construction(nmbNodes, nmbEdges);
        char *oldstr, *newstr;

        for (NodeID i = 0; i < nmbNodes; ++i) {
                NodeID node = G.new_node();
                G.setPartitionIndex(node, 0);
                G.setNodeWeight(node, 1);

                //Fetch current line
                in.getline(line, MAXLINE);

                //******************************************
                // insert the edges given in this line
                // Parse the Line
                // *****************************************
                oldstr = line;
                newstr = NULL;
                NodeID source = i;
                NodeID target = 0;

                for (;;) {
                        target = (NodeID) strtol(oldstr, &newstr, 10);

                        if (target == 0) {
                                break;
                        }

                        oldstr = newstr; target -= 1; // -1 since there are no nodes with id 0 in the file

                        ASSERT_NEQ(source, target);
                        ASSERT_LEQ(target, nmbNodes);

                        EdgeID e = G.new_edge(source, target);
                        G.setEdgeWeight(e, 1);
                }
        }

        G.finish_construction();
        delete[] line;
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


