/******************************************************************************
 * graphchecker.cpp 
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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <unordered_set>

using namespace std;

// this program implements the functions to check the metis graph 
// format
int main(int argn, char **argv)
{

        if( argn != 2 ) {
                std::cout <<  "Usage: graphchecker FILE"  << std::endl;
                exit(0);
        }

        std::string line;
        std::string filename(argv[1]);

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening " << filename << std::endl;
                return 1;
        }

        std::cout <<  "*******************************************************************************"  << std::endl;
        std::cout <<  "KaHIP -- graph format checker."  << std::endl;
        std::cout <<  "Output will be given using the IDs from file, i.e. the IDs are starting from 1."  << std::endl;
        std::cout <<  "*******************************************************************************"  << std::endl;

        std::getline(in,line);
        //skip comments
        while( line[0] == '%' ) {
                std::getline(in, line);
        }

        long nmbNodes;
        long nmbEdges;
        long ew = 0;

        std::stringstream ss(line);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;


        std::vector< long > node_starts;
        node_starts.reserve(nmbNodes + 1);
        node_starts.push_back(0);

        std::vector< long > adjacent_nodes;
        adjacent_nodes.reserve(nmbEdges * 2);

        std::vector< long > graph_edgeweights;
        graph_edgeweights.reserve(nmbEdges * 2);
        long node_weight;
        long total_nodeweight = 0;
        long total_edgeweight = 0;

        long node_degree = 0;
        long node_counter = 0;
        long edge_counter = 0;

        bool node_weights = false;
        bool edge_weights = false;
        if( ew == 11 ) {
                node_weights = true;
                edge_weights = true;
        } else if ( ew == 10) {
                node_weights = true;
        } else if ( ew == 1 ) {
                edge_weights = true;
        }

        while( std::getline(in, line) ) {
                
                if( node_counter > nmbNodes ) {
                        std::cout <<  "There are more nodes in the file than specified in the first line of the file."  << std::endl;
                        std::cout <<  "You specified " <<  nmbNodes << " nodes." << std::endl;
                        std::cout <<  node_counter  << std::endl;
                        std::cout <<  "*******************************************************************************"  << std::endl;
                        exit(0);
                }

                if (line[0] == '%') { // a comment in the file
                        continue;
                }
                node_degree = 0;

                std::stringstream ss(line);
                if( node_weights ) {
                        ss >> node_weight;
                        if( node_weight < 0 ) {
                            std::cout <<  "The node " <<  node_counter+1 << " has weight < 0."  << std::endl;
                            std::cout <<  "See line " << node_counter+2 << " of your file."  << std::endl;
                            std::cout <<  "*******************************************************************************"  << std::endl;
                            exit(0);
                        }
                        total_nodeweight += node_weight;
                        if(total_nodeweight > (long)std::numeric_limits<unsigned int>::max()) {
                                std::cout <<  "The sum of the node weights exeeds 32 bits. Currently not supported."  << std::endl;
                                std::cout <<  "Please scale weights of the graph."  << std::endl;
                                std::cout <<  "*******************************************************************************"  << std::endl;
                                exit(0);
                        }
                }

                long target;
                while( ss >> target ) {
                        node_degree++;
                        if( target > nmbNodes || target <= 0 ) {
                                 std::cout <<  "Node " << node_counter+1 << " has an edge to a node greater than the number of nodes specified in the file or smaller or equal to zero, i.e. it has target " <<  target << " and the number of nodes specified was " <<  nmbNodes << std::endl;
                                 std::cout <<  "See line " << node_counter+2 << " of your file."  << std::endl;
                                 std::cout <<  "*******************************************************************************"  << std::endl;
                                 exit(0);
                        }
                        adjacent_nodes.push_back(target - 1);
                        long edge_weight = 1;
                        if( edge_weights ) {
                                if( ss.eof() ) {
                                        std::cout <<  "Something is wrong."  << std::endl;
                                        std::cout <<  "See line " << node_counter+2 << " of your file."  << std::endl;
                                        std::cout <<  "There is not the right amount of numbers in line " << node_counter+2 << " of the file. " << std::endl;
                                        if( node_weights ) {
                                                std::cout <<  "That means either the node weight is missing, "
                                                          <<  "or there is an edge without a weight specified, " 
                                                          <<  "or there are no edge weights at all despite the specification " 
                                                          <<  ew  << " in the first line of the file."<< std::endl;
                                                std::cout <<  "*******************************************************************************"  << std::endl;
                                        } else {
                                                std::cout <<  "That there is may be an edge without an edge weight specified "
                                                          <<  "or there are no edge weights at all despite the specification " 
                                                          <<  ew  << " in the first line of the file."<< std::endl;
                                                std::cout <<  "*******************************************************************************"  << std::endl;

                                        }
                                        exit(0);

                                }
                                ss >> edge_weight;
                                total_edgeweight += edge_weight;

                                if(total_edgeweight > (long)std::numeric_limits<unsigned int>::max()) {
                                        std::cout <<  "The sum of the edge weights exeeds 32 bits. Currently not supported."  << std::endl;
                                        std::cout <<  "Please scale weights of the graph."  << std::endl;
                                        std::cout <<  "*******************************************************************************"  << std::endl;
                                        exit(0);
                                }

                                if( edge_weight <= 0 ) {
                                        std::cout <<  "The edge starting from node " <<  (node_counter+1) << " and ending in node " << target  
                                                  <<  " has weight <= 0. " << std::endl;
                                        std::cout <<  "See line " << node_counter+2 << " of your file."  << std::endl;
                                        std::cout <<  "*******************************************************************************"  << std::endl;
                                        exit(0);
                                }
                        }
                        graph_edgeweights.push_back(edge_weight);
                        edge_counter++;
                }
                node_counter++;
                node_starts.push_back(node_starts.back() + node_degree);
                if(in.eof()) {
                        break;
                }
        }
        std::cout <<  "IO done. Now checking the graph .... "  << std::endl;


        // check node counter
        if( node_counter != nmbNodes ) {
                std::cout <<  "The number of nodes specified in the beginning of the file "
                          <<  "does not match the number of nodes that are in the file."  << std::endl;
                std::cout <<  "You specified " <<  nmbNodes <<  " but there are " <<  node_counter  << std::endl;
                std::cout <<  "*******************************************************************************"  << std::endl;
                exit(0);
        }

        // check edge counter
        if( edge_counter != 2*nmbEdges || node_starts.back() != 2*nmbEdges) {
                std::cout <<  "The number of edges specified in the beginning of the file " 
                          <<  "does not match the number of edges that are in the file."  << std::endl;
                std::cout <<  "You specified " <<  2*nmbEdges <<  " but there are " <<  edge_counter << std::endl;
                std::cout <<  "*******************************************************************************"  << std::endl;
                exit(0);
        }


        // check if there are parallel edges
        for( long node = 0; node < nmbNodes; node++) {
                std::unordered_set< long > seen_adjacent_nodes;
                for( long e = node_starts[node]; e < node_starts[node + 1]; e++) {
                        long target = adjacent_nodes[e];
                        if( !seen_adjacent_nodes.insert(target).second ) {
                                std::cout <<  "The file contains parallel edges."  << std::endl;
                                std::cout <<  "In line " <<  node+2 << " of the file " <<  (target+1) << " is listed twice."   << std::endl;
                                std::cout <<  "*******************************************************************************"  << std::endl;
                                exit(0);
                        }

                        // check for self loops
                        if( target == node ) {
                                std::cout <<  "The file contains a graph with self-loops."  << std::endl;
                                std::cout <<  "In line " <<  node+2 << " of the file (node=" 
                                          << node+1 << ") the target " <<  (target+1) << " is listed."   << std::endl;
                                std::cout <<  "*******************************************************************************"  << std::endl;
                                exit(0);
                        }
                }
        }

        // if for backward edges and also check if the weights match
        for( long node = 0; node < nmbNodes; node++) {
                for( long e = node_starts[node]; e < node_starts[node + 1]; e++) {
                        long target = adjacent_nodes[e];
                        long forward_weight = graph_edgeweights[e];

                        bool found = false;
                        for( long e_bar = node_starts[target]; e_bar < node_starts[target + 1]; e_bar++) {
                                if( adjacent_nodes[e_bar] == node) {
                                        found = true;
                                        if( graph_edgeweights[e_bar] != forward_weight) {
                                                std::cout <<  "The file does not contain valid edge weights. " 
                                                          <<  "The weights of the forward edges must be equal " 
                                                          <<  "to the weight of the backward edges. "<< std::endl;
                                                std::cout <<  "Node " <<  node+1 << " does contain an arc to node " 
                                                          << target+1 << " with weight " << forward_weight 
                                                          <<  " but the weight of the backward edge (" <<  target+1 << "," << node+1 << ") " 
                                                          <<  " is  " << graph_edgeweights[e_bar] << std::endl;
                                                std::cout <<  "You can find the backward edge in line " << target+2 << " of the file."  << std::endl;
                                                std::cout <<  "You can find the forward edge in line " << node+2 << " of the file."  << std::endl;
                                                std::cout <<  "*******************************************************************************"  << std::endl;
                                                exit(0);
                                        } 
                                        break;
                                }
                        }
                        if(!found) {
                                std::cout <<  "The file does not contain all forward and backward edges. "  << std::endl;
                                std::cout <<  "Node " <<  node+1 << " (line " << node+2 
                                          << ") does contain an arc to node " << target+1 << " but there is no edge ("
                                          <<  target+1 << "," << node+1 << ") in the file. "<<  std::endl;
                                std::cout <<  "Please insert this edge in line " << target+2 << " of the file." << std::endl;
                                std::cout <<  "*******************************************************************************"  << std::endl;
                                exit(0);
                        }
                }
        }

        std::cout <<  "The graph format seems correct."  << std::endl;
        std::cout <<  "*******************************************************************************"  << std::endl;


        return 0;
}

