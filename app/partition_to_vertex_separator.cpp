/******************************************************************************
 * partition_to_vertex_separator.cpp 
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
#include <iostream>
#include <sstream>
#include <string.h> 

#include <argtable2.h>
#include <regex.h>

#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "tools/timer.h"
#include "tools/quality_metrics.h"
#include "tools/macros_assertions.h"
#include "tools/random_functions.h"
#include "partition/partition_config.h"
#include "partition/graph_partitioner.h"
#include "partition/uncoarsening/separator/vertex_separator_algorithm.h"

#include "parse_parameters.h"
#include <math.h>

using namespace std;
int main(int argn, char **argv) {

        PartitionConfig partition_config;
        string graph_filename;

        bool is_graph_weighted = false;
        bool suppress_output   = false;
        bool recursive         = false;
       
        int ret_code = parse_parameters(argn, argv, 
                                        partition_config, graph_filename, 
                                        is_graph_weighted, suppress_output, 
                                        recursive); 

        if(ret_code) {
                return 0;
        }

        streambuf* backup = cout.rdbuf();
        ofstream ofs;
        ofs.open("/dev/null");

        if(suppress_output) {
               cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.LogDump(stdout);

        graph_access G;     
        timer t;
        graph_io::readGraphWeighted(G, graph_filename);
        cout << "io time: " << t.elapsed()  << endl;
       
        G.set_partition_count(partition_config.k); 
 
        balance_configuration bc;
        bc.configurate_balance( partition_config, G);

        if(partition_config.input_partition != "") {
                cout <<  "reading input partition" << endl;
                graph_io::readPartition(G, partition_config.input_partition);
                partition_config.graph_allready_partitioned = true;
        } else {
                cout <<  "please specify an input partition"  << endl;
                exit(0);
        }

        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        cout <<  "graph has " <<  G.number_of_nodes() <<  " nodes and " <<  G.number_of_edges() <<  " edges"  << endl;
        // ***************************** perform conversion ***************************************       
        t.restart();
        quality_metrics qm;
        complete_boundary boundary(&G);
        boundary.build();

        vertex_separator_algorithm vsa;
        vector<NodeID> separator;
        vsa.compute_vertex_separator(partition_config, G, boundary, separator);
        cout << "separatorsize " << separator.size() << endl;
        cout << "cut \t\t"       << qm.edge_cut(G)   << endl;
        // ********************************** done ************************************************       
        
        ofs.close();
        cout.rdbuf(backup);
        cout <<  "time spent for conversion " << t.elapsed()  << endl;
       
        // write the partition to the disc 
        for( unsigned i = 0; i < separator.size(); i++) {
                G.setPartitionIndex(separator[i], partition_config.k);
        }

         // write the partition to the disc 
        std::stringstream filename;
        if(!partition_config.filename_output.compare("")) {
                // no output filename given
                filename << "tmpseparator" << partition_config.k;
        } else {
                filename << partition_config.filename_output;
        }

        graph_io::writePartition(G, filename.str());
}
