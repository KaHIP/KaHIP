/******************************************************************************
 * graph_io.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPHIO_H_
#define GRAPHIO_H_

#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "definitions.h"
#include "data_structure/graph_access.h"

class graph_io {
        public:
                graph_io();
                virtual ~graph_io () ;

                static
                int readGraphWeighted(graph_access & G, const std::string & filename);

                static
                int writeGraphWeighted(graph_access & G, const std::string & filename);

                static
                int writeGraph(graph_access & G, const std::string & filename);

                static
                int readPartition(graph_access& G, const std::string & filename);

                static
                void writePartition(graph_access& G, const std::string & filename);

                template<typename vectortype>
                static void writeVector(std::vector<vectortype> & vec, const std::string & filename);

                template<typename vectortype>
                static void readVector(std::vector<vectortype> & vec, const std::string & filename);


};

template<typename vectortype>
void graph_io::writeVector(std::vector<vectortype> & vec, const std::string & filename) {
        std::ofstream f(filename.c_str());
        for( unsigned i = 0; i < vec.size(); ++i) {
                f << vec[i] <<  std::endl;
        }

        f.close();
}

template<typename vectortype>
void graph_io::readVector(std::vector<vectortype> & vec, const std::string & filename) {

        std::string line;

        // open file for reading
        std::ifstream in(filename.c_str());
        if (!in) {
                std::cerr << "Error opening vectorfile" << filename << std::endl;
                return;
        }

        unsigned pos = 0;
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

#endif /*GRAPHIO_H_*/
