/******************************************************************************
 * matching.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "matching.h"

matching::matching() {

}

matching::~matching() {

}

void matching::print_matching(FILE * out, Matching & edge_matching) {
        for (NodeID n = 0; n < edge_matching.size(); n++) {
                fprintf(out, "%d:%d\n", n, edge_matching[n]);
        }        
}


