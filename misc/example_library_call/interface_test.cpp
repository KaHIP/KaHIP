/******************************************************************************
 * kaffpa.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#include <iostream>
#include <sstream>

#include "kaHIP_interface.h"


int main(int argn, char **argv) {

        std::cout <<  "partitioning graph from the manual"  << std::endl;

        int n            = 5;
        int* xadj        = new int[6];
        xadj[0] = 0; xadj[1] = 2; xadj[2] = 5; xadj[3] = 7; xadj[4] = 9; xadj[5] = 12;

        int* adjncy      = new int[12];
        adjncy[0]  = 1; adjncy[1]  = 4; adjncy[2]  = 0; adjncy[3]  = 2; adjncy[4]  = 4; adjncy[5]  = 1; 
        adjncy[6]  = 3; adjncy[7]  = 2; adjncy[8]  = 4; adjncy[9]  = 0; adjncy[10] = 1; adjncy[11] = 3; 
        
        double imbalance = 0.03;
        int* part        = new int[n];
        int edge_cut     = 0;
        int nparts       = 2;
        int* vwgt        = NULL;
        int* adjcwgt     = NULL;

        //void kaffpa(int* n, int* vwgt, int* xadj, 
                   //int* adjcwgt, int* adjncy, int* nparts, 
                   //double* imbalance,  bool suppress_output, int seed, int mode, 
                   //int* edgecut, int* part);

        kaffpa(&n, vwgt, xadj, adjcwgt, adjncy, &nparts, &imbalance, false, 0, ECO, & edge_cut, part);

        std::cout <<  "edge cut " <<  edge_cut  << std::endl;

        //void process_mapping(int* n, int* vwgt, int* xadj, 
                   //int* adjcwgt, int* adjncy, 
                   //int* hierarchy_parameter,  int* distance_parameter, int hierarchy_depth, 
                   //int mode_partitioning, int mode_mapping,
                   //double* imbalance,  
                   //bool suppress_output, int seed,
                   //int* edgecut, int* qap, int* part); 

        int* hierarchy_parameter = new int[2];
        int* distance_parameter = new int[2];
        hierarchy_parameter[0] = 2;
        hierarchy_parameter[1] = 2;
        distance_parameter[0] = 1;
        distance_parameter[1] = 100;
        int qap = 0;

        process_mapping(&n, vwgt, xadj, adjcwgt, adjncy, hierarchy_parameter, distance_parameter, 2, STRONG, MAPMODE_MULTISECTION, &imbalance, false, 0, & edge_cut, & qap, part);

        std::cout <<  "edge cut " <<  edge_cut  << std::endl;
        std::cout <<  "qap " <<  qap << std::endl;
                
}
