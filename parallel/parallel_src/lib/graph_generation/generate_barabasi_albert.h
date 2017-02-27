/******************************************************************************
 * generate_barabasi_albert.h
 *
 * Source of the Parallel Partitioning Program
 ******************************************************************************
 * Copyright (C) 2014 Christian Schulz <christian.schulz@kit.edu>
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


#ifndef GENERATE_BERTASSIALBERT_PLHS3WMW
#define GENERATE_BERTASSIALBERT_PLHS3WMW

#include <x86intrin.h>
#include <algorithm>
#include <iomanip>
#include <unordered_map>
#include <sstream>
#include "communication/mpi_tools.h"
#include "data_structure/parallel_graph_access.h"
#include "data_structure/linear_probing_hashmap_ll.h"
#include "io/parallel_graph_io.h"
#include "macros_assertions.h"
#include "util_rand.h"
#include "partition_config.h"
#include "helpers.h"
#include "build_from_edge_list.h"

#if __GNUC__
  #if __x86_64__ || __ppc64__
    #define ENV64BIT
  #else
    #define ENV32BIT
  #endif
#endif

#define UPPER_MASK 0xffff0000 /* most significant w-r bits */
#define LOWER_MASK 0x0000ffff /* least significant r bits */

#define SEEDA 28475421
#define SEEDB 52150599

using namespace distributed_graph;

class generate_barabasialbert {
public:
        generate_barabasialbert() ;
        virtual ~generate_barabasialbert() ;

        void generate_k_deghist(PPartitionConfig & config, parallel_graph_access & G, bool n_given = false, bool create_graph = false) {
                timer t;
                ULONG n = 0;
                if( n_given ) {
                        n = config.n;
                } else {
                        n = (ULONG) pow(2, config.log_num_verts);
                }

                int d      = (int) config.barabasi_albert_mindegree;
                ULONG twod = 2*d;

                PEID rank, size;
                MPI_Comm_rank( MPI_COMM_WORLD, &rank);
                MPI_Comm_size( MPI_COMM_WORLD, &size);
                
                ULONG from  = rank * ceil(n / (double)size);
                ULONG to    = (rank+1) * ceil(n / (double)size) - 1;
                to = std::min(to, n-1); // boundary case

                ULONG * degrees = new ULONG [config.k_deg];
                ULONG * global_degrees = new ULONG [config.k_deg];
                for( ULONG i = 0; i < config.k_deg; i++) {
                        degrees[i] = 0;
                        global_degrees[i] = 0;
                }

                MPI_Barrier(MPI_COMM_WORLD);
                
                t.restart();

                for( ULONG v = from ; v <= to; v++) {
                        for( int i = 0; i < d; i++) {
                                ULONG r = 2*(v*d + i)+1;
                                do {
                                        //compute hash h(r)
#ifdef ENV64BIT 
                                        ULONG hash  = _mm_crc32_u64(SEEDA, r);
                                        hash        = hash << 32; 
                                        hash        += _mm_crc32_u64(SEEDB, r);
#else
                                        ULONG hash  = 0;
                                        hash += _mm_crc32_u32(SEEDA, r & UPPER_MASK );
                                        hash += _mm_crc32_u32(SEEDA, r & LOWER_MASK);
                                        hash  = hash << 32; 
                                        hash += _mm_crc32_u32(SEEDB, r & UPPER_MASK );
                                        hash += _mm_crc32_u32(SEEDB, r & LOWER_MASK);
 #endif

                                        r = hash % r;
                                } while( r % 2 == 1 );
                                if( v < config.k_deg ) {
                                        degrees[v]++;
                                }
                                ULONG w = r / twod;
                                if( w < config.k_deg ) {
                                        degrees[w]++;
                                } 
                        }
                }

                MPI_Barrier(MPI_COMM_WORLD);
                
                if( rank == ROOT ) {
                        double elapsed = t.elapsed();
                        std::cout <<  "generation took " <<  elapsed  << std::endl;
                        std::cout <<  "total number of edges generated " <<  n*d  << std::endl;

                        long double terra_bytes = 128*n*d;
                        terra_bytes /= (8);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);

                        std::cout <<  "memory of graph in tera bytes  " << std::setprecision(40) << terra_bytes  << std::endl;
                        std::cout <<  "edges per second " << std::setprecision(40) <<  (n*d)/elapsed  << std::endl;
                        std::cout <<  "edges per (second*proc) "<< std::setprecision(40) <<  (n*d)/(size*elapsed)  << std::endl;
                        std::cout <<  "Medges/s " <<  std::setprecision(40) << (n*d)/ (elapsed) * 1.e-6<< std::endl;
                        std::cout <<  "num pes used " <<  size  << std::endl;
                        std::cout <<  "n is set to " <<  n << std::endl;
                        std::cout <<  "d is set to " <<  d << std::endl;
                }

                MPI_Reduce(degrees, global_degrees, config.k_deg, MPI_UNSIGNED_LONG_LONG, MPI_SUM, ROOT, MPI_COMM_WORLD); 
                MPI_Barrier(MPI_COMM_WORLD);
                
                if( rank == ROOT ) {
                        std::cout <<  "computing degree sequence took " <<  t.elapsed()  << std::endl;
                }

                if( rank == ROOT ) {
                        std::stringstream filename;
                        filename << "degrees_"<< config.k_deg<<"_logn_" << config.log_num_verts << "_d_" << d << "_s0n1_pes_" << size;
                        std::ofstream f(filename.str().c_str());
                        for( ULONG i = 0; i < config.k_deg+1; i++) {
                                if( global_degrees[i] != 0 ) {
                                        f <<   i <<  " " <<  global_degrees[i] << std::endl;

                                }
                        }
                        f.close();
                }

                MPI_Barrier(MPI_COMM_WORLD);
        }

        void generate_withgraph(PPartitionConfig & config, parallel_graph_access & G, bool n_given = false) {
                timer t;
                ULONG n = 0;
                if( n_given ) {
                        n = config.n;
                } else {
                        n = (ULONG) pow(2, config.log_num_verts);
                }
 
                int d      = (int) config.barabasi_albert_mindegree;
                ULONG twod = 2*d;

                PEID rank, size;
                MPI_Comm_rank( MPI_COMM_WORLD, &rank);
                MPI_Comm_size( MPI_COMM_WORLD, &size);
                
                ULONG from  = rank * ceil(n / (double)size);
                ULONG to    = (rank+1) * ceil(n / (double)size) - 1;
                to = std::min(to, n-1); // boundary case

                ULONG local_no_nodes = to - from + 1;
                std::cout <<  "rank " <<  rank <<  " generates nodes from " <<  from <<  " to " <<  to  <<  ". Total amount " <<  local_no_nodes << "." << std::endl;

                std::vector< source_target_pair > edge_list;
                MPI_Barrier(MPI_COMM_WORLD);
                
                t.restart();

                for( ULONG v = from ; v <= to; v++) {
                        for( int i = 0; i < d; i++) {
                                ULONG r = 2*(v*d + i)+1;
                                do {
                                        //compute hash h(r)
#ifdef ENV64BIT 
                                        ULONG hash  = _mm_crc32_u64(SEEDA, r);
                                        hash        = hash << 32; 
                                        hash        += _mm_crc32_u64(SEEDB, r);
#else
                                        ULONG hash  = 0;
                                        hash += _mm_crc32_u32(SEEDA, r & UPPER_MASK);
                                        hash += _mm_crc32_u32(SEEDA, r & LOWER_MASK);
                                        hash  = hash << 32; 
                                        hash += _mm_crc32_u32(SEEDB, r & UPPER_MASK);
                                        hash += _mm_crc32_u32(SEEDB, r & LOWER_MASK);
 #endif
                                        r = hash % r;
                                } while( r % 2 == 1 );
                                
                                source_target_pair stp;
                                stp.source = v;
                                stp.target = r / twod;

                                source_target_pair bstp;
                                bstp.source = r / twod;
                                bstp.target = v;

                                edge_list.push_back(stp);
                                edge_list.push_back(bstp);
                        }
                }

                MPI_Barrier(MPI_COMM_WORLD);
       
                if( rank == ROOT ) {
                        double elapsed = t.elapsed();
                        std::cout <<  "generation took " <<  elapsed  << std::endl;
                        std::cout <<  "total number of edges generated " <<  n*d  << std::endl;

                        long double terra_bytes = 128*n*d;
                        terra_bytes /= (8);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);

                        std::cout <<  "memory of graph in tera bytes  " << std::setprecision(40) << terra_bytes  << std::endl;

                        std::cout <<  "edges per second " << std::setprecision(40) <<  (n*d)/elapsed  << std::endl;
                        std::cout <<  "edges per (second*proc) "<< std::setprecision(40) <<  (n*d)/(size*elapsed)  << std::endl;
                        std::cout <<  "Medges/s " <<  std::setprecision(40) << (n*d)/ (elapsed) * 1.e-6<< std::endl;
                        std::cout <<  "num pes used " <<  size  << std::endl;
                        std::cout <<  "n is set to " <<  n << std::endl;
                        std::cout <<  "d is set to " <<  d << std::endl;
                        std::cout <<  "now remapping "  << std::endl;
                }

                // remap vertex ids so that the edges are fairly well distributed
                NodeID divisor = ceil( n / (double)size );
                NodeID val = size * ceil( n / (double)size );
                for( source_target_pair stp : edge_list ) {
                        stp.source = (stp.source * divisor) % val + ceil( stp.source / (double)size);
                        stp.target = (stp.target * divisor) % val + ceil( stp.target/ (double)size);
                }

                build_from_edge_list().build_graph_from_edge_list( config, G, edge_list);
                parallel_graph_io().writeGraphParallelSimple( G, "outputba");
                MPI_Barrier(MPI_COMM_WORLD);
                if( rank == ROOT ) std::cout <<  "generation done"  << std::endl;
        }

        void generate(PPartitionConfig & config, parallel_graph_access & G, bool n_given = false, bool create_graph = false) {
                timer t;
                ULONG n = 0;
                if( n_given ) {
                        n = config.n;
                } else {
                        n = (ULONG) pow(2, config.log_num_verts);
                }
 
                int d      = (int) config.barabasi_albert_mindegree;
                ULONG twod = 2*d;

                PEID rank, size;
                MPI_Comm_rank( MPI_COMM_WORLD, &rank);
                MPI_Comm_size( MPI_COMM_WORLD, &size);
                
                ULONG from  = rank * ceil(n / (double)size);
                ULONG to    = (rank+1) * ceil(n / (double)size) - 1;
                to = std::min(to, n-1); // boundary case

                ULONG local_no_nodes = to - from + 1;
                std::cout <<  "rank " <<  rank <<  " generates nodes from " <<  from <<  " to " <<  to  <<  ". Total amount " <<  local_no_nodes << "." << std::endl;

                ULONG * edge_heads = new ULONG [2*local_no_nodes*d];
                MPI_Barrier(MPI_COMM_WORLD);
                
                t.restart();

                for( ULONG v = from ; v <= to; v++) {
                        for( int i = 0; i < d; i++) {
                                ULONG r = 2*(v*d + i)+1;
                                do {
                                        //compute hash h(r)
#ifdef ENV64BIT 
                                        ULONG hash  = _mm_crc32_u64(SEEDA, r);
                                        hash        = hash << 32; 
                                        hash        += _mm_crc32_u64(SEEDB, r);
#else
                                        ULONG hash  = 0;
                                        hash += _mm_crc32_u32(SEEDA, r & UPPER_MASK);
                                        hash += _mm_crc32_u32(SEEDA, r & LOWER_MASK);
                                        hash  = hash << 32; 
                                        hash += _mm_crc32_u32(SEEDB, r & UPPER_MASK);
                                        hash += _mm_crc32_u32(SEEDB, r & LOWER_MASK);
 #endif
                                        r = hash % r;
                                } while( r % 2 == 1 );
                                
                                ULONG pos = 2*(d*(v-from) + i);
                                edge_heads[pos]   = v; 
                                edge_heads[pos+1] = r / twod; // E union {v, r div 2d}
                        }
                }

                MPI_Barrier(MPI_COMM_WORLD);
                
                if( rank == ROOT ) {
                        double elapsed = t.elapsed();
                        std::cout <<  "generation took " <<  elapsed  << std::endl;
                        std::cout <<  "total number of edges generated " <<  n*d  << std::endl;

                        long double terra_bytes = 128*n*d;
                        terra_bytes /= (8);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);

                        std::cout <<  "memory of graph in tera bytes  " << std::setprecision(40) << terra_bytes  << std::endl;

                        std::cout <<  "edges per second " << std::setprecision(40) <<  (n*d)/elapsed  << std::endl;
                        std::cout <<  "edges per (second*proc) "<< std::setprecision(40) <<  (n*d)/(size*elapsed)  << std::endl;
                        std::cout <<  "Medges/s " <<  std::setprecision(40) << (n*d)/ (elapsed) * 1.e-6<< std::endl;
                        std::cout <<  "num pes used " <<  size  << std::endl;
                        std::cout <<  "n is set to " <<  n << std::endl;
                        std::cout <<  "d is set to " <<  d << std::endl;
                }
                if(config.compute_degree_sequence_ba) {
                        if( rank == ROOT ) std::cout <<  "now computeing degree sequence"  << std::endl;
                        compute_degree_sequence_fast_lcnt( local_no_nodes, n, d, edge_heads, true);
                } else {
                        delete[] edge_heads;            
                }
        }

        void generate_32bit(PPartitionConfig & config, parallel_graph_access & G, bool n_given = false, bool create_graph = false) {
                timer t;
                UINT n = 0;
                if( n_given ) {
                        n = config.n;
                } else {
                        n = (UINT) pow(2, config.log_num_verts);
                }

                int d      = (int) config.barabasi_albert_mindegree;
                UINT twod = 2*d;

                PEID rank, size;
                MPI_Comm_rank( MPI_COMM_WORLD, &rank);
                MPI_Comm_size( MPI_COMM_WORLD, &size);
 
                ULONG total = (ULONG)n*(ULONG)d;
                if(total > std::numeric_limits< UINT >::max() ) {
                        if(rank == ROOT) std::cout <<  "total number of edges has to be smaller than max int"  << std::endl;
                        MPI_Finalize();
                } 

               
                UINT from  = rank * ceil(n / (double)size);
                UINT to    = (rank+1) * ceil(n / (double)size) - 1;
                to = std::min(to, n-1); // boundary case

                UINT local_no_nodes = to - from + 1;
                std::cout <<  "rank " <<  rank <<  " generates nodes from " <<  from <<  " to " <<  to  <<  ". Total amount " <<  local_no_nodes << "." << std::endl;

                UINT * edge_heads = new UINT [2*local_no_nodes*d];
                MPI_Barrier(MPI_COMM_WORLD);
                
                t.restart();

                for( UINT v = from ; v <= to; v++) {
                        for( int i = 0; i < d; i++) {
                                UINT r = 2*(((UINT)v)*d + i)+1; 
                                do {
                                        //compute hash h(r)
                                        UINT hash  = _mm_crc32_u32(SEEDA, r);
                                        r = hash % r;
                                } while( r % 2 == 1 );
                                
                                UINT pos = 2*(d*(v-from) + i);
                                edge_heads[pos]   = v; 
                                edge_heads[pos+1] = r / twod; // E union {v, r div 2d}
                        }
                }

                MPI_Barrier(MPI_COMM_WORLD);
                
                if( rank == ROOT ) {
                        double elapsed = t.elapsed();
                        std::cout <<  "generation took " <<  elapsed  << std::endl;
                        std::cout <<  "total number of edges generated " <<  (ULONG)n*(ULONG)d  << std::endl;

                        long double terra_bytes = 64*(ULONG)n*(ULONG)d;
                        terra_bytes /= (8);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);
                        terra_bytes /= (1024);

                        std::cout <<  "memory of graph in tera bytes  " << std::setprecision(40) << terra_bytes  << std::endl;

                        std::cout <<  "edges per second " << std::setprecision(40) <<  ((ULONG)n*(ULONG)d)/elapsed  << std::endl;
                        std::cout <<  "edges per (second*proc) "<< std::setprecision(40) <<  ((ULONG)n*(ULONG)d)/(size*elapsed)  << std::endl;
                        std::cout <<  "Medges/s " <<  std::setprecision(40) << ((ULONG)n*(ULONG)d)/ (elapsed) * 1.e-6<< std::endl;
                        std::cout <<  "num pes used " <<  size  << std::endl;
                        std::cout <<  "n is set to " << (ULONG) n << std::endl;
                        std::cout <<  "d is set to " <<  (ULONG)d << std::endl;
                }
                delete[] edge_heads;            
        }


        void compute_degree_sequence( ULONG local_no_nodes, 
                                      ULONG global_no_nodes,
                                      unsigned int d,
                                      ULONG * edge_heads,
                                      bool print_output=false) {

                timer t; t.restart();

                PEID rank, size;
                MPI_Comm_rank( MPI_COMM_WORLD, &rank);
                MPI_Comm_size( MPI_COMM_WORLD, &size);

                // default initializor is 0
                std::unordered_map< ULONG, ULONG > local_degrees;
                for( ULONG pos = 0; pos < 2*local_no_nodes*d; pos+=2) {
                        local_degrees[edge_heads[pos]]++;
                        if( edge_heads[pos] != edge_heads[pos+1] ) {
                                local_degrees[edge_heads[pos+1]]++;
                        }
                }

                std::vector< std::vector< ULONG > >  messages; messages.resize(size);
                for( std::pair< ULONG, ULONG > p : local_degrees) {
                        PEID peID = p.first % size;
                        if( peID != rank ) {
                                messages[ peID ].push_back( p.first  );
                                messages[ peID ].push_back( p.second );
                        }
                }

                //// now flood the network
                for( PEID peID = 0; peID < size; peID++) {
                        if( peID != rank ) {
                                if( messages[peID].size() == 0 ){
                                        // length 1 encode no message
                                        messages[peID].push_back(0);
                                }

                                MPI_Request rq; 
                                MPI_Isend( &messages[peID][0], messages[peID].size(), MPI_UNSIGNED_LONG_LONG, peID, peID, MPI_COMM_WORLD, &rq);
                        }
                }

                PEID counter = 0;
                while( counter < size - 1) {
                        // wait for incomming message of an adjacent processor
                        MPI_Status st;
                        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

                        int message_length;
                        MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                        std::vector<ULONG> incmessage; incmessage.resize(message_length);

                        MPI_Status rst;
                        MPI_Recv( &incmessage[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, rank, MPI_COMM_WORLD, &rst); 
                        
                        counter++;

                        // now integrate the changes
                        if(message_length == 1) continue; // nothing to do

                        for( int i = 0; i < message_length-1; i+=2) {
                                local_degrees[incmessage[i]] += incmessage[i+1];
                        }
                }

                
                // now finding largest overall degree
                ULONG local_max_degree = 0;
                for( auto & p : local_degrees) {
                        if( p.second > local_max_degree ) {
                                local_max_degree = p.second;
                        }
                }
                ULONG global_max_degree = 0;
                MPI_Allreduce(&local_max_degree, &global_max_degree, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD); 
                ULONG* local_histogram  = new ULONG[global_max_degree+1];
                ULONG* global_histogram = new ULONG[global_max_degree+1];
                for( ULONG i = 0; i < global_max_degree+1; i++) {
                        local_histogram[i]  = 0;
                        global_histogram[i] = 0;
                }

                for( auto & p : local_degrees) {
                        PEID peID = p.first % size;
                        if( peID == rank ) { local_histogram[ p.second ]++; }

                }

                MPI_Allreduce(local_histogram, global_histogram, global_max_degree+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD); 
                MPI_Barrier(MPI_COMM_WORLD);
                
                if( rank == ROOT ) {
                        std::cout <<  "took " <<  t.elapsed()  << std::endl;
                }

                if( rank == ROOT && print_output ) {
                        for( ULONG i = 0; i < global_max_degree+1; i++) {
                                if( global_histogram[i] != 0 ) {
                                        std::cout <<  "deg " <<  i <<  " cnt " <<  global_histogram[i] << std::endl;
                                }
                        }
                }
        }

        // using real alltoall
        // using linear probing hash tables (two of them)
        // trying to free memory as soon as it is not needed anymore
        void compute_degree_sequence_fast_lcnt( ULONG local_no_nodes, 
                                                ULONG global_no_nodes,
                                                unsigned int d,
                                                ULONG * edge_heads,
                                                bool print_output = false) {

                timer t; t.restart();

                PEID rank, size;
                MPI_Comm_rank( MPI_COMM_WORLD, &rank);
                MPI_Comm_size( MPI_COMM_WORLD, &size);
                
                ULONG sbktsize[size];
                ULONG rbktsize[size];
                for( PEID i = 0; i < size; i++) {
                        sbktsize[i] = 0;
                        rbktsize[i] = 0;
                }

                // default initializor is 0
                linear_probing_hashmap_ll* internal_lph = new linear_probing_hashmap_ll;
                linear_probing_hashmap_ll& lph = *internal_lph;
                lph.init(1.5*local_no_nodes*d, 2);

                for( ULONG pos = 0; pos < 2*local_no_nodes*d; pos+=2) {
                        lph[edge_heads[pos]]++;
                        if( edge_heads[pos] != edge_heads[pos+1] ) {
                                lph[edge_heads[pos+1]]++;
                        }
                }
                delete[] edge_heads; // not needed anymore

                std::vector< KeyValuePair > * internal_hmap = lph.internal_access();
                for( KeyValuePair & p : (*internal_hmap) ) {
                        if(p.key != NOT_CONTAINED_LL ) {
                                PEID peID = p.key % size;
                                sbktsize[ peID ] += 2;
                        }
                }
                MPI_Alltoall(sbktsize, 1, MPI_UNSIGNED_LONG_LONG, rbktsize, 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);

                ULONG sdispl[size+1];
                ULONG rdispl[size+1];
                ULONG posptr[size+1];

                sdispl[0]=0;
                rdispl[0]=0;
                posptr[0]=0;

                // now compute displacements
                for( PEID i = 1; i <= size; i++) {
                        sdispl[i] = sdispl[i-1] + sbktsize[i-1];
                        rdispl[i] = rdispl[i-1] + rbktsize[i-1];
                        posptr[i] = sdispl[i];
                }

                ULONG* sdata = new ULONG[sdispl[size]];
                ULONG* rdata = new ULONG[rdispl[size]];

                for( KeyValuePair & p : (*internal_hmap) ) {
                        if(p.key != NOT_CONTAINED_LL ) {
                                PEID peID = p.key % size;
                                ULONG cur_pos    = posptr[peID]; 
                                sdata[cur_pos]   = p.key;
                                sdata[cur_pos+1] = p.value;
                                posptr[peID]    += 2;
                        }
                }

                delete internal_lph; // save some memory

                mpi_tools().alltoallv(sdata, sbktsize, sdispl, MPI_UNSIGNED_LONG_LONG, 
                                      rdata, rbktsize, rdispl, MPI_UNSIGNED_LONG_LONG);

                delete[] sdata; // not needed anymore

                linear_probing_hashmap_ll* internal_ndph = new linear_probing_hashmap_ll;
                linear_probing_hashmap_ll& ndph = *internal_ndph;
                ndph.init(rdispl[size], 2);

                ULONG pos = 0;
                for( PEID peID = 0; peID < (PEID)size; peID ++) {
                        for( ULONG j = 0; j < rbktsize[peID]; j+=2, pos+=2) {
                                ndph[rdata[pos]]+= rdata[pos+1];
                        }
                }

                delete[] rdata; //not needed anymore

                internal_hmap = ndph.internal_access();
                ULONG local_max_degree = 0;
                for( KeyValuePair & p : (*internal_hmap) ) {
                        if(p.key != NOT_CONTAINED_LL ) {
                                if( p.value > local_max_degree ) {
                                        local_max_degree = p.value;
                                }
                        }
                }

                // now finding largest overall degree
                ULONG global_max_degree = 0;
                MPI_Allreduce(&local_max_degree, &global_max_degree, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD); 
                ULONG* local_histogram  = new ULONG[global_max_degree+1];
                for( ULONG i = 0; i < global_max_degree+1; i++) { local_histogram[i]  = 0; }

                for( KeyValuePair & p : (*internal_hmap) ) {
                        if(p.key != NOT_CONTAINED_LL ) {
                                local_histogram[ p.value ]++; 
                        }
                }

                // free the hash table since it is not needed anymore
                delete internal_ndph;
                ULONG* global_histogram = new ULONG[global_max_degree+1];
                for( ULONG i = 0; i < global_max_degree+1; i++) { global_histogram[i] = 0; }

                MPI_Allreduce(local_histogram, global_histogram, global_max_degree+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD); 
                MPI_Barrier(MPI_COMM_WORLD);
                
                if( rank == ROOT ) {
                        std::cout <<  "computing degree sequence took " <<  t.elapsed()  << std::endl;
                }

                if( rank == ROOT && print_output ) {
                        std::stringstream filename;
                        filename << "degree_seq_logn_" << log2(global_no_nodes) << "_d_" << d << "_s0n1_pes_" << size;
                        std::ofstream f(filename.str().c_str());
                        for( ULONG i = 0; i < global_max_degree+1; i++) {
                                if( global_histogram[i] != 0 ) {
                                        f <<   i <<  " " <<  global_histogram[i] << std::endl;

                                }
                        }
                        f.close();
                }
                delete[] local_histogram;
                delete[] global_histogram;
        }

};


#endif /* end of include guard: GENERATE_RGG_PLHS3WMW */
