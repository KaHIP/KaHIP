/******************************************************************************
 * exchanger.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <mpi.h>

#include "exchanger.h"
#include "tools/quality_metrics.h"
#include "tools/random_functions.h"

exchanger::exchanger(MPI_Comm communicator) {
        m_prev_best_objective = std::numeric_limits<EdgeWeight>::max();

        m_communicator = communicator;

        int rank, comm_size;
        MPI_Comm_rank( m_communicator, &rank);
        MPI_Comm_size( m_communicator, &comm_size);

        m_cur_num_pushes = 0;
        if(comm_size > 2) m_max_num_pushes = ceil(log2(comm_size));
        else              m_max_num_pushes = 1;

        std::cout <<  "max num pushes " <<  m_max_num_pushes  << std::endl;

        m_allready_send_to.resize(comm_size);

        for( unsigned i = 0; i < m_allready_send_to.size(); i++) {
                m_allready_send_to[i] = false;
        }

        m_allready_send_to[rank] = true;
}

exchanger::~exchanger() {
        MPI_Barrier( m_communicator );
        int rank;
        MPI_Comm_rank( m_communicator, &rank);
        
        int flag; MPI_Status st;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_communicator, &flag, &st);
        
        while(flag) {
                int message_length;
                MPI_Get_count(&st, MPI_INT, &message_length);
                 
                int* partition_map = new int[message_length];
                MPI_Status rst;
                MPI_Recv( partition_map, message_length, MPI_INT, st.MPI_SOURCE, rank, m_communicator, &rst); 
                
                delete[] partition_map;
                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_communicator, &flag, &st);
        }

        MPI_Barrier( m_communicator );
        for( unsigned i = 0; i < m_request_pointers.size(); i++) {
                MPI_Cancel( m_request_pointers[i] );
        }
        
        for( unsigned i = 0; i < m_request_pointers.size(); i++) {
                MPI_Status st;
                MPI_Wait( m_request_pointers[i], & st );
                delete[] m_partition_map_buffers[i];
                delete   m_request_pointers[i];
        }
                
}

void exchanger::diversify_population( PartitionConfig & config, graph_access & G,  population & island, bool replace ) {
       
        int rank, comm_size;
        MPI_Comm_rank( m_communicator, &rank);
        MPI_Comm_size( m_communicator, &comm_size);
        
        std::vector<unsigned> permutation(comm_size, 0);

        if( rank == ROOT ) {
                random_functions::circular_permutation(permutation); 
        }

        MPI_Bcast(&permutation[0], comm_size, MPI_INT, ROOT, m_communicator);

        int from = 0;
        int to   = permutation[rank];
        for( unsigned i = 0; i < permutation.size(); i++) {
                if( permutation[i] == (unsigned)rank ) {
                        from = (int)i;
                        break;
                }
        }

        Individuum in;
        Individuum out;

        if(config.mh_diversify_best) {
                island.get_best_individuum(in);
        } else {
                island.get_random_individuum(in);
        }
        exchange_individum( config, G, from, rank, to, in, out);

        if( replace ) {
                island.replace( in, out );
        } else {
                island.insert( G, out );
        }

}

void exchanger::quick_start( PartitionConfig & config, graph_access & G, population & island ) {
        int comm_size;
        MPI_Comm_size( m_communicator, &comm_size);
        
        unsigned no_of_individuals = ceil(config.mh_pool_size / (double)comm_size) - 1;

        std::cout <<  "creating " <<  no_of_individuals << std::endl;

        for(unsigned i = 0; i < no_of_individuals; i++) {
                PartitionConfig copy            = config;
                copy.combine                    = false;
                copy.graph_allready_partitioned = false;

                Individuum ind;
                island.createIndividuum(config, G, ind, true);
                island.insert(G, ind);
        }

        int reps = config.mh_pool_size - no_of_individuals;
        if(reps < 0) reps = 0;

        PartitionConfig div_config   = config;
        div_config.mh_diversify_best = false;
        for( unsigned i = 0; i < (unsigned) reps; i++) {
                diversify_population( div_config , G, island, false); 
        }
}


void exchanger::exchange_individum( const PartitionConfig & config,  graph_access & G, 
                                    int & from, int & rank, int & to, 
                                    Individuum & in, Individuum & out) {
        //recv. edge cut, partition_map, cut_edges from "from"
        //send in to "to"

        int* partition_map = new int[G.number_of_nodes()];
        out.partition_map  = partition_map;
        out.cut_edges      = new std::vector<EdgeID>();

        MPI_Status st;
        MPI_Sendrecv( in.partition_map , G.number_of_nodes(), MPI_INT, to, 0, 
                      out.partition_map, G.number_of_nodes(), MPI_INT, from, 0, m_communicator, &st); 

        //recompute cut edges and edge cut locally
        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] != partition_map[target]) {
                                out.cut_edges->push_back(e);
                        }
                } endfor
        } endfor

        out.objective = m_qm.objective(config, G, partition_map);
}


//extended push protocol -- see paper for details
void exchanger::push_best( PartitionConfig & config, graph_access & G, population & island ) {
        int rank, size;
        MPI_Comm_rank( m_communicator, &rank);
        MPI_Comm_size( m_communicator, &size);

        Individuum best_ind;
        island.get_best_individuum(best_ind);

        if( best_ind.objective < m_prev_best_objective) {
                m_prev_best_objective = best_ind.objective;
                for( unsigned i = 0; i < m_allready_send_to.size(); i++) {
                        m_allready_send_to[i] = false;
                }
                
                m_allready_send_to[rank] = true;
                m_cur_num_pushes         = 0;

                std::cout << "rank " <<  rank 
                          << ": pool improved *************************************** " 
                          <<  best_ind.objective << std::endl;
        }

        bool something_todo = false;
        for( unsigned i = 0; i < m_allready_send_to.size(); i++) {
                if(!m_allready_send_to[i]) {
                      something_todo = true;
                      break;
                }
        }

        if( m_cur_num_pushes > m_max_num_pushes ) {
                something_todo = false;
        }

        if(something_todo) {
                int* partition_map = new int[G.number_of_nodes()];
                forall_nodes(G, node) {
                        partition_map[node] = G.getPartitionIndex(node);
                } endfor

                int target = rank;
                while( m_allready_send_to[target] ) {
                        //while (target == rank) { // m_allready_send_to[rank] always true
                                target = random_functions::nextInt(0, size-1);
                        //}
                }

                MPI_Request* rq = new MPI_Request;
                MPI_Isend( partition_map, G.number_of_nodes(), MPI_INT, target, target, m_communicator, rq);
                
                m_cur_num_pushes++;

                m_request_pointers.push_back( rq );
                m_partition_map_buffers.push_back( partition_map );

                m_allready_send_to[target] = true;
        }

        for( unsigned i = 0; i < m_request_pointers.size(); i++) {
                int finished = 0;
                MPI_Status st;
                MPI_Test( m_request_pointers[i], &finished, &st);

                if(finished) {
                        std::swap(m_request_pointers[i], m_request_pointers[m_request_pointers.size()-1]);
                        std::swap(m_partition_map_buffers[i], m_partition_map_buffers[m_request_pointers.size()-1]);

                        delete[] m_partition_map_buffers[m_partition_map_buffers.size() - 1];
                        delete   m_request_pointers[m_request_pointers.size() - 1];

                        m_partition_map_buffers.pop_back();
                        m_request_pointers.pop_back();
                }
        }
}

void exchanger::recv_incoming( PartitionConfig & config, graph_access & G, population & island ) {
        int rank;
        MPI_Comm_rank( m_communicator, &rank);
        
        int flag; MPI_Status st;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_communicator, &flag, &st);
        
        while(flag) {
                Individuum out;
                int* partition_map = new int[G.number_of_nodes()];
                out.partition_map  = partition_map;
                out.cut_edges      = new std::vector<EdgeID>();

                MPI_Status rst;
                MPI_Recv( out.partition_map, G.number_of_nodes(), MPI_INT, st.MPI_SOURCE, rank, m_communicator, &rst); 
                
                //recompute cut edges and edge cut locally
                forall_nodes(G, node) {
                        forall_out_edges(G, e, node) {
                                NodeID target = G.getEdgeTarget(e);
                                if(partition_map[node] != partition_map[target]) {
                                        out.cut_edges->push_back(e);
                                }
                        } endfor
                } endfor

                out.objective = m_qm.objective(config, G, partition_map);
                island.insert( G, out );

                if( (unsigned)out.objective < (unsigned)m_prev_best_objective) {
                        m_prev_best_objective = out.objective;
                        std::cout << "rank " <<  rank 
                                  <<   ": pool improved (inc) **************************************** " 
                                  <<  out.objective << std::endl;

                        for( unsigned i = 0; i < m_allready_send_to.size(); i++) {
                                m_allready_send_to[i] = false;
                        }

                        m_allready_send_to[rank] = true;
                        m_cur_num_pushes         = 0;
                }

                m_allready_send_to[st.MPI_SOURCE] = true; // we dont need to send it back - saves us P * 1 messages of length n

                MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_communicator, &flag, &st);
        }
}

