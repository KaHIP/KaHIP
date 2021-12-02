/******************************************************************************
 * distributed_evolutionary_partitioning.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "communication/mpi_tools.h"
#include "distributed_evolutionary_partitioning.h"
#include "kaHIP_interface.h"
#include "parallel_contraction_projection/parallel_projection.h"
#include "io/parallel_graph_io.h"
#include "tools/distributed_quality_metrics.h"

distributed_evolutionary_partitioning::distributed_evolutionary_partitioning() {
                
}

distributed_evolutionary_partitioning::~distributed_evolutionary_partitioning() {
                
}

void distributed_evolutionary_partitioning::perform_partitioning( MPI_Comm communicator, PPartitionConfig & config, 
                                                                  parallel_graph_access & Q) {

        mpi_tools mpitools;
        parallel_graph_access Q_bar;
        distributed_quality_metrics dqm;
        mpitools.collect_parallel_graph_to_local_graph( communicator, config, Q, Q_bar);
        mpitools.distribute_local_graph( communicator, config, Q_bar);



        int n       = Q_bar.number_of_local_nodes();
        int nparts  = config.k;    // k-way partitioning.
        
        int* xadj   = Q_bar.UNSAFE_metis_style_xadj_array();
        int* adjncy = Q_bar.UNSAFE_metis_style_adjncy_array();

        int* vwgt   = Q_bar.UNSAFE_metis_style_vwgt_array();
        int* adjwgt = Q_bar.UNSAFE_metis_style_adjwgt_array();

        int* partition_map = new int[n];

        PEID rank; MPI_Comm_rank( communicator, &rank);
        
        EdgeWeight prev_cut              = 0;
        NodeWeight prev_max_block_weight = 0;

        if( config.vcycle && rank == 0) {
                forall_local_nodes(Q_bar, node) {
                        partition_map[node] = Q_bar.getSecondPartitionIndex(node); 
                } endfor

                prev_cut = dqm.local_edge_cut(Q_bar, partition_map, communicator);
                prev_max_block_weight = dqm.local_max_block_weight(config, Q_bar, partition_map, communicator);
                //std::cout <<  "prev cut "  <<  prev_cut << std::endl;
                //std::cout <<  "prev max block "  <<  prev_max_block_weight << std::endl;
        }

        if( config.vcycle ) {
                MPI_Bcast(partition_map, n, MPI_INT, ROOT, communicator);
                MPI_Bcast(&prev_cut, 1 , MPI_LONG, ROOT, communicator);
                MPI_Bcast(&prev_max_block_weight, 1 , MPI_LONG, ROOT, communicator);
        }

        double inbalance       = config.inbalance/100.0;
        int edgecut            = 0;
        double balance         = 0;
        bool graph_partitioned = config.vcycle;

        int mode = 0; 
 
        switch( config.initial_partitioning_algorithm ) {
                case KAFFPAESTRONG:
                        mode = STRONG;
                        break;
                case KAFFPAEECO:
                        mode = ECO;
                        break;
                case KAFFPAEFAST:
                        mode = FAST;
                        break;
                case KAFFPAEULTRAFASTSNW:
                        mode = ULTRAFASTSOCIAL;
                        break;
                case KAFFPAEFASTSNW:
                        mode = FASTSOCIAL;
                        break;
                case KAFFPAEECOSNW:
                        mode = ECOSOCIAL;
                        break;
                case KAFFPAESTRONGSNW:
                        mode = STRONGSOCIAL;
                        break;
                default: 
                        mode = FASTSOCIAL;
                        break;
        }

        if(config.vcycle) {
                forall_local_nodes(Q_bar, node) {
                        Q_bar.setNodeLabel(node, partition_map[node]);
                } endfor
        }

        timer t;
        
#ifdef NOOUTPUT
        std::streambuf* backup = std::cout.rdbuf();
        std::ofstream ofs;
        ofs.open("/dev/null");
        std::cout.rdbuf(ofs.rdbuf()); 
#endif

#ifdef DETERMINISTIC_PARHIP
        kaffpaE(&n, 
                vwgt, 
                xadj, 
                adjwgt, 
                adjncy, 
                &nparts, 
                &inbalance, 
                false,  // supress output
                graph_partitioned,
                0, // time limit set to zero, so only the initial population is created
                config.seed, 
                mode,
                communicator, 
                &edgecut, 
                &balance, 
                partition_map);
#else
        kaffpaE(&n, 
                vwgt, 
                xadj, 
                adjwgt, 
                adjncy, 
                &nparts, 
                &inbalance, 
                false,  // supress output
                graph_partitioned,
                config.evolutionary_time_limit, // time limit
                config.seed, 
                mode,
                communicator, 
                &edgecut, 
                &balance, 
                partition_map);
#endif

        

#ifdef NOOUTPUT
        ofs.close();
        std::cout.rdbuf(backup);
#endif

        if( rank == (int)ROOT) {
                PRINT(std::cout <<  "partitioner call took " <<  t.elapsed() << std::endl;);
        }

#ifndef NOOUTPUT
        if( rank == (int)ROOT) {
                std::cout <<  "log>cut computed by IP algorithm " <<  edgecut  << std::endl;
                std::cout <<  "log>balance computed by IP algorithm "<<  balance << std::endl;
        }
#endif

        if( !config.vcycle ) {
                forall_local_nodes(Q_bar, node) {
                        Q_bar.setNodeLabel(node, partition_map[node]);
                } endfor
        } else {
                NodeWeight cur_max_block_weight = dqm.local_max_block_weight(config, Q_bar, partition_map, communicator);

                //balance and cut improved
                bool accept = (cur_max_block_weight <= prev_max_block_weight || balance <= 1.03) && (EdgeWeight)edgecut <= prev_cut;
                // or we previously have not been feasible and now are feasible
                accept = accept || (prev_max_block_weight >= config.upper_bound_partition && cur_max_block_weight <= config.upper_bound_partition);
                         
                if( accept ) {
                        if( rank == (int)ROOT) {
                                PRINT(std::cout <<  "log>update criterion reached, updating partition"  << std::endl;)

                        }
                        forall_local_nodes(Q_bar, node) {
                                Q_bar.setNodeLabel(node, partition_map[node]);
                        } endfor
                } else {
                        if( rank == (int)ROOT) {
                                 PRINT(std::cout <<  "update criterion not reached, not updating partition"  << std::endl;)
                        }
                }
        }

        parallel_projection parallel_project_init;
        parallel_project_init.initial_assignment( Q, Q_bar );

#ifndef NOOUTPUT
        edgecut = dqm.edge_cut(Q, communicator); 
        balance = dqm.balance(config, Q, communicator);
        if( rank == (int)ROOT) {
                std::cout <<  "log>cur edge cut " <<  edgecut  << std::endl;
                std::cout <<  "log>cur balance  " <<  balance << std::endl;
        }
#endif

        delete[] xadj;
        delete[] adjncy;
        delete[] vwgt;
        delete[] adjwgt;
	delete[] partition_map;
}

