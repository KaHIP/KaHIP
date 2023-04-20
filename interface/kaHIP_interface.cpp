/******************************************************************************
 * kaffpa_interface.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#include <iostream>
#include <sstream>

#ifdef USEMETIS
        #include "metis.h"
#endif

#include "kaHIP_interface.h"
#include "../lib/data_structure/graph_access.h"
#include "../lib/io/graph_io.h"
#include "../lib/node_ordering/nested_dissection.h"
#include "../lib/tools/timer.h"
#include "../lib/tools/quality_metrics.h"
#include "../lib/tools/macros_assertions.h"
#include "../lib/tools/random_functions.h"
//#include "../lib/parallel_mh/parallel_mh_async.h"
#include "../lib/partition/uncoarsening/separator/area_bfs.h"
#include "../lib/partition/partition_config.h"
#include "../lib/partition/graph_partitioner.h"
#include "../lib/partition/uncoarsening/separator/vertex_separator_algorithm.h"
#include "../lib/partition/uncoarsening/refinement/cycle_improvements/cycle_refinement.h"
#include "../app/configuration.h"
#include "../app/balance_configuration.h"
#include "../lib/data_structure/matrix/normal_matrix.h"
#include "../lib/data_structure/matrix/online_distance_matrix.h"
#include "../lib/mapping/mapping_algorithms.h"


using namespace std;

void internal_kaffpa_set_configuration( configuration & cfg,
                                 PartitionConfig & partition_config,
                                 int mode) {
        switch( mode ) {
                case FAST: 
                        cfg.fast(partition_config);
                        break;
                case ECO: 
                        cfg.eco(partition_config);
                        break;
                case STRONG: 
                        cfg.strong(partition_config);
                        break;
                case FASTSOCIAL: 
                        cfg.fastsocial(partition_config);
                        break;
                case ECOSOCIAL: 
                        cfg.ecosocial(partition_config);
                        break;
                case STRONGSOCIAL: 
                        cfg.strongsocial(partition_config);
                        break;
                default: 
                        cfg.eco(partition_config);
                        break;
        }
}

void internal_build_graph( PartitionConfig & partition_config, 
                           int* n, 
                           int* vwgt, 
                           int* xadj, 
                           int* adjcwgt, 
                           int* adjncy,
                           graph_access & G) {
        G.build_from_metis(*n, xadj, adjncy); 
        G.set_partition_count(partition_config.k); 
 
        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);
       
        if(vwgt != NULL) {
                forall_nodes(G, node) {
                        G.setNodeWeight(node, vwgt[node]);
                } endfor
        }

        if(adjcwgt != NULL) {
                forall_edges(G, e) {
                        G.setEdgeWeight(e, adjcwgt[e]);
                } endfor 
        }

        balance_configuration bc;
        bc.configurate_balance( partition_config, G);
}

void internal_kaffpa_call(PartitionConfig & partition_config, 
                          bool suppress_output, 
                          int* n, 
                          int* vwgt, 
                          int* xadj, 
                          int* adjcwgt, 
                          int* adjncy, 
                          int* nparts, 
                          double* imbalance, 
                          bool perfectly_balance,
                          int* edgecut, 
                          int* part) {

        //streambuf* backup = cout.rdbuf();
        //ofstream ofs;
        //ofs.open("/dev/null");
        //if(suppress_output) {
               //cout.rdbuf(ofs.rdbuf()); 
        //}

        partition_config.imbalance = 100*(*imbalance);
        partition_config.kaffpa_perfectly_balance = perfectly_balance;
        graph_access G;     
        internal_build_graph( partition_config, n, vwgt, xadj, adjcwgt, adjncy, G);

        graph_partitioner partitioner;
        partitioner.perform_partitioning(partition_config, G);

        if( partition_config.kaffpa_perfectly_balance ) {
                double epsilon                         = partition_config.imbalance/100.0;
                partition_config.upper_bound_partition = (1+epsilon)*ceil(partition_config.largest_graph_weight/(double)partition_config.k);

                complete_boundary boundary(&G);
                boundary.build();

                cycle_refinement cr;
                cr.perform_refinement(partition_config, G, boundary);
        }


        forall_nodes(G, node) {
                part[node] = G.getPartitionIndex(node);
        } endfor

        quality_metrics qm;
        *edgecut = qm.edge_cut(G);

        //ofs.close();
        //cout.rdbuf(backup);
}

void kaffpa(int* n, 
                   int* vwgt, 
                   int* xadj, 
                   int* adjcwgt, 
                   int* adjncy, 
                   int* nparts, 
                   double* imbalance, 
                   bool suppress_output, 
                   int seed,
                   int mode,
                   int* edgecut, 
                   int* part) {
        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = *nparts;

        internal_kaffpa_set_configuration(cfg, partition_config, mode);

        partition_config.seed = seed;
        internal_kaffpa_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, false, edgecut, part);
}

void kaffpa_balance(int* n, 
                   int* vwgt, 
                   int* xadj, 
                   int* adjcwgt, 
                   int* adjncy, 
                   int* nparts, 
                   double* imbalance, 
                   bool perfectly_balance, 
                   bool suppress_output, 
                   int seed, 
                   int mode, 
                   int* edgecut, 
                   int* part) {
        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = *nparts;

        internal_kaffpa_set_configuration(cfg, partition_config, mode);

        partition_config.seed = seed;
        internal_kaffpa_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, perfectly_balance, edgecut, part);
}

void kaffpa_balance_NE(int* n, 
                   int* vwgt, 
                   int* xadj, 
                   int* adjcwgt, 
                   int* adjncy, 
                   int* nparts, 
                   double* imbalance, 
                   bool suppress_output, 
                   int seed,
                   int mode,
                   int* edgecut, 
                   int* part) {
        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = *nparts;

        internal_kaffpa_set_configuration(cfg, partition_config, mode);

        partition_config.seed = seed;
        partition_config.balance_edges = true;
        internal_kaffpa_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, false, edgecut, part);
}

void internal_nodeseparator_call(PartitionConfig & partition_config, 
                          bool suppress_output, 
                          int* n, 
                          int* vwgt, 
                          int* xadj, 
                          int* adjcwgt, 
                          int* adjncy, 
                          int* nparts, 
                          double* imbalance, 
                          int mode,
                          int* num_nodeseparator_vertices, 
                          int** separator) {

        //first perform std partitioning using KaFFPa
        streambuf* backup = cout.rdbuf();
        ofstream ofs;
        ofs.open("/dev/null");
        if(suppress_output) {
               cout.rdbuf(ofs.rdbuf()); 
        }

        partition_config.k         = *nparts;
        partition_config.imbalance = 100*(*imbalance);
        graph_access G;     
        internal_build_graph( partition_config, n, vwgt, xadj, adjcwgt, adjncy, G);
        graph_partitioner partitioner;

        area_bfs::m_deepth.resize(G.number_of_nodes());
        forall_nodes(G, node) {
                area_bfs::m_deepth[node] = 0;
        } endfor

        if( partition_config.k > 2 ) {
                partitioner.perform_partitioning(partition_config, G);

                // now compute a node separator from the partition of the graph
                complete_boundary boundary(&G);
                boundary.build();

                vertex_separator_algorithm vsa;
                std::vector<NodeID> internal_separator;
                vsa.compute_vertex_separator(partition_config, G, boundary, internal_separator);

                // copy to output variables
                *num_nodeseparator_vertices =  internal_separator.size();
                *separator = new int[*num_nodeseparator_vertices];
                for( unsigned int i = 0; i < internal_separator.size(); i++) {
                        (*separator)[i] = internal_separator[i];
                }
        } else {
                
                configuration cfg;
                switch( mode ) {
                        case FAST: 
                                cfg.fast_separator(partition_config);
                                break;
                        case ECO: 
                                cfg.eco_separator(partition_config);
                                break;
                        case STRONG: 
                                cfg.strong_separator(partition_config);
                                break;
                        case FASTSOCIAL: 
                                cfg.fast_separator(partition_config);
                                //cfg.fastsocial_separator(partition_config);
                                break;
                        case ECOSOCIAL: 
                                cfg.eco_separator(partition_config);
                                //cfg.ecosocial_separator(partition_config);
                                break;
                        case STRONGSOCIAL: 
                                //cfg.strongsocial_separator(partition_config);
                                cfg.strong_separator(partition_config);
                                break;
                        default: 
                                cfg.strong_separator(partition_config);
                                break;
                }       
                partition_config.mode_node_separators = true;
                partitioner.perform_partitioning(partition_config, G);
                NodeWeight ns_size = 0;
                forall_nodes(G, node) {
                        if(G.getPartitionIndex(node) == G.getSeparatorBlock()) {
                                ns_size++;
                        }
                } endfor
                *num_nodeseparator_vertices = ns_size;
                *separator = new int[*num_nodeseparator_vertices];
                unsigned int i = 0;
                forall_nodes(G, node) {
                        if(G.getPartitionIndex(node) == G.getSeparatorBlock()) {
                                (*separator)[i] = node;
                                i++;
                        }
                } endfor
        }

        ofs.close();
        cout.rdbuf(backup);
}


void node_separator(int* n, 
                    int* vwgt, 
                    int* xadj, 
                    int* adjcwgt, 
                    int* adjncy, 
                    int* nparts, 
                    double* imbalance, 
                    bool suppress_output, 
                    int seed,
                    int mode,
                    int* num_separator_vertices, 
                    int** separator) {
        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = *nparts;

        switch( mode ) {
                case FAST: 
                        cfg.fast(partition_config);
                        break;
                case ECO: 
                        cfg.eco(partition_config);
                        break;
                case STRONG: 
                        cfg.strong(partition_config);
                        break;
                case FASTSOCIAL: 
                        cfg.fastsocial(partition_config);
                        break;
                case ECOSOCIAL: 
                        cfg.ecosocial(partition_config);
                        break;
                case STRONGSOCIAL: 
                        cfg.strongsocial(partition_config);
                        break;
                default: 
                        cfg.eco(partition_config);
                        break;
        }
        partition_config.seed = seed;

        internal_nodeseparator_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy, nparts, imbalance, mode, num_separator_vertices, separator);
}

void reduced_nd(int* n,
                int* xadj,
                int* adjncy,
                bool suppress_output,
                int seed,
                int mode,
                int* ordering) {
        std::streambuf* backup = std::cout.rdbuf();
        if(suppress_output) {
                std::cout.rdbuf(nullptr);
        }

        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = 2;
        partition_config.dissection_rec_limit = 120;
        partition_config.max_simplicial_degree = 12;
        partition_config.disable_reductions = false;
        partition_config.convergence_factor = 1;
        partition_config.reduction_order = {simplicial_nodes, degree_2_nodes};


        partition_config.seed = seed;
        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);

        switch( mode ) {
                case FAST: 
                        cfg.fast_separator(partition_config);
                        break;
                case ECO: 
                        cfg.eco_separator(partition_config);
                        break;
                case STRONG: 
                        cfg.strong_separator(partition_config);
                        break;
                case FASTSOCIAL: 
                        cfg.fastsocial_separator(partition_config);
                        break;
                case ECOSOCIAL: 
                        cfg.ecosocial_separator(partition_config);
                        break;
                case STRONGSOCIAL: 
                        cfg.strongsocial_separator(partition_config);
                        break;
                default: 
                        cfg.eco_separator(partition_config);
                        break;
        }

        partition_config.seed = seed;

        graph_access G;     
        internal_build_graph( partition_config, n, nullptr, xadj, nullptr, adjncy, G);
        
        partition_config.imbalance = 20;// 20 percent
        balance_configuration bc;
        bc.configurate_balance(partition_config, G);
        
        nested_dissection dissection(&G);
        dissection.perform_nested_dissection(partition_config);

        for (int i = 0; i < *n; ++i) {
                ordering[i] = dissection.ordering()[i];
        }

        // Restore cout output stream
        std::cout.rdbuf(backup);
}

#ifdef USEMETIS
void reduced_nd_fast(int* n,
                      int* xadj,
                      int* adjncy,
                      bool suppress_output,
                      int seed,
                      int* ordering) {
        std::streambuf* backup = std::cout.rdbuf();
        if(suppress_output) {
                std::cout.rdbuf(nullptr);
        }

        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = 2;
        partition_config.dissection_rec_limit = 120;
        partition_config.max_simplicial_degree = 12;
        partition_config.disable_reductions = false;
        partition_config.convergence_factor = 1;
        partition_config.reduction_order = {simplicial_nodes, degree_2_nodes};
        
        partition_config.seed = seed;
        srand(partition_config.seed);
        random_functions::setSeed(partition_config.seed);
        partition_config.seed = seed;
       
        graph_access input_graph;
        internal_build_graph( partition_config, n, nullptr, xadj, nullptr, adjncy, input_graph);
        
        // 'active_graph' is the graph to use after reductions have been applied.
        // If no reductions have been applied, 'active_graph' points to 'input_graph'.
        // Otherwise, it points to 'reduction_stack.back()->get_reduced_graph()'.
        graph_access *active_graph;
        std::vector<std::unique_ptr<Reduction>> reduction_stack;
        bool used_reductions = apply_reductions(partition_config, input_graph, reduction_stack);
        if (used_reductions) {
                active_graph = &reduction_stack.back()->get_reduced_graph();
        } else {
                active_graph = &input_graph;
        }

        idx_t num_nodes = active_graph->number_of_nodes();
        // convert the graph into metis-style
        idx_t* m_xadj = new idx_t[num_nodes + 1];
        forall_nodes((*active_graph), node) {
                m_xadj[node] = (idx_t)active_graph->get_first_edge(node);
        } endfor
        m_xadj[num_nodes] = (idx_t)active_graph->number_of_edges();
        idx_t* m_adjncy = new idx_t[active_graph->number_of_edges()];
        forall_edges((*active_graph), edge) {
                m_adjncy[edge] = (idx_t)active_graph->getEdgeTarget(edge);
        } endfor

        idx_t* m_perm = new idx_t[active_graph->number_of_nodes()];
        idx_t* m_iperm = new idx_t[active_graph->number_of_nodes()];      // inverse ordering. This is the one we are interested in.
        idx_t* metis_options = new idx_t[METIS_NOPTIONS];

        // Perform nested dissection with Metis
        if (num_nodes > 0) {
                METIS_SetDefaultOptions(metis_options);
                metis_options[METIS_OPTION_SEED] = seed;
                METIS_NodeND(&num_nodes, m_xadj, m_adjncy, nullptr, metis_options, m_perm, m_iperm);
        }

        std::vector<NodeID> reduced_labels(active_graph->number_of_nodes(), 0);
        for (size_t i = 0; i < active_graph->number_of_nodes(); ++i) {
                reduced_labels[i] = m_iperm[i];
        }

        // Map ordering of reduced graph to input graph
        std::vector<NodeID> final_labels;
        if (used_reductions) {
                map_ordering(reduction_stack, reduced_labels, final_labels);
        } else {
                final_labels = reduced_labels;
        }

        for (int i = 0; i < *n; ++i) {
                ordering[i] = final_labels[i];
        }

        // Restore cout output stream
        std::cout.rdbuf(backup);

        // Delete temporary graph
        delete[] m_xadj;
        delete[] m_adjncy;
        delete[] m_perm;
        delete[] m_iperm;
        delete[] metis_options;

}
#endif

void internal_processmapping_call(PartitionConfig & partition_config, 
                          bool suppress_output, 
                          int* n, 
                          int* vwgt, 
                          int* xadj, 
                          int* adjcwgt, 
                          int* adjncy, 
                          int mode_mapping,
                          double* imbalance, 
                          int* edgecut, 
                          int* qap,
                          int* part) {

        //streambuf* backup = cout.rdbuf();
        //ofstream ofs;
        //ofs.open("/dev/null");
        //if(suppress_output) {
               //cout.rdbuf(ofs.rdbuf()); 
        //}

        partition_config.imbalance = 100*(*imbalance);
        graph_access G;     
        internal_build_graph( partition_config, n, vwgt, xadj, adjcwgt, adjncy, G);

        graph_partitioner partitioner;
        if( mode_mapping == MAPMODE_BISECTION ) {
                partitioner.perform_partitioning(partition_config, G);
        } else {
                partitioner.perform_partitioning_krec_hierarchy(partition_config, G);
        }

        forall_nodes(G, node) {
                part[node] = G.getPartitionIndex(node);
        } endfor

        quality_metrics qm;
        *edgecut = qm.edge_cut(G);

        int internal_qap = 0;
        //check if k is a power of 2 
        bool power_of_two = (partition_config.k & (partition_config.k-1)) == 0;
        std::vector< NodeID > perm_rank(partition_config.k);
        graph_access C;
        complete_boundary boundary(&G);
        boundary.build();
        boundary.getUnderlyingQuotientGraph(C);

        forall_nodes(C, node) {
                C.setNodeWeight(node, 1);
        } endfor

        if(!power_of_two ) {
                mapping_algorithms ma;
                if( partition_config.distance_construction_algorithm != DIST_CONST_HIERARCHY_ONLINE) {
                        normal_matrix D(partition_config.k, partition_config.k);
                        ma.construct_a_mapping(partition_config, C, D, perm_rank);
                        internal_qap = qm.total_qap(C, D, perm_rank );
                } else {
                        online_distance_matrix D(partition_config.k, partition_config.k);
                        D.setPartitionConfig(partition_config);
                        ma.construct_a_mapping(partition_config, C, D, perm_rank);
                        internal_qap = qm.total_qap(C, D, perm_rank );
                }
        } else {
                for( unsigned i = 0; i < perm_rank.size(); i++) {
                        perm_rank[i] = i;
                }

                online_distance_matrix D(partition_config.k, partition_config.k);
                D.setPartitionConfig(partition_config);
                internal_qap = qm.total_qap(C, D, perm_rank );
        }

        forall_nodes(G, node) {
                G.setPartitionIndex(node, perm_rank[G.getPartitionIndex(node)]);
        } endfor

        *qap = (int)internal_qap;

        //ofs.close();
        //cout.rdbuf(backup);
}

void process_mapping(int* n, int* vwgt, int* xadj, 
                   int* adjcwgt, int* adjncy, 
                   int* hierarchy_parameter,  int* distance_parameter, int hierarchy_depth, 
                   int mode_partitioning, int mode_mapping,
                   double* imbalance,  
                   bool suppress_output, int seed,
                   int* edgecut, int* qap, int* part) {

        configuration cfg;
        PartitionConfig partition_config;
        partition_config.k = 1;

        switch( mode_partitioning ) {
                case FAST: 
                        cfg.fast(partition_config);
                        break;
                case ECO: 
                        cfg.eco(partition_config);
                        break;
                case STRONG: 
                        cfg.strong(partition_config);
                        break;
                case FASTSOCIAL: 
                        cfg.fastsocial(partition_config);
                        break;
                case ECOSOCIAL: 
                        cfg.ecosocial(partition_config);
                        break;
                case STRONGSOCIAL: 
                        cfg.strongsocial(partition_config);
                        break;
                default: 
                        cfg.eco(partition_config);
                        break;
        }

        partition_config.group_sizes.clear();
        partition_config.distances.clear();
        for( int i = 0; i < hierarchy_depth; i++) {
                partition_config.group_sizes.push_back(hierarchy_parameter[i]);
                partition_config.distances.push_back(distance_parameter[i]);
        }

        partition_config.seed = seed;
        internal_processmapping_call(partition_config, suppress_output, n, vwgt, xadj, adjcwgt, adjncy,  mode_mapping, imbalance, edgecut, qap, part);

};


