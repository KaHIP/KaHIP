/******************************************************************************
 * population.cpp 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <sstream>

#include "diversifyer.h"
#include "galinier_combine/gal_combine.h"
#include "graph_partitioner.h"
#include "population.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"
#include "uncoarsening/refinement/cycle_improvements/cycle_refinement.h"

population::population( MPI_Comm communicator, const PartitionConfig & partition_config ) {
        m_population_size    = partition_config.mh_pool_size;
        m_no_partition_calls = 0;
        m_num_NCs            = partition_config.mh_num_ncs_to_compute;
        m_num_NCs_computed   = 0;
        m_num_ENCs           = 0;
        m_time_stamp         = 0;
        m_communicator       = communicator;
        m_global_timer.restart();
}

population::~population() {
        for( unsigned i = 0; i < m_internal_population.size(); i++) {
                delete[] (m_internal_population[i].partition_map);
                delete m_internal_population[i].cut_edges;
        }         
}

void population::set_pool_size(int size) {
        m_population_size = size;
}

void population::createIndividuum(const PartitionConfig & config, 
                                  graph_access & G, 
			          Individuum & ind, bool output) {

        PartitionConfig copy = config;
        graph_partitioner partitioner;
        quality_metrics qm;

        std::ofstream ofs;
        std::streambuf* backup = std::cout.rdbuf();
        ofs.open("/dev/null");
        std::cout.rdbuf(ofs.rdbuf()); 

        timer t; t.restart();

        if(config.buffoon) { // graph is weighted -> no negative cycle detection yet
                partitioner.perform_partitioning(copy, G);
                ofs.close();
                std::cout.rdbuf(backup);
        } else {
                if(config.kabapE) {
                        double real_epsilon        = config.imbalance/100.0;
                        double lb 		   = real_epsilon+0.005;
                        double ub 	           = real_epsilon+config.kabaE_internal_bal;
                        double epsilon             = random_functions::nextDouble(lb,ub);
                        copy.upper_bound_partition = (1+epsilon)*ceil(config.largest_graph_weight/(double)config.k);

                        partitioner.perform_partitioning(copy, G);

                        ofs.close();
                        std::cout.rdbuf(backup);

                        complete_boundary boundary(&G);
                        boundary.build();

                        copy = config;

                        diversifyer df;
                        df.diversify_kaba(copy);

                        cycle_refinement cr;
                        cr.perform_refinement(copy, G, boundary);
                } else {
                        partitioner.perform_partitioning(copy, G);
                        ofs.close();
                        std::cout.rdbuf(backup);
                }
        }

        int* partition_map = new int[G.number_of_nodes()];

        forall_nodes(G, node) {
                partition_map[node] = G.getPartitionIndex(node);
        } endfor

        ind.objective     = qm.objective(config, G, partition_map);
        ind.partition_map = partition_map;
        ind.cut_edges     = new std::vector<EdgeID>();

        forall_nodes(G, node) {
                forall_out_edges(G, e, node) {
                        NodeID target = G.getEdgeTarget(e);
                        if(partition_map[node] != partition_map[target]) {
                                ind.cut_edges->push_back(e);
                        }
                } endfor
        } endfor

        if(output) {
                 m_filebuffer_string <<  m_global_timer.elapsed() <<  " " <<  ind.cut_edges->size()/2 <<  std::endl;
                 m_time_stamp++;
        }
}

void population::insert(graph_access & G, Individuum & ind) {

        m_no_partition_calls++;
        if(m_internal_population.size() < m_population_size) {
                m_internal_population.push_back(ind);
        } else {
                EdgeWeight worst_objective = 0;
                for( unsigned i = 0; i < m_internal_population.size(); i++) {
                        if(m_internal_population[i].objective > worst_objective) {
                                worst_objective = m_internal_population[i].objective;
                        }
                }         
                if(ind.objective > worst_objective ) {
                        delete[] (ind.partition_map);
                        delete ind.cut_edges;
                        return; // do nothing
                }
                //else measure similarity
                unsigned max_similarity = std::numeric_limits<unsigned>::max();
                unsigned max_similarity_idx = 0;
                for( unsigned i = 0; i < m_internal_population.size(); i++) {
                        if(m_internal_population[i].objective >= ind.objective) {
                                //now measure
				int diff_size = m_internal_population[i].cut_edges->size() + ind.cut_edges->size();
                                std::vector<EdgeID> output_diff(diff_size,std::numeric_limits<EdgeID>::max());

                                set_symmetric_difference(m_internal_population[i].cut_edges->begin(),
                                                         m_internal_population[i].cut_edges->end(),
                                                         ind.cut_edges->begin(),
                                                         ind.cut_edges->end(),
                                                         output_diff.begin());

                                unsigned similarity = 0;
                                for( unsigned j = 0; j < output_diff.size(); j++) {
                                        if(output_diff[j] < std::numeric_limits<EdgeID>::max()) {
                                                similarity++;
                                        } else {
                                                break;
                                        }
                                }

                                if( similarity < max_similarity) {
                                        max_similarity     = similarity;
                                        max_similarity_idx = i;
                                }
                        }
                }         

                delete[] (m_internal_population[max_similarity_idx].partition_map);
                delete m_internal_population[max_similarity_idx].cut_edges;

                m_internal_population[max_similarity_idx] = ind;
        }
}

void population::replace(Individuum & in, Individuum & out) {
        //first find it:
        for( unsigned i = 0; i < m_internal_population.size(); i++) {
                if(m_internal_population[i].partition_map == in.partition_map) {
                        //found it
                        delete[] (m_internal_population[i].partition_map);
                        delete m_internal_population[i].cut_edges;

                        m_internal_population[i] = out;
                        break;
                }
        }
}

void population::combine(const PartitionConfig & partition_config, 
                         graph_access & G, 
                         Individuum & first_ind, 
                         Individuum & second_ind, 
                         Individuum & output_ind) {

        PartitionConfig config = partition_config;
        G.resizeSecondPartitionIndex(G.number_of_nodes());
        if( first_ind.objective < second_ind.objective ) {
                forall_nodes(G, node) {
                        G.setPartitionIndex(node, first_ind.partition_map[node]);
                        G.setSecondPartitionIndex(node, second_ind.partition_map[node]);

                } endfor
        } else {
                forall_nodes(G, node) {
                        G.setPartitionIndex(node, second_ind.partition_map[node]);
                        G.setSecondPartitionIndex(node, first_ind.partition_map[node]);
                } endfor
        }

        config.combine                     = true;
        config.graph_allready_partitioned  = true;
        config.no_new_initial_partitioning = true;

	bool coin = false;
	if( partition_config.mh_enable_gal_combine ) {
		coin = random_functions::nextBool();	
	}

	if( coin ) {
		gal_combine combine_operator;
		combine_operator.perform_gal_combine( config, G);
		int* partition_map = new int[G.number_of_nodes()];

		forall_nodes(G, node) {
			partition_map[node] = G.getPartitionIndex(node);
		} endfor

		quality_metrics qm;
		output_ind.objective     = qm.objective(config, G, partition_map);
		output_ind.partition_map = partition_map;
		output_ind.cut_edges     = new std::vector<EdgeID>();

		forall_nodes(G, node) {
			forall_out_edges(G, e, node) {
				NodeID target = G.getEdgeTarget(e);
				if(partition_map[node] != partition_map[target]) {
					output_ind.cut_edges->push_back(e);
				}
			} endfor
		} endfor
	} else {
	        createIndividuum(config, G, output_ind, true);
	}
        std::cout <<  "objective mh " <<  output_ind.objective << std::endl;
}

void population::combine_cross(const PartitionConfig & partition_config, 
                graph_access & G, 
                Individuum & first_ind, 
                Individuum & output_ind) {

        PartitionConfig config = partition_config;
        G.resizeSecondPartitionIndex(G.number_of_nodes());

        int lowerbound = config.k / 4;
        lowerbound     = std::max(2, lowerbound);
        int kfactor    = random_functions::nextInt(lowerbound,4*config.k);
        kfactor = std::min( kfactor, (int)G.number_of_nodes());

        if( config.mh_cross_combine_original_k ) {
                MPI_Bcast(&kfactor, 1, MPI_INT, 0, m_communicator);
        }

        unsigned larger_imbalance = random_functions::nextInt(config.epsilon,25);
        double epsilon = larger_imbalance/100.0;

        
        PartitionConfig cross_config                      = config;
        cross_config.k                                    = kfactor;
        cross_config.kaffpa_perfectly_balanced_refinement = false;
        cross_config.upper_bound_partition                = (1+epsilon)*ceil(partition_config.largest_graph_weight/(double)partition_config.k);
        cross_config.refinement_scheduling_algorithm      = REFINEMENT_SCHEDULING_ACTIVE_BLOCKS;
        cross_config.combine                              = false;
        cross_config.graph_allready_partitioned           = false;

	std::ofstream ofs;
	std::streambuf* backup = std::cout.rdbuf();
        ofs.open("/dev/null");
        std::cout.rdbuf(ofs.rdbuf()); 

        graph_partitioner partitioner;
        partitioner.perform_partitioning(cross_config, G);

        ofs.close();
        std::cout.rdbuf(backup);

        forall_nodes(G, node) {
                G.setSecondPartitionIndex(node, G.getPartitionIndex(node));
                G.setPartitionIndex(node, first_ind.partition_map[node]);
        } endfor

        config.combine                     = true;
        config.graph_allready_partitioned  = true;
        config.no_new_initial_partitioning = true;

        createIndividuum(config, G, output_ind, true);
        std::cout << "objective cross " << output_ind.objective
                  << " k "              << kfactor
                  << " imbal "          << larger_imbalance
                  << " impro "          << (first_ind.objective - output_ind.objective) << std::endl;

}

void population::mutate_random( const PartitionConfig & partition_config, graph_access & G, Individuum & first_ind) {
        int number = random_functions::nextInt(0,5);

        PartitionConfig config            = partition_config;
        config.combine                    = false;
        config.graph_allready_partitioned = true;
        get_random_individuum(first_ind);

        if(number < 5) {
                forall_nodes(G, node) {
                        G.setPartitionIndex(node, first_ind.partition_map[node]);
                } endfor

                config.no_new_initial_partitioning = true;
                createIndividuum( config, G, first_ind, true);

        } else {
                forall_nodes(G, node) {
                        G.setPartitionIndex(node, first_ind.partition_map[node]);
                } endfor

                config.graph_allready_partitioned  = false;
                createIndividuum( config, G, first_ind, true);
        }
}

void population::extinction( ) {
        for( unsigned i = 0; i < m_internal_population.size(); i++) {
                delete[] m_internal_population[i].partition_map;
                delete m_internal_population[i].cut_edges; 
        }

        m_internal_population.clear();
        m_internal_population.resize(0);
}

void population::get_two_random_individuals(Individuum & first, Individuum & second) {
        int first_idx = random_functions::nextInt(0, m_internal_population.size()-1);
        first = m_internal_population[first_idx];

        int second_idx = random_functions::nextInt(0, m_internal_population.size()-1);
        while( first_idx == second_idx ) {
                second_idx = random_functions::nextInt(0, m_internal_population.size()-1);
        }

        second = m_internal_population[second_idx];
}

void population::get_one_individual_tournament(Individuum & first) {
        Individuum one, two;
        get_two_random_individuals(one, two);
        first  =  one.objective < two.objective ? one : two;
}

void population::get_two_individuals_tournament(Individuum & first, Individuum & second) {
        Individuum one, two;
        get_two_random_individuals(one, two);
        first  =  one.objective < two.objective? one : two;

        get_two_random_individuals(one, two);
        second =  one.objective < two.objective ? one : two;

        if( first.objective == second.objective) {
                second = one.objective >= two.objective? one : two;
        }
}

void population::get_random_individuum(Individuum & ind) {
        int idx = random_functions::nextInt(0, m_internal_population.size()-1);
        ind     = m_internal_population[idx];
}

void population::get_best_individuum(Individuum & ind) {
        EdgeWeight min_objective = std::numeric_limits<EdgeWeight>::max();
        unsigned idx             = 0;

        for( unsigned i = 0; i < m_internal_population.size(); i++) {
                if((EdgeWeight)m_internal_population[i].objective < min_objective) {
                        min_objective = m_internal_population[i].objective;
                        idx           = i;
                }
        }

        ind = m_internal_population[idx];
}

bool population::is_full() {
        return m_internal_population.size() == m_population_size;
}

void population::apply_fittest( graph_access & G, EdgeWeight & objective ) {
        EdgeWeight min_objective = std::numeric_limits<EdgeWeight>::max();
	double best_balance      = std::numeric_limits<EdgeWeight>::max();
        unsigned idx             = 0;

	quality_metrics qm;
        for( unsigned i = 0; i < m_internal_population.size(); i++) {
		forall_nodes(G, node) {
			G.setPartitionIndex(node, m_internal_population[i].partition_map[node]);
		} endfor
		double cur_balance = qm.balance(G);
                if((EdgeWeight)m_internal_population[i].objective < min_objective 
	          || ((EdgeWeight)m_internal_population[i].objective == min_objective && cur_balance < best_balance)) {
                        min_objective = m_internal_population[i].objective;
                        idx           = i;
			best_balance  = cur_balance;
                }
        }

        forall_nodes(G, node) {
                G.setPartitionIndex(node, m_internal_population[idx].partition_map[node]);
        } endfor

        objective = min_objective;
}

void population::print() {
        int rank;
        MPI_Comm_rank( m_communicator, &rank);
        
        std::cout <<  "rank " <<  rank << " fingerprint ";

        for( unsigned i = 0; i < m_internal_population.size(); i++) {
                std::cout <<  m_internal_population[i].objective << " ";
        }         

        std::cout <<  std::endl;
}

void population::write_log(std::string & filename) {
        std::ofstream f(filename.c_str());
        f << m_filebuffer_string.str();
        f.close();
}

