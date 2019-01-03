/******************************************************************************
 * population.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef POPULATION_AEFH46G6
#define POPULATION_AEFH46G6

#include <sstream>
#include <mpi.h>

#include "data_structure/graph_access.h"
#include "partition_config.h"
#include "timer.h"

struct Individuum {
        int* partition_map;
        EdgeWeight objective;
        std::vector<EdgeID>* cut_edges; //sorted
};

struct ENC {
        std::vector<NodeID> vertices;
};

class population {
        public:
                population( MPI_Comm comm, const PartitionConfig & config );
                virtual ~population();

                void createIndividuum(const PartitionConfig & config, 
                                      graph_access & G, 
				      Individuum & ind, 
				      bool output); 

                void combine(const PartitionConfig & config, 
                             graph_access & G, 
                             Individuum & first_ind, 
                             Individuum & second_ind, 
                             Individuum & output_ind); 

                void combine_cross(const PartitionConfig & partition_config, 
				   graph_access & G, 
				   Individuum & first_ind, 
				   Individuum & output_ind);

                void mutate_random(const PartitionConfig & partition_config, 
                                   graph_access & G, 
                                   Individuum & first_ind);

                void insert(graph_access & G, Individuum & ind);

                void set_pool_size(int size);

                void extinction();

                void get_two_random_individuals(Individuum & first, Individuum & second);
               
                void get_one_individual_tournament(Individuum & first); 

                void get_two_individuals_tournament(Individuum & first, Individuum & second);

                void replace(Individuum & in, Individuum & out);

                void get_random_individuum(Individuum & ind);

                void get_best_individuum(Individuum & ind);

                bool is_full(); 

                void apply_fittest( graph_access & G, EdgeWeight & objective);

                unsigned size() { return m_internal_population.size(); }
                
                void print();

                void write_log(std::string & filename);


        private:

                unsigned                m_no_partition_calls;
                unsigned 		m_population_size;
                std::vector<Individuum> m_internal_population;
                std::vector< std::vector< unsigned int > > m_vertex_ENCs;
                std::vector< ENC > m_ENCs;

                int m_num_NCs;
                int m_num_NCs_computed;
                int m_num_ENCs;
                int m_time_stamp;

                MPI_Comm m_communicator;

                std::stringstream m_filebuffer_string;
                timer   	  m_global_timer;
};


#endif /* end of include guard: POPULATION_AEFH46G6 */
