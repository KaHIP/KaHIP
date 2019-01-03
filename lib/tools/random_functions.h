/******************************************************************************
 * random_functions.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef RANDOM_FUNCTIONS_RMEPKWYT
#define RANDOM_FUNCTIONS_RMEPKWYT

#include <iostream>
#include <random>
#include <vector>

#include "definitions.h"
#include "partition_config.h"

typedef std::mt19937 MersenneTwister;

class random_functions {
        public:
                random_functions();
                virtual ~random_functions();

                template<typename sometype>
                        static void circular_permutation(std::vector<sometype> & vec) {
                                if(vec.size() < 2) return;
                                for( unsigned int i = 0; i < vec.size(); i++) {
                                        vec[i] = i;
                                }                   

                                unsigned int size = vec.size();
                                std::uniform_int_distribution<unsigned int> A(0,size-1);
                                std::uniform_int_distribution<unsigned int> B(0,size-1);

                                for( unsigned int i = 0; i < size; i++) {
                                        unsigned int posA = A(m_mt);
                                        unsigned int posB = B(m_mt);

                                        while(posB == posA) {
                                                posB = B(m_mt);
                                        }

                                        if( posA != vec[posB] && posB != vec[posA]) {
                                                std::swap(vec[posA], vec[posB]); 
                                        }
                                } 

                         }

                template<typename sometype>
                        static void permutate_vector_fast(std::vector<sometype> & vec, bool init) {
                                if(init) {
                                        for( unsigned int i = 0; i < vec.size(); i++) {
                                                vec[i] = i;
                                        }                   
                                }
                                
                                if(vec.size() < 10) return;
                                        
                                int distance = 20; 
                                std::uniform_int_distribution<unsigned int> A(0, distance);
                                unsigned int size = vec.size()-4;
                                for( unsigned int i = 0; i < size; i++) {
                                        unsigned int posA = i;
                                        unsigned int posB = (posA + A(m_mt))%size;
                                        std::swap(vec[posA], vec[posB]);
                                        std::swap(vec[posA+1], vec[posB+1]); 
                                        std::swap(vec[posA+2], vec[posB+2]); 
                                        std::swap(vec[posA+3], vec[posB+3]); 
                                }               
                        }

                static void permutate_vector_good(std::vector<std::pair< NodeID, NodeID >> & vec) {

                        unsigned int size = vec.size();
                        if(size < 4) return;

                        std::uniform_int_distribution<unsigned int> A(0,size - 4);
                        std::uniform_int_distribution<unsigned int> B(0,size - 4);

                        for( unsigned int i = 0; i < size; i++) {
                                unsigned int posA = A(m_mt);
                                unsigned int posB = B(m_mt);
                                std::swap(vec[posA], vec[posB]); 
                                std::swap(vec[posA+1], vec[posB+1]); 
                                std::swap(vec[posA+2], vec[posB+2]); 
                                std::swap(vec[posA+3], vec[posB+3]); 

                        } 
                }

                template<typename sometype>
                        static void permutate_vector_good(std::vector<sometype> & vec, bool init) {
                                if(init) {
                                        for( unsigned int i = 0; i < vec.size(); i++) {
                                                vec[i] = (sometype)i;
                                        }                   
                                }

                                if(vec.size() < 10) { 
                                        permutate_vector_good_small(vec);
                                        return;
                                }
                                unsigned int size = vec.size();
                                std::uniform_int_distribution<unsigned int> A(0,size - 4);
                                std::uniform_int_distribution<unsigned int> B(0,size - 4);

                                for( unsigned int i = 0; i < size; i++) {
                                        unsigned int posA = A(m_mt);
                                        unsigned int posB = B(m_mt);
                                        std::swap(vec[posA], vec[posB]); 
                                        std::swap(vec[posA+1], vec[posB+1]); 
                                        std::swap(vec[posA+2], vec[posB+2]); 
                                        std::swap(vec[posA+3], vec[posB+3]); 

                                } 
                        }

                template<typename sometype>
                        static void permutate_vector_good_small(std::vector<sometype> & vec) {
                                if(vec.size() < 2) return;
                                unsigned int size = vec.size();
                                std::uniform_int_distribution<unsigned int> A(0,size-1);
                                std::uniform_int_distribution<unsigned int> B(0,size-1);

                                for( unsigned int i = 0; i < size; i++) {
                                        unsigned int posA = A(m_mt);
                                        unsigned int posB = B(m_mt);
                                        std::swap(vec[posA], vec[posB]); 
                                } 
                        }

                template<typename sometype>
                        static void permutate_entries(const PartitionConfig & partition_config, 
                                                      std::vector<sometype> & vec, 
                                                      bool init) {
                                if(init) {
                                        for( unsigned int i = 0; i < vec.size(); i++) {
                                                vec[i] = i;
                                        }                   
                                }

                                switch(partition_config.permutation_quality) {
                                        case PERMUTATION_QUALITY_NONE: break;
                                        case PERMUTATION_QUALITY_FAST: permutate_vector_fast(vec, false); break;
                                        case PERMUTATION_QUALITY_GOOD: permutate_vector_good(vec, false); break;
                                }      

                        }

                static bool nextBool() {
                        std::uniform_int_distribution<unsigned int> A(0,1);
                        return (bool) A(m_mt); 
                }


                //including lb and rb
                static unsigned nextInt(unsigned int lb, unsigned int rb) {
                        std::uniform_int_distribution<unsigned int> A(lb,rb);
                        return A(m_mt); 
                }

                static double nextDouble(double lb, double rb) {
                        double rnbr   = (double) rand() / (double) RAND_MAX; // rnd in 0,1
                        double length = rb - lb;
                        rnbr         *= length;
                        rnbr         += lb;

                        return rnbr; 
                }

                static void setSeed(int seed) {
                        m_seed = seed;
                        srand(seed);
                        m_mt.seed(m_seed);
                }

        private:
                static int m_seed;
                static MersenneTwister m_mt;
};

#endif /* end of include guard: RANDOM_FUNCTIONS_RMEPKWYT */
