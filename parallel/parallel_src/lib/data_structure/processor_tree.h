#ifndef PROCESSOR_TREE_H
#define PROCESSOR_TREE_H

using namespace std;

#include <vector>
//#include <iostream>
#include <bitset>

#include "parallel_graph_access.h"

class processor_tree
{
public:


        processor_tree();

	processor_tree(const vector<int> &distances, const vector<int> &descendants );
  
	~processor_tree(){};

	/**	@brief Function to calculate communication costs between two application vertices (based on labels)
		@param[in] label_size[int]: size of the bit-label representation of the application graph
		@param[in] x, y[int]: application vertices ( = their label)
		@brief[out] their communication cost based on the distances in the hierarchy of processor tree   
	*/
  
	inline int getDistance_xy(int label_size, int x, int y) const  {
		int labelDiff = x ^ y;
		if(!labelDiff)
			return 0;
		int count_leading_zeros = __builtin_clzll(labelDiff); // index of highest bit
		int total_n_bits = 8*sizeof(unsigned long long int) - 1;
		int idx = total_n_bits - count_leading_zeros;
		assert(idx < label_size); // label of no more than label_size bits
		int j = label_size - idx - 1;
		if(j >= traversalDistances.size())
			return 0;
		return traversalDistances[j];
	}
       

	/**	@brief Function to calculate communication costs between two PUs (based on labels)
		@param[in] x, y[int]: PU, nodes of the processor tree ( = their label)
		@brief[out] their distances in the hierarchy of processor tree   
	*/

        inline int getDistance_PxPy(int x, int y) const {
               assert((x <= numPUs) and (y <= numPUs) );
                int k = 0;
                unsigned long long int xor_x_y = (*compact_bin_id)[x] ^ (*compact_bin_id)[y];
                if (!xor_x_y) return 0;
                int count_leading_zeros = __builtin_clzll(xor_x_y);
                int total_n_bits = 8*sizeof(unsigned long long int);
                int clz = total_n_bits - count_leading_zeros -1;
                if (clz >= 0) {
                        k = (int)floor(clz / bit_sec_len);
                        return traversalDistances[k];
                } else  {
                        return 0;
                }       

		
 	}	
  


  
	inline vector<int> get_traversalDistances() const {return traversalDistances;};
	inline vector<int> get_traversalDescendants() const  {return traversalDescendants;};
	inline unsigned int get_numOfLevels() const {return numOfLevels;};
	inline unsigned int get_numPUs() const {return numPUs;};
	void print() const;
	void print_allPairDistances() const;
	int printDistance_PxPy(int x, int y) const;
	void create_predecessorMatrix(parallel_graph_access & P, vector< vector<int>> & predecessorMatrix) const;
	void print_predecessorMatrix(parallel_graph_access & P, vector< vector<int>> & predecessorMatrix) const;
	
	void create_procGraph(parallel_graph_access & cg, MPI_Comm communicator) const;
	void create_parallelprocGraph(parallel_graph_access & cg, MPI_Comm communicator) const;
	void create_commGraph(parallel_graph_access & C, parallel_graph_access & cg, MPI_Comm communicator) const;
	void setCompactRep();

	
private:

	
	unsigned int numOfLevels;
	unsigned int numPUs = 1;
	vector<int> traversalDistances;
	vector<int> traversalDescendants;
	std::vector<unsigned int>  * compact_bin_id;
	int bit_sec_len;
	
	/* parallel_graph_access & P; */
	/* vector< vector<int> > predecessorMatrix; */
	

	

};


#endif/* PROCESSOR_ΤREE_H */


