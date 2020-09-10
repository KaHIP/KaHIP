#ifndef PROCESSOR_TREE_H
#define PROCESSOR_TREE_H

using namespace std;

#include <vector>

class processor_tree
{
public:

	/** @brief Default constructor. **/
	processor_tree(){};

	/**	@brief Constructor to create a processor tree based on leaves.
		@param[in] distances[vector]: each vector element corresponds to the communication costs on each level from the leaves
		@param[in] descendants[vector]: each vector element corresponds to the number of descendants on each level
		@brief convention that ascending element positions (0,1,2, ...) corresponds to levels in the tree from higher to lower. 
	*/
	processor_tree(const vector<int> &distances, const vector<int> descendants ) {
		assert( distances.size() == descendants.size());
		traversalDistances = distances;
		traversalDescendants = descendants;
		numOfLevels = distances.size();
		for( unsigned int i = 0; i < get_numOfLevels(); i++)
			numPUs *= traversalDescendants[i];
	}
  
	~processor_tree(){};

	/**	@brief Function to calculate communication costs between two application vertices (based on labels)
		@param[in] label_size[int]: size of the bit-label representation of the application graph
		@param[in] x, y[int]: application vertices ( = their label)
		@brief[out] their communication cost based on the distances in the hierarchy of processor tree   
	*/
  
	inline int getDistance_xy(int label_size, int x, int y) {
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


	/* inline int getDistance_xy(unsigned int x, unsigned int y) { */
	/*   //now depending on x and y, generate distance */
	/*   int k = 0; */
	/*   unsigned long long int xor_x_y = x ^ y; */
	/*   int count_leading_zeros = __builtin_clzll(xor_x_y); */
	/*   int total_n_bits = 8*sizeof(unsigned long long int); */
	/*   int clz = total_n_bits - count_leading_zeros -1; */
	/*   /\* if (clz >= 0) { *\/ */
	/*   /\*   k = (int)floor(clz / config.bit_sec_len);                       *\/ */
	/*   /\*   return config.distances[k]; *\/ */
	/*   /\* } else  { *\/ */
	/*   /\*   return 0; *\/ */
	/*   /\* }        *\/ */
	/* }; */
	

	/**	@brief Function to calculate communication costs between two PUs (based on labels)
		@param[in] x, y[int]: PU, nodes of the processor tree ( = their label)
		@brief[out] their distances in the hierarchy of processor tree   
	*/
  
	inline int getDistance_PxPy(int x, int y) {
		assert((x <= numPUs) and (y <= numPUs) );
		int labelDiff = x ^ y;
		if(!labelDiff)
			return 0;
		
		int count_leading_zeros = __builtin_clzll(labelDiff); // index of highest bit
		int total_n_bits = 8*sizeof(unsigned long long int) - 1;
		int idx = total_n_bits - count_leading_zeros;
		assert(idx <= numOfLevels); // index of no more than number of levels
		if(idx >= traversalDistances.size())
			return 0;
		return traversalDistances[idx];
	}


  
	inline vector<int> get_traversalDistances() {return traversalDistances;};
	inline vector<int> get_traversalDescendants() {return traversalDescendants;};
	inline unsigned int get_numOfLevels() {return numOfLevels;};
	inline unsigned int get_numPUs() {return numPUs;};

	void print() {
		std::cout << " ===== Printing Tree Information ===== " << std::endl;
		for( unsigned int i = 0; i < get_numOfLevels(); i++) {
			std::cout << "Level ==" << i << "== Distance : "
				  <<   traversalDistances[i] << " Descedants : "
				  <<   traversalDescendants[i]  << std::endl;
		}
		std::cout << "Total number of processors = " << get_numPUs() << std::endl;
		std::cout << " ===================================== " << std::endl;
	}


	void print_allPairDistances() {
		std::cout << " ========== Distance Matrix ==========" << std::endl;
		std::cout << " ===================================== " << std::endl;
		for( unsigned int i = 0; i < get_numPUs(); i++) {
			for( unsigned int j = 0; j < get_numPUs(); j++) {
				std::cout << getDistance_PxPy(i, j) << "  "; //<< std::endl;
			}
			std::cout  << std::endl;
		}
		std::cout << " ===================================== " << std::endl;
	}


  
  
private:
	unsigned int numOfLevels;
	unsigned int numPUs = 1;
	// Q: make distances double?
	vector<int> traversalDistances;
	vector<int> traversalDescendants;

};


#endif/* PROCESSOR_ΤREE_H */
