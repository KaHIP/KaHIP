/******************************************************************************
 * compare_degrees.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef COMPARE_DEGREES_750FUZ7Z
#define COMPARE_DEGREES_750FUZ7Z

class compare_degrees : public std::binary_function<EdgeWeight, EdgeWeight, bool> {
        public:
                compare_degrees(std::vector<EdgeWeight> * degrees) : m_node_degrees(degrees) {};
                virtual ~compare_degrees() {};

                bool operator() (const EdgeWeight left, const EdgeWeight right ) {
                        return (*m_node_degrees)[left] < (*m_node_degrees)[right];
                }

        private:
                std::vector<EdgeWeight> * m_node_degrees;
};


#endif /* end of include guard: COMPARE_DEGREES_750FUZ7Z */
