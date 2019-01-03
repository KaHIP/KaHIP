/******************************************************************************
 * union_find.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef UNION_FIND_H
#define UNION_FIND_H

#include <vector>

// A simple Union-Find datastructure implementation.
// This is sometimes also caled "disjoint sets datastructure.
class union_find 
{
        public:
                union_find(unsigned n) : m_parent(n), m_rank(n), m_n(n) {
                        for( unsigned i = 0; i < m_parent.size(); i++) {
                                m_parent[i] = i;
                                m_rank[i]   = 0;
                        }
                };
                inline void Union(unsigned lhs, unsigned rhs)
                {
                        int set_lhs = Find(lhs);
                        int set_rhs = Find(rhs);
                        if( set_lhs != set_rhs ) {
                                if( m_rank[set_lhs] < m_rank[set_rhs]) {
                                        m_parent[set_lhs] = set_rhs;
                                } else {
                                        m_parent[set_rhs] = set_lhs;
                                        if( m_rank[set_lhs] == m_rank[set_rhs] ) m_rank[set_lhs]++;
                                }
                                --m_n;
                        }
                };

                inline unsigned Find(unsigned element)
                {
                        if( m_parent[element] != element ) {
                                unsigned retValue = Find( m_parent[element] );  
                                m_parent[element] = retValue; // path compression
                                return retValue;
                        }
                        return element;
                };

                // Returns:
                //   The total number of sets.
                inline unsigned n() const
                { return m_n; };

        private:
                std::vector< unsigned > m_parent;
                std::vector< unsigned > m_rank;

                // Number of elements in UF data structure.
                unsigned m_n;
};



#endif // ifndef UNION_FIND_H

