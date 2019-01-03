/******************************************************************************
 * augmented_Qgraph.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef AUGMENTED_QUOTIENT_GRAPH_E5ZEJUBV
#define AUGMENTED_QUOTIENT_GRAPH_E5ZEJUBV

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "definitions.h"
#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "uncoarsening/refinement/quotient_graph_refinement/2way_fm_refinement/two_way_fm.h"

struct pairwise_local_search { // a single two way local search
        std::vector<Gain>        gains;
        std::vector<NodeID>      vertex_movements;
        std::vector<PartitionID> block_movements;
        std::vector<int>         load_difference;

        pairwise_local_search() {
        }

};

// bp |-> set pairwise local searches
struct set_pairwise_local_searches { // this is the structure which is associated with a directed qgraph edge
        std::vector<pairwise_local_search> local_searches;
        std::vector<int>                   search_to_use;
        std::vector<Gain>                  search_gain;
        std::vector<int>                   search_num_moves; // number of moves to perform in the associated local search

        set_pairwise_local_searches() {
        }
};

struct simple_move {
        PartitionID from;
        PartitionID to;
        NodeID      node;
};

struct block_pair_difference {
        PartitionID lhs;
        PartitionID rhs;
        int weight_difference;
 
        block_pair_difference() {
        }
};

typedef std::unordered_map<const boundary_pair, set_pairwise_local_searches, hash_boundary_pair_directed, compare_boundary_pair_directed> augmented_Qgraph_internal;

class augmented_Qgraph {
public:
        augmented_Qgraph();
        virtual ~augmented_Qgraph();

        //prepare data structures
        void prepare( PartitionConfig & config, graph_access & G, graph_access & G_bar, unsigned & steps);

        // get information on gains coresponding to vertex differences
        bool exists_vmovements_of_diff( boundary_pair & bp, unsigned  & diff); 
        Gain get_gain_of_vmovements( boundary_pair & bp, unsigned  & diff); 

        // get actual movement coresponding to vertex differences
        void get_associated_vertices(boundary_pair & pair, unsigned & load_diff, std::vector<NodeID> & vertices_of_move);
        void get_associated_blocks  (boundary_pair & pair, unsigned & load_diff, std::vector<PartitionID> & blocks_of_move);

        // commit a pairwise local search
        void commit_pairwise_local_search( boundary_pair & pair, pairwise_local_search & pls);

        //query wether to local searches are conflicted
        bool check_conflict( const  PartitionConfig & config, 
                             PartitionID & lhs, PartitionID & rhs, 
                             unsigned forward_load_diff, unsigned backward_load_diff);

        int get_max_vertex_weight_difference() { return m_max_vertex_weight_difference; };

private:
        augmented_Qgraph_internal m_aqg;
        int m_max_vertex_weight_difference;
};



inline
void augmented_Qgraph::get_associated_vertices(boundary_pair & pair, 
                                               unsigned & load_diff, 
                                               std::vector<NodeID> & vertices_of_move) {
        unsigned internal_idx     = load_diff - 1;
        unsigned search_to_use    = m_aqg[pair].search_to_use[internal_idx]; 
        unsigned search_num_moves = m_aqg[pair].search_num_moves[internal_idx];
        for( unsigned j = 0; j <=  search_num_moves; j++) {
                vertices_of_move.push_back(m_aqg[pair].local_searches[search_to_use].vertex_movements[j]);
        }
}

inline
void augmented_Qgraph::get_associated_blocks(boundary_pair & pair, 
                                             unsigned & load_diff, 
                                             std::vector<PartitionID> & blocks_of_move) {
        unsigned internal_idx     = load_diff - 1;
        unsigned search_to_use    = m_aqg[pair].search_to_use[internal_idx];
        unsigned search_num_moves = m_aqg[pair].search_num_moves[internal_idx];
        for( unsigned j = 0; j <=  search_num_moves; j++) {
                blocks_of_move.push_back(m_aqg[pair].local_searches[search_to_use].block_movements[j]);
        }
}

inline
bool augmented_Qgraph::exists_vmovements_of_diff( boundary_pair & bp, unsigned & diff) {
        unsigned internal_idx = diff - 1;
        if( m_aqg[bp].local_searches.size() > 0 ) {
            if(m_aqg[bp].search_to_use.size() > internal_idx) {
                    if(m_aqg[bp].search_to_use[internal_idx] != -1) {
                            return true;
                    }
            }
        }

        return false;
}

inline 
Gain augmented_Qgraph::get_gain_of_vmovements( boundary_pair & bp, unsigned & diff) {
        unsigned internal_idx = diff - 1;
        return m_aqg[bp].search_gain[internal_idx];
}

inline 
void augmented_Qgraph::prepare( PartitionConfig & config, graph_access & G, graph_access & G_bar, unsigned & steps) {
        m_max_vertex_weight_difference = 0;
        forall_nodes(G_bar, lhs) {
                forall_out_edges(G_bar, e, lhs) {
                        EdgeID rhs = G_bar.getEdgeTarget(e);

                        //find the right edge in the augmented quotient graph
                        boundary_pair bp;
                        bp.k   = config.k;
                        bp.lhs = lhs;
                        bp.rhs = rhs;

                        if( m_aqg[bp].local_searches.size() == 0 ) continue;

                        //estimate the maximum load difference from lhs to rhs on this edge
                        int max_difference = std::numeric_limits<int>::min();
                        for( unsigned i = 0; i < m_aqg[bp].local_searches.size(); i++) {
                                unsigned local_search_size = m_aqg[bp].local_searches[i].vertex_movements.size();
                                m_aqg[bp].local_searches[i].load_difference.resize( local_search_size );

                                int cur_difference  = 0;
                                for( unsigned j = 0; j < local_search_size; j++) {
                                        NodeID node = m_aqg[bp].local_searches[i].vertex_movements[j];
                                        if( G.getPartitionIndex(node) == lhs ) { // hence it will be moved to the other side
                                                cur_difference += G.getNodeWeight(node);
                                        } else {
                                                cur_difference -= G.getNodeWeight(node);
                                        }
                                        m_aqg[bp].local_searches[i].load_difference[j] = cur_difference;

                                        if( cur_difference > max_difference ) {
                                                max_difference = cur_difference;
                                        }
                                }

                        }

                        if(max_difference <= 0) continue;
                        if(max_difference > m_max_vertex_weight_difference) {
                                m_max_vertex_weight_difference = max_difference;
                        }

                        // init arrays
                        int UNDEF = -1;
                        m_aqg[bp].search_to_use.resize(max_difference);
                        m_aqg[bp].search_gain.resize(max_difference);
                        m_aqg[bp].search_num_moves.resize(max_difference);
                        for( int i = 0; i < max_difference; i++) {
                                m_aqg[bp].search_to_use[i]    = UNDEF;
                                m_aqg[bp].search_gain[i]      = std::numeric_limits<Gain>::min();
                                m_aqg[bp].search_num_moves[i] = UNDEF;
                        }

                        ////now create the local search to use array
                        // for each local search
                        for( unsigned i = 0; i < m_aqg[bp].local_searches.size(); i++) {
                                for( unsigned j = 0; j < m_aqg[bp].local_searches[i].vertex_movements.size(); j++) {
                                        int load_diff = m_aqg[bp].local_searches[i].load_difference[j];
                                        if( load_diff <= 0 ) continue; 

                                        unsigned internal_idx  = load_diff - 1;
                                        if( m_aqg[bp].search_gain[internal_idx] < m_aqg[bp].local_searches[i].gains[j]) {
                                                m_aqg[bp].search_num_moves[internal_idx] = j;       
                                                m_aqg[bp].search_gain[internal_idx]      = m_aqg[bp].local_searches[i].gains[j];     
                                                m_aqg[bp].search_to_use[internal_idx]    = i;       
                                        }
                                }
                        }
                } endfor
        } endfor
}

inline
void augmented_Qgraph::commit_pairwise_local_search( boundary_pair & pair, pairwise_local_search & pls) {
        m_aqg[pair].local_searches.push_back(pls);
}

inline
bool augmented_Qgraph::check_conflict( const  PartitionConfig & config, 
                                       PartitionID & lhs, PartitionID & rhs, 
                                       unsigned forward_load_diff, unsigned backward_load_diff) {
        // look at the first vertices to be moved of the 
        // corr. local searches ... 
        // if they are equal then this two cycle contains a conflict
        boundary_pair bp;
        bp.k   = config.k;
        bp.lhs = lhs;
        bp.rhs = rhs;

        unsigned internal_idx          = forward_load_diff - 1;
        unsigned local_searches_to_use = m_aqg[bp].search_to_use[internal_idx];
        NodeID node = m_aqg[bp].local_searches[local_searches_to_use].vertex_movements[0];

        bp.k                  = config.k;
        bp.lhs                = rhs;
        bp.rhs                = lhs;
        internal_idx          = backward_load_diff - 1;
        local_searches_to_use = m_aqg[bp].search_to_use[internal_idx];

        return (node == m_aqg[bp].local_searches[local_searches_to_use].vertex_movements[0]);
}

#endif /* end of include guard: AUGMENTED_QUOTIENT_GRAPH_E5ZEJUBV */
