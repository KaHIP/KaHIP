/******************************************************************************
 * complete_boundary.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef COMPLETE_BOUNDARY_URZZFDEI
#define COMPLETE_BOUNDARY_URZZFDEI

#include <execinfo.h>
#include <unordered_map>
#include <utility>

#include "boundary_lookup.h"
#include "data_structure/graph_access.h"
#include "partial_boundary.h"
#include "partition_config.h"

struct block_informations {
        NodeWeight block_weight;
        NodeID block_no_nodes;
};

typedef std::vector<boundary_pair> QuotientGraphEdges;

class complete_boundary {
        public:
                complete_boundary(graph_access * G );
                virtual ~complete_boundary();

                void build();
                void build_from_coarser(complete_boundary * coarser_boundary, NodeID coarser_no_nodes, CoarseMapping * cmapping);

                inline void insert(NodeID node, PartitionID insert_node_into, boundary_pair * pair);
                inline bool contains(NodeID node, PartitionID partition, boundary_pair * pair);
                inline void deleteNode(NodeID node, PartitionID partition, boundary_pair * pair);
                void postMovedBoundaryNodeUpdates(NodeID target, boundary_pair * pair, 
                                                  bool update_edge_cuts, bool update_all_boundaries);
                void balance_singletons(const PartitionConfig & config, graph_access & G); 

                inline NodeID size(PartitionID partition, boundary_pair * pair);

                inline NodeWeight getBlockWeight(PartitionID partition);
                inline NodeWeight getBlockNoNodes(PartitionID partition);
                inline EdgeWeight getEdgeCut(boundary_pair * pair);
                inline EdgeWeight getEdgeCut(PartitionID lhs, PartitionID rhs);

                inline void setBlockWeight(PartitionID partition, NodeWeight weight);
                inline void setBlockNoNodes(PartitionID partition, NodeID no_nodes);
                inline void setEdgeCut(boundary_pair * pair, EdgeWeight edge_cut);

                inline void getQuotientGraphEdges(QuotientGraphEdges & qgraph_edges);
                inline PartialBoundary&  getDirectedBoundary(PartitionID partition, PartitionID lhs, PartitionID rhs);

                inline void setup_start_nodes(graph_access & G, PartitionID partition, 
                                              boundary_pair & bp, boundary_starting_nodes & start_nodes); 

                inline void setup_start_nodes_around_blocks(graph_access & G, PartitionID & lhs, PartitionID & rhs,
                                                            boundary_starting_nodes & start_nodes);

                inline void setup_start_nodes_all(graph_access & G, boundary_starting_nodes & start_nodes);

                inline void get_max_norm();
                inline void getUnderlyingQuotientGraph( graph_access & qgraph );
                inline void getNeighbors(PartitionID & block, std::vector<PartitionID> & neighbors);

        private:
                //updates lazy values that the access functions need
                inline void update_lazy_values(boundary_pair * pair);
                
                //lazy members to avoid hashtable loop ups
                PartialBoundary*   m_pb_lhs_lazy;
                PartialBoundary*   m_pb_rhs_lazy;
                PartitionID        m_lazy_lhs;
                PartitionID        m_lazy_rhs;
                boundary_pair*     m_last_pair;
                size_t             m_last_key; 
                hash_boundary_pair m_hbp;

                graph_access * m_graph_ref;
                //implicit quotient graph structure
                //
                block_pairs m_pairs;
                std::vector<block_informations> m_block_infos;

                //explicit quotient graph structure / may be outdated!
                graph_access Q;
                std::vector< NodeID > m_singletons;

                //////////////////////////////////////////////////////////////
                ///////// Data Structure Invariants 
                //////////////////////////////////////////////////////////////
#ifndef NDEBUG
        public:
                bool assert_bnodes_in_boundaries();
                bool assert_boundaries_are_bnodes();
#endif 
};



inline void complete_boundary::build() {
        graph_access & G = *m_graph_ref;

        for(PartitionID block = 0; block < G.get_partition_count(); block++) {
                m_block_infos[block].block_weight   = 0;
                m_block_infos[block].block_no_nodes = 0;
        }

        forall_nodes(G, n) {
                PartitionID source_partition = G.getPartitionIndex(n);
                m_block_infos[source_partition].block_weight   += G.getNodeWeight(n);
                m_block_infos[source_partition].block_no_nodes += 1;

                if(G.getNodeDegree(n) == 0) {
                        m_singletons.push_back(n);
                }

                forall_out_edges(G, e, n) {
                        NodeID targetID              = G.getEdgeTarget(e);
                        PartitionID target_partition = G.getPartitionIndex(targetID);
                        bool is_cut_edge             = (source_partition != target_partition);

                        if(is_cut_edge) {
                                boundary_pair bp;
                                bp.k   = m_graph_ref->get_partition_count();
                                bp.lhs = source_partition;
                                bp.rhs = target_partition;
                                update_lazy_values(&bp);
                                m_pairs[bp].edge_cut += G.getEdgeWeight(e);    
                                insert(n, source_partition, &bp);
                        }
                } endfor
        } endfor

        block_pairs::iterator iter; 
        for(iter = m_pairs.begin(); iter != m_pairs.end(); iter++ ) { 
                data_boundary_pair& value = iter->second;
                value.edge_cut /= 2;
        }

}

inline void complete_boundary::build_from_coarser(complete_boundary * coarser_boundary, 
                                                  NodeID coarser_no_nodes, 
                                                  CoarseMapping * cmapping) {

        graph_access & G = *m_graph_ref;

        std::vector<bool> coarse_is_border_node(coarser_no_nodes, false);
        QuotientGraphEdges coarser_qgraph_edges;
        coarser_boundary->getQuotientGraphEdges(coarser_qgraph_edges);

        for(unsigned int i = 0; i < coarser_qgraph_edges.size(); i++) {
                PartitionID lhs        = coarser_qgraph_edges[i].lhs;
                PartitionID rhs        = coarser_qgraph_edges[i].rhs;
                PartialBoundary& lhs_b = coarser_boundary->getDirectedBoundary(lhs, lhs, rhs);
                PartialBoundary& rhs_b = coarser_boundary->getDirectedBoundary(rhs, lhs, rhs);
        
                forall_boundary_nodes(lhs_b, n) {
                        coarse_is_border_node[n] = true;
                } endfor 
                
                forall_boundary_nodes(rhs_b, n) {
                        coarse_is_border_node[n] = true;
                } endfor
                
        }

        for(PartitionID block = 0; block < G.get_partition_count(); block++) {
                m_block_infos[block].block_weight   = 0;
                m_block_infos[block].block_no_nodes = 0;
        }

        forall_nodes(G, n) {
                PartitionID source_partition = G.getPartitionIndex(n);
                m_block_infos[source_partition].block_no_nodes += 1;

                if( G.getNodeDegree(n) == 0 ) {
                        m_singletons.push_back(n);
                }

                NodeID coarse_node = (*cmapping)[n];
                if(!coarse_is_border_node[coarse_node]) continue;
                
                forall_out_edges(G, e, n) {
                        NodeID targetID = G.getEdgeTarget(e);
                        PartitionID target_partition = G.getPartitionIndex(targetID);
                        bool is_cut_edge             = (source_partition != target_partition);

                        if(is_cut_edge) {
                                boundary_pair bp;
                                bp.k   = m_graph_ref->get_partition_count();
                                bp.lhs = source_partition;
                                bp.rhs = target_partition;
                                update_lazy_values(&bp);
                                m_pairs[bp].edge_cut += G.getEdgeWeight(e);    
                                insert(n, source_partition, &bp);
                        }
                } endfor
        } endfor

        for(PartitionID p = 0; p < G.get_partition_count(); p++) {
                setBlockWeight(p, coarser_boundary->getBlockWeight(p));
        }

        block_pairs::iterator iter; 
        for(iter = m_pairs.begin(); iter != m_pairs.end(); iter++ ) { 
                data_boundary_pair& value = iter->second;
                value.edge_cut /= 2;
        }
}

inline void complete_boundary::insert(NodeID node, PartitionID insert_node_into, boundary_pair * pair) {
        update_lazy_values(pair);
        ASSERT_TRUE((m_lazy_lhs == pair->lhs && m_lazy_rhs == pair->rhs) 
                 || (m_lazy_lhs == pair->rhs && m_lazy_rhs == pair->lhs));

        if(insert_node_into == m_lazy_lhs) {
                ASSERT_EQ(m_graph_ref->getPartitionIndex(node),m_lazy_lhs);
                m_pb_lhs_lazy->insert(node);
        } else {
                ASSERT_EQ(m_graph_ref->getPartitionIndex(node),m_lazy_rhs);
                m_pb_rhs_lazy->insert(node);
        }    
}

inline bool complete_boundary::contains(NodeID node, PartitionID partition, boundary_pair * pair){
        update_lazy_values(pair);
        if(partition == m_lazy_lhs) {
                ASSERT_EQ(m_graph_ref->getPartitionIndex(node),m_lazy_lhs);
                return m_pb_lhs_lazy->contains(node);
        } else {
                ASSERT_EQ(m_graph_ref->getPartitionIndex(node),m_lazy_rhs);
                return m_pb_rhs_lazy->contains(node);
        }    
}

inline void complete_boundary::deleteNode(NodeID node, PartitionID partition, boundary_pair * pair) {
        update_lazy_values(pair);
        if(partition == m_lazy_lhs) {
                m_pb_lhs_lazy->deleteNode(node);
        } else {
                m_pb_rhs_lazy->deleteNode(node);
        }    
}

inline NodeID complete_boundary::size(PartitionID partition, boundary_pair * pair){
        update_lazy_values(pair);
        if(partition == m_lazy_lhs) {
                return m_pb_lhs_lazy->size();
        } else {
                return m_pb_rhs_lazy->size();
        }    
}

inline NodeWeight complete_boundary::getBlockWeight(PartitionID partition){
        return m_block_infos[partition].block_weight;
}

inline NodeWeight complete_boundary::getBlockNoNodes(PartitionID partition){
        return m_block_infos[partition].block_no_nodes;
}

inline EdgeWeight complete_boundary::getEdgeCut(boundary_pair * pair){
        update_lazy_values(pair);
        return m_pairs[*pair].edge_cut;
}

inline EdgeWeight complete_boundary::getEdgeCut(PartitionID lhs, PartitionID rhs) {
        boundary_pair bp;
        bp.k   = m_graph_ref->get_partition_count();
        bp.lhs = lhs;
        bp.rhs = rhs;

        return getEdgeCut(&bp);
}

inline void complete_boundary::setBlockWeight(PartitionID partition, NodeWeight weight){
        m_block_infos[partition].block_weight = weight;
}

inline void complete_boundary::setBlockNoNodes(PartitionID partition, NodeID no_nodes){
        m_block_infos[partition].block_no_nodes = no_nodes;
}

inline void complete_boundary::setEdgeCut(boundary_pair * pair, EdgeWeight edge_cut){
        update_lazy_values(pair);
        m_pairs[*pair].edge_cut = edge_cut;
}

inline void complete_boundary::getQuotientGraphEdges(QuotientGraphEdges & qgraph_edges) {
        //the quotient graph is stored implicitly in the pairs hashtable
        block_pairs::iterator iter; 
        for(iter = m_pairs.begin(); iter != m_pairs.end(); iter++ ) { 
                boundary_pair key = iter->first;
                qgraph_edges.push_back(key);
        }
}

inline PartialBoundary& complete_boundary::getDirectedBoundary(PartitionID partition, PartitionID lhs, PartitionID rhs) {
        boundary_pair bp;
        bp.k   = m_graph_ref->get_partition_count();
        bp.lhs = lhs;
        bp.rhs = rhs;

        update_lazy_values(&bp);
        if(partition == m_lazy_lhs) {
                return *m_pb_lhs_lazy;
        } else {
                return *m_pb_rhs_lazy;
        }    
}

inline void complete_boundary::update_lazy_values(boundary_pair * pair) {
        ASSERT_NEQ(pair->lhs, pair->rhs);
        
        boundary_pair & bp = *pair;
        size_t key = m_hbp(bp); 
        if(key != m_last_key) {
                data_boundary_pair & dbp = m_pairs[*pair]; 
                if(!dbp.initialized) {
                        m_pairs[*pair].lhs = pair->lhs;
                        m_pairs[*pair].rhs = pair->rhs;
                        dbp.initialized = true;
                }

                m_pb_lhs_lazy = &dbp.pb_lhs;
                m_pb_rhs_lazy = &dbp.pb_rhs;
                m_lazy_lhs    = dbp.lhs;
                m_lazy_rhs    = dbp.rhs;
                m_last_pair   = pair;
                m_last_key    = key;
        }
}
void complete_boundary::setup_start_nodes(graph_access & G, 
                PartitionID partition, 
                boundary_pair & bp, 
                boundary_starting_nodes & start_nodes) {

        start_nodes.resize(size(partition, &bp));
        NodeID cur_idx = 0;

        PartitionID lhs         = bp.lhs;
        PartitionID rhs         = bp.rhs;
        PartialBoundary & lhs_b = getDirectedBoundary(partition, lhs, rhs);

        forall_boundary_nodes(lhs_b, cur_bnd_node) {
                ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), partition);
                start_nodes[cur_idx++] = cur_bnd_node;
        } endfor
}

inline void complete_boundary::get_max_norm() {
         QuotientGraphEdges qgraph_edges;
         getQuotientGraphEdges(qgraph_edges);
         double max = 0;
         for( unsigned i = 0; i < qgraph_edges.size(); i++) {
                 boundary_pair & pair = qgraph_edges[i];

                 if( m_pairs[pair].edge_cut > max ) {
                         max = m_pairs[pair].edge_cut;
                 }
         }
        
         std::cout <<  "max norm is " <<  max  << std::endl;
}

inline void complete_boundary::getUnderlyingQuotientGraph( graph_access & Q_bar ) {
         basicGraph * graphref = new basicGraph; 
         
         if(Q_bar.graphref != NULL) {
                delete Q_bar.graphref;
         }
         Q_bar.graphref = graphref;
 
         std::vector< std::vector< std::pair<PartitionID, EdgeWeight> > >  building_tool;
         building_tool.resize(m_block_infos.size());

         block_pairs::iterator iter; 
         for(iter = m_pairs.begin(); iter != m_pairs.end(); iter++ ) { 
                 boundary_pair cur_pair = iter->first;

                 std::pair<PartitionID, EdgeWeight> qedge_lhs;
                 qedge_lhs.first  = cur_pair.rhs;
                 qedge_lhs.second = m_pairs[cur_pair].edge_cut;
                 building_tool[cur_pair.lhs].push_back(qedge_lhs);

                 std::pair<PartitionID, EdgeWeight> qedge_rhs;
                 qedge_rhs.first  = cur_pair.lhs;
                 qedge_rhs.second = m_pairs[cur_pair].edge_cut;
                 building_tool[cur_pair.rhs].push_back(qedge_rhs);
         }

         Q_bar.start_construction(building_tool.size(), 2*m_pairs.size());
         
         for( unsigned p = 0; p < building_tool.size(); p++) {
                 NodeID node = Q_bar.new_node();
                 Q_bar.setNodeWeight(node,  m_block_infos[p].block_weight);

                 for( unsigned j = 0; j < building_tool[p].size(); j++) {
                         EdgeID e = Q_bar.new_edge(node, building_tool[p][j].first);
                         Q_bar.setEdgeWeight(e, building_tool[p][j].second);
                 }
         }

         Q_bar.finish_construction();
}

inline void complete_boundary::getNeighbors(PartitionID & block, std::vector<PartitionID> & neighbors) {
        //lazy
        if(Q.graphref == NULL) {
                getUnderlyingQuotientGraph(Q); 
                // note that the quotient graph structure currently does not get updated
        }

        forall_out_edges(Q, e, block) {
                ASSERT_TRUE(Q.getEdgeTarget(e) != block);
                neighbors.push_back(Q.getEdgeTarget(e));
        } endfor
}

void complete_boundary::setup_start_nodes_around_blocks(graph_access & G, 
                                                        PartitionID & lhs, PartitionID & rhs, 
                                                        boundary_starting_nodes & start_nodes) {

        std::vector<PartitionID> lhs_neighbors;
        getNeighbors(lhs, lhs_neighbors);

        std::vector<PartitionID> rhs_neighbors;
        getNeighbors(rhs, rhs_neighbors);

        std::unordered_map<NodeID, bool> allready_contained;
        for( unsigned i = 0; i < lhs_neighbors.size(); i++) {
                PartitionID neighbor = lhs_neighbors[i];
                PartialBoundary & partial_boundary_lhs = getDirectedBoundary(lhs, lhs, neighbor);
                forall_boundary_nodes(partial_boundary_lhs, cur_bnd_node) {
                        ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), lhs);
                        if(allready_contained.find(cur_bnd_node) == allready_contained.end() ) { 
                                start_nodes.push_back(cur_bnd_node);
                                allready_contained[cur_bnd_node] = true;
                        }
                } endfor

                PartialBoundary & partial_boundary_neighbor = getDirectedBoundary(neighbor, lhs, neighbor);
                forall_boundary_nodes(partial_boundary_neighbor, cur_bnd_node) {
                        ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), neighbor);
                        if(allready_contained.find(cur_bnd_node) == allready_contained.end()) { 
                                start_nodes.push_back(cur_bnd_node);
                                allready_contained[cur_bnd_node] = true;
                        }
                } endfor
        }

        for( unsigned i = 0; i < rhs_neighbors.size(); i++) {
                PartitionID neighbor = rhs_neighbors[i];
                PartialBoundary & partial_boundary_rhs = getDirectedBoundary(rhs, rhs, neighbor);
                forall_boundary_nodes(partial_boundary_rhs, cur_bnd_node) {
                        ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), rhs);
                        if(allready_contained.find(cur_bnd_node) == allready_contained.end() ) { 
                                start_nodes.push_back(cur_bnd_node);
                                allready_contained[cur_bnd_node] = true;
                        }
                } endfor

                PartialBoundary & partial_boundary_neighbor = getDirectedBoundary(neighbor, rhs, neighbor);
                forall_boundary_nodes(partial_boundary_neighbor, cur_bnd_node) {
                        ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), neighbor);
                        if(allready_contained.find(cur_bnd_node) == allready_contained.end()) { 
                                start_nodes.push_back(cur_bnd_node);
                                allready_contained[cur_bnd_node] = true;
                        }
                } endfor
        }
}


void complete_boundary::setup_start_nodes_all(graph_access & G, boundary_starting_nodes & start_nodes) {
        QuotientGraphEdges quotient_graph_edges;
        getQuotientGraphEdges(quotient_graph_edges);

        std::unordered_map<NodeID, bool> allready_contained;
        
        for( unsigned i = 0; i < quotient_graph_edges.size(); i++) {
                boundary_pair & ret_value = quotient_graph_edges[i];
                PartitionID lhs = ret_value.lhs; 
                PartitionID rhs = ret_value.rhs;

                PartialBoundary & partial_boundary_lhs = getDirectedBoundary(lhs, lhs, rhs);
                forall_boundary_nodes(partial_boundary_lhs, cur_bnd_node) {
                        ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), lhs);
                        if(allready_contained.find(cur_bnd_node) == allready_contained.end() ) { 
                                start_nodes.push_back(cur_bnd_node);
                                allready_contained[cur_bnd_node] = true;
                        }
                } endfor

                PartialBoundary & partial_boundary_rhs = getDirectedBoundary(rhs, lhs, rhs);
                forall_boundary_nodes(partial_boundary_rhs, cur_bnd_node) {
                        ASSERT_EQ(G.getPartitionIndex(cur_bnd_node), rhs);
                        if(allready_contained.find(cur_bnd_node) == allready_contained.end()) { 
                                start_nodes.push_back(cur_bnd_node);
                                allready_contained[cur_bnd_node] = true;
                        }
                } endfor
        }
}


#ifndef NDEBUG
inline bool complete_boundary::assert_bnodes_in_boundaries() {
        PartitionID k = m_graph_ref->get_partition_count();

        for(PartitionID lhs = 0; lhs < k; lhs++) {
                for(PartitionID rhs = 0; rhs < k; rhs++) {
                        if(rhs == lhs || lhs > rhs) continue;

                        boundary_pair bp;
                        bp.k = m_graph_ref->get_partition_count();
                        bp.lhs = lhs;
                        bp.rhs = rhs;
                        graph_access & G = *m_graph_ref;

                        NodeWeight lhs_part_weight = 0;
                        NodeWeight rhs_part_weight = 0;

                        NodeID lhs_no_nodes = 0;
                        NodeID rhs_no_nodes = 0;

                        EdgeWeight edge_cut = 0;
                        forall_nodes(G, n) {
                                PartitionID source_partition = G.getPartitionIndex(n);
                                if(source_partition == lhs){
                                        lhs_part_weight += G.getNodeWeight(n);
                                        lhs_no_nodes++;
                                } else if(source_partition == rhs){
                                        rhs_part_weight += G.getNodeWeight(n); 
                                        rhs_no_nodes++;
                                }

                                forall_out_edges(G, e, n) {
                                        NodeID targetID = G.getEdgeTarget(e);
                                        PartitionID target_partition = G.getPartitionIndex(targetID);
                                        bool is_cut_edge =  (source_partition == lhs &&  target_partition == rhs) 
                                                         || (source_partition == rhs &&  target_partition == lhs);

                                        if(is_cut_edge) {
                                                edge_cut += G.getEdgeWeight(e);
                                                ASSERT_TRUE(contains(n, source_partition,&bp));
                                        }

                                } endfor
                        } endfor

                        ASSERT_EQ(m_block_infos[lhs].block_weight, lhs_part_weight);
                        ASSERT_EQ(m_block_infos[rhs].block_weight, rhs_part_weight);
                        ASSERT_EQ(m_block_infos[lhs].block_no_nodes, lhs_no_nodes);
                        ASSERT_EQ(m_block_infos[rhs].block_no_nodes, rhs_no_nodes);
                        ASSERT_EQ(m_pairs[bp].edge_cut,edge_cut/2);        
                }
        }

        return true;
}

inline bool complete_boundary::assert_boundaries_are_bnodes() {
        graph_access & G = *m_graph_ref;
        forall_nodes(G, n) {
                 PartitionID partition = G.getPartitionIndex(n);
                 forall_out_edges(G, e, n) {
                         NodeID target = G.getEdgeTarget(e);
                         PartitionID targets_partition = G.getPartitionIndex(target);

                         if(partition != targets_partition) {
                                boundary_pair bp;
                                bp.k   = G.get_partition_count();
                                bp.lhs = partition;
                                bp.rhs = targets_partition;

                                ASSERT_TRUE(contains(n, partition, &bp));
                                ASSERT_TRUE(contains(target, targets_partition, &bp));

                         }
                 } endfor
                 
         } endfor
         QuotientGraphEdges qgraph_edges;
         getQuotientGraphEdges(qgraph_edges);
         for( unsigned i = 0; i < qgraph_edges.size(); i++) {
                 boundary_pair & pair = qgraph_edges[i];
                 ASSERT_NEQ(pair.lhs, pair.rhs);
         }
        
        return true;
}
#endif // #ifndef NDEBUG

#endif /* end of include guard: COMPLETE_BOUNDARY_URZZFDEI */
