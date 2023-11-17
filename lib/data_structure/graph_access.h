/******************************************************************************
 * graph_access.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef GRAPH_ACCESS_EFRXO4X2
#define GRAPH_ACCESS_EFRXO4X2

#include <bitset>
#include <cassert>
#include <iostream>
#include <vector>

#include "definitions.h"

struct Node {
    EdgeID firstEdge;
    NodeWeight weight;
};

struct Edge {
    NodeID target;
    EdgeWeight weight;
};

struct refinementNode {
    PartitionID partitionIndex; 
    //Count queueIndex;
};

struct coarseningEdge {
    EdgeRatingType rating;
};

class graph_access;

//construction etc. is encapsulated in basicGraph / access to properties etc. is encapsulated in graph_access
class basicGraph {
    friend class graph_access;

public:
    basicGraph() : m_building_graph(false) {
    }

private:
    //methods only to be used by friend class
    EdgeID number_of_edges() {
        return m_edges.size();
    }

    NodeID number_of_nodes() {
        return m_nodes.size()-1;
    }

    inline EdgeID get_first_edge(const NodeID & node) {
        return m_nodes[node].firstEdge;
    }

    inline EdgeID get_first_invalid_edge(const NodeID & node) {
        return m_nodes[node+1].firstEdge;
    }

    // construction of the graph
    void start_construction(NodeID n, EdgeID m) {
        m_building_graph = true;
        node             = 0;
        e                = 0;
        m_last_source    = -1;

        //resizes property arrays
        m_nodes.resize(n+1);
        m_refinement_node_props.resize(n+1);
        m_edges.resize(m);
        m_coarsening_edge_props.resize(m);

        m_contraction_offset.resize(n+1, 0);

        m_nodes[node].firstEdge = e;
    }

    // Add a new edge from node 'source' to node 'target'.
    // If an edge with source = n has been added, adding
    // edges with source < n will lead to a broken graph.
    EdgeID new_edge(NodeID source, NodeID target) {
        ASSERT_TRUE(m_building_graph);
        ASSERT_TRUE(e < m_edges.size());
       
        m_edges[e].target = target;
        EdgeID e_bar = e;
        ++e;

        ASSERT_TRUE(source+1 < m_nodes.size());
        m_nodes[source+1].firstEdge = e;

        //fill isolated sources at the end
        if ((NodeID)(m_last_source+1) < source) {
            for (NodeID i = source; i>(NodeID)(m_last_source+1); i--) {
                m_nodes[i].firstEdge = m_nodes[m_last_source+1].firstEdge;
            }
        }
        m_last_source = source;
        return e_bar;
    }

    NodeID new_node() {
        ASSERT_TRUE(m_building_graph);
        return node++;
    }

    void finish_construction() {
        // inert dummy node
        m_nodes.resize(node+1);
        m_refinement_node_props.resize(node+1);

        m_contraction_offset.resize(node+1);

        m_edges.resize(e);
        m_coarsening_edge_props.resize(e);

        m_building_graph = false;

        //fill isolated sources at the end
        if ((unsigned int)(m_last_source) != node-1) {
                //in that case at least the last node was an isolated node
                for (NodeID i = node; i>(unsigned int)(m_last_source+1); i--) {
                        m_nodes[i].firstEdge = m_nodes[m_last_source+1].firstEdge;
                }
        }
    }

    // %%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // split properties for coarsening and uncoarsening
    std::vector<Node> m_nodes;
    std::vector<Edge> m_edges;
    
    std::vector<refinementNode> m_refinement_node_props;
    std::vector<coarseningEdge> m_coarsening_edge_props;

    // Offsets for computing sizes of reachable sets for contracted nodes
    std::vector<NodeWeight> m_contraction_offset;
        
    // construction properties
    bool m_building_graph;
    int m_last_source;
    NodeID node; //current node that is constructed
    EdgeID e;    //current edge that is constructed
};

//makros - graph access
#define forall_edges(G,e) { for(EdgeID e = 0, end = G.number_of_edges(); e < end; ++e) {
#define forall_nodes(G,n) { for(NodeID n = 0, end = G.number_of_nodes(); n < end; ++n) {
#define forall_out_edges(G,e,n) { for(EdgeID e = G.get_first_edge(n), end = G.get_first_invalid_edge(n); e < end; ++e) {
#define forall_out_edges_starting_at(G,e,n,e_bar) { for(EdgeID e = e_bar, end = G.get_first_invalid_edge(n); e < end; ++e) {
#define forall_blocks(G,p) { for (PartitionID p = 0, end = G.get_partition_count(); p < end; p++) {
#define endfor }}


class complete_boundary;

class graph_access {
        friend class complete_boundary;
        public:
                graph_access() { m_max_degree_computed = false; m_max_degree = 0; graphref = new basicGraph(); m_separator_block_ID = 2;}
                virtual ~graph_access(){ delete graphref; };

                graph_access(const graph_access&) = delete;

                /* ============================================================= */
                /* build methods */
                /* ============================================================= */
                void start_construction(NodeID nodes, EdgeID edges);
                NodeID new_node();
                EdgeID new_edge(NodeID source, NodeID target);
                void finish_construction();

                /* ============================================================= */
                /* graph access methods */
                /* ============================================================= */
                NodeID number_of_nodes();
                EdgeID number_of_edges();

                EdgeID get_first_edge(NodeID node);
                EdgeID get_first_invalid_edge(NodeID node);

                PartitionID get_partition_count(); 
                void set_partition_count(PartitionID count); 

                PartitionID getSeparatorBlock();
                void setSeparatorBlock(PartitionID id);

                PartitionID getPartitionIndex(NodeID node);
                void setPartitionIndex(NodeID node, PartitionID id);

                PartitionID getSecondPartitionIndex(NodeID node);
                void setSecondPartitionIndex(NodeID node, PartitionID id);

                //to be called if combine in meta heuristic is used
                void resizeSecondPartitionIndex(unsigned no_nodes);

                NodeWeight getNodeWeight(NodeID node);
                void setNodeWeight(NodeID node, NodeWeight weight);

                EdgeWeight getNodeDegree(NodeID node);
                EdgeWeight getWeightedNodeDegree(NodeID node);
                EdgeWeight getMaxDegree();

                EdgeWeight getEdgeWeight(EdgeID edge);
                void setEdgeWeight(EdgeID edge, EdgeWeight weight);

                NodeID getEdgeTarget(EdgeID edge);

                EdgeRatingType getEdgeRating(EdgeID edge);
                void setEdgeRating(EdgeID edge, EdgeRatingType rating);

                // access the contraction offset of a node
                NodeWeight get_contraction_offset(NodeID node) const;
                void set_contraction_offset(NodeID node, NodeWeight offset);

                int* UNSAFE_metis_style_xadj_array();
                int* UNSAFE_metis_style_adjncy_array();

                int* UNSAFE_metis_style_vwgt_array();
                int* UNSAFE_metis_style_adjwgt_array();

                int build_from_metis(int n, int* xadj, int* adjncy);
                int build_from_metis_weighted(int n, int* xadj, int* adjncy, int * vwgt, int* adjwgt);

                //void set_node_queue_index(NodeID node, Count queue_index); 
                //Count get_node_queue_index(NodeID node);

                void copy(graph_access & Gcopy);
        private:
                basicGraph * graphref;     
                bool         m_max_degree_computed;
                unsigned int m_partition_count;
                EdgeWeight   m_max_degree;
                PartitionID  m_separator_block_ID;
                std::vector<PartitionID> m_second_partition_index;
};


inline NodeWeight graph_access::get_contraction_offset(NodeID node) const {
        return graphref->m_contraction_offset[node];
}

inline void graph_access::set_contraction_offset(NodeID node, NodeWeight offset) {
        graphref->m_contraction_offset[node] = offset;
}



/* graph build methods */
inline void graph_access::start_construction(NodeID nodes, EdgeID edges) {
        graphref->start_construction(nodes, edges);
}

inline NodeID graph_access::new_node() {
        return graphref->new_node();
}

inline EdgeID graph_access::new_edge(NodeID source, NodeID target) {
        return graphref->new_edge(source, target);
}

inline void graph_access::finish_construction() {
        graphref->finish_construction();
}

/* graph access methods */
inline NodeID graph_access::number_of_nodes() {
        return graphref->number_of_nodes();
}

inline EdgeID graph_access::number_of_edges() {
        return graphref->number_of_edges();
}

inline void graph_access::resizeSecondPartitionIndex(unsigned no_nodes) {
        m_second_partition_index.resize(no_nodes);
}

inline EdgeID graph_access::get_first_edge(NodeID node) {
#ifdef NDEBUG
        return graphref->m_nodes[node].firstEdge;
#else
        return graphref->m_nodes.at(node).firstEdge;
#endif
}

inline EdgeID graph_access::get_first_invalid_edge(NodeID node) {
        return graphref->m_nodes[node+1].firstEdge;
}

inline PartitionID graph_access::get_partition_count() {
        return m_partition_count;
}

inline PartitionID graph_access::getSecondPartitionIndex(NodeID node) {
#ifdef NDEBUG
        return m_second_partition_index[node];
#else
        return m_second_partition_index.at(node);
#endif
}

inline void graph_access::setSecondPartitionIndex(NodeID node, PartitionID id) {
#ifdef NDEBUG
        m_second_partition_index[node] = id;
#else
        m_second_partition_index.at(node) = id;
#endif
}


inline PartitionID graph_access::getSeparatorBlock() {
        return m_separator_block_ID;
}

inline void graph_access::setSeparatorBlock(PartitionID id) {
        m_separator_block_ID = id;
}

inline PartitionID graph_access::getPartitionIndex(NodeID node) {
#ifdef NDEBUG
        return graphref->m_refinement_node_props[node].partitionIndex;
#else
        return graphref->m_refinement_node_props.at(node).partitionIndex;
#endif
}

inline void graph_access::setPartitionIndex(NodeID node, PartitionID id) {
#ifdef NDEBUG
        graphref->m_refinement_node_props[node].partitionIndex = id;
#else
        graphref->m_refinement_node_props.at(node).partitionIndex = id;
#endif
}

inline NodeWeight graph_access::getNodeWeight(NodeID node){
#ifdef NDEBUG
        return graphref->m_nodes[node].weight;        
#else
        return graphref->m_nodes.at(node).weight;        
#endif
}

inline void graph_access::setNodeWeight(NodeID node, NodeWeight weight){
#ifdef NDEBUG
        graphref->m_nodes[node].weight = weight;        
#else
        graphref->m_nodes.at(node).weight = weight;        
#endif
}

inline EdgeWeight graph_access::getEdgeWeight(EdgeID edge){
#ifdef NDEBUG
        return graphref->m_edges[edge].weight;        
#else
        return graphref->m_edges.at(edge).weight;        
#endif
}

inline void graph_access::setEdgeWeight(EdgeID edge, EdgeWeight weight){
#ifdef NDEBUG
        graphref->m_edges[edge].weight = weight;        
#else
        graphref->m_edges.at(edge).weight = weight;        
#endif
}

inline NodeID graph_access::getEdgeTarget(EdgeID edge){
#ifdef NDEBUG
        return graphref->m_edges[edge].target;        
#else
        return graphref->m_edges.at(edge).target;        
#endif
}

inline EdgeRatingType graph_access::getEdgeRating(EdgeID edge) {
#ifdef NDEBUG
        return graphref->m_coarsening_edge_props[edge].rating;        
#else
        return graphref->m_coarsening_edge_props.at(edge).rating;        
#endif
}

inline void graph_access::setEdgeRating(EdgeID edge, EdgeRatingType rating){
#ifdef NDEBUG
        graphref->m_coarsening_edge_props[edge].rating = rating;
#else
        graphref->m_coarsening_edge_props.at(edge).rating = rating;
#endif
}

inline EdgeWeight graph_access::getNodeDegree(NodeID node) {
        return graphref->m_nodes[node+1].firstEdge-graphref->m_nodes[node].firstEdge;
}

inline EdgeWeight graph_access::getWeightedNodeDegree(NodeID node) {
	EdgeWeight degree = 0;
	for( unsigned e = graphref->m_nodes[node].firstEdge; e < graphref->m_nodes[node+1].firstEdge; ++e) {
		degree += getEdgeWeight(e);
	}
        return degree;
}

inline EdgeWeight graph_access::getMaxDegree() {
        if(!m_max_degree_computed) {
                //compute it
                basicGraph& ref = *graphref;
                forall_nodes(ref, node) {
                        EdgeWeight cur_degree = 0;
                        forall_out_edges(ref, e, node) {
                                cur_degree += getEdgeWeight(e);
                        } endfor
                        if(cur_degree > m_max_degree) {
                                m_max_degree = cur_degree;
                        }
                } endfor
                m_max_degree_computed = true;
        }

        return m_max_degree;
}

inline int* graph_access::UNSAFE_metis_style_xadj_array() {
        int* xadj      = new int[graphref->number_of_nodes()+1];
        basicGraph& ref = *graphref;

        forall_nodes(ref, n) {
                xadj[n] = graphref->m_nodes[n].firstEdge;
        } endfor
        xadj[graphref->number_of_nodes()] = graphref->m_nodes[graphref->number_of_nodes()].firstEdge;
        return xadj;
}


inline int* graph_access::UNSAFE_metis_style_adjncy_array() {
        int* adjncy    = new int[graphref->number_of_edges()];
        basicGraph& ref = *graphref;
        forall_edges(ref, e) {
                adjncy[e] = graphref->m_edges[e].target;
        } endfor 

        return adjncy;
}


inline int* graph_access::UNSAFE_metis_style_vwgt_array() {
        int* vwgt      = new int[graphref->number_of_nodes()];
        basicGraph& ref = *graphref;

        forall_nodes(ref, n) {
                vwgt[n] = (int)graphref->m_nodes[n].weight;
        } endfor
        return vwgt;
}

inline int* graph_access::UNSAFE_metis_style_adjwgt_array() {
        int* adjwgt    = new int[graphref->number_of_edges()];
        basicGraph& ref = *graphref;

        forall_edges(ref, e) {
                adjwgt[e] = (int)graphref->m_edges[e].weight;
        } endfor 

        return adjwgt;
}

inline void graph_access::set_partition_count(PartitionID count) {
        m_partition_count = count;
}

inline int graph_access::build_from_metis(int n, int* xadj, int* adjncy) {
        if(graphref != NULL) {
                delete graphref;
        }
        graphref = new basicGraph();
        start_construction(n, xadj[n]);

        for( unsigned i = 0; i < (unsigned)n; i++) {
                NodeID node = new_node();
                setNodeWeight(node, 1);
                setPartitionIndex(node, 0);

                for( unsigned e = xadj[i]; e < (unsigned)xadj[i+1]; e++) {
                        EdgeID e_bar = new_edge(node, adjncy[e]);
                        setEdgeWeight(e_bar, 1);
                }

        }
        
        finish_construction();
        return 0;
}

inline int graph_access::build_from_metis_weighted(int n, int* xadj, int* adjncy, int * vwgt, int* adjwgt) {
        if(graphref != NULL) {
                delete graphref;
        }
        graphref = new basicGraph();
        start_construction(n, xadj[n]);

        for( unsigned i = 0; i < (unsigned)n; i++) {
                NodeID node = new_node();
                setNodeWeight(node, vwgt[i]);
                setPartitionIndex(node, 0);

                for( unsigned e = xadj[i]; e < (unsigned)xadj[i+1]; e++) {
                        EdgeID e_bar = new_edge(node, adjncy[e]);
                        setEdgeWeight(e_bar, adjwgt[e]);
                }
        }
        
        finish_construction();
        return 0;
}

inline void graph_access::copy(graph_access & G_bar) {
        G_bar.start_construction(number_of_nodes(), number_of_edges());

        basicGraph& ref = *graphref;
        forall_nodes(ref, node) {
                NodeID shadow_node = G_bar.new_node();
                G_bar.setNodeWeight(shadow_node, getNodeWeight(node));
                forall_out_edges(ref, e, node) {
                        NodeID target                   = getEdgeTarget(e);
                        EdgeID shadow_edge              = G_bar.new_edge(shadow_node, target);
                        G_bar.setEdgeWeight(shadow_edge, getEdgeWeight(e));
                } endfor
        } endfor

        G_bar.finish_construction();
}

#endif /* end of include guard: GRAPH_ACCESS_EFRXO4X2 */
