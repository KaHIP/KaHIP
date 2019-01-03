/******************************************************************************
 * path_set.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PATH_SET_80E9CQT1
#define PATH_SET_80E9CQT1

#include "data_structure/graph_access.h"
#include "macros_assertions.h"
#include "partition_config.h"
#include "path.h"

class path_set {
        public:

                path_set( graph_access * G, const PartitionConfig * config );
                virtual ~path_set();

                //returns the path that v lies on iff v is an endpoint
                const path& get_path(const NodeID & v) const;

                //returns the number of paths in the set
                PathID path_count() const;

                // add the edge with given id to the path set if it is applicable
                // returns true iff the edge was applicable
                bool add_if_applicable(const NodeID & source, const EdgeID & e); 

                //**********
                //Navigation
                //**********

                //returns the if of vertex next to v on the path
                NodeID next_vertex( const NodeID & v ) const;       

                //returns the if of vertex previous to v on the path
                NodeID prev_vertex( const NodeID & v ) const;      

                //returns the id of the edge to the next vertex on the path
                EdgeID edge_to_next(const NodeID & v) const;

                //returns the id of the edge to the previous vertex on the path
                EdgeID edge_to_prev(const NodeID & v) const;
        private:
                graph_access * pG;

                const PartitionConfig * config;


                // Number of Paths
                PathID m_no_of_paths;

                // for every vertex v, vertex_to_path[v] is the id of the path
                std::vector<PathID> m_vertex_to_path;

                // for every path id p, paths[p] is the path for this id
                std::vector<path> m_paths;

                // for every vertex v, next[v] is the id of the vertex that is next on its path.
                // for the head v of a path, next[v] == v
                std::vector<NodeID> m_next;


                // for every vertex v, prev[v] is the id of the vertex that is previouson its path.
                // for the tail v of a path, prev[v] == v
                std::vector<NodeID> m_prev;

                // for every vertex v, next_edge[v] is the id of the vertex that is used to
                // connect the vertex v to the next vertex in the path.
                // if next[v] == v the next_edge[v] = UNDEFINED_EDGE
                std::vector<EdgeID> m_next_edge;

                // for every vertex v, prev_edge[v] is the id of the vertex that is used to
                // connect the vertex v to the previous vertex in the path.
                // if prev[v] == v the prev_edge[v] = UNDEFINED_EDGE
                std::vector<EdgeID> m_prev_edge;

                inline bool is_endpoint(const NodeID & v) const {
                        return (m_next[v] == v or m_prev[v] == v);
                } 
};


inline const path& path_set::get_path(const NodeID & v) const {
        PathID path_id = m_vertex_to_path[v];
        ASSERT_TRUE(path_id < m_vertex_to_path.size());
        return m_paths[path_id];
}

inline PathID path_set::path_count() const {
        return m_no_of_paths;
}

inline NodeID path_set::next_vertex( const NodeID & v ) const {
        return m_next[v];
}       

inline NodeID path_set::prev_vertex( const NodeID & v ) const {
        return m_prev[v];
}      

inline EdgeID path_set::edge_to_next(const NodeID & v) const {
        return m_next_edge[v];
}

inline EdgeID path_set::edge_to_prev(const NodeID & v) const {
        return m_prev_edge[v];
}

inline bool path_set::add_if_applicable(const NodeID & source, const EdgeID & e) {
        graph_access & G = *pG;

        NodeID target = G.getEdgeTarget(e);

        if(config->graph_allready_partitioned && !config->gpa_grow_paths_between_blocks) {
                // in this case we only grow paths inside blocks
                if(G.getPartitionIndex(source) != G.getPartitionIndex(target))
                        return false;
        
                if(config->combine) {
                        if(G.getSecondPartitionIndex(source) != G.getSecondPartitionIndex(target)) {
                                return false;
                        }
                }
                         
        }

        PathID sourcePathID = m_vertex_to_path[source]; 
        PathID targetPathID = m_vertex_to_path[target]; 

        ASSERT_NEQ(source, target);

        path & source_path = m_paths[sourcePathID];
        path & target_path = m_paths[targetPathID];

        if(not is_endpoint(source) or not is_endpoint(target)) {
                // both vertices must be endpoints. otherwise, the edge is not applicable
                return false;
        }

        ASSERT_TRUE(source_path.is_active());
        ASSERT_TRUE(target_path.is_active());
        
         if(source_path.is_cycle() or target_path.is_cycle()) {
                // if one of the paths is a cycle then it is not applicable 
                return false;
        }

        if(sourcePathID != targetPathID) {
                // then we wont close a cycle, and we will join the paths
                // else case handles cycles
                m_no_of_paths--;
                source_path.set_length(source_path.get_length() + target_path.get_length() + 1);

                // first we update the path data structure
                // handle four cases / see implementation details for a picture
                if(source_path.get_head() == source && target_path.get_head() == target) {
                        m_vertex_to_path[target_path.get_tail()] = sourcePathID;
                        source_path.set_head(target_path.get_tail()); 
                } else if(source_path.get_head() == source && target_path.get_tail() == target) {
                        m_vertex_to_path[target_path.get_head()] = sourcePathID;                       
                        source_path.set_head(target_path.get_head()); 
                } else if(source_path.get_tail() == source && target_path.get_head() == target) {
                        m_vertex_to_path[target_path.get_tail()] = sourcePathID;
                        source_path.set_tail(target_path.get_tail()); 
                } else if(source_path.get_tail() == source && target_path.get_tail() == target) {
                        m_vertex_to_path[target_path.get_head()] = sourcePathID;                       
                        source_path.set_tail(target_path.get_head()); 
                }

                //update the double linked list so that we can navigate through the paths later
                //first handle the source node
                if(m_next[source] == source) {
                        ASSERT_TRUE(edge_to_next(source) == UNDEFINED_EDGE);
                        m_next[source]      = target;
                        m_next_edge[source] = e;
                } else {
                        ASSERT_TRUE(edge_to_prev(source) == UNDEFINED_EDGE);
                        m_prev[source]      = target;
                        m_prev_edge[source] = e;
                }

                //then handle the target node
                if(m_next[target] == target) {
                        ASSERT_TRUE(edge_to_next(target) == UNDEFINED_EDGE);
                        m_next[target]      = source;
                        m_next_edge[target] = e;
                } else {
                        ASSERT_TRUE(edge_to_prev(target) == UNDEFINED_EDGE);
                        m_prev[target]      = source;
                        m_prev_edge[target] = e;
                }

                target_path.set_active(false);

                return true;
        } else if(sourcePathID == targetPathID && source_path.get_length() % 2 == 1) {

                //first we update the path data structure
                source_path.set_length(source_path.get_length()+1);
                
                //close the cycle by updateing the doubly linked list
                if(m_next[source_path.get_head()] == source_path.get_head()) {
                        m_next[source_path.get_head()]      = source_path.get_tail();
                        m_next_edge[source_path.get_head()] = e;
                } else {
                        m_prev[source_path.get_head()]      = source_path.get_tail();
                        m_prev_edge[source_path.get_head()] = e;
                }

                if(m_next[source_path.get_tail()] == source_path.get_tail()) {
                        m_next[source_path.get_tail()]      = source_path.get_head();
                        m_next_edge[source_path.get_tail()] = e;
                } else {
                        m_prev[source_path.get_tail()]      = source_path.get_head();
                        m_prev_edge[source_path.get_tail()] = e;
                }

                source_path.set_tail(source_path.get_head());
                return true;
        }       
        return false;
} 


#endif /* end of include guard: PATH_SET_80E9CQT1 */
