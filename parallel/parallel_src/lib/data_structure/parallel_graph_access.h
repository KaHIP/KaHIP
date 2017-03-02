/******************************************************************************
 * parallel_graph_access.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef PARALLEL_GRAPH_ACCESS_X6O9MRS8
#define PARALLEL_GRAPH_ACCESS_X6O9MRS8


#include <mpi.h>
#include <unordered_map>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>

#include "data_structure/balance_management.h"
#include "definitions.h"
#include "partition_config.h"
#include "tools/timer.h"

struct Node {
    EdgeID firstEdge;
};
struct NodeData {
    NodeID     label;
    PartitionID block; // a given partition of the graph (for v-cycles)
    NodeWeight weight; // save a little bit of memory
    bool       is_interface_node; // save a little bit of memory
};

//struct NodeData {
    //NodeID     label;
    //PartitionID block:15; // a given partition of the graph (for v-cycles)
    //NodeWeight weight:47; // save a little bit of memory
    //bool       is_interface_node:1; // save a little bit of memory
//};

//struct AdditionalNonLocalNodeData {
    //PEID   peID:15; // save a little bit of memory
    //NodeID globalID:48;
//};

struct AdditionalNonLocalNodeData {
    PEID   peID; // save a little bit of memory
    NodeID globalID;
};

struct Edge {
    NodeID     local_target;
    EdgeWeight weight;
};

//makros - graph access
#define forall_local_nodes(G,n) { for(NodeID n = 0, end = G.number_of_local_nodes(); n < end; ++n) {
#define forall_ghost_nodes(G,n) { for(NodeID node = G.number_of_local_nodes()+1, end = G.number_of_local_nodes()+1+G.number_of_ghost_nodes(); node < end; ++node) { n = node;
#define forall_local_edges(G,e) { for(EdgeID e = 0, end = G.number_of_local_edges(); e < end; ++e) {
#define forall_out_edges(G,e,n) { for(EdgeID e = G.get_first_edge(n), end = G.get_first_invalid_edge(n); e < end; ++e) {
#define endfor }}

class parallel_graph_access;

//handle communication of data associated with ghost nodes
class ghost_node_communication {
public:
        ghost_node_communication(MPI_Comm communicator) : m_iteration_counter(0), m_first_send(true) {
                m_communicator     = communicator;

                MPI_Comm_rank( m_communicator, &m_rank);
                MPI_Comm_size( m_communicator, &m_size);
                
                m_PE_packed.resize(m_size); 
                m_adjacent_processors.resize(m_size); 
                for( PEID peID = 0; peID < (PEID) m_PE_packed.size(); peID++) {
                        m_PE_packed[ peID ]           = false;
                        m_adjacent_processors[ peID ] = false;
                }

                m_send_buffers_A.resize(m_size);
                m_send_buffers_B.resize(m_size);
                m_send_buffers_ptr = & m_send_buffers_A;
                m_send_iteration   = 1;
                m_recv_iteration   = 1;

                m_send_tag         = 100*m_size;
                m_recv_tag         = 100*m_size;

        };

        virtual ~ghost_node_communication() {};

        inline 
        void setGraphReference( parallel_graph_access * G ) {
                m_G = G;
        }; 

        inline 
        void init( ) {
                m_num_adjacent = 0;
                for( PEID peID = 0; peID < (PEID)m_adjacent_processors.size(); peID++) {
                        if( m_adjacent_processors[peID] ) {
                                m_num_adjacent++;
                        }
                }
        }; 


        inline 
        void add_adjacent_processor( PEID peID) {
                m_adjacent_processors[peID] = true;
        }; 

        inline 
        void set_skip_limit( ULONG skip_limit ) {
                m_skip_limit = skip_limit;
        }

        inline 
        void set_desired_rounds( ULONG desired_rounds) {
                m_desired_rounds = desired_rounds;
        }

        inline 
        void update_ghost_node_data( bool check_iteration_counter );

        inline 
        void update_ghost_node_data_finish();

        inline 
        void update_ghost_node_data_global();

        inline
        void addLabel(NodeID node, NodeID label);

        inline
        bool is_adjacent_PE(PEID peID) {
                return m_adjacent_processors[peID];
        }

        inline
        PEID getNumberOfAdjacentPEs() {
                PEID counter = 0;
                for( PEID peID = 0; peID < (PEID)m_adjacent_processors.size(); peID++) {
                        if( m_adjacent_processors[peID] ) counter++;
                }
                return counter;
        }

private:

        inline 
        void receive_messages_of_neighbors();

        parallel_graph_access * m_G;
        PEID m_size;
        PEID m_rank;
        NodeID m_iteration_counter; // this counter is used to manage the communication rounds
        ULONG m_skip_limit; 
        bool m_first_send;

        ULONG m_send_iteration;
        ULONG m_recv_iteration;

        ULONG m_send_tag;
        ULONG m_recv_tag;

        ULONG m_desired_rounds;

        // store the number of adjacent processors ( a block is a neighbor iff there is an edge between the subgraphs )
        PEID m_num_adjacent; 

        std::vector< bool >                   m_PE_packed;
        std::vector< std::vector< NodeID > >  m_send_buffers_A; // buffers to send messages
        std::vector< std::vector< NodeID > >  m_send_buffers_B; // buffers to send messages
        std::vector< std::vector< NodeID > >* m_send_buffers_ptr; // pointer to current buffers to send messages
        std::vector< bool >                   m_adjacent_processors; // buffers to send messages
        std::vector< MPI_Request* >           m_isend_requests;

        MPI_Comm m_communicator;
};


class parallel_graph_access {
public:

        friend class ghost_node_communication;

        parallel_graph_access( ) : m_num_local_nodes(0), 
                                     from(0), 
                                     to(0),
                                     m_num_ghost_nodes(0), m_max_node_degree(0), m_bm(NULL)  { 
                                             m_communicator = MPI_COMM_WORLD;
                                             MPI_Comm_rank( m_communicator, &rank);
                                             MPI_Comm_size( m_communicator, &size);

                                             m_gnc = new ghost_node_communication(m_communicator);
                                             m_gnc->setGraphReference(this);
                                     };

        parallel_graph_access( MPI_Comm communicator );

        virtual ~parallel_graph_access();

        /* ============================================================= */
        /* build methods */
        /* ============================================================= */
        void start_construction(NodeID n, EdgeID m, NodeID global_n, NodeID global_m, bool update_comm_rounds = true) {
                m_building_graph             = true;
                node                         = 0;
                e                            = 0;
                m_last_source                = -1;
                m_num_nodes                  = n+1;
                m_num_local_nodes            = n;
                m_global_n                   = global_n;
                m_global_m                   = global_m;
                m_ghost_adddata_array_offset = n+1;
                m_bm                         = NULL;
		m_cur_degree                 = 0;

                //resizes property arrays
                m_nodes.resize(n+1);
                m_nodes_data.resize(n+1);
                m_edges.resize(m);

                m_nodes[node].firstEdge = e;
                m_divisor = ceil(global_n / (double)size);
                // every PE has to make same amount communication iterations 
                // we use ceil an check afterwards wether everyone has done the right 
                // amount of communication rounds 
                if( update_comm_rounds ) {
                        m_comm_rounds = std::max(m_comm_rounds, 8ULL);
                        m_gnc->set_desired_rounds(m_comm_rounds); 
                        m_gnc->set_skip_limit(ceil(n/(double)m_comm_rounds)); 
                }
        };

        void set_range(NodeID l, NodeID r) {
                from        = l;
                to          = r;
        };

        NodeID get_from_range() {
                return from;
        };

        NodeID get_to_range() {
                return to;
        };

        void set_range_array(std::vector< NodeID > & vertex_dist) {
                m_range_array = vertex_dist;
        };

        PEID get_PEID_from_range_array(NodeID node) {
                // TODO optimize with binary search
                for( PEID peID = 1; peID < (PEID)m_range_array.size(); peID++) {
                        if( node < m_range_array[peID] ) {
                                return (peID-1);
                        }
                }
                return -1;
        };

        NodeID new_node() {
		m_cur_degree = 0;
                ASSERT_TRUE(m_building_graph);
                return node++;
        };

        EdgeID new_edge(NodeID source, NodeID target) {
                ASSERT_TRUE(m_building_graph);
                ASSERT_TRUE(e < m_edges.size());

                // build ghost nodes on the fly
                if( from <= target && target <= to) {
                        m_edges[e].local_target = target - from; 
                } else {
                        m_nodes_data[source].is_interface_node = true;

                        // check wether this is already a ghost node
                        if(m_global_to_local_id.find(target) != m_global_to_local_id.end()) {
                                // this node is already a ghost node
                                m_edges[e].local_target = m_global_to_local_id[target]; 
                        } else {
                                // we need to create a new ghost node
                                m_global_to_local_id[target] = m_num_nodes++;
                                m_edges[e].local_target      = m_global_to_local_id[target]; 

                                //create the ghost node in the array
                                Node dummy;
                                dummy.firstEdge = 0;
                                m_nodes.push_back(dummy);

                                NodeData dummy_data; 
                                dummy_data.label             = target;
                                dummy_data.block             = 0;
                                dummy_data.is_interface_node = false;
                                dummy_data.weight            = 1;
                                m_nodes_data.push_back(dummy_data);

                                // add addtional data
                                AdditionalNonLocalNodeData add_data;
                                //has to be changed once we implement better load balancing 
                                //add_data.peID     = target / m_divisor; 
                                add_data.peID     = get_PEID_from_range_array(target);
                                add_data.globalID = target;

                                m_add_non_local_node_data.push_back(add_data);
                                m_gnc->add_adjacent_processor(add_data.peID);
                        }
                }

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
		m_cur_degree++;

		if( m_cur_degree > m_max_node_degree ) {
			m_max_node_degree = m_cur_degree;
		}
                return e_bar;
        };


        void finish_construction() {
                m_edges.resize(e);
                m_building_graph = false;

                //fill isolated sources at the end
                if ((NodeID)(m_last_source) != node-1) {
                        //in that case at least the last node was an isolated node
                        for (NodeID i = node; i>(NodeID)(m_last_source+1); i--) {
                                m_nodes[i].firstEdge = m_nodes[m_last_source+1].firstEdge;
                        }
                }

                m_gnc->init();
        };

	NodeID get_max_degree() {
		return m_max_node_degree;
	}
        /* ============================================================= */
        /* methods handeling balance */
        /* ============================================================= */
        void init_balance_management( PPartitionConfig & config );
        void update_block_weights();

        // if a ghost node changes its block we update the fuzzy block weight
        void update_non_contained_block_balance( PartitionID from, PartitionID to, NodeWeight node_weight);

        NodeWeight getBlockSize( PartitionID block );
        void setBlockSize( PartitionID block, NodeWeight block_size );

        /* ============================================================= */
        /* parallel graph access methods */
        /* ============================================================= */
        NodeID number_of_local_nodes() {return m_num_local_nodes;};
        NodeID number_of_ghost_nodes() {return m_nodes.size() - m_num_local_nodes - 1;};
        NodeID number_of_global_nodes() {return m_global_n;};
        EdgeID number_of_local_edges() {return m_edges.size();};
        EdgeID number_of_global_edges() {return m_global_m;};
        void set_number_of_global_edges( EdgeID global_edges) {m_global_m = global_edges;};

        void allocate_node_to_cnode() {
                m_nodes_to_cnode.resize( m_nodes.size() );
        }

        void setCNode( NodeID node, NodeID cnode) {
                m_nodes_to_cnode[ node ] = cnode;
        }

        NodeID getCNode( NodeID node ) {
                return m_nodes_to_cnode[node];
        };

        EdgeID get_first_edge(NodeID node);
        EdgeID get_first_invalid_edge(NodeID node);

        NodeID getNodeLabel(NodeID node); 
        void setNodeLabel(NodeID node, NodeID label); 


        NodeID getSecondPartitionIndex(NodeID node); 
        void setSecondPartitionIndex(NodeID node, NodeID label); 

        NodeWeight getNodeWeight(NodeID node); 
        void setNodeWeight(NodeID node, NodeWeight weight); 

        EdgeID getNodeDegree(NodeID node);
        EdgeID getNodeNumGhostNodes(NodeID node);

        bool is_interface_node(NodeID node);
        bool is_local_node(NodeID node);
        bool is_local_node_from_global_id(NodeID node);

        bool is_adjacent_PE(PEID peID) {
                return m_gnc->is_adjacent_PE(peID);
        };

        PEID getNumberOfAdjacentPEs() {
                return m_gnc->getNumberOfAdjacentPEs();
        }

        MPI_Comm getCommunicator() {
                return m_communicator;
        }

        EdgeWeight getEdgeWeight(EdgeID e); 
        void setEdgeWeight(EdgeID e, EdgeWeight weight); 

        NodeID getEdgeTarget(EdgeID e);

        //methods for non-local / ghost nodes only
        //these methods are usally called to communicate data
        PEID getTargetPE(NodeID node);

        //input is a global id 
        //output is the local id
        NodeID getLocalID(NodeID node) {
                if( from <= node && node <= to ) {
                        return node - from;
                } else {
                        return m_global_to_local_id[node];
                }
        };

        //methods for local nodes only
        PEID getGlobalID(NodeID node);

        // these functions should only be called if the graph completely resides on the current PE
        // they are for the call to KaFFPaE
        int* UNSAFE_metis_style_xadj_array();
        int* UNSAFE_metis_style_adjncy_array();

        int* UNSAFE_metis_style_vwgt_array();
        int* UNSAFE_metis_style_adjwgt_array();

        int build_from_metis(int n, int* xadj, int* adjncy);
        int build_from_metis_weighted(int n, int* xadj, int* adjncy, int * vwgt, int* adjwgt);

        /* ============================================================= */
        /* inter process communication routines  */
        /* ============================================================= */
        void update_ghost_node_data( bool check_iteration_counter = true );
        void update_ghost_node_data_finish();
        void update_ghost_node_data_global();

        static void set_comm_rounds(ULONG comm_rounds); 
        static void set_comm_rounds_up(ULONG comm_rounds); 

        /* ============================================================= */
        /* info  */
        /* ============================================================= */
        void printMemoryUsage(std::ostream& out) const {
#ifndef NOOUTPUT
                out << "** approx. local memory usage on hard disk per node [MB (bytes per node)] **" << std::endl;

                unsigned int memoryTotal = 0;
                memoryTotal += printMemoryUsage(out, "nodes", (m_nodes.size()-1) * (sizeof(Node)+sizeof(NodeData)+sizeof(NodeID)+sizeof(AdditionalNonLocalNodeData)));
                memoryTotal += printMemoryUsage(out, "edges", (m_edges.size()-1) * sizeof(Edge));

                printMemoryUsage(out, "TOTAL", memoryTotal);
                out << std::endl;
#endif
        }

        /** Prints the memory usage of one particular data structure of this UpdateableGraph. */
        unsigned int printMemoryUsage(std::ostream& out, const std::string descr, const unsigned int mem) const {
#ifndef NOOUTPUT
                unsigned int megaBytes = (unsigned int)((mem / (double)(1024*1024)) + 0.5);
                unsigned int bytesPerNode = (unsigned int)((mem / (double)(m_nodes.size()-1)) + 0.5);
                out << "   " << descr << ": " << megaBytes << " (" << bytesPerNode << ")" << std::endl;
#endif
                return mem;
        }

        void reinit();
        /* ============================================================= */
        /* parallel graph data structure  */
        /* ============================================================= */
private:
        // the graph representation itself
        // local and ghost nodes in one array, 
        // local nodes are stored in the beginning
        // ghost nodes in the end of the array
        std::vector<Node>                       m_nodes; 
        std::vector<NodeData>                   m_nodes_data;
        std::vector<Edge>                       m_edges;

        //Ghost Node Stuff
        std::vector<AdditionalNonLocalNodeData> m_add_non_local_node_data;

        // NodeID to CNode for ghost nodes and local nodes
        std::vector<NodeID>                     m_nodes_to_cnode; 

        // stores the ranges for which a processor is responsible for
        // m_range_array[i]= starting position of PE i
        std::vector<NodeID>                     m_range_array; 

        std::unordered_map<NodeID, NodeID> m_global_to_local_id;

        NodeID m_ghost_adddata_array_offset; // node id of ghost node - offset to get the position in add data  
        NodeID m_divisor; // needed to compute the target id of a ghost node
        NodeID m_num_local_nodes; // store the number of local / non-ghost nodes
        NodeID from; // each process stores nodes [from. to]
        NodeID to; 
        
        // construction properties
        bool   m_building_graph;
        NodeID m_last_source;
        NodeID m_num_ghost_nodes;
        NodeID node; //current node that is constructed
        EdgeID e;    //current edge that is constructed
        NodeID m_num_nodes; 

        NodeID m_global_n; // global number of nodes
        NodeID m_global_m; // global number of edges
        static ULONG m_comm_rounds; // global number of edges
        static ULONG m_comm_rounds_up; // global number of edges

	NodeID m_max_node_degree;
	NodeID m_cur_degree;

        PEID size;
        PEID rank;

        ghost_node_communication* m_gnc;
        balance_management* m_bm;
        MPI_Comm m_communicator;

};

typedef parallel_graph_access complete_graph_access; // this is just a naming convention for a graph that is completely local

inline
NodeWeight parallel_graph_access::getBlockSize( PartitionID block ) {
        return m_bm->getBlockSize(block);
}

inline 
void parallel_graph_access::setBlockSize( PartitionID block, NodeWeight block_size ) {
        m_bm->setBlockSize(block, block_size);
}

inline EdgeID parallel_graph_access::get_first_edge(NodeID node) {
#ifdef NDEBUG
        return m_nodes[node].firstEdge;
#else
        return m_nodes.at(node).firstEdge;
#endif
}

inline EdgeID parallel_graph_access::get_first_invalid_edge(NodeID node) {
        return m_nodes[node+1].firstEdge;
}

inline EdgeID parallel_graph_access::getNodeLabel(NodeID node) {
#ifdef NDEBUG
        return m_nodes_data[node].label;
#else
        return m_nodes_data.at(node).label;
#endif
}

inline void parallel_graph_access::setNodeLabel(NodeID node, NodeID label) {
        if( m_nodes_data[node].label != label && is_interface_node(node)) {
                m_gnc->addLabel(node, label);
        }
#ifdef NDEBUG
        m_nodes_data[node].label = label;
#else
        m_nodes_data.at(node).label = label;
#endif
}

inline
NodeID parallel_graph_access::getSecondPartitionIndex(NodeID node) {
#ifdef NDEBUG
        return m_nodes_data[node].block;
#else
        return m_nodes_data.at(node).block;
#endif
}

inline
void parallel_graph_access::setSecondPartitionIndex(NodeID node, NodeID block) {
#ifdef NDEBUG
        m_nodes_data[node].block = block;
#else
        m_nodes_data.at(node).block = block;
#endif
}

inline void parallel_graph_access::setNodeWeight(NodeID node, NodeWeight weight) {
#ifdef NDEBUG
        m_nodes_data[node].weight = weight;
#else
        m_nodes_data.at(node).weight = weight;
#endif
}

inline NodeWeight parallel_graph_access::getNodeWeight(NodeID node) {
#ifdef NDEBUG
        return m_nodes_data[node].weight;
#else
        return m_nodes_data.at(node).weight;
#endif
}
inline EdgeID parallel_graph_access::getNodeDegree(NodeID node) {
        return m_nodes[node+1].firstEdge-m_nodes[node].firstEdge;
}

inline EdgeID parallel_graph_access::getNodeNumGhostNodes(NodeID node) {
        return m_nodes[node+1].firstEdge-m_nodes[node].firstEdge;
}
inline bool parallel_graph_access::is_interface_node(NodeID node) {
#ifdef NDEBUG
        return m_nodes_data[node].is_interface_node;
#else
        return m_nodes_data.at(node).is_interface_node;
#endif
}

inline bool parallel_graph_access::is_local_node_from_global_id(NodeID node) {
        return from <= node && node <= to;
}

inline bool parallel_graph_access::is_local_node(NodeID node) {
        return (node < m_num_local_nodes);
}

inline EdgeWeight parallel_graph_access::getEdgeWeight(EdgeID e) {
#ifdef NDEBUG
        return m_edges[e].weight;
#else
        return m_edges.at(e).weight;
#endif
}

inline void parallel_graph_access::setEdgeWeight(EdgeID e, EdgeWeight weight) {
#ifdef NDEBUG
        m_edges[e].weight = weight;
#else
        m_edges.at(e).weight = weight;
#endif
}

inline NodeID parallel_graph_access::getEdgeTarget(EdgeID e){
#ifdef NDEBUG
        return m_edges[e].local_target;        
#else
        return m_edges.at(e).local_target;        
#endif
}

//function should only be called for ghost nodes
inline PEID parallel_graph_access::getTargetPE(NodeID node) {
#ifdef NDEBUG
        return m_add_non_local_node_data[node-m_ghost_adddata_array_offset].peID;
#else
        ASSERT_GEQ(node, m_ghost_adddata_array_offset);
        return m_add_non_local_node_data.at(node-m_ghost_adddata_array_offset).peID;
#endif
}

//function should only be called for local nodes
inline PEID parallel_graph_access::getGlobalID(NodeID node) {
        if( is_local_node(node) ) {
                return from + node;
        } else {
#ifdef NDEBUG
                return m_add_non_local_node_data[node-m_ghost_adddata_array_offset].globalID;
#else
                ASSERT_GEQ(node, m_ghost_adddata_array_offset);
                return m_add_non_local_node_data.at(node-m_ghost_adddata_array_offset).globalID;
#endif
        }
}



// this function should only be called if the graph is completely stored on the root PE
inline int* parallel_graph_access::UNSAFE_metis_style_xadj_array() {
        int * xadj      = new int[number_of_local_nodes()+1];
        forall_local_nodes((*this), node) {
                xadj[node] = m_nodes[node].firstEdge;
        } endfor

        xadj[number_of_local_nodes()] = m_nodes[number_of_local_nodes()].firstEdge;

        return xadj;
}


// this function should only be called if the graph is completely stored on the root PE
inline int* parallel_graph_access::UNSAFE_metis_style_adjncy_array() {
        int * adjncy    = new int[number_of_local_edges()];
        forall_local_edges((*this), e) {
                adjncy[e] = m_edges[e].local_target;
        } endfor 

        return adjncy;
}


inline int* parallel_graph_access::UNSAFE_metis_style_vwgt_array() {
        int * vwgt      = new int[number_of_local_nodes()];
        forall_local_nodes((*this), node) {
                vwgt[node] = m_nodes_data[node].weight;
        } endfor

        return vwgt;
}

inline int* parallel_graph_access::UNSAFE_metis_style_adjwgt_array() {
        int * adjwgt    = new int[number_of_local_edges()];
        forall_local_edges((*this), e) {
                adjwgt[e] = m_edges[e].weight;
        } endfor 

        return adjwgt;
}

inline int parallel_graph_access::build_from_metis(int n, int* xadj, int* adjncy) {
        start_construction(n, xadj[n], n, xadj[n]);
        set_range(0,n);

        for( unsigned i = 0; i < (unsigned)n; i++) {
                NodeID node = new_node();
                setNodeWeight(node, 1);
                setNodeLabel(node, 0);

                for( unsigned e = xadj[i]; e < (unsigned)xadj[i+1]; e++) {
                        EdgeID e_bar = new_edge(node, adjncy[e]);
                        setEdgeWeight(e_bar, 1);
                }

        }

        finish_construction();
        return 0;
}

inline int parallel_graph_access::build_from_metis_weighted(int n, int* xadj, int* adjncy, int * vwgt, int* adjwgt) {
        start_construction(n, xadj[n], n, xadj[n]);
        set_range(0,n);

        for( unsigned i = 0; i < (unsigned)n; i++) {
                NodeID node = new_node();
                setNodeWeight(node, vwgt[i]);
                setNodeLabel(node, 0);

                for( unsigned e = xadj[i]; e < (unsigned)xadj[i+1]; e++) {
                        EdgeID e_bar = new_edge(node, adjncy[e]);
                        setEdgeWeight(e_bar, adjwgt[e]);
                }
        }

        finish_construction();
        return 0;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%% Handle Communication of Ghost Node Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inline 
void ghost_node_communication::addLabel(NodeID node, NodeID label) {
        forall_out_edges((*m_G), e, node) {
                NodeID target = m_G->getEdgeTarget(e);
                if( !m_G->is_local_node(target)  ) {
                        PEID peID = m_G->getTargetPE(target);
                        if( !m_PE_packed[peID] ) { // make sure a node is sent at most once
                                (*m_send_buffers_ptr)[peID].push_back(m_G->getGlobalID(node));
                                (*m_send_buffers_ptr)[peID].push_back(label);
                                m_PE_packed[peID] = true;
                        }
                }
        } endfor
        forall_out_edges((*m_G), e, node) {
                NodeID target = m_G->getEdgeTarget(e);
                if( !m_G->is_local_node(target)  ) {
                        m_PE_packed[m_G->getTargetPE(target)] = false;
                }
        } endfor
}

// we want to interleave computation and communication
// check_iteration_counter default is true
inline void ghost_node_communication::update_ghost_node_data( bool check_iteration_counter ) {
        if( check_iteration_counter ) {
                if( ++m_iteration_counter <= m_skip_limit || m_size == 1 ) return; 
        }
        m_iteration_counter = 0;
        m_send_iteration++;
        m_send_tag++;

        //send all neighbors their packages using Isends
        //a neighbor that does not receive something gets a specific token
        if( m_first_send ) {
                for( PEID peID = 0; peID < m_size; peID++) {
                        if( m_adjacent_processors[peID] ) {
                                //now we have to send a message
                                if( (*m_send_buffers_ptr)[peID].size() == 0 ){
                                        // length 1 encode no message
                                        (*m_send_buffers_ptr)[peID].push_back(0);
                                }

                                MPI_Request * request = new MPI_Request();
                                MPI_Isend( &(*m_send_buffers_ptr)[peID][0], 
                                           (*m_send_buffers_ptr)[peID].size(), 
                                           MPI_UNSIGNED_LONG_LONG, 
                                           peID, m_send_tag, m_communicator, request);
                                
                                m_isend_requests.push_back( request );
                        }
                }

                m_first_send = false;
                if( m_send_buffers_ptr == & m_send_buffers_A) {
                        m_send_buffers_ptr = & m_send_buffers_B;
                } else {
                        m_send_buffers_ptr = & m_send_buffers_A;
                }
                return; // compute a little bit more
        }

        m_G->update_block_weights();

        //receive incomming
        receive_messages_of_neighbors();
       
        if( m_send_buffers_ptr == & m_send_buffers_A) {
                for( int i = 0; i < m_size; i++) {
                        m_send_buffers_B[i].clear();
                }
        } else {
                for( int i = 0; i < m_size; i++) {
                        m_send_buffers_A[i].clear();
                }
        }

        for( PEID peID = 0; peID < (PEID)(*m_send_buffers_ptr).size(); peID++) {
                if( m_adjacent_processors[peID] ) {
                        //now we have to send a message
                        if( (*m_send_buffers_ptr)[peID].size() == 0 ){
                                // length 1 encode no message
                                (*m_send_buffers_ptr)[peID].push_back(0);
                        }

                        MPI_Request * request = new MPI_Request();
                        MPI_Isend( &(*m_send_buffers_ptr)[peID][0], 
                                   (*m_send_buffers_ptr)[peID].size(), 
                                   MPI_UNSIGNED_LONG_LONG, 
                                   peID, m_send_tag, m_communicator, request);
                        
                        m_isend_requests.push_back( request );
                }
        }

        // switch send buffers
        if( m_send_buffers_ptr == & m_send_buffers_A) {
                m_send_buffers_ptr = & m_send_buffers_B;
        } else {
                m_send_buffers_ptr = & m_send_buffers_A;
        }

}

inline 
void ghost_node_communication::receive_messages_of_neighbors() {
        PEID counter = 0;
        m_recv_iteration++;
        m_recv_tag++;
        while( counter < m_num_adjacent ) {
                // wait for incomming message of an adjacent processor
                MPI_Status st;
                MPI_Probe(MPI_ANY_SOURCE, m_recv_tag, m_communicator,  &st);

                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);

                std::vector<NodeID> message; message.resize(message_length);
                MPI_Status rst;
                MPI_Recv( &message[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, m_recv_tag, m_communicator, &rst); 

                counter++;

                // now integrate the changes
                if(message_length == 1) continue; // nothing to do

                for( int i = 0; i < message_length-1; i+=2) {
                        NodeID global_id = message[i];
                        NodeID label     = message[i+1];

                        NodeID local_id = m_G->m_global_to_local_id[global_id];
                        m_G->update_non_contained_block_balance(m_G->getNodeLabel(local_id), label, m_G->getNodeWeight(local_id));
                        m_G->setNodeLabel(local_id, label);
                }
        }

        // wait for previous iteration to finish
        for( unsigned i = 0; i < m_isend_requests.size(); i++) {
                MPI_Status st;
                MPI_Wait( m_isend_requests[i], &st);
                delete m_isend_requests[i];
        }
        m_isend_requests.clear();

}

inline void ghost_node_communication::update_ghost_node_data_finish() {
        while( m_send_iteration < m_desired_rounds) {
                // we have to do another send
                update_ghost_node_data( false ); // flush the lokal buffers to our neighbors
        }
        
        while( m_recv_iteration < m_desired_rounds) {
                receive_messages_of_neighbors(); // last receive
        }
        
        m_first_send = true; // last send
        update_ghost_node_data(false);
        m_G->update_block_weights();
        receive_messages_of_neighbors();

        m_send_iteration = 0;
        m_recv_iteration = 0;

        m_send_tag = 100*m_size-1;
        m_recv_tag = 100*m_size-1;
        m_first_send     = true;

        for( int i = 0; i < m_size; i++) {
                m_send_buffers_B[i].clear();
        }

        for( int i = 0; i < m_size; i++) {
                m_send_buffers_A[i].clear();
        }

        MPI_Barrier(m_communicator);
        
}

inline void ghost_node_communication::update_ghost_node_data_global() {
        std::vector< std::vector< NodeID > > send_buffers; // buffers to send messages
        send_buffers.resize(m_size);
        forall_local_nodes((*m_G), node) {
                forall_out_edges((*m_G), e, node) {
                        NodeID target = m_G->getEdgeTarget(e);
                        if( !m_G->is_local_node(target)  ) {
                                PEID peID = m_G->getTargetPE(target);
                                if( !m_PE_packed[peID] ) { // make sure a node is sent at most once
                                        send_buffers[peID].push_back(m_G->getGlobalID(node));
                                        send_buffers[peID].push_back(m_G->getNodeLabel(node));
                                        m_PE_packed[peID] = true;
                                }
                        }
                } endfor
                forall_out_edges((*m_G), e, node) {
                        NodeID target = m_G->getEdgeTarget(e);
                        if( !m_G->is_local_node(target)  ) {
                                m_PE_packed[m_G->getTargetPE(target)] = false;
                        }
                } endfor
        } endfor

        //send all neighbors their packages using Isends
        //a neighbor that does not receive something gets a specific token
        for( PEID peID = 0; peID < (PEID)send_buffers.size(); peID++) {
                if( m_adjacent_processors[peID] ) {
                        //now we have to send a message
                        if( send_buffers[peID].size() == 0 ){
                                // length 1 encode no message
                                send_buffers[peID].push_back(0);
                        }

                        MPI_Request rq; 
                        MPI_Isend( &send_buffers[peID][0], 
                                    send_buffers[peID].size(), MPI_UNSIGNED_LONG_LONG, peID, peID+3*m_size, m_communicator, &rq);
                }
        }

        //receive incomming
        PEID counter = 0;
        while( counter < m_num_adjacent ) {
                // wait for incomming message of an adjacent processor
                MPI_Status st; unsigned int tag = m_rank+3*m_size;
                MPI_Probe(MPI_ANY_SOURCE, tag, m_communicator,  &st);

                int message_length;
                MPI_Get_count(&st, MPI_UNSIGNED_LONG_LONG, &message_length);
                std::vector<NodeID> message; message.resize(message_length);

                MPI_Status rst;
                MPI_Recv( &message[0], message_length, MPI_UNSIGNED_LONG_LONG, st.MPI_SOURCE, tag, m_communicator, &rst); 
                counter++;

                // now integrate the changes
                if(message_length == 1) continue; // nothing to do

                for( int i = 0; i < message_length-1; i+=2) {
                        NodeID global_id = message[i];
                        NodeID label     = message[i+1];

                        m_G->setNodeLabel( m_G->m_global_to_local_id[global_id], label);
                }
        }

        MPI_Barrier(m_communicator);
}

#endif /* end of include guard: PARALLEL_GRAPH_ACCESS_X6O9MRS8 */
