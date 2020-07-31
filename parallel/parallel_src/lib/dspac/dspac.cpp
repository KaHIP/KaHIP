/******************************************************************************
 * dspac.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Author: Daniel Seemaier <daniel.seemaier@student.kit.edu>
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include "dspac.h"

dspac::dspac(parallel_graph_access &graph, MPI_Comm comm, EdgeWeight infinity)
        : m_comm(comm), m_infinity(infinity), m_input_graph(graph) {
}

void dspac::construct(parallel_graph_access &split_graph) {
    assert(assert_adjacency_lists_sorted());
    MPI_Barrier(m_comm);
    internal_construct(split_graph);
    MPI_Barrier(m_comm);
    assert(assert_sanity_checks(split_graph));
}

void dspac::internal_construct(parallel_graph_access &split_graph) {
    int size, rank;
    MPI_Comm_size(m_comm, &size);
    MPI_Comm_rank(m_comm, &rank);
    const NodeID n = m_input_graph.number_of_global_nodes();
#ifndef NDEBUG
    const NodeID m = m_input_graph.number_of_global_edges();
#endif

    timer construction_timer;

    // we construct the split nodes from..(to - 1) on this node
    auto edge_range_array = m_input_graph.get_edge_range_array();
    assert(assert_edge_range_array_ok(edge_range_array));
    const std::size_t from = edge_range_array[rank]; // inclusive
    const std::size_t to = edge_range_array[rank + 1]; // exclusive
    assert(to <= m_input_graph.number_of_global_edges());

    // we need the number of vertices of degree 1 or 2 to calculate the dimension of the split graph
    NodeID local_number_of_deg_1_or_2_vertices = 0;
    for (NodeID v = 0; v < m_input_graph.number_of_local_nodes(); ++v) {
        EdgeID deg = m_input_graph.getNodeDegree(v);
        if (deg == 1 || deg == 2) {
            ++local_number_of_deg_1_or_2_vertices;
        }
    }

    NodeID global_number_of_deg_1_or_2_vertices = 0;
    MPI_Allreduce(&local_number_of_deg_1_or_2_vertices, &global_number_of_deg_1_or_2_vertices, 1,
                  MPI_UNSIGNED_LONG_LONG, MPI_SUM, m_comm);
    if (rank == 0) {
        std::cout << "[dspac::internal_construct()] Up to MPI_Allreduce() took "
                  << construction_timer.elapsed() << std::endl;
        construction_timer.restart();
    }

    // calculate split graph dimensions
    const NodeID local_number_of_split_nodes = m_input_graph.number_of_local_edges();
    const EdgeID global_number_of_split_nodes = m_input_graph.number_of_global_edges();

    assert(3 * m_input_graph.number_of_local_edges() >= 2 * local_number_of_deg_1_or_2_vertices);
    const EdgeID local_number_of_split_edges = 3 * m_input_graph.number_of_local_edges()
            - 2 * local_number_of_deg_1_or_2_vertices;

    assert(3 * m_input_graph.number_of_global_edges() >= 2 * global_number_of_deg_1_or_2_vertices);
    const NodeID global_number_of_split_edges = 3 * m_input_graph.number_of_global_edges()
            - 2 * global_number_of_deg_1_or_2_vertices;

    // this array stores the distribution of nodes across PEs, namely PE i stores nodes
    // node_range_array[i]..node_range_array[i + 1]-1
    // the default loader sets a wrong value for node_range_array[size] though, so we need to fix that for our purposes
    // here
    auto node_range_array = m_input_graph.get_range_array();
    node_range_array[size] = m_input_graph.number_of_global_nodes();
    assert(assert_node_range_array_ok(node_range_array));

    std::vector<std::vector<NodeID>> first_split_node_on(size);

    // first, reserve memory for adjacent PEs
    for (PEID pe = 0; pe < size; ++pe) {
        if (m_input_graph.is_adjacent_PE(pe) || pe == rank) {
            first_split_node_on[pe].resize(m_input_graph.number_of_local_nodes());
        }
    }

    // then fill the reserved memory
    for (NodeID v = 0; v < m_input_graph.number_of_local_nodes(); ++v) {
        PEID current_pe = -1;
        for (EdgeID e = m_input_graph.get_first_edge(v); e < m_input_graph.get_first_invalid_edge(v); ++e) {
            NodeID u = m_input_graph.getEdgeTarget(e);
            PEID pe = m_input_graph.is_local_node(u) ? static_cast<PEID>(rank) : m_input_graph.getTargetPE(u);

            if (pe != current_pe) {
                assert(first_split_node_on[pe].size() == m_input_graph.number_of_local_nodes());
                first_split_node_on[pe][v] = from + e;
                current_pe = pe;
            }
        }
    }

    if (rank == 0) {
        std::cout << "[dspac::internal_construct()] Preparation of first_split_node_on[] took "
                  << construction_timer.elapsed() << std::endl;
        construction_timer.restart();
    }

    // once created, this array has the following semantic: say we have an edge vu in the original graph where
    // v is on our PE and u is on any PE
    // when we create the split graph, when need to connect one split vertex of v and one of u to represent the vu edge
    // in the split graph
    // so when we connect the split vertices of v, we use first_split_node[globalId(u)] as the split node id of u and
    // then increment it by one, so that when need a split vertex of u again on this PE, we use the next one and so on
    std::vector<NodeID> first_split_node(n); // contains global node ids

    // receive the messages from adjacent PEs and place them at the right position in first_split_node: the messages
    // from PE i should be placed starting at node_range_array[i]
    std::vector<MPI_Request *> requests;

    // send the messages to adjacent PEs
    for (PEID pe = 0; pe < size; ++pe) {
        if (m_input_graph.is_adjacent_PE(pe)) {
            assert(rank != pe);

            NodeID *buf = &first_split_node_on[pe][0];
            const std::size_t count = first_split_node_on[pe].size();

            assert(count == node_range_array[rank + 1] - node_range_array[rank]);
            assert(count < std::numeric_limits<int>::max());

            MPI_Request *request = new MPI_Request;
            MPI_Isend(buf, static_cast<int>(count), MPI_UNSIGNED_LONG_LONG, pe, 0, m_comm, request);
            requests.push_back(request);
        }
    }

    // copy own data from first_split_node_on to first_split_node
    assert(first_split_node_on[rank].size() == m_input_graph.number_of_local_nodes());
    assert(first_split_node_on[rank].size() == node_range_array[rank + 1] - node_range_array[rank]);
    assert(first_split_node.data() + node_range_array[rank] + m_input_graph.number_of_local_nodes()
           <= (&first_split_node[n - 1]) + 1);
    std::copy(first_split_node_on[rank].begin(), first_split_node_on[rank].end(),
              first_split_node.begin() + node_range_array[rank]);

    // receive messages from adjacent neighbors
    for (PEID pe = 0; pe < size; ++pe) {
        if (m_input_graph.is_adjacent_PE(pe)) {
            assert(rank != pe);

            NodeID *buf = &first_split_node[node_range_array[pe]];
            const NodeID count = node_range_array[pe + 1] - node_range_array[pe];

            assert(node_range_array[pe] + count <= first_split_node.size());
            assert(count < std::numeric_limits<int>::max());

            MPI_Recv(buf, static_cast<int>(count), MPI_UNSIGNED_LONG_LONG, pe, 0, m_comm, MPI_STATUS_IGNORE);
        }
    }

    // wait for own messages to be received
    for (MPI_Request *request : requests) {
        MPI_Wait(request, MPI_STATUS_IGNORE);
        delete request;
    }

    if (rank == 0) {
        std::cout << "[dspac::internal_construct()] first_split_node[] communication took "
                  << construction_timer.elapsed() << std::endl;
        construction_timer.restart();
    }

    // we no longer need first_split_node_on from now on since it's copied to first_split_node on each PE
    first_split_node_on.clear();

    // now we construct the split graph
    split_graph.start_construction(local_number_of_split_nodes, local_number_of_split_edges,
                                   global_number_of_split_nodes, global_number_of_split_edges);
    split_graph.set_range_array(edge_range_array);
    split_graph.set_range(from, to - 1);

    NodeID nodes_created = 0;
    EdgeID edges_created = 0;

    for (NodeID v = 0; v < m_input_graph.number_of_local_nodes(); ++v) {
        EdgeID deg = m_input_graph.getNodeDegree(v);
        if (deg == 0) { // explicitly skip isolated nodes
            continue;
        }

        for (EdgeID e = m_input_graph.get_first_edge(v); e < m_input_graph.get_first_invalid_edge(v); ++e) {
            NodeID u = m_input_graph.getEdgeTarget(e);
            NodeID global_u = m_input_graph.getGlobalID(u);

            // create the split node
            ++nodes_created;
            NodeID split_node = split_graph.new_node();
            assert(split_node == e);

            split_graph.setNodeWeight(split_node, 1);
            split_graph.setNodeLabel(split_node, from + split_node);
            split_graph.setSecondPartitionIndex(split_node, 0);

            // create dominant edge
            ++edges_created;
            assert(global_u < first_split_node.size());
            NodeID target_node = first_split_node[global_u];
            EdgeID dominant_edge = split_graph.new_edge(split_node, target_node);
            ++first_split_node[global_u];
            split_graph.setEdgeWeight(dominant_edge, m_infinity);

            // create auxiliary edges
            bool first = (e == m_input_graph.get_first_edge(v));
            bool last = (e + 1 == m_input_graph.get_first_invalid_edge(v));

            if (deg == 2) {
                // degree 2: we create a path with a single edge in the split graph
                ++edges_created;
                int target_offset = first ? 1 : -1;
                assert(0 <= split_node + target_offset && split_node + target_offset < local_number_of_split_nodes);
                EdgeID auxiliary_edge = split_graph.new_edge(split_node, from + split_node + target_offset);
                split_graph.setEdgeWeight(auxiliary_edge, 1);
            } else if (deg > 2) {
                // degree > 2: we create a cycle with all split nodes, thus we need a edge to the previous and one to
                // the next node in the cycle
                ++edges_created;
                int next_offset = last ? -(static_cast<int>(deg) - 1) : 1;
                NodeID global_next = from + split_node + next_offset;
                assert(from == split_graph.get_from_range());
                assert(split_graph.get_from_range() <= global_next && global_next <= split_graph.get_to_range());

                EdgeID next_auxiliary_edge = split_graph.new_edge(split_node, global_next);
                split_graph.setEdgeWeight(next_auxiliary_edge, 1);

                ++edges_created;
                int prev_offset = first ? static_cast<int>(deg) - 1 : -1;
                NodeID global_prev = from + split_node + prev_offset;
                assert(split_graph.get_from_range() <= global_prev && global_prev <= split_graph.get_to_range());
                EdgeID prev_auxiliary_edge = split_graph.new_edge(split_node, global_prev);
                split_graph.setEdgeWeight(prev_auxiliary_edge, 1);
            } else {
                assert(deg == 1);
                // nothing to do for leaves
            }
        }
    }

    if (rank == 0) {
        std::cout << "[dspac::internal_construct()] Local construction took "
                  << construction_timer.elapsed() << std::endl;
        construction_timer.restart();
    }

    assert(nodes_created == local_number_of_split_nodes);
    assert(edges_created == local_number_of_split_edges);
    split_graph.finish_construction();
}

/**
 * assert()'s some sanity checks on the split graph.
 * @return Pointless bool so that the method call can be used as expression.
 */
bool dspac::assert_sanity_checks(parallel_graph_access &split_graph) {
#ifndef NDEBUG
    assert(split_graph.number_of_local_nodes() == m_input_graph.number_of_local_edges());
    for (NodeID v = 0; v < split_graph.number_of_local_nodes(); ++v) {
        // isolated vertices should be removed for now
        assert(0 < split_graph.getNodeDegree(v) && split_graph.getNodeDegree(v) <= 3);

        // make sure that the edge weights are correct, i.e. auxiliary edges have edge weight 1 and
        // dominant edges have edge weight m_infinity
        EdgeID firstEdge = split_graph.get_first_edge(v);
        switch (split_graph.getNodeDegree(v)) {
            case 3: // fall through intended
                assert(split_graph.getEdgeWeight(firstEdge + 2) == 1);

            case 2:
                assert(split_graph.getEdgeWeight(firstEdge + 1) == 1);

            case 1:
                assert(split_graph.getEdgeWeight(firstEdge) == m_infinity);
                break;

            default:
                assert(false);
        }
    }

    // this part checks that the auxiliary edges are connected to the right nodes
    for (NodeID v = 0; v < m_input_graph.number_of_local_nodes(); ++v) {
        EdgeID deg = m_input_graph.getNodeDegree(v);
        if (deg == 0) { // explicitly skip isolated nodes
            continue;
        }

        for (EdgeID e = m_input_graph.get_first_edge(v); e < m_input_graph.get_first_invalid_edge(v); ++e) {
            bool first = (e == m_input_graph.get_first_edge(v));
            bool last = (e + 1 == m_input_graph.get_first_invalid_edge(v));

            if (deg == 1) {
                if (split_graph.get_first_edge(e) + 1 < split_graph.number_of_local_edges()) {
                    // degree 1 node --> no auxiliary edges --> next edge must be a dominant edge of another node
                    assert(split_graph.getEdgeWeight(split_graph.get_first_edge(e) + 1) == m_infinity);
                }
            } else if (deg == 2) {
                if (first) {
                    // first split node --> auxiliary edge must target the second split node
                    assert(split_graph.getEdgeTarget(split_graph.get_first_edge(e) + 1) == e + 1);
                } else if (last) {
                    // second split node --> auxiliary edge must target the first split node
                    assert(split_graph.getEdgeTarget(split_graph.get_first_edge(e) + 1) == e - 1);
                } else {
                    assert(false);
                }

                // a dominant edge must follow a single auxiliary edge
                if (split_graph.get_first_edge(e) + 2 < split_graph.number_of_local_edges()) {
                    assert(split_graph.getEdgeWeight(split_graph.get_first_edge(e) + 2) == m_infinity);
                }
            } else if (deg > 2) {
                if (first) {
                    assert(split_graph.getEdgeTarget(split_graph.get_first_edge(e) + 1) == e + 1);
                    assert(split_graph.getEdgeTarget(split_graph.get_first_edge(e) + 2) == e + (deg - 1));
                } else if (last) {
                    assert(split_graph.getEdgeTarget(split_graph.get_first_edge(e) + 1) == e - (deg - 1));
                    assert(split_graph.getEdgeTarget(split_graph.get_first_edge(e) + 2) == e - 1);
                } else {
                    assert(split_graph.getEdgeTarget(split_graph.get_first_edge(e) + 1) == e + 1);
                    assert(split_graph.getEdgeTarget(split_graph.get_first_edge(e) + 2) == e - 1);
                }

                // a dominant edge must follow after two auxiliary edges
                if (split_graph.get_first_edge(e) + 3 < split_graph.number_of_local_edges()) {
                    assert(split_graph.getEdgeWeight(split_graph.get_first_edge(e) + 3) == m_infinity);
                }
            } else {
                assert(false);
            }
        }
    }
#endif
    return true;
}

/**
 * assert()'s that the adjacency lists of the input graph are sorted.
 * @return Pointless bool so that the method call can be used as expression.
 */
bool dspac::assert_adjacency_lists_sorted() {
#ifndef NDEBUG
    for (NodeID v = 0; v < m_input_graph.number_of_local_nodes(); ++v) {
        if (m_input_graph.getNodeDegree(v) == 0) {
            continue;
        }

        NodeID local_first_neighbor = m_input_graph.getEdgeTarget(m_input_graph.get_first_edge(v));
        auto global_first_neighbor = static_cast<NodeID>(m_input_graph.getGlobalID(local_first_neighbor));
        NodeID cur = global_first_neighbor;
        for (EdgeID e = m_input_graph.get_first_edge(v); e < m_input_graph.get_first_invalid_edge(v); ++e) {
            NodeID u = m_input_graph.getEdgeTarget(e);
            assert(cur <= m_input_graph.getGlobalID(u));
        }
    }
#endif
    return true;
}

bool dspac::assert_edge_range_array_ok(const std::vector<NodeID> &edge_range_array) {
    int size, rank;
    MPI_Comm_size(m_comm, &size);
    MPI_Comm_rank(m_comm, &rank);
    assert(edge_range_array.size() == size + 1);
    assert(edge_range_array[0] == 0);
    assert(edge_range_array[size] == m_input_graph.number_of_global_edges());
    assert(m_input_graph.number_of_local_edges() == edge_range_array[rank + 1] - edge_range_array[rank]);
    for (std::size_t pe = 0; pe < (size_t)size; ++pe)
        assert(edge_range_array[pe] <= edge_range_array[pe + 1]);
    return true;
}

bool dspac::assert_node_range_array_ok(const std::vector<NodeID> &node_range_array) {
    int size, rank;
    MPI_Comm_size(m_comm, &size);
    MPI_Comm_rank(m_comm, &rank);
    assert(node_range_array.size() == size + 1);
    assert(node_range_array[0] == 0);
    assert(node_range_array[size] == m_input_graph.number_of_global_nodes());
    assert(m_input_graph.number_of_local_nodes() == node_range_array[rank + 1] - node_range_array[rank]);
    return true;
}

std::vector<PartitionID> dspac::project_partition(parallel_graph_access &split_graph) {
    std::vector<PartitionID> edge_partition(m_input_graph.number_of_local_edges());

    for (NodeID v = 0; v < m_input_graph.number_of_local_nodes(); ++v) {
        for (EdgeID e = m_input_graph.get_first_edge(v); e < m_input_graph.get_first_invalid_edge(v); ++e) {
            edge_partition[e] = split_graph.getNodeLabel(e);
        }
    }

    MPI_Barrier(m_comm);
    return edge_partition;
}

EdgeWeight dspac::calculate_vertex_cut(PartitionID k, const std::vector<PartitionID> &edge_partition) {
    EdgeWeight local_cost = 0;

    std::vector<bool> counted(k);
    for (NodeID v = 0; v < m_input_graph.number_of_local_nodes(); ++v) {
        if (m_input_graph.getNodeDegree(v) == 0) {
            continue;
        }

        for (EdgeID e = m_input_graph.get_first_edge(v); e < m_input_graph.get_first_invalid_edge(v); ++e) {
            PartitionID p = edge_partition[e];
            if (!counted[p]) {
                counted[p] = true;
                ++local_cost;
            }
        }

        counted.clear();
        counted.resize(k);

        assert(local_cost > 0);
        --local_cost;
    }

    EdgeWeight global_cost;
    MPI_Reduce(&local_cost, &global_cost, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, m_comm);
    return global_cost;
}

void dspac::fix_cut_dominant_edges(parallel_graph_access &split_graph) {
    for (NodeID v = 0; v < split_graph.number_of_local_nodes(); ++v) {
        EdgeID e_vu = split_graph.get_first_edge(v);
        NodeID u = split_graph.getEdgeTarget(e_vu);

        PartitionID part_v = split_graph.getNodeLabel(v);
        PartitionID part_u = split_graph.getNodeLabel(u);
        if (part_v != part_u) {
            NodeWeight part_v_size = split_graph.getBlockSize(part_v);
            NodeWeight part_u_size = split_graph.getBlockSize(part_u);

            if (part_v_size < part_u_size) {
                split_graph.setNodeLabel(u, part_v);
                split_graph.setBlockSize(part_v, part_v_size + 1);
                split_graph.setBlockSize(part_u, part_u_size - 1);
            } else {
                split_graph.setNodeLabel(v, part_u);
                split_graph.setBlockSize(part_v, part_v_size - 1);
                split_graph.setBlockSize(part_u, part_u_size + 1);
            }
        }
    }
    split_graph.update_block_weights();
}
