/******************************************************************************
 * spac.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Author: Daniel Seemaier <daniel.seemaier@student.kit.edu>
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/
#include "spac.h"

spac::spac(graph_access &input_graph, EdgeWeight infinity)
        : m_input_graph(input_graph),  m_split_graph(), m_infinity(infinity), m_reverse_edge() {
}

graph_access &spac::construct_split_graph() {
    find_reverse_edges();

    // find number of degree 1 and 2 vertices
    NodeID number_of_deg1_vertices = 0;
    NodeID number_of_deg2_vertices = 0;
    for (NodeID v = 0; v < m_input_graph.number_of_nodes(); ++v) {
        //NodeID deg = m_input_graph.getNodeDegree(v);

        if (m_input_graph.getNodeDegree(v) == 1) {
            ++number_of_deg1_vertices;
        } else if (m_input_graph.getNodeDegree(v) == 2) {
            ++number_of_deg2_vertices;
        }
    }

    // calculation of the dimension of the split graph
    NodeID split_n = m_input_graph.number_of_edges(); // two times the number of edges
    EdgeID split_m = 3 * m_input_graph.number_of_edges() - 2 * (number_of_deg1_vertices + number_of_deg2_vertices);
    m_split_graph.start_construction(split_n, split_m);

    for (NodeID u = 0; u < m_input_graph.number_of_nodes(); ++u) {
        NodeID deg = m_input_graph.getNodeDegree(u);

        for (EdgeID e = m_input_graph.get_first_edge(u); e < m_input_graph.get_first_invalid_edge(u); ++e) {
            EdgeID nth_edge_at_u = (e - m_input_graph.get_first_edge(u));

            NodeID split_node = m_split_graph.new_node();
            assert(e == split_node);
            m_split_graph.setNodeWeight(split_node, 1);

            EdgeID dominant_edge = m_split_graph.new_edge(split_node, m_reverse_edge[e]);
            m_split_graph.setEdgeWeight(dominant_edge, m_infinity);

            if (deg == 2) {
                EdgeID auxiliary_edge;
                if (nth_edge_at_u == 0) { // split split vertex
                    auxiliary_edge = m_split_graph.new_edge(split_node, split_node + 1);
                } else { // second split vertex
                    assert(e + 1 == m_input_graph.get_first_invalid_edge(u));
                    assert(split_node > 0);

                    auxiliary_edge = m_split_graph.new_edge(split_node, split_node - 1);
                }
                m_split_graph.setEdgeWeight(auxiliary_edge, 1);
            } else if (deg > 2) {
                // calculate offsets between split_node and the next / prev node in the cycle
                int prev_offset = -1;
                int next_offset = 1;
                if (nth_edge_at_u == 0) {
                    prev_offset = deg - 1;
                } else if (nth_edge_at_u == deg - 1) {
                    next_offset = -(deg - 1);
                }

                m_split_graph.setEdgeWeight(m_split_graph.new_edge(split_node, split_node + prev_offset), 1);
                m_split_graph.setEdgeWeight(m_split_graph.new_edge(split_node, split_node + next_offset), 1);
            } else {
                // nothing to do for leaves
                assert(deg == 1);
            }
        }
    }

#ifndef NDEBUG
    for (NodeID u = 0; u < m_split_graph.number_of_nodes(); ++u) {
        assert(m_split_graph.getNodeDegree(u) > 0);
        assert(m_split_graph.getNodeDegree(u) <= 3);

        NodeID v = m_split_graph.getEdgeTarget(m_split_graph.get_first_edge(u));
        assert(u == m_split_graph.getEdgeTarget(m_split_graph.get_first_edge(v)));
    }
#endif

    m_split_graph.finish_construction();
    return m_split_graph;
}

void spac::fix_cut_dominant_edges() {
    EdgeID number_of_bad_edges = 0;

    // check if there are bad edges, i.e. dominant edges with endpoints in different blocks
    for (NodeID u = 0; u < m_split_graph.number_of_nodes(); ++u) {
        for (EdgeID e = m_split_graph.get_first_edge(u); e < m_split_graph.get_first_invalid_edge(u); ++e) {
            NodeID v = m_split_graph.getEdgeTarget(e);

            PartitionID u_part = m_split_graph.getPartitionIndex(u);
            PartitionID v_part = m_split_graph.getPartitionIndex(v);

            if (m_split_graph.getEdgeWeight(e) > 1 && u_part != v_part) {
                ++number_of_bad_edges;
            }
        }
    }

    // if there are bad edges, fix them
    if (number_of_bad_edges > 0) {
        std::vector<NodeID> partition_sizes(m_split_graph.get_partition_count());
        for (NodeID u = 0; u < m_split_graph.number_of_nodes(); ++u) {
            ++partition_sizes[m_split_graph.getPartitionIndex(u)];
        }

        for (NodeID u = 0; u < m_split_graph.number_of_nodes(); ++u) {
            for (EdgeID e = m_split_graph.get_first_edge(u); e < m_split_graph.get_first_invalid_edge(u); ++e) {
                NodeID v = m_split_graph.getEdgeTarget(e);

                PartitionID u_part = m_split_graph.getPartitionIndex(u);
                PartitionID v_part = m_split_graph.getPartitionIndex(v);

                // move one endpoint to the smaller partition
                if (m_split_graph.getEdgeWeight(e) > 1 && u_part != v_part) {
                    if (partition_sizes[u_part] < partition_sizes[v_part]) {
                        ++partition_sizes[u_part];
                        --partition_sizes[v_part];
                        m_split_graph.setPartitionIndex(v, u_part);
                    } else {
                        --partition_sizes[u_part];
                        ++partition_sizes[v_part];
                        m_split_graph.setPartitionIndex(u, v_part);
                    }
                }
            }
        }

    }
}

std::vector<PartitionID> spac::project_partition() {
    std::vector<PartitionID> edge_partition(m_input_graph.number_of_edges());

    for (NodeID u = 0; u < m_input_graph.number_of_nodes(); ++u) {
        for (EdgeID e = m_input_graph.get_first_edge(u); e < m_input_graph.get_first_invalid_edge(u); ++e) {
            edge_partition[e] = m_split_graph.getPartitionIndex(e);
        }
    }

    return edge_partition;
}

unsigned spac::calculate_vertex_cut(const std::vector<PartitionID> &edge_partition) {
    unsigned cost = 0;

    for (NodeID u = 0; u < m_input_graph.number_of_nodes(); ++u) {
        if (m_input_graph.getNodeDegree(u) == 0) continue;
        std::vector<bool> counted(m_split_graph.get_partition_count());

        for (EdgeID e = m_input_graph.get_first_edge(u); e < m_input_graph.get_first_invalid_edge(u); ++e) {
            PartitionID part = edge_partition[e];
            if (!counted[part]) {
                counted[part] = true;
                ++cost;
            }
        }

        --cost;
    }

    return cost;
}

void spac::find_reverse_edges() {
    // reverse edges were already calculated, thus nothing to do
    if (!m_reverse_edge.empty()) {
        return;
    }

    const EdgeID guard = std::numeric_limits<EdgeID>::max();
    m_reverse_edge.resize(m_input_graph.number_of_edges(), guard);

    for (NodeID u = 0; u < m_input_graph.number_of_nodes(); ++u) {
        for (EdgeID e_uv = m_input_graph.get_first_edge(u); e_uv < m_input_graph.get_first_invalid_edge(u); ++e_uv) {
            NodeID v = m_input_graph.getEdgeTarget(e_uv);

            if (u < v) {
                assert(m_reverse_edge[e_uv] == guard);
            } else {
                assert(m_reverse_edge[e_uv] != guard);
                continue;
            }

#ifndef NDEBUG
            bool found_reverse_edge = false;
#endif
            for (EdgeID e_vu = m_input_graph.get_first_edge(v); e_vu < m_input_graph.get_first_invalid_edge(v); ++e_vu) {
                if (m_input_graph.getEdgeTarget(e_vu) == u) {
                    m_reverse_edge[e_uv] = e_vu;
                    m_reverse_edge[e_vu] = e_uv;
#ifndef NDEBUG
                    found_reverse_edge = true;
#endif
                    break;
                }
            }

            assert(found_reverse_edge);
            assert(m_reverse_edge[e_uv] != guard);
        }
    }
}
