/******************************************************************************
 * edge_balanced_graph_io.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#include <algorithm>
#include <array>

#include "edge_balanced_graph_io.h"

static constexpr ULONG FILE_TYPE_VERSION = 3;

static constexpr ULONG HEADER_SIZE = 3;

static ULONG calculateFromNode(std::ifstream &in, ULONG numberOfNodes, ULONG numberOfEdges, int rank, int size);

static ULONG calculateToNode(std::ifstream &in, ULONG numberOfNodes, ULONG numberOfEdges, int rank, int size);

static ULONG readNumberOfEdgesInRange(std::ifstream &in, ULONG numberOfNodes, ULONG from, ULONG to);

static ULONG readFirstEdge(std::ifstream &in, ULONG numberOfNodes, ULONG node);

static ULONG readFirstInvalidEdge(std::ifstream &in, ULONG numberOfNodes, ULONG node);

static ULONG adjacencyListOffsetToEdgeID(ULONG numberOfNodes, ULONG offset);

void edge_balanced_graph_io::read_binary_graph_edge_balanced(parallel_graph_access &G, const std::string &filename,
                                                             const PPartitionConfig &config) {
    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    read_binary_graph_edge_balanced(G, filename, config, rank, size);
}

void edge_balanced_graph_io::read_binary_graph_edge_balanced(parallel_graph_access &G, const std::string &filename,
                                                             const PPartitionConfig &config, int rank, int size) {
    // header[0] = version number
    // header[1] = number of nodes
    // header[2] = number of edges
    std::array<ULONG, HEADER_SIZE> header{0};

    // read header on ROOT; if it fails, throw, otherwise broadcast this to other PEs
    int success = 0;
    if (rank == ROOT) {
        std::ifstream headerIn(filename, std::ios::binary | std::ios::in);
        if (headerIn) {
            success = 1;
            headerIn.read((char *) (&header[0]), 3 * sizeof(ULONG));
        }
        headerIn.close();
    }

    MPI_Bcast(&success, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    if (success != 1) {
        throw std::ios_base::failure("unable to read graph file");
    }

    MPI_Bcast(&header[0], header.size(), MPI_UNSIGNED_LONG_LONG, ROOT, MPI_COMM_WORLD);
    ULONG version = header[0];
    ULONG n = header[1];
    ULONG m = header[2];

    if (rank == ROOT) {
        std::cout << "n=" << n << ", m=" << m << std::endl;
    }

    if (version != FILE_TYPE_VERSION) {
        throw std::ios_base::failure("wrong file type version");
    }

    /*
     * next, we determine the number of vertices on each PE such that the number of edges are almost evenly
     * distributed
     */
    std::ifstream in(filename, std::ios::binary | std::ios::in);
    in.exceptions(std::ios_base::failbit | std::ios_base::badbit);
    ULONG from = calculateFromNode(in, n, m, rank, size); // inclusive!
    ULONG to = calculateToNode(in, n, m, rank, size); // inclusive!
    ULONG numberOfLocalNodes = to - from + 1;
    ULONG numberOfLocalEdges = readNumberOfEdgesInRange(in, n, from, to);

    std::cout << "peID=" << rank << ": from=" << from << ", to=" << to << ", numberOfLocalNodes="
              << numberOfLocalNodes << ", numberOfLocalEdges=" << numberOfLocalEdges << std::endl;

    in.close();

    // to construct the vertex range array, send 'from' to all other PEs
    std::vector<NodeID> nodeRanges(static_cast<std::size_t>(size + 1));
    MPI_Allgather(&from, 1, MPI_UNSIGNED_LONG_LONG, &nodeRanges[0], 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD);
    nodeRanges[size] = n;

    /*
     * the last part, i.e. loading the graph, is the same as in parallel_graph_io::readGraphBinary()
     */
    PEID windowSize = std::min(size, config.binary_io_window_size);
    PEID lowPE = 0;
    PEID highPE = windowSize;

    while (lowPE < size) {
        if (rank >= lowPE && rank < highPE) {
            std::ifstream in(filename, std::ios::binary | std::ios::in);
            in.exceptions(std::ios_base::failbit | std::ios_base::badbit);

            // extract splitters for split graph construction
            std::vector <EdgeID> edgeRanges(static_cast<std::size_t>(size + 1));
            for (PEID pe = 0; pe < size; ++pe) {
                ULONG peFrom = nodeRanges[pe]; // inclusive
                ULONG peTo = nodeRanges[pe + 1]; // exclusive
                ULONG peLocalNodes = peTo - peFrom;

                ULONG peStartPos = (HEADER_SIZE + peFrom) * sizeof(ULONG);
                NodeID peFirstNodeOffset, peLastNodeOffset;

                in.seekg(peStartPos);
                in.read((char *) &peFirstNodeOffset, sizeof(ULONG));
                in.seekg(peStartPos + peLocalNodes * sizeof(ULONG));
                in.read((char *) &peLastNodeOffset, sizeof(ULONG));

                EdgeID peLocalEdges = (peLastNodeOffset - peFirstNodeOffset) / sizeof(ULONG);
                edgeRanges[pe + 1] = edgeRanges[pe] + peLocalEdges;
            }

            // load and construction, just like parallel_graph_io::readGraphBinary()
            ULONG startPos = (HEADER_SIZE + from) * sizeof(ULONG);
            NodeID *vertexOffsets = new NodeID[numberOfLocalNodes + 1];
            in.seekg(startPos);
            in.read((char *) vertexOffsets, static_cast<std::size_t>((numberOfLocalNodes + 1) * sizeof(ULONG)));

            ULONG edgeStartPos = vertexOffsets[0];
            EdgeID numReads = vertexOffsets[numberOfLocalNodes] - vertexOffsets[0];
            EdgeID numEdgesToRead = numReads / sizeof(ULONG);
            EdgeID *edges = new EdgeID[numEdgesToRead];
            in.seekg(edgeStartPos);
            in.read((char *) edges, static_cast<std::streamsize>(numEdgesToRead * sizeof(ULONG)));

            G.start_construction(numberOfLocalNodes, numberOfLocalEdges, n, m);
            G.set_range(from, to);
            G.set_range_array(nodeRanges);
            G.set_edge_range_array(edgeRanges);

            ULONG pos = 0;
            for (NodeID i = 0; i < numberOfLocalNodes; ++i) {
                NodeID node = G.new_node();
                G.setNodeWeight(node, 1);
                G.setNodeLabel(node, from + node);
                G.setSecondPartitionIndex(node, 0);

                NodeID degree = (vertexOffsets[i + 1] - vertexOffsets[i]) / sizeof(ULONG);
                std::sort(edges + pos, edges + pos + degree);
                for (ULONG j = 0; j < degree; ++j, ++pos) {
                    NodeID target = edges[pos];
                    EdgeID edge = G.new_edge(node, target);
                    G.setEdgeWeight(edge, 1);
                }
            }

            G.finish_construction();
            in.close();

            delete[] edges;
            delete[] vertexOffsets;
        }

        lowPE += windowSize;
        highPE += windowSize;
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

/**
 * Calculates the first node that should be on the given PE such that edges are roughly balanced across all PEs.
 * To achieve this, we use binary search to find the first node such that enough edges are incident to nodes before
 * that one.
 */
static ULONG calculateFromNode(std::ifstream &in, ULONG numberOfNodes, ULONG numberOfEdges, int rank, int size) {
    if (rank == 0) {
        return 0;
    }
    if (rank == size) {
        return numberOfNodes;
    }

    // calculate the number of edges that should come before the first edge on this PE
    ULONG chunk = numberOfEdges / size;
    ULONG remainder = numberOfEdges % size;
    ULONG target = rank * chunk + std::min(static_cast<ULONG>(rank), remainder);

    // find the first node that is incident to an edge greater than target
    // this should be the first edge on this PE
    // a.first = node id
    // a.second = first edge id
    std::pair<ULONG, ULONG> a{0, 0};
    std::pair<ULONG, ULONG> b{numberOfNodes - 1, numberOfEdges - 1};

    while (b.first - a.first > 1) {
        std::pair <ULONG, ULONG> mid;
        mid.first = (a.first + b.first) / 2;
        mid.second = readFirstEdge(in, numberOfNodes, mid.first);

        if (mid.second < target) {
            a = mid;
        } else {
            b = mid;
        }

        assert(b.first >= a.first);
    }

    assert(a.second <= target && target <= b.second);
    assert(b.first < numberOfNodes);
    return b.first;
}

/**
 * Same as calculateFromNode(), but calculates that last node that should be on the given PE.
 */
static ULONG calculateToNode(std::ifstream &in, ULONG numberOfNodes, ULONG numberOfEdges, int rank, int size) {
    return calculateFromNode(in, numberOfNodes, numberOfEdges, rank + 1, size) - 1;
}

/**
 * Determines the number of edges that are incident to nodes in [from, to].
 */
static ULONG readNumberOfEdgesInRange(std::ifstream &in, ULONG numberOfNodes, ULONG from, ULONG to) {
    if (from == to + 1) {
        return 0;
    }
    ULONG firstEdge = readFirstEdge(in, numberOfNodes, from);
    ULONG firstInvalidEdge = readFirstInvalidEdge(in, numberOfNodes, to);
    return firstInvalidEdge - firstEdge;

}

/**
 * Reads G.get_first_edge(NodeID) from a binary graph file.
 */
static ULONG readFirstEdge(std::ifstream &in, ULONG numberOfNodes, ULONG node) {
    assert(node <= numberOfNodes);

    ULONG pos = (HEADER_SIZE + node) * sizeof(ULONG);
    in.seekg(pos);

    ULONG entry = 0;
    in.read((char *) (&entry), sizeof(ULONG));
    return adjacencyListOffsetToEdgeID(numberOfNodes, entry);
}

/**
 * Translates a vertex offset read from a binary graph file to an EdgeID.
 */
static ULONG adjacencyListOffsetToEdgeID(ULONG numberOfNodes, ULONG offset) {
    return (offset / sizeof(ULONG)) - HEADER_SIZE - (numberOfNodes + 1);
}

/**
 * Reads G.get_first_invalid_edge(NodeID) from a binary graph file.
 */
static ULONG readFirstInvalidEdge(std::ifstream &in, ULONG numberOfNodes, ULONG node) {
    return readFirstEdge(in, numberOfNodes, node + 1);
}
