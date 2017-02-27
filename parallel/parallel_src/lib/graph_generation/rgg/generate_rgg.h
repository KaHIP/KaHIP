// File:   generate_rgg.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>

#ifndef GENERATE_RGG_H
#define GENERATE_RGG_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <tr1/random>
#include <tr1/tuple>
#include <tr1/random>
#include <vector>

#include "graph_types.h"
#include "util.h"

namespace distributed_graph {

// Forward-Declarations.
//
class StaticDistributedGraph;


// Interface.
//

// The kind of random geometric graph to generate.
enum RggKind {
  RGG_INVALID,
  RGG_UNIT_SQUARE,    // unit square graph
  RGG_UNIT_DISK       // unit disk graph
};


// Configuration for the RGG generator.
struct RandomGeometricGraphGeneratorConfig {
  // The kind of graph to create.
  RggKind kind;

  // The number of points to create.
  long n;

  // The distance below which to connect two points.
  double distance;

  RandomGeometricGraphGeneratorConfig() :
    kind(RGG_INVALID), n(0), distance(0.0)
  {}

  RandomGeometricGraphGeneratorConfig(RggKind arg_kind, long arg_n,
                                      double arg_distance) :
    kind(arg_kind), n(arg_n), distance(arg_distance)
  {}
};


// This class implements a generator for random geometric graph.  It
// generates graphs according to a RandomGeometricGraphGeneratorConfig.
//
// Internally, the generator works with integers as coordinates and then scales
// them down to [0, 1] x [0, 1], the unit disk and so on.
class RandomGeometricGraphGenerator {
public:
  // TODO(manuel): Rename to Quadruple *consistently*
  // Representation of points with identifiers.
  typedef std::tr1::tuple<long, long, long, long> Triple;

  RandomGeometricGraphGenerator() :
    comm_(NULL), rank_(0), process_count_(0), sqrt_process_count_(0),
    int_distance_(0)
  {}

  // Initialize the generator with a configuration.
  //
  // Args:
  //
  //  comm    MPI communicator to use.
  //  config  Generator configuration to use.
  void Init(const MPI::Intracomm *comm,
            const RandomGeometricGraphGeneratorConfig &config);

  // Run the generator with the given random number generator, writing the
  // output to graph.
  //
  // Args:
  //
  //  mt      Mersenne Twister to use for random number generation.
  //  graph   Distributed Graph to create.
  void Run(std::tr1::mt19937 *mt,
           StaticDistributedGraph **graph);

private:

  // Generate points locally.
  //
  // points will be cleared.
  //
  // Args:
  //
  //  n_local Number of local points to create.
  //  mt      Mersenne Twister to use for random number generation.
  //  points  Vector of points-with-id triple.
  void GeneratePointsLocally(const long n_local,
                             std::tr1::mt19937 *mt,
                             std::vector<Triple> *points);

  // Redistribute points to the owners.
  //
  // Args:
  //
  //  points            Vector with points-with-id to redistribute.
  void RedistributePointsToOwners(std::vector<Triple> *points);

  // Sort points, assign ids and build splitters.
  //
  // Args:
  //
  //  vertex_splitters  Vector to create point id splitters in
  //  points            The points to work on.
  //  coordinates       Coordinates to write out to.
  void SortAndAssignIds(std::vector<long> *vertex_splitters,
                        std::vector<Triple> *points,
                        VertexCoordinate **coordinates);

  // Exchange ghost points.  After calling this method, each process has all
  // points it needs to know.
  //
  // Args:
  //
  //  points  Vector with the point set to augment.
  void ExchangeGhostPoints(std::vector<Triple> *points);

  // Build edges from point vector.  Note that points must include the ghost
  // points.
  void BuildEdges(const std::vector<long> &vertex_splitters,
                  std::vector<Triple> *points,
                  std::vector<Edge> *edges);

  // Convert process coordinates to rank.
  void process_xy_to_process_rank(long x, long y, long *rank)
  {
    ASSERT_GEQ( x, 0 );
    ASSERT_LT( x, sqrt_process_count_ );
    ASSERT_GEQ( y, 0 );
    ASSERT_LT( y, sqrt_process_count_ );
    *rank = y * sqrt_process_count_ + x;
  }

  // Convert process rank to process coordinates.
  void process_rank_to_process_xy(long rank, long *x, long *y)
  {
    ASSERT_GEQ( rank, 0 );
    ASSERT_LT( rank, process_count_ );
    *x = rank_ % sqrt_process_count_;
    *y = rank_ / sqrt_process_count_;
  }

  // Convert point coordinates to coordinates of owning point.
  void point_xy_to_process_xy(long x, long y, long *px, long *py)
  {
    ASSERT_GEQ( x, 0 );
    ASSERT_GEQ( y, 0 );
    *px = LocateInSplitters(coordinate_splitters_.begin(), coordinate_splitters_.end(), x);
    *py = LocateInSplitters(coordinate_splitters_.begin(), coordinate_splitters_.end(), y);
    ASSERT_GEQ( *px, 0 );
    ASSERT_LT( *px, sqrt_process_count_ );
    ASSERT_GEQ( *py, 0 );
    ASSERT_LT( *py, sqrt_process_count_ );
  }

  // Convert point coordinates to rank of owning process.
  void point_xy_to_process_rank(long x, long y, long *rank)
  {
    ASSERT_GEQ( x, 0 );
    ASSERT_GEQ( y, 0 );
    long px, py;
    point_xy_to_process_xy(x, y, &px, &py);
    *rank = py * sqrt_process_count_ + px;
  }

  // Returns true iff (x, y) is a valid process coordinate.
  bool valid_process_coordinate(long x, long y)
  {
    if (x < 0) return false;
    if (x >= sqrt_process_count_) return false;
    if (y < 0) return false;
    if (y >= sqrt_process_count_) return false;
    return true;
  }

  bool distance_leq(const Triple &a, const Triple &b, long d);

  // MPI Intracommunicator to use for communication.
  const MPI::Intracomm *comm_;

  // The generator's configuration.
  RandomGeometricGraphGeneratorConfig config_;

  std::vector<long > vertex_splitters_;

  // Rank of this process.
  long rank_;

  // Process count.
  long process_count_;

  // Square root of process count.
  long sqrt_process_count_;

  // Splitters for the coordinates.
  std::vector<long> coordinate_splitters_;

  // Scaled-up integer value of the distance to connect points if closer.
  long int_distance_;

  // largest value for coordinate entries.
  long max_;

  DISALLOW_COPY_AND_ASSIGN(RandomGeometricGraphGenerator);
};

}   // namespace distributed_graph

#endif  // ifndef GENERATE_RGG_H


