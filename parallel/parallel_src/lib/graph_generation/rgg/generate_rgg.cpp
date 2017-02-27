// File:   generate_rgg.cpp
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>

#include "generate_rgg.h"  // This file's header.

#include <tr1/unordered_map>
#include <limits>

//#include "distributed_coarsening.h"
#include "protocol/all_to_all.h"
#include "static_distributed_graph.h"
#include "static_local_graph.h"
#include "static_local_graph_properties.h"
#include "timer.h"

namespace distributed_graph {

// Private Interface
//

// Rasters a vector of points (with ids) and allows the fast access of points
// withing a given raster.
class PointRaster {
public:

  typedef std::tr1::tuple<long, long, long, long> PointWithId;

  PointRaster( const MPI::Intracomm *comm) { 
        comm_ = comm;
  }

  // Initialize raster with a vector of points.
  //
  // The points in the vector will be rearranged (sorted by raster).
  //
  // Args:
  //
  //  idx_begin             Identifier of first point owned by this process.
  //  idx_end               Identifier of first point owned by next process.
  //  rank
  //  process_x
  //  process_y
  //  coordinate_splitters  Splitter values for the coordinates.
  //  points                Points to allow access to.
  void Init(long idx_begin, long idx_end, long rank, long process_x, long process_y,
            long radius,
            const std::vector<long> &coordinate_splitters,
            std::vector<PointWithId> *points)
  {
    idx_begin_ = idx_begin;
    idx_end_ = idx_end;
    rank_ = rank;
    process_x_ = process_x;
    process_y_ = process_y;
    radius_ = radius;
    coordinate_splitters_ = &coordinate_splitters;
    points_ = points;

    BuildRaster();
  }

  // Returns the point with the given index.
  const PointWithId &point(long idx) const
  {
    return (*points_)[idx];
  }

  // Returns count of rasters.
  long raster_count() const
  {
    long x_min = (*coordinate_splitters_)[process_x_];
    long x_max = (*coordinate_splitters_)[process_x_ + 1];
    long raster_x = ceil(1.0 * (x_max - x_min) / radius_) + 3;
    long y_min = (*coordinate_splitters_)[process_y_];
    long y_max = (*coordinate_splitters_)[process_y_ + 1];
    long raster_y = ceil(1.0 * (y_max - y_min) / radius_) + 3;
    long raster_count = raster_x * raster_y;
    return raster_count;
  }

  // Returns index of the first point in the raster with index i.
  long raster_begin(long i) const
  {
    return offsets_[i];
  }

  // Returns index of raster that has a location delta of (dx, dy) from the
  // raster with index i.
  //
  // Returns -1 if such a raster does not exist.
  long raster_with_offset(long i, long dx, long dy) const
  {
    long x_min = (*coordinate_splitters_)[process_x_];
    long x_max = (*coordinate_splitters_)[process_x_ + 1];
    long raster_x = ceil(1.0 * (x_max - x_min) / radius_) + 3;
    long y_min = (*coordinate_splitters_)[process_y_];
    long y_max = (*coordinate_splitters_)[process_y_ + 1];
    long raster_y = ceil(1.0 * (y_max - y_min) / radius_) + 3;
    long x = i % raster_x;
    long y = i / raster_x;
    x += dx;
    y += dy;
    if (x < 0 or x >= raster_x) return -1;
    if (y < 0 or y >= raster_y) return -1;
    return y * raster_x + x;
  }

  // Returns true if the given point is a ghost point.
  bool is_ghost_point(int i) const
  {
    long idx = std::tr1::get<2>((*points_)[i]);
    return (idx < idx_begin_) or (idx >= idx_end_);
  }

private:

  void point_xy_to_raster_xy(long x, long y, long *rx, long *ry) const
  {
    // We create a raster entry left of the origin and two right of the
    // right border.  Create one above the origin and two below the lower
    // border.
    long x_min = (*coordinate_splitters_)[process_x_];
    long y_min = (*coordinate_splitters_)[process_y_];
    long x_off = x_min - radius_;
    long y_off = y_min - radius_;
    *rx = (x - x_off) / radius_;
    *ry = (y - y_off) / radius_;
  }

  void raster_xy_to_raster_id(long x, long y, long *i) const
  {
    long x_min = (*coordinate_splitters_)[process_x_];
    long x_max = (*coordinate_splitters_)[process_x_ + 1];
    long raster_x = ceil(1.0 * (x_max - x_min) / radius_) + 3;
    *i = x + y * raster_x;
    //std::cout <<  *i << " " <<  x <<  " " <<  y << std::endl;
  }

  void point_xy_to_raster_id(long x, long y, long *i) const
  {
    long rx, ry;
    point_xy_to_raster_xy(x, y, &rx, &ry);
    raster_xy_to_raster_id(rx, ry, i);
  }

  struct Comparator
  {
    inline bool operator()(const PointWithId &left, const PointWithId &right) const
    { return std::tr1::get<3>(left) < std::tr1::get<3>(right); }
  };

  void BuildRaster()
  {
#ifndef NDEBUG
    for (long i = 0, iend = points_->size(); i < iend; ++i) {
      long x = std::tr1::get<0>((*points_)[i]);
      long y = std::tr1::get<1>((*points_)[i]);
    }
#endif

    long x_min = (*coordinate_splitters_)[process_x_];
    long x_max = (*coordinate_splitters_)[process_x_ + 1];
    long raster_x = ceil(1.0 * (x_max - x_min) / radius_) + 3;
    long y_min = (*coordinate_splitters_)[process_y_];
    long y_max = (*coordinate_splitters_)[process_y_ + 1];
    long raster_y = ceil(1.0 * (y_max - y_min) / radius_) + 3;
    long raster_count = raster_x * raster_y;

    if (rank_ == 0) {
      printf("  Sorting points by raster element\n");
      printf("    Building map\n");
    }
    ClockTimer timer; timer.Start();
    // Sort the points by the raster entry they fall into.
    //printf("raster x = %d, raster y = %d\n", raster_x, raster_y);
    //printf("radius == %d\n", radius_);
    for (long i = 0, iend = points_->size(); i < iend; ++i) {
      long x = std::tr1::get<0>((*points_)[i]);
      long y = std::tr1::get<1>((*points_)[i]);
      long rid;
      point_xy_to_raster_id(x, y, &rid);
      std::tr1::get<3>((*points_)[i]) = rid;
      //std::cout <<  "x " <<  x <<  " y " <<  y  << std::endl;
    }
    comm_->Barrier();
    if (rank_ == 0)
      printf("    Took %.2f s\n", timer.elapsed());
    timer.Restart();
    if (rank_ == 0)
      printf("  Actually sorting\n");
    Comparator cmp;
    std::sort(points_->begin(), points_->end(), cmp);
    /*
    for (int i = 0, iend = points_->size(); i < iend; ++i) {
      int x = std::tr1::get<0>((*points_)[i]);
      int y = std::tr1::get<1>((*points_)[i]);
      int rid;
      point_xy_to_raster_id(x, y, &rid);
      printf("%d %d\n", rid, entry_map[std::tr1::get<2>((*points_)[i])]);
    }
    */
    comm_->Barrier();
    if (rank_ == 0)
      printf("    Took %.2f s\n", timer.elapsed());

    if (rank_ == 0)
      printf("  Building offsets\n");
    timer.Restart();
    // Build offsets_ array.
    offsets_.clear();
    offsets_.reserve(raster_count + 1);
    long old;
    point_xy_to_raster_id(std::tr1::get<0>((*points_)[0]),
                          std::tr1::get<1>((*points_)[0]), &old);
    //std::cout <<  "old " << old  << std::endl;
    for (long i = 0; i <= old; ++i) {
      //printf("i = %d, push 0\n", i);
      offsets_.push_back(0);
    }
    //std::cout <<  points_->size() << std::endl;
    for (long i = 1, iend = points_->size(); i < iend; ++i) {
      long idx;
      point_xy_to_raster_id(std::tr1::get<0>((*points_)[i]), 
                            std::tr1::get<1>((*points_)[i]), &idx);
      //printf("idx = %d\n", idx);
      //std::cout <<  rank_ << " " <<idx  << std::endl;
      while (offsets_.size() <= (unsigned long)idx) {
      //for (int j = old; j < idx; ++j) {
        //printf("j = %lu, push i = %d\n", offsets_.size(), i);
        offsets_.push_back(i);
        old = idx;
      }
    }
    //std::cout <<  "B"  << std::endl;
    while (offsets_.size() <= (unsigned long)raster_count) {
      //printf("i = %lu, push points_->size() = %lu\n", offsets_.size(), points_->size());
      offsets_.push_back(points_->size());
    }

    //std::cout <<  "C"  << std::endl;
    comm_->Barrier();
    if (rank_ == 0)
      printf("    Took %.2f s\n", timer.elapsed());
  }

  // See documentation of Init().
  long idx_begin_;
  long idx_end_;
  long rank_;
  long process_x_;
  long process_y_;
  long radius_;
  const std::vector<long> *coordinate_splitters_;

  // Offsets of the rasters.
  std::vector<long> offsets_;

  // The rastered vector of points.
  std::vector<PointWithId> *points_;
  const MPI::Intracomm *comm_;

  DISALLOW_COPY_AND_ASSIGN(PointRaster);
};


// Implementation
//

void RandomGeometricGraphGenerator::Init(
    const MPI::Intracomm *comm,
    const RandomGeometricGraphGeneratorConfig &config)
{
  comm_ = comm;
  config_ = config;

  rank_ = comm_->Get_rank();
  process_count_ = comm_->Get_size();
  sqrt_process_count_ = sqrt(process_count_);
  // The number of processes must be a square.
  if (sqrt_process_count_ * sqrt_process_count_ != process_count_)
    AbortProgram("Communicator size is not a square number!\n");

  max_ = std::numeric_limits<int>::max() / 2;
  int_distance_ = round(config_.distance * max_);

  // Compute coordinate splitters.
  VertexId *splitters;
  ComputeSplitters(max_, sqrt_process_count_, &splitters);
  std::vector<long> splitters2(splitters, splitters + sqrt_process_count_ + 1);
  coordinate_splitters_.swap(splitters2);
  delete[] splitters;
}


void RandomGeometricGraphGenerator::Run(
    std::tr1::mt19937 *mt,
    StaticDistributedGraph **graph)
{
  // MPI information.
  long process_count = comm_->Get_size();
  long rank = comm_->Get_rank();
  // The set o fpoints {(x, y, id)}.
  std::vector<Triple> points;
  // Build splitters to get number of wanted local points.
  VertexId *vertex_splitters;
  ComputeSplitters(config_.n, process_count, &vertex_splitters);
  long n_local = vertex_splitters[rank + 1] - vertex_splitters[rank];
  delete[] vertex_splitters;

  // Generate the points, distribute them to their real owners, kick out
  // duplicates and fix the local point sets such that they have the number
  // of vertices they should.
  if (rank == 0) printf("Configuration: n = %lu, radius = %f, max = %lu, int radius = %lu\n", config_.n, config_.distance, max_, int_distance_);
  comm_->Barrier();
  if (rank == 0) printf("Generating points locally.\n");
  GeneratePointsLocally(n_local, mt, &points);
  comm_->Barrier();
  if (rank == 0) printf("Redistributing points to owners.\n");
  RedistributePointsToOwners(&points);
  VertexCoordinate *coordinates;
  comm_->Barrier();
  if (rank == 0) printf("Sort, assign ids, write coordinates.\n");
  SortAndAssignIds(&vertex_splitters_, &points, &coordinates);
  comm_->Barrier();
  if (rank == 0) printf("Exchange ghost points.\n");
  ExchangeGhostPoints(&points);
  std::vector<Edge> edges;
  comm_->Barrier();
  if (rank == 0) printf("Building edges.\n");
  BuildEdges(vertex_splitters_, &points, &edges);
  if (rank == 0) printf("Sorting edges.\n");
  ClockTimer timer;  timer.Start();
  //std::cout <<  edges.size()  << std::endl;
  std::sort(edges.begin(), edges.end());
  comm_->Barrier();
  if (rank_ == 0)
    printf("  Took %.2f\n", timer.elapsed());
  if (rank == 0) printf("Building graph.\n");


  // Create the graph from the edges and check it.
  long n = vertex_splitters_[process_count];
  long local_n = vertex_splitters_[rank + 1] - vertex_splitters_[rank];
  long local_m = edges.size();
  long vertex_offset = vertex_splitters_[rank];
  EdgeId *edge_splitters = new EdgeId[process_count_ + 1];
  edge_splitters[0] = 0;
  comm_->Allgather(&local_m, 1, MPI::LONG,
                   edge_splitters + 1, 1, MPI::LONG);
  std::partial_sum(edge_splitters, edge_splitters + process_count_ + 1,
                   edge_splitters);
  long m = edge_splitters[process_count];
  //std::cout <<  m  << std::endl;
  long edge_offset = edge_splitters[rank_];
  StaticLocalGraph *local_graph = new StaticLocalGraph();
  local_graph->Init(vertex_offset, edge_offset, local_n, local_m,
                    &edges[0]);
  StaticLocalGraphProperties *properties = new StaticLocalGraphProperties();
  properties->InitEmpty(vertex_offset, edge_offset, local_n, local_m);
  properties->InitCoordinates(coordinates);
  *graph = new StaticDistributedGraph();
  vertex_splitters = new VertexId[process_count_ + 1];
  std::copy(vertex_splitters_.begin(), vertex_splitters_.end(),
            vertex_splitters);
  (*graph)->Init(comm_, n, m, vertex_splitters, edge_splitters,
                 local_graph, properties);

}


void RandomGeometricGraphGenerator::GeneratePointsLocally(
    const long n_local, std::tr1::mt19937 *mt, std::vector<Triple> *points)
{
  // Prepare points and reserve memory for more than 10% than expected.
  points->clear();
  if (config_.kind == RGG_UNIT_SQUARE) {
    points->reserve(1.1 * n_local);
  } else if (config_.kind == RGG_UNIT_DISK) {
    AbortProgram("RGG_UNIT_DISK unsupported as of now!\n");
    points->reserve(1.1 * 4.0 / 3.141 * n_local);
  } else {
    AbortProgram("Invalid random geometric graph kind.\n");
  }

  // Generate points.
  std::tr1::uniform_int<int> dist(0, max_ - 1);
  for (long i = 0; i < n_local; ++i) {
    long x = dist(*mt);
    long y = dist(*mt);
    //std::cout << "a " <<   x   << " "  <<  y   << std::endl;
    points->push_back(std::tr1::make_tuple(x, y, 0, 0));
  }
}


void RandomGeometricGraphGenerator::RedistributePointsToOwners(
    std::vector<Triple> *points)
{
  // Send data.
  AllToAllHelper<long> all2all;
  all2all.Init(comm_);
  for (long i = 0, iend = points->size(); i < iend; ++i) {
    long pid;
    point_xy_to_process_rank(std::tr1::get<0>((*points)[i]),
                             std::tr1::get<1>((*points)[i]), &pid);
    all2all.set_send_count(pid, all2all.send_count(pid) + 2);
  }
  all2all.ExchangeSendCounts();
  for (long i = 0, iend = points->size(); i < iend; ++i) {
    long pid; 
    point_xy_to_process_rank(std::tr1::get<0>((*points)[i]),
                             std::tr1::get<1>((*points)[i]), &pid);
    //std::cout <<  pid  << std::endl;
    all2all.PushForProcess(std::tr1::get<0>((*points)[i]), pid);
    all2all.PushForProcess(std::tr1::get<1>((*points)[i]), pid);
  }
  all2all.ExchangeData();

  // Copy received data into points.
  points->clear();
  for (long i = 0, iend = all2all.receive_displacement(process_count_);
       i < iend; i += 2) {
    points->push_back(std::tr1::make_tuple(
            all2all.received_element(i), all2all.received_element(i + 1), 0, 0));
    //std::cout <<  std::tr1::get<0>((*points)[i])  << std::endl;
  }
}


void RandomGeometricGraphGenerator::SortAndAssignIds(
    std::vector<long> *vertex_splitters,
    std::vector<Triple> *points,
    VertexCoordinate **coordinates)
{
  if (config_.kind != RGG_UNIT_SQUARE) {
    AbortProgram("Only RGG_UNIT_SQURE is supported at the moment!\n");
  }

  // Build vertex splitters.
  vertex_splitters->resize(process_count_ + 1);
  (*vertex_splitters)[0] = 0;
  long local_n = points->size();
  comm_->Allgather(&local_n, 1, MPI::LONG,
                   &(*vertex_splitters)[1], 1, MPI::LONG);
  std::partial_sum(vertex_splitters->begin(), vertex_splitters->end(),
                   vertex_splitters->begin());

  // Sort vertices lexicographically.  This sucks in terms of locality.
  // However, I could not find a spatial sorting routine except CGAL's.
  // TODO(manuel): Use spatial sorting here, rest should work as always.
  std::sort(points->begin(), points->end());

  // Assign identifiers to points.
 *coordinates = new VertexCoordinate[3 * local_n];
  long offset = (*vertex_splitters)[rank_];
  //if (rank_ == 0)
  //  printf("sort and assign ids\n");
  for (long i = 0, iend = points->size(); i < iend; ++i) {
    std::tr1::get<2>((*points)[i]) = offset + i;
    (*coordinates)[i * 3 + 0] = 1.0 * std::tr1::get<0>((*points)[i]) / max_;
    (*coordinates)[i * 3 + 1] = 1.0 * std::tr1::get<1>((*points)[i]) / max_;
    (*coordinates)[i * 3 + 2] = 0.0;
    //printf("  rank %d -- %2d (%f, %f) sai\n", rank_, std::tr1::get<2>((*points)[i]),
    //       (*coordinates)[i * 3], (*coordinates)[i * 3 + 1]);
  }
  comm_->Barrier();
}


void RandomGeometricGraphGenerator::ExchangeGhostPoints(
    std::vector<Triple> *points)
{
  // NOTE(manuel):  I decided not to abstract the case distinction away
  // for now because it is clearer this way.

  // This process' coordinates.
  long px, py;
  process_rank_to_process_xy(rank_, &px, &py);
  // Exchange the data.
  long deltas[4][2] = {{0, -1}, {0, 1}, {-1, 0}, {1, 0}};
  for (long dir = 0; dir < 4; ++dir) {
    // Compute process coordinates.
    long dx = deltas[dir][0];
    long dy = deltas[dir][1];
    // Neighbour x and y coordinates.
    long nx = px + dx;
    long ny = py + dy;
    long neighbour_rank = 0;
    if (valid_process_coordinate(nx, ny)) {
      process_xy_to_process_rank(nx, ny, &neighbour_rank);
    }
    AllToAllHelper<long> all2all;
    all2all.Init(comm_);
    if (valid_process_coordinate(nx, ny)) {
      if (dx == 0) {
        if (dy == -1) {
          long bottom_of_upper = coordinate_splitters_[py] - 1;
          for (long i = 0, iend = points->size(); i < iend; ++i) {
            if (abs(std::tr1::get<1>((*points)[i]) - bottom_of_upper) <= int_distance_) {
              all2all.set_send_count(neighbour_rank, all2all.send_count(neighbour_rank) + 3);
            }
          }
        } else {  // dy == 1
          long top_of_lower = coordinate_splitters_[ny];
          for (long i = 0, iend = points->size(); i < iend; ++i) {
            if (abs(std::tr1::get<1>((*points)[i]) - top_of_lower) <= int_distance_) {
              all2all.set_send_count(neighbour_rank, all2all.send_count(neighbour_rank) + 3);
            }
          }
        }
      } else {
        if (dx == -1) {
          long right_of_left_hand = coordinate_splitters_[px] - 1;
          for (long i = 0, iend = points->size(); i < iend; ++i) {
            if (abs(std::tr1::get<0>((*points)[i]) - right_of_left_hand) <= int_distance_) {
              all2all.set_send_count(neighbour_rank, all2all.send_count(neighbour_rank) + 3);
            }
          }
        } else {  // dx == 1
          long left_of_right_hand = coordinate_splitters_[nx];
          for (long i = 0, iend = points->size(); i < iend; ++i) {
            if (abs(std::tr1::get<0>((*points)[i]) - left_of_right_hand) <= int_distance_) {
              all2all.set_send_count(neighbour_rank, all2all.send_count(neighbour_rank) + 3);
            }
          }
        }
      }
    }
    all2all.ExchangeSendCounts();
    if (valid_process_coordinate(nx, ny)) {
      if (dx == 0) {
        if (dy == -1) {
          long bottom_of_upper = coordinate_splitters_[py] - 1;
          for (long i = 0, iend = points->size(); i < iend; ++i) {
            if (abs(std::tr1::get<1>((*points)[i]) - bottom_of_upper) <= int_distance_) {
              all2all.PushForProcess(std::tr1::get<0>((*points)[i]), neighbour_rank);
              all2all.PushForProcess(std::tr1::get<1>((*points)[i]), neighbour_rank);
              all2all.PushForProcess(std::tr1::get<2>((*points)[i]), neighbour_rank);
            }
          }
        } else {  // dy == 1
          long top_of_lower = coordinate_splitters_[ny];
          for (long i = 0, iend = points->size(); i < iend; ++i) {
            if (abs(std::tr1::get<1>((*points)[i]) - top_of_lower) <= int_distance_) {
              all2all.PushForProcess(std::tr1::get<0>((*points)[i]), neighbour_rank);
              all2all.PushForProcess(std::tr1::get<1>((*points)[i]), neighbour_rank);
              all2all.PushForProcess(std::tr1::get<2>((*points)[i]), neighbour_rank);
            }
          }
        }
      } else {
        if (dx == -1) {
          long right_of_left_hand = coordinate_splitters_[px] - 1;
          for (long i = 0, iend = points->size(); i < iend; ++i) {
            if (abs(std::tr1::get<0>((*points)[i]) - right_of_left_hand) <= int_distance_) {
              all2all.PushForProcess(std::tr1::get<0>((*points)[i]), neighbour_rank);
              all2all.PushForProcess(std::tr1::get<1>((*points)[i]), neighbour_rank);
              all2all.PushForProcess(std::tr1::get<2>((*points)[i]), neighbour_rank);
            }
          }
        } else {  // dx == 1
          long left_of_right_hand = coordinate_splitters_[nx];
          for (long i = 0, iend = points->size(); i < iend; ++i) {
            if (abs(std::tr1::get<0>((*points)[i]) - left_of_right_hand) <= int_distance_) {
              all2all.PushForProcess(std::tr1::get<0>((*points)[i]), neighbour_rank);
              all2all.PushForProcess(std::tr1::get<1>((*points)[i]), neighbour_rank);
              all2all.PushForProcess(std::tr1::get<2>((*points)[i]), neighbour_rank);
            }
          }
        }
      }
    }
    all2all.ExchangeData();

    // Append points to points
    for (long i = 0, iend = all2all.receive_displacement(process_count_);
         i < iend; i += 3) {
      long x = all2all.received_element(i);
      long y = all2all.received_element(i + 1);
      long id = all2all.received_element(i + 2);
      points->push_back(std::tr1::make_tuple(x, y, id, 0));
    }
  }

  // Sort points and remove duplicates.
  std::sort(points->begin(), points->end());
  std::vector<Triple>::iterator end = std::unique(points->begin(), points->end());
  points->resize(end - points->begin());


  /*
  printf("%d After receiving ghosts\n", rank_);
  for (int i = 0, iend = points->size(); i < iend; ++i) {
    printf("  rank %d -- i=%2d, id=%2d (%f, %f) rg\n", rank_,
           i, 
           std::tr1::get<2>((*points)[i]),
           1.0 * std::tr1::get<0>((*points)[i]) / max_,
           1.0 * std::tr1::get<1>((*points)[i]) / max_); 
  }
  */
}


bool RandomGeometricGraphGenerator::distance_leq(
    const Triple &a, const Triple &b, long d)
{
  double dx = 1.0 * (std::tr1::get<0>(a) - std::tr1::get<0>(b)) / max_;
  double dy = 1.0 * (std::tr1::get<1>(a) - std::tr1::get<1>(b)) / max_;
  double dd = 1.0 * d / max_;
  //printf("P_%d @ (%d, %d), P_%d @ (%d, %d), d = %d, < = %d\n",
        //std::tr1::get<2>(a), std::tr1::get<0>(a), std::tr1::get<1>(a),
        //std::tr1::get<2>(b), std::tr1::get<0>(b), std::tr1::get<1>(b), d,
        //dx*dx+dy*dy<dd*dd);
  return dx*dx+dy*dy<dd*dd;
}


void RandomGeometricGraphGenerator::BuildEdges(
    const std::vector<long> &vertex_splitters,
    std::vector<Triple> *points,
    std::vector<Edge> *edges)
{
  long process_x, process_y;
  process_rank_to_process_xy(rank_, &process_x, &process_y);
  if (rank_ == 0) printf("Building raster.\n");
  ClockTimer timer; timer.Start();
  PointRaster raster(comm_);
 //printf("Rank %d from %d to %d\n", rank_,
        //vertex_splitters[rank_], vertex_splitters[rank_ + 1]);
  raster.Init(vertex_splitters[rank_], vertex_splitters[rank_ + 1], rank_,
              process_x, process_y, int_distance_, coordinate_splitters_, points);
  comm_->Barrier();
  timer.Stop();
  if (rank_ == 0)
    printf("  Took %.2f s\n", timer.elapsed());
  if (rank_ == 0) printf("  raster count == %lu\n", raster.raster_count());
  long global_n = points->size();
  comm_->Allreduce(MPI::IN_PLACE, &global_n, 1, MPI::LONG, MPI::SUM);
  double expected_average_degree = 3.141 * config_.distance * config_.distance * global_n;
  long reserved = expected_average_degree * points->size();
  edges->reserve(1.1 * reserved);

  //if (rank_ == 0)
  //  printf("splitters %d, %d, %d, %d, %d\n", vertex_splitters_[0], vertex_splitters_[1],
  //         vertex_splitters_[2], vertex_splitters_[3], vertex_splitters_[4]);

  timer.Restart(); timer.Start();
  if (rank_ == 0) {
    printf("Building edges.\n");
    printf("  Expected average degree = %f, n = %lu\n", expected_average_degree, (unsigned long) points->size());
    printf("  Reserving memory for %lu\n", static_cast<long>(1.1 * reserved));
  }
  for (long i = 0, iend = raster.raster_count(); i < iend; ++i) {
    // i is the index of the current raster
    for (long dx = -1; dx <= 1; ++dx) {
      for (long dy = -1; dy <= 1; ++dy) {
        // j is the index of the raster with (dx, dy) off i.
        long j = raster.raster_with_offset(i, dx, dy);
        if (j == -1) continue;
        //if (raster.raster_begin(i) != raster.raster_begin(i + 1)) {
        //  if (raster.raster_begin(j) != raster.raster_begin(j + 1)) {
        //    printf("rank %d -- %d <-> %d\n", rank_, i, j);
        //  }
        //}
        for (long a = raster.raster_begin(i), aend = raster.raster_begin(i + 1); a < aend; ++a) {
          const Triple &p1 = raster.point(a);
          if (raster.is_ghost_point(a)) {
            //printf("rank %d -- ghost point: %d\n", rank_, std::tr1::get<2>(p1));
            continue;
          }
          for (long b = raster.raster_begin(j), bend = raster.raster_begin(j + 1); b < bend; ++b) {
            const Triple &p2 = raster.point(b);
            //std::cout <<  std::tr1::get<0>(p1) << std::endl;
            //std::cout <<  std::tr1::get<1>(p1) << std::endl;
            //std::cout <<  std::tr1::get<0>(p2) << std::endl;
            //std::cout <<  std::tr1::get<1>(p2) << std::endl;
            if (std::tr1::get<2>(p1) == std::tr1::get<2>(p2)) continue;
            if (distance_leq(p1, p2, int_distance_)) {
              //printf("rank %d -- (%d, %d)\n", rank_, std::tr1::get<2>(p1), std::tr1::get<2>(p2));
              edges->push_back(Edge(std::tr1::get<2>(p1), std::tr1::get<2>(p2)));
            }
          }
        }
      }
    }
  }
  comm_->Barrier();
  timer.Stop();
  if (rank_ == 0)
    printf("  Took %.2f s\n", timer.elapsed());
}


};  // namespace distributed_graph

