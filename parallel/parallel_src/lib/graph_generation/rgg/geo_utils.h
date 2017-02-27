// File:   geo_utils.h
// Author: Manuel Holtgrewe <holtgrewe@ira.uka.de>
//
// Geometric helper code.
//
// Classes:
//   BoundingBox

#ifndef GEO_UTILS_H
#define GEO_UTILS_H

#include <algorithm>

#include <mpi.h>

#include "instrumentation.h"
#include "macros_common.h"

namespace distributed_graph {

// A helper that simplifies the computation of a bounding box of a number of
// points incrementally.
//
// TODO(manuel): Move out functions with > 8 lines into .cpp file?
class BoundingBox
{
public:
  // Obligatory default constructor.
  BoundingBox() :
    initialized_(false)
  {}

  // Register a new point, possibly expand bounding box.
  void AddPoint(double x, double y, double z)
  {
    if (not initialized_) {
      initialized_ = true;
      x_min_ = x;
      x_max_ = x;
      y_min_ = y;
      y_max_ = y;
      z_min_ = z;
      z_max_ = z;
    } else {
      x_min_ = std::min(x_min_, x);
      x_max_ = std::max(x_max_, x);
      y_min_ = std::min(y_min_, y);
      y_max_ = std::max(y_max_, y);
      z_min_ = std::min(z_min_, z);
      z_max_ = std::max(z_max_, z);
    }
  }

  // Get the current bounds of the point set.
  void GetBounds(double *x_min, double *x_max, double *y_min, double *y_max,
                 double *z_min, double *z_max) const
  {
    *x_min = x_min_;
    *x_max = x_max_;
    *y_min = y_min_;
    *y_max = y_max_;
    *z_min = z_min_;
    *z_max = z_max_;
  }

  // Set bounds of the bounding box.
  void SetBounds(double x_min, double x_max, double y_min, double y_max,
                 double z_min, double z_max)
  {
    x_min_ = x_min;
    x_max_ = x_max;
    y_min_ = y_min;
    y_max_ = y_max;
    z_min_ = z_min;
    z_max_ = z_max;
  }

  // Return length in the given dimension, d = 0 is x, d = 1 is y, d = 2 is z.
  double dimension_length(long d) const
  { 
    if (d == 0)
      return x_max_ - x_min_;
    else if (d == 1)
      return y_max_ - y_min_;
    else if (d == 2)
      return z_max_ - z_min_;
    else
      ASSERT_TRUE( false );
    return 0.0;
  }

  // Returns hte longest dimension.
  long longest_dimension() const
  {
    long result = 0;
    if (dimension_length(1) > dimension_length(0))
      result = 1;
    if (dimension_length(2) > dimension_length(result))
      result= 2;
    return result;
  }

  // Compute the bounding box of all BoundingBox in the given communicator.
  void ComputeGlobalBoundingBox(const MPI::Intracomm *comm)
  {
    ENTER_SECTION(comm, COMPUTE_GLOBAL_BOUNDING_BOX);
    // Take the min entries times -1 so we can use MPI::MAX below.
    double bounds[6] = {-x_min_, x_max_, -y_min_, y_max_, -z_min_, z_max_};
    // Convert the local bounds to global ones using the in place feature.
    comm->Allreduce(MPI::IN_PLACE, bounds, 6, MPI::DOUBLE, MPI::MAX);
    // Again, multiply the {x,y,z}_min values by -1 to work with them again.
    bounds[0] *= -1;
    bounds[2] *= -1;
    bounds[4] *= -1;
    SetBounds(bounds[0], bounds[1], bounds[2], bounds[3],
              bounds[4], bounds[5]);
    LEAVE_SECTION(comm, COMPUTE_GLOBAL_BOUNDING_BOX);
  }

private:
  // true iff at least one point has been added.
  bool initialized_;

  // The current bounds in each direction.
  double x_min_;
  double x_max_;
  double y_min_;
  double y_max_;
  double z_min_;
  double z_max_;

  DISALLOW_ASSIGN(BoundingBox);
};

}  // namespace distributed_graph

#endif // ifndef GEO_UTILS_H

