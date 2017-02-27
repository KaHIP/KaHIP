/******************************************************************************
 * statistics.h
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

#ifndef STATISTICS_H
#define STATISTICS_H

#include <algorithm>
#include <cassert>
#include <cmath>    // log(), exp()
#include <cstdio>
#include <vector>

// Wrapper around a vector with methods to calculate statistic values.
// 
// Note that we allow copy and assignment and use the default ones from the
// compiler.  This class only contains an STL vector so this is fairly safe.
template <typename T>
class Statistics
{
public:

  Statistics()
  {}

  /* Dump the statistics. */
  void Dump(FILE *out) const
  {
    fprintf(out, "{ n = %lu, mean = %f, median = %f, stddev = %f"
        ", min = %f, max = %f}\n",
        values_.size(), static_cast<double>(Mean()), Median(), StdDev(),
        static_cast<double>(Min()), static_cast<double>(Max()));
  }

  /* Add a value to the list of values. 
   */
  void AddValue(const T value)
  {
    values_.push_back(value);
  }

  /* Returns the arithmetic mean of all values.
   */
  double Mean() const
  {
    double sum = 0;
    typedef typename std::vector<T>::const_iterator iterator;
    for (iterator it = values_.begin(), itend = values_.end();
         it != itend; ++it) {
      sum += *it;
    }
    unsigned n = values_.size() > 0 ? values_.size() : 1;
    return sum / n;
  }

  /* Returns the geometric mean of all values.
   */
  double GeometricMean() const
  {
    double sum = 0;
    typedef typename std::vector<T>::const_iterator iterator;
    for (iterator it = values_.begin(), itend = values_.end();
         it != itend; ++it) {
      assert( *it >= 0 );
      sum += log(*it);
    }
    unsigned n = values_.size() > 0 ? values_.size() : 1;
    return exp(sum / n);
  }

  /* Returns the median.  If there is an even number of values then return the
   * arithmetic mean of the two "middle" elements. If there is no value then
   * 0 is returned.
   */
  double Median() const
  {
    if (values_.size() == 0) return 0;
    if (values_.size() == 1) return values_[0];

    std::vector<T> values(values_);

    std::sort(values.begin(), values.end());
    if (values.size() % 2 != 0) {
      return values[values.size() / 2 + 1];
    } else {
      double x1 = static_cast<double>(values[values.size() / 2]) / 2;
      double x2 = static_cast<double>(values[values.size() / 2 + 1]) / 2;
      return x1 + x2;
    }
  }

  /* Return the standard deviation of the value array.  If there are less than
   * two values then -1 is returned.
   */
  double StdDev() const
  {
    if (values_.size() < 2) return -1;

    const double m = Mean();
    double sum = 0;

    typedef typename std::vector<T>::const_iterator iterator;
    for (iterator it = values_.begin(), itend = values_.end();
         it != itend; ++it) {
      sum += (*it - m) * (*it - m);
    }

    return sqrt((1.0 / values_.size()) * sum);
  }

  /* Reset the object. */
  void Clear()
  { values_.clear(); }

  /* Return maximum.  Return 0 if there is no value. */
  T Max() const
  {
    if (values_.size() == 0) return 0;
    T max = values_[0];
    for (int i = 0, iend = values_.size(); i != iend; ++i) {
      if (max < values_[i])
        max = values_[i];
    }
    return max;
  }

  /* Return minimum.  Return 0 if there is no value. */
  T Min() const
  {
    if (values_.size() == 0) return 0;
    T min = values_[0];
    for (int i = 0, iend = values_.size(); i != iend; ++i) {
      if (min > values_[i])
        min = values_[i];
    }
    return min;
  }

private:
  // The values we have accumulated so far.
  std::vector<T> values_;
};


#endif // ifndef STATISTICS_H

