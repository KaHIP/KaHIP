/******************************************************************************
 * timer.h
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


#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <cstdio>

#include "macros_common.h"

#define RUN_TIMED(command, timer) \
  (timer).Restart();\
  command; \
  (timer).Stop(); \
  printf("%s took %f s\n", STR(command), timer.elapsed());


// TODO(manuel): Using a double to accumulate time is not that great. Maybe restrict to 1/x-th of a second and use integers?

// Small timer class, based on <ctime>'s clock().
//
// Examples of equivalent usage:
//
//  clock_timer timer1(true):
//  // do something
//  double t1 = timer1.elapsed();
//
//  clock_timer timer2;
//  // do something
//  timer2.Start();
//  // do something
//  double t2 = timer2.Stop();
class ClockTimer
{
public:
  
  // Star the timer if autostart.
  explicit
  ClockTimer(bool autostart = false) : start_(0), end_(0), running_(0),
    accumulated_time_(0)
  {
    if (autostart) Start();
  }

  // Restart the timer.
  void Restart()
  { Start(); }

  // Start or restart the timer.
  void Start()
  {
    running_ = true;
    start_ = clock();
  }

  // Stop the timer, returning the elapsed time.  If the timer is not running,
  // the last elapsed time will be returned.
  double Stop()
  {
    if (not running_) return elapsed();
    end_ = clock();
    running_ = false;
    double t = elapsed();
    accumulated_time_ += t;
    return t;
  }

  // Return whether the timer is running at the moment.
  bool running()
  { return running_; }

  // Returns the elapsed time in seconds.
  //
  // If the timer is still running then it returns the time until
  // now.  Otherwise, it returns the time between starting and stopping
  // the timer.
  double elapsed() const
  {
    if (running_)
      return static_cast<double>(clock() - start_) / CLOCKS_PER_SEC;
    else
      return static_cast<double>(end_ - start_) / CLOCKS_PER_SEC;
  }

  // Returns the total time elapsed between start/stop pairs so far.  The
  // currently elapsing time is not considered.
  double accumulated() const
  { return accumulated_time_; }

private:
  // Point of time where the timer was started/stopped.
  clock_t start_, end_;

  // Flag whether timer is running.
  bool running_;

  // Accumulator for the total time elapsed so far.
  double accumulated_time_;
};

#endif // ifndef TIMER_H 

