/******************************************************************************
 * timer.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#ifndef TIMER_9KPDEP
#define TIMER_9KPDEP

#include <sys/time.h>
#include <sys/resource.h> 
#include <unistd.h> 

class timer {
        public:
                timer() {
                        m_start = timestamp(); 
                } 

                void restart() { 
                        m_start = timestamp(); 
                } 

                double elapsed() { 
                        return timestamp()-m_start;
                }

        private:

                /** Returns a timestamp ('now') in seconds (incl. a fractional part). */
                inline double timestamp() {
                        struct timeval tp;
                        gettimeofday(&tp, NULL);
                        return double(tp.tv_sec) + tp.tv_usec / 1000000.;
                }

                double m_start;
}; 

#endif /* end of include guard: TIMER_9KPDEP */
