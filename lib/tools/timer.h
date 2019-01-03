/******************************************************************************
 * timer.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
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
