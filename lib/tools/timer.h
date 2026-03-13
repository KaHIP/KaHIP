/******************************************************************************
 * timer.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef TIMER_9KPDEP
#define TIMER_9KPDEP

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <stdint.h>
#else
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

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

#ifdef _WIN32
                inline int gettimeofday(struct timeval* tp, struct timezone* /*tzp*/) {
                        static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);
                        SYSTEMTIME  system_time;
                        FILETIME    file_time;
                        uint64_t    time;
                        GetSystemTime(&system_time);
                        SystemTimeToFileTime(&system_time, &file_time);
                        time = ((uint64_t)file_time.dwLowDateTime);
                        time += ((uint64_t)file_time.dwHighDateTime) << 32;
                        tp->tv_sec = (long)((time - EPOCH) / 10000000L);
                        tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
                        return 0;
                }
#endif

                /** Returns a timestamp ('now') in seconds (incl. a fractional part). */
                inline double timestamp() {
                        struct timeval tp;
                        gettimeofday(&tp, NULL);
                        return double(tp.tv_sec) + tp.tv_usec / 1000000.;
                }

                double m_start;
}; 

#endif /* end of include guard: TIMER_9KPDEP */
