/******************************************************************************
 * aae_timer.h
 *
 * Cheap per-component wall-clock timers. Top-level phases are exclusive: a
 * phase scope only counts when no other phase scope is open, so work nested
 * inside initial partitioning (which recursively coarsens and refines) is
 * attributed to "initial". Detail scopes count only outside initial
 * partitioning. print() writes "aae_time <name> <seconds>" lines to stdout.
 *****************************************************************************/

#ifndef AAE_TIMER_H
#define AAE_TIMER_H

#include <chrono>
#include <cstdio>

namespace aae {

enum Phase { IO = 0, COARSEN, INITIAL, REFINE, POST,
             D_MATCH, D_CONTRACT, D_KWAY, D_QGRAPH, D_LPREF, D_FLOW, D_MULTITRY, NPHASES };

struct timer_state {
        double acc[NPHASES];
        int depth;
        int ip_depth;
};

inline timer_state & state() {
        static timer_state s = {{0}, 0, 0};
        return s;
}

inline double now() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

class phase_scope {
public:
        explicit phase_scope(Phase p) : m_p(p) {
                timer_state & s = state();
                m_on = s.depth == 0;
                ++s.depth;
                if(p == INITIAL) ++s.ip_depth;
                if(m_on) m_t = now();
        }
        ~phase_scope() {
                timer_state & s = state();
                --s.depth;
                if(m_p == INITIAL) --s.ip_depth;
                if(m_on) s.acc[m_p] += now() - m_t;
        }
private:
        Phase m_p;
        bool m_on;
        double m_t;
};

class detail_scope {
public:
        explicit detail_scope(Phase p) : m_p(p) {
                m_on = state().ip_depth == 0;
                if(m_on) m_t = now();
        }
        ~detail_scope() {
                if(m_on) state().acc[m_p] += now() - m_t;
        }
private:
        Phase m_p;
        bool m_on;
        double m_t;
};

inline void print() {
        static const char * names[NPHASES] = {"io", "coarsen", "initial", "refine", "post",
                                              "d_match", "d_contract", "d_kway", "d_qgraph", "d_lpref",
                                              "d_flow", "d_multitry"};
        timer_state & s = state();
        for(int i = 0; i < NPHASES; i++) {
                std::printf("aae_time %s %.6f\n", names[i], s.acc[i]);
        }
        std::fflush(stdout);
}

}

#endif
