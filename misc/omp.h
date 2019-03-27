#pragma once

//workaround dummy header
template<typename T>
void omp_set_num_threads(T) {}

inline int omp_get_thread_num() {
        return 1;
}

inline int omp_get_max_threads() {
        return 1;
}
