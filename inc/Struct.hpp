#ifndef STRUCT_H
#define STRUCT_H

#include <vector>

struct Instance {
    int n_jobs;                    // liczba zadań (n)
    int m_machines;                // liczba maszyn (m)
    // times[job][machine] = czas przetwarzania zadania na maszynie
    std::vector<std::vector<int>> times;
};

#endif // STRUCT_H