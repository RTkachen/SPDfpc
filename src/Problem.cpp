#include "Problem.hpp"
Problem::Problem(const std::vector<Instance>& instances)
    : instances_(instances) {}

    int Problem::computeCmax(const std::vector<int>& sequence, int instanceIndex) const {
        const Instance& inst = instances_[instanceIndex];
        int k = static_cast<int>(sequence.size());    // tylko tyle zadań bierzemy pod uwagę
        int m = inst.m_machines;
        std::vector<std::vector<int>> C(k+1, std::vector<int>(m+1, 0));
    
        for (int j = 1; j <= k; ++j) {
            for (int i = 1; i <= m; ++i) {
                int job = sequence[j-1];
                int p   = inst.times[job][i-1];
                C[j][i] = std::max(C[j][i-1], C[j-1][i]) + p;
            }
        }
        return C[k][m];
    }
    
    std::pair<std::vector<int>,int> Problem::bruteForce(int instanceIndex) const {
        const Instance& inst = instances_[instanceIndex];
        int n = inst.n_jobs;
    
        // guard przed zbyt dużym n
        if (n > 12) {
            return {{}, 0};
        }
    
        std::vector<int> perm(n);
        std::iota(perm.begin(), perm.end(), 0);
        std::vector<int> bestPerm = perm;
        int bestCmax = computeCmax(perm, instanceIndex);
    
        do {
            int c = computeCmax(perm, instanceIndex);
            if (c < bestCmax) {
                bestCmax = c;
                bestPerm = perm;
            }
        } while (std::next_permutation(perm.begin(), perm.end()));
    
        return {bestPerm, bestCmax};
    }
    
    std::pair<std::vector<int>, int> Problem::neh(int instanceIndex) const {
        const Instance& inst = instances_[instanceIndex];
        int n = inst.n_jobs;
        int m = inst.m_machines;
    
        // 1. Oblicz sumy czasów dla każdego zadania i posortuj malejąco
        std::vector<std::pair<int,int>> sums(n);
        for (int j = 0; j < n; ++j) {
            int s = 0;
            for (int i = 0; i < m; ++i)
                s += inst.times[j][i];
            sums[j] = { -s, j };
        }
        std::sort(sums.begin(), sums.end());
    
        // 2. Buduj sekwencję stopniowo
        std::vector<int> sequence;
        for (auto [_, job] : sums) {
            std::vector<int> bestSeq;
            int bestC = std::numeric_limits<int>::max();
    
            // wstaw job we wszystkie możliwe pozycje
            for (int pos = 0; pos <= static_cast<int>(sequence.size()); ++pos) {
                auto temp = sequence;
                temp.insert(temp.begin() + pos, job);
                int c = computeCmax(temp, instanceIndex);
                if (c < bestC) {
                    bestC = c;
                    bestSeq = std::move(temp);
                }
            }
            sequence = std::move(bestSeq);
        }
    
        int finalC = computeCmax(sequence, instanceIndex);
        return {sequence, finalC};
    }
    
    void Problem::printSolution(const std::vector<int>& sequence, int cmax) const {
        std::cout << "Kolejnosc zadan: ";
        for (int j : sequence)
            std::cout << j << " ";
        std::cout << "\nCmax = " << cmax << "\n";
    }