#ifndef PROBLEM_H
#define PROBLEM_H

#include <vector>
#include <iostream>
#include <algorithm>  // do next_permutation
#include <limits>     // dla std::numeric_limits
#include <numeric>
#include <utility>
#include <bits/stdc++.h>
#include "Struct.hpp"
class Problem {
    public:
        // Konstruktor przyjmujący wektor instancji
        explicit Problem(const std::vector<Instance>& instances);
    
        // Pełny przegląd zupełny dla instancji o danym indeksie
        // Zwraca parę <najlepsza_perm, najlepszy_Cmax>
        std::pair<std::vector<int>, int> bruteForce(int instanceIndex) const;
    
        // Algorytm NEH dla instancji o danym indeksie
        // Zwraca parę <najlepsza_perm, najlepszy_Cmax>
        std::pair<std::vector<int>, int> neh(int instanceIndex) const;
    
        // Oblicza Cmax dla sekwencji zadań w danej instancji
        int computeCmax(const std::vector<int>& sequence, int instanceIndex) const;
    
        // Wyświetla rozwiązanie: kolejność i Cmax
        void printSolution(const std::vector<int>& sequence, int cmax) const;

        std::pair<std::vector<int>, int> quickNEH(int instanceIndex) const;

        std::pair<std::vector<int>, int> johnson(int instanceIndex) const;

        std::pair<std::vector<int>, int> branchAndBound(int instanceIndex) const;

        std::pair<std::vector<int>, int> simulatedAnnealing(int instanceIndex) const;
    
    private:
        std::vector<Instance> instances_;  // wszystkie instancje
        std::vector<int> updateDPBackward(const std::vector<int>& currentDP, const std::vector<int>& processingTimes) const;
        std::vector<int> updateDP(const std::vector<int>& currentDP, const std::vector<int>& processingTimes) const;
        std::vector<std::vector<int>> computePrefixDP(const std::vector<int>& seq, int instanceIndex) const;
        std::vector<std::vector<int>> computeSuffixDP_Exact(const std::vector<int>& seq, int instanceIndex) const;
        //std::vector<std::vector<int>> computeSuffixDP(const std::vector<int>& seq, int instanceIndex) const;
        int quickEvaluateInsertion(const std::vector<int>& sequence,
            int insertPos,
            int job,
            int instanceIndex,
            const std::vector<std::vector<int>>& prefixDP,
            const std::vector<std::vector<int>>& suffixDP) const;
        void bnbRecursive(std::vector<int>& partial,
                      std::vector<bool>& used,
                      int instanceIndex,
                      int& bestCmax,
                      std::vector<int>& bestSequence) const; 
    };

#endif // PROBLEM_H