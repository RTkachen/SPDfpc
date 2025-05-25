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
    
    private:
        std::vector<Instance> instances_;  // wszystkie instancje
    };

#endif // PROBLEM_H