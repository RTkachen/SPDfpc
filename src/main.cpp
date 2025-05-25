#include "Zadanie.hpp"
#include "Problem.hpp"

// Funkcja pomocnicza do mierzenia czasu działania
/*void measureExecutionTime(Problem& problem, const std::function<void()>& func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end   = std::chrono::high_resolution_clock::now();

    auto us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "="<< us << "*0,000001"<< '\n' ;
    std::cout << problem.maxSum() << '\n';
}*/

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Użycie: " << argv[0] << " <ścieżka do pliku>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    Zadanie zadanie(filename);
    const auto& instVec = zadanie.getInstances();
    int total = static_cast<int>(instVec.size());
    Problem problem(instVec);

    for (int idx = 0; idx < total; ++idx) {
        std::cout << "\n--- Instancja " << idx+1 << " z " << total << " ---" << std::endl;

        // Pełny przegląd zupełny
        auto [seq, cmax] = problem.bruteForce(idx);
        if (!seq.empty()) {
            problem.printSolution(seq, cmax);
        }

        // Heurystyka NEH
        auto [nehSeq, nehCmax] = problem.neh(idx);
        problem.printSolution(nehSeq, nehCmax);
    }
 
    return 0;
}