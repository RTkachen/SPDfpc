#include "Zadanie.hpp"
#include "Problem.hpp"

// -----------------------------------------------------------------------------
// Funkcja measureExecutionTime
// -----------------------------------------------------------------------------
// Mierzy czas wykonania przekazanej lambdy (func) i wypisuje wynik
// w sekundach w notacji naukowej (np. 1.234e-03 s) wraz z etykietą label.
void measureExecutionTime(const std::function<void()>& func, const std::string& label) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end   = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << label << " execution time: " 
              << std::scientific << elapsed.count() << " s" << std::endl;
}

// Szablonowa funkcja pomocnicza do pomiaru czasu wykonania
template <typename F>
double runAlgorithm(F func) {
    auto start = std::chrono::high_resolution_clock::now();
    func();
    auto end   = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Użycie: " << argv[0] << " <ścieżka do pliku>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    // Wczytanie zadania oraz instancji z pliku wejściowego
    Zadanie zadanie(filename);
    const auto& instVec = zadanie.getInstances();
    int total = static_cast<int>(instVec.size());
    Problem problem(instVec);

    // Dla każdej instancji wywołujemy wszystkie algorytmy
    for (int idx = 0; idx < total; ++idx) {
        std::cout << "\nInstancja " << idx + 1 << std::endl;
        
        std::pair<std::vector<int>, int> bfResult, johnsonResult, bnbResult, quickResult, nehResult, saResult;
        //double bfTime = runAlgorithm([&]() { bfResult = problem.bruteForce(idx); });
        //double johnsonTime = runAlgorithm([&]() { johnsonResult = problem.johnson(idx); });
        //double bnbTime = runAlgorithm([&]() { bnbResult = problem.branchAndBound(idx); });
        //double quickTime = runAlgorithm([&]() { quickResult = problem.quickNEH(idx); });
        double nehTime = runAlgorithm([&]() { nehResult = problem.neh(idx); });
        double saTime = runAlgorithm([&]() { saResult = problem.simulatedAnnealing(idx); });

        // Wypisujemy wyniki – format: <algorytm>: <czas>;<Cmax>
        // Używamy std::scientific do zapisu czasu w notacji naukowej.
        //std::cout <<  std::scientific << bfTime << ";" << bfResult.second << std::endl;
        //std::cout << std::scientific << johnsonTime << ";" << johnsonResult.second << std::endl;
        //std::cout <<  std::scientific << bnbTime << ";" << bnbResult.second << std::endl;
        //std::cout <<  std::scientific << quickTime << ";" << quickResult.second << std::endl;
        std::cout <<  std::scientific << nehTime << ";" << nehResult.second << std::endl;
        std::cout <<  std::scientific << saTime << ";" << saResult.second << std::endl;
    }

    return 0;
}
