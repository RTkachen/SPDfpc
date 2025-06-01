#include "Problem.hpp"
Problem::Problem(const std::vector<Instance>& instances)
    : instances_(instances) {}

// Funkcja computeCmax
// --------------------
// Cel: Oblicza maksymalny czas ukończenia (Cmax) dla zadanej sekwencji zadań na danej instancji problemu.
// Opis: Dla danej kolejności zadań (sequence) funkcja krok po kroku oblicza czasy zakończenia na każdej maszynie, 
//       bazując na regule: moment ukończenia zadania = max(zakończenie poprzedniego zadania na tej maszynie,
//       zakończenie tego samego zadania na poprzedniej maszynie) + czas przetwarzania.
// Parametry:
//   - sequence: wektor numerów zadań (indeksy) określający kolejność ich wykonania,
//   - instanceIndex: indeks instancji w tablicy instancji, z której pobierane są dane (liczba maszyn oraz czasy przetwarzania).
// Zwraca: Obliczony Cmax (całkowity czas ukończenia ostatniego zadania na ostatniej maszynie).
int Problem::computeCmax(const std::vector<int>& sequence, int instanceIndex) const {
    // Pobranie instancji problemu na podstawie indeksu
    const Instance& inst = instances_[instanceIndex];

    // Liczba zadań, które należy uwzględnić (na podstawie długości podanej sekwencji)
    int k = static_cast<int>(sequence.size());
    // Liczba maszyn dostępnych w instancji
    int m = inst.m_machines;
    
    // Inicjalizacja macierzy C: tablica (k+1) x (m+1) przechowująca czasy zakończenia zadań.
    // Pierwszy wiersz i pierwsza kolumna pozostają zerowe, co upraszcza obliczenia (brzegowe warunki).
    std::vector<std::vector<int>> C(k+1, std::vector<int>(m+1, 0));
    
    // Iteracja po zadaniach (indeks j) i maszynach (indeks i)
    for (int j = 1; j <= k; ++j) {
        for (int i = 1; i <= m; ++i) {
            // Pobranie numeru zadania z sekwencji (uwzględnienie różnicy w indeksowaniu: sekwencja 0-indeksowana)
            int job = sequence[j-1];
            // Pobranie czasu przetwarzania zadania 'job' na maszynie 'i'
            int p   = inst.times[job][i-1];
            // Obliczenie czasu zakończenia:
            // Wybieramy maksymalny czas zakończenia między:
            // - zakończeniem tego samego zadania na poprzedniej maszynie (C[j][i-1]),
            // - zakończeniem poprzedniego zadania na tej samej maszynie (C[j-1][i])
            // Następnie dodajemy czas przetwarzania p.
            C[j][i] = std::max(C[j][i-1], C[j-1][i]) + p;
        }
    }
    // Zwracamy Cmax, czyli czas ukończenia ostatniego zadania na ostatniej maszynie.
    return C[k][m];
}

// Funkcja bruteForce
// --------------------
// Cel: Znalezienie optymalnej (minimalnego Cmax) sekwencji zadań poprzez pełne przeszukanie przestrzeni permutacji.
// Opis: Funkcja generuje wszystkie możliwe permutacje zadań (o ile liczba zadań nie przekracza 12),
//       oblicza Cmax dla każdej sekwencji i wybiera tę o najmniejszym Cmax.
// Parametry:
//   - instanceIndex: indeks instancji problemu, dla której szukamy optymalnego harmonogramu.
// Zwraca: Parę (bestPerm, bestCmax) gdzie bestPerm to optymalna sekwencja zadań, a bestCmax to odpowiadający jej Cmax.
std::pair<std::vector<int>,int> Problem::bruteForce(int instanceIndex) const {
    // Pobranie instancji problemu
    const Instance& inst = instances_[instanceIndex];
    int n = inst.n_jobs;
    
    // Guard przed zbyt dużą liczbą zadań
    if (n > 12) {
        // Zwracamy pustą sekwencję, gdy n > 12 (brak możliwości obliczenia z uwagi na zbyt duży koszt obliczeniowy)
        return {{}, 0};
    }
    
    // Inicjalizacja początkowej sekwencji zadań: {0, 1, 2, ..., n-1}
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::vector<int> bestPerm = perm;
    // Obliczenie Cmax dla początkowej sekwencji
    int bestCmax = computeCmax(perm, instanceIndex);
    
    // Przegląd wszystkich możliwych permutacji zadań
    do {
        // Obliczenie Cmax dla bieżącej permutacji
        int c = computeCmax(perm, instanceIndex);
        // Jeśli znajdziemy lepsze rozwiązanie (mniejsze Cmax), aktualizujemy najlepszy wynik i sekwencję
        if (c < bestCmax) {
            bestCmax = c;
            bestPerm = perm;
        }
    } while (std::next_permutation(perm.begin(), perm.end()));
    
    // Zwracamy najlepszą sekwencję oraz odpowiadające jej Cmax
    return {bestPerm, bestCmax};
}

// Funkcja neh
// -----------
// Cel: Znalezienie przybliżonej (heurystycznej) sekwencji zadań minimalizującej Cmax za pomocą algorytmu NEH.
// Opis: Algorytm NEH składa się z trzech głównych etapów:
//   1. Wyznaczenie sum czasów przetwarzania dla każdego zadania i posortowanie zadań malejąco według tej sumy.
//   2. Stopniowe budowanie sekwencji: dla każdego zadania próbujemy wstawić je na każdą możliwą pozycję w aktualnej sekwencji
//      i wybieramy tę pozycję, która daje najmniejszy Cmax.
//   3. Zwrócenie ostatecznej sekwencji i obliczonego Cmax.
// Parametry:
//   - instanceIndex: indeks instancji problemu.
// Zwraca: Parę (sequence, finalC) gdzie sequence to wyznaczona sekwencja zadań, a finalC to wynikowy Cmax.
std::pair<std::vector<int>, int> Problem::neh(int instanceIndex) const {
    // Pobranie instancji problemu i ustalenie liczby zadań oraz maszyn
    const Instance& inst = instances_[instanceIndex];
    int n = inst.n_jobs;
    int m = inst.m_machines;
    
    // 1. Obliczenie sum czasów przetwarzania dla każdego zadania.
    // Sumy traktujemy jako ujemne wartości, co pozwala na sortowanie w kolejności malejącej.
    std::vector<std::pair<int,int>> sums(n);
    for (int j = 0; j < n; ++j) {
        int s = 0;
        for (int i = 0; i < m; ++i)
            s += inst.times[j][i];
        // Przechowujemy ujemną sumę (aby późniejsze sortowanie dawało kolejność malejącą) oraz indeks zadania.
        sums[j] = { -s, j };
    }
    // Sortowanie zadań według obliczonej sumy czasów (malejąco po użyciu ujemnych wartości)
    std::sort(sums.begin(), sums.end());
    
    // 2. Budowanie sekwencji zadań poprzez stopniowe wstawianie.
    std::vector<int> sequence;
    // Iteracja dla każdego zadania w posortowanej kolejności
    for (auto [_, job] : sums) {
        std::vector<int> bestSeq; // Tymczasowa najlepsza sekwencja dla bieżącego wstawiania
        int bestC = std::numeric_limits<int>::max(); // Początkowo ustawiamy najlepszy Cmax na bardzo wysoką wartość
        
        // Próba wstawienia zadania 'job' we wszystkie możliwe pozycje w dotychczasowej sekwencji
        for (int pos = 0; pos <= static_cast<int>(sequence.size()); ++pos) {
            auto temp = sequence;                // Kopia aktualnej sekwencji
            temp.insert(temp.begin() + pos, job);  // Wstawienie zadania na pozycję 'pos'
            int c = computeCmax(temp, instanceIndex); // Obliczenie Cmax dla tej tymczasowej sekwencji
            // Jeśli uzyskany Cmax jest lepszy niż dotychczasowy najlepszy wynik, aktualizujemy najlepszą sekwencję
            if (c < bestC) {
                bestC = c;
                bestSeq = std::move(temp);
            }
        }
        // Aktualizacja głównej sekwencji: dodajemy zadanie w miejscu zapewniającym najniższy Cmax
        sequence = std::move(bestSeq);
    }
    
    // Obliczenie końcowego Cmax dla uzyskanej sekwencji
    int finalC = computeCmax(sequence, instanceIndex);
    // Zwracamy ostateczną sekwencję oraz odpowiadający jej Cmax
    return {sequence, finalC};
}

std::pair<std::vector<int>, int> Problem::johnson(int instanceIndex) const {
    // Pobieramy instancję problemu
    const Instance& inst = instances_[instanceIndex];
    int n = inst.n_jobs;
    int m = inst.m_machines;

    // Algorytm Johnsona jest optymalny dla 2 maszyn. Jeśli nie mamy dokładnie 2 maszyn,
    // informujemy użytkownika i przerywamy wykonanie.
    if (m != 2) {
        std::cerr << "Johnson algorithm jest dostępny tylko dla 2 maszyn!" << std::endl;
        return { {}, 0 };
    }

    // Dzielimy zadania na dwie grupy:
    // Group A: zadania, dla których czas na maszynie 1 (p1) <= czas na maszynie 2 (p2)
    // Group B: zadania, dla których p1 > p2
    std::vector<int> groupA;
    std::vector<int> groupB;
    for (int j = 0; j < n; ++j) {
        int p1 = inst.times[j][0];
        int p2 = inst.times[j][1];
        if (p1 <= p2)
            groupA.push_back(j);
        else
            groupB.push_back(j);
    }

    // Sortujemy grupę A rosnąco według p1
    std::sort(groupA.begin(), groupA.end(), [&](int a, int b) {
        return inst.times[a][0] < inst.times[b][0];
    });

    // Sortujemy grupę B malejąco według p2
    std::sort(groupB.begin(), groupB.end(), [&](int a, int b) {
        return inst.times[a][1] > inst.times[b][1];
    });

    // Łączymy obie grupy: najpierw A, potem B
    std::vector<int> sequence;
    sequence.reserve(n);
    sequence.insert(sequence.end(), groupA.begin(), groupA.end());
    sequence.insert(sequence.end(), groupB.begin(), groupB.end());

    // Obliczamy ostateczne C_max dla otrzymanej sekwencji, korzystając z
    // istniejącej funkcji computeCmax (kompatybilnej z naszym modelem danych).
    int finalCmax = computeCmax(sequence, instanceIndex);

    return { sequence, finalCmax };
}

// Funkcja printSolution
// ----------------------
// Cel: Wypisanie na standardowe wyjście (konsolę) sekwencji zadań oraz odpowiadającej jej wartości Cmax.
// Opis: Przekazana sekwencja zadań jest iterowana i każdy numer zadania jest wypisywany.
//       Następnie wypisywany jest całkowity czas ukończenia (Cmax).
void Problem::printSolution(const std::vector<int>& sequence, int cmax) const {
   // std::cout << "Kolejnosc zadan: ";
    //for (int j : sequence)
        //std::cout << j << " ";
    std::cout << "\nCmax = " << cmax << "\n";
}




// Aktualizuje bieżący stan (wsteczny) DP dla przetworzenia zadania reprezentowanego przez processingTimes.
// currentDP – wektor długości m+1, gdzie indeks m odpowiada ostatniej maszynie.
// Zakładamy, że newDP[m] = currentDP[m] i robimy iterację od i = m-1 do 0.
std::vector<int> Problem::updateDPBackward(const std::vector<int>& currentDP, const std::vector<int>& processingTimes) const {
    int m = processingTimes.size();
    std::vector<int> newDP(m + 1, 0);
    newDP[m] = currentDP[m];
    for (int i = m - 1; i >= 0; --i) {
        newDP[i] = std::max(newDP[i + 1], currentDP[i]) + processingTimes[i];
    }
    return newDP;
}


// Funkcja pomocnicza updateDP
// -------------------------------------
// Aktualizuje bieżący stan DP (wektor currentDP o długości m+1)
// o przetworzenie zadania reprezentowanego przez wektor processingTimes (długość m).
// Schemat: newDP[i] = max(newDP[i-1], currentDP[i]) + processingTimes[i-1],
// dla i = 1..m, przy czym newDP[0] = currentDP[0] (w praktyce zakładamy, że jest to 0).
std::vector<int> Problem::updateDP(const std::vector<int>& currentDP, const std::vector<int>& processingTimes) const {
    int m = processingTimes.size();
    std::vector<int> newDP(m + 1, 0);
    newDP[0] = currentDP[0];
    for (int i = 1; i <= m; ++i) {
        newDP[i] = std::max(newDP[i - 1], currentDP[i]) + processingTimes[i - 1];
    }
    return newDP;
}

// Funkcja computePrefixDP
// -----------------------------------
// Dla danej sekwencji (vector<int> seq) i instancji o danym indeksie,
// zwraca vector DPState, czyli dla każdego 0 ≤ i ≤ L (L = seq.size())
// prefixDP[i] to wektor DP (długości m+1) osiągnięty po przetworzeniu pierwszych i zadań.
std::vector<std::vector<int>> Problem::computePrefixDP(const std::vector<int>& seq, int instanceIndex) const {
    const Instance& inst = instances_[instanceIndex];
    int m = inst.m_machines;
    int L = seq.size();
    std::vector<std::vector<int>> prefixDP(L + 1, std::vector<int>(m + 1, 0));
    // Stan początkowy – dla 0 zadań (wszystkie 0)
    for (int j = 0; j <= m; ++j)
        prefixDP[0][j] = 0;
    for (int i = 1; i <= L; ++i) {
        int job = seq[i - 1];
        prefixDP[i] = updateDP(prefixDP[i - 1], inst.times[job]);
    }
    return prefixDP;
}

// Funkcja computeSuffixDP_Exact
// Dla danej sekwencji seq i instancji o indeksie instanceIndex, oblicza DP dla suffixu (od i-tego zadania do końca)
// w sposób jawny, używając backward update.
std::vector<std::vector<int>> Problem::computeSuffixDP_Exact(const std::vector<int>& seq, int instanceIndex) const {
    const Instance& inst = instances_[instanceIndex];
    int m = inst.m_machines;
    int L = seq.size();
    std::vector<std::vector<int>> suffixDP(L + 1, std::vector<int>(m + 1, 0));
    
    // Stan początkowy dla pustego suffixu: wszystkie 0.
    suffixDP[L] = std::vector<int>(m + 1, 0);
    
    // Iterujemy wstecz (od ostatniego zadania w sekwencji do pierwszego)
    for (int i = L - 1; i >= 0; --i) {
        int job = seq[i];
        suffixDP[i] = updateDPBackward(suffixDP[i + 1], inst.times[job]);
    }
    return suffixDP;
}

// --------------------------------------------------------------------------
// Funkcja quickEvaluateInsertion – wersja bidirectional
// ------------------------------------------
// Dla danej sekwencji (już ustalonej) oraz potencjalnego zadania 'job'
// oraz pozycji insertPos, wykorzystujemy prekomputowane stany DP z prefixu i suffixu.
// Estymujemy Cmax jako: 
//    estimated = max_{i=0}^{m} { candidateDP[i] + suffixDP[insertPos][i] },
// gdzie candidateDP = updateDP(prefixDP[insertPos], inst.times[job])
int Problem::quickEvaluateInsertion(const std::vector<int>& sequence,
    int insertPos,
    int job,
    int instanceIndex,
    const std::vector<std::vector<int>>& prefixDP,
    const std::vector<std::vector<int>>& suffixDP) const {
const Instance& inst = instances_[instanceIndex];
int m = inst.m_machines;

// Uzyskujemy stan DP dla prefixu już zcache'owany: prefixDP[insertPos]
std::vector<int> candidateDP = updateDP(prefixDP[insertPos], inst.times[job]);
int estimated = 0;
// Łączymy: dla każdego poziomu maszyny wybieramy sumę: candidateDP[i] + suffixDP[insertPos][i]
for (int i = 0; i <= m; ++i) {
estimated = std::max(estimated, candidateDP[i] + suffixDP[insertPos][i]);
}
return estimated;
}

// --------------------------------------------------------------------------
// Funkcja quickNEH – wykorzystująca zarówno DP forward (prefixDP) jak i backward (suffixDP)
// ------------------------------------------
std::pair<std::vector<int>, int> Problem::quickNEH(int instanceIndex) const {
    const Instance& inst = instances_[instanceIndex];
    int n = inst.n_jobs;
    int m = inst.m_machines;
    
    // Sortowanie zadań malejąco według sumy czasów (sztuczka z ujemnymi sumami)
    std::vector<std::pair<int, int>> sums(n);
    for (int j = 0; j < n; ++j) {
        int totalTime = 0;
        for (int i = 0; i < m; ++i) {
            totalTime += inst.times[j][i];
        }
        sums[j] = { -totalTime, j };
    }
    std::sort(sums.begin(), sums.end());
    
    // Inicjujemy sekwencję – pierwszy element
    std::vector<int> sequence;
    sequence.push_back(sums[0].second);
    
    // Iteracyjnie dodajemy kolejne zadania wg uporządkowania z 'sums'
    for (size_t idx = 1; idx < sums.size(); ++idx) {
        int job = sums[idx].second;
        int bestCandidate = std::numeric_limits<int>::max();
        int bestPos = 0;
        int L = sequence.size();
        // Obliczamy prefixDP dla bieżącej sekwencji
        std::vector<std::vector<int>> prefixDP = computePrefixDP(sequence, instanceIndex);
        // Obliczamy suffixDP dokładnie – używamy computeSuffixDP_Exact
        std::vector<std::vector<int>> suffixDP = computeSuffixDP_Exact(sequence, instanceIndex);
        
        // Sprawdzamy wszystkie możliwe pozycje wstawienia
        for (int pos = 0; pos <= L; ++pos) {
            int candidateCmax = quickEvaluateInsertion(sequence, pos, job, instanceIndex, prefixDP, suffixDP);
            if (candidateCmax < bestCandidate) {
                bestCandidate = candidateCmax;
                bestPos = pos;
            }
        }
        // Wstawiamy kandydujące zadanie w znalezionej najlepszej pozycji
        sequence.insert(sequence.begin() + bestPos, job);
    }
    
    // Obliczamy końcowy Cmax dla całej sekwencji metodą klasyczną
    int finalCmax = computeCmax(sequence, instanceIndex);
    return { sequence, finalCmax };
}

// Rekurencyjna funkcja branch and bound
void Problem::bnbRecursive(std::vector<int>& partial,
    std::vector<bool>& used,
    int instanceIndex,
    int& bestCmax,
    std::vector<int>& bestSequence) const {
const Instance& inst = instances_[instanceIndex];
int n = inst.n_jobs;
int m = inst.m_machines;

// Jeżeli mamy już pełną sekwencję, porównujemy wartość C_max
if (partial.size() == static_cast<size_t>(n)) {
int cmax = computeCmax(partial, instanceIndex);
if (cmax < bestCmax) {
bestCmax = cmax;
bestSequence = partial;
}
return;
}

// Obliczenie bieżącego "kosztu" dla częściowej sekwencji:
// jeżeli partial jest pusta, zakładamy 0; inaczej wykorzystujemy computeCmax.
int currentCost = partial.empty() ? 0 : computeCmax(partial, instanceIndex);

// Obliczamy LB dla maszyny 1: suma czasów na maszynie 1 dla partial plus dla nieprzypisanych zadań
int partialMachine1 = 0;
for (int job : partial) {
partialMachine1 += inst.times[job][0];
}
int unscheduledMachine1 = 0;
for (int j = 0; j < n; ++j) {
if (!used[j]) {
unscheduledMachine1 += inst.times[j][0];
}
}
int LB1 = partialMachine1 + unscheduledMachine1;

// Obliczamy LB dla ostatniej maszyny (maszyny m)
// currentCost = ukończenie częściowej sekwencji na maszynie m
int unscheduledMachineM = 0;
for (int j = 0; j < n; ++j) {
if (!used[j]) {
unscheduledMachineM += inst.times[j][m-1];
}
}
int LB2 = currentCost + unscheduledMachineM;

// Prostym dolnym ogranicznikiem (LB) przyjmujemy maksimum z: currentCost, LB1 i LB2.
int LB = std::max({currentCost, LB1, LB2});

// Jeśli LB >= najlepszy znaleziony dotychczas bestCmax, przenosimy (przycinamy drzewo)
if (LB >= bestCmax)
return;

// W przeciwnym razie, dla każdego nieużytego zadania rozszerzamy częśćową sekwencję.
for (int j = 0; j < n; ++j) {
if (!used[j]) {
used[j] = true;
partial.push_back(j);
bnbRecursive(partial, used, instanceIndex, bestCmax, bestSequence);
partial.pop_back();
used[j] = false;
}
}
}

// Metoda branchAndBound – rozpoczyna rekurencję
std::pair<std::vector<int>, int> Problem::branchAndBound(int instanceIndex) const {
const Instance& inst = instances_[instanceIndex];
int n = inst.n_jobs;

int bestCmax = std::numeric_limits<int>::max();
std::vector<int> bestSequence;
std::vector<int> partial;
std::vector<bool> used(n, false);

bnbRecursive(partial, used, instanceIndex, bestCmax, bestSequence);

return { bestSequence, bestCmax };
}

// Algorytm symulowanego wyżarzania dla problemu flow shop.
// Zwraca parę: <najlepsza sekwencja, C_max dla tej sekwencji>.
std::pair<std::vector<int>, int> Problem::simulatedAnnealing(int instanceIndex) const {
    const Instance& inst = instances_[instanceIndex];
    int n = inst.n_jobs;
    
    // Początkowe rozwiązanie uzyskujemy z heurystyki NEH (możesz zmodyfikować – np. losowe lub NEH).
    std::pair<std::vector<int>, int> initSolution = neh(instanceIndex);
    std::vector<int> currentSequence = initSolution.first;
    int currentCost = initSolution.second;
    
    // Zachowujemy najlepsze znalezione rozwiązanie.
    std::vector<int> bestSequence = currentSequence;
    int bestCost = currentCost;
    
    // Ustawienia parametrów symulowanego wyżarzania.
    double T = 10000.0;      // temperatura początkowa – można dobrać eksperymentalnie
    double Tmin = 0.001;    // temperatura końcowa
    double alpha = 0.95;    // współczynnik chłodzenia (temperatura jest mnożona przez alpha)
    int iterPerTemp = 1000;  // liczba prób (iteracji) dla danej temperatury
    
    // Zainicjujemy generator liczb losowych.
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> distReal(0.0, 1.0);
    std::uniform_int_distribution<int> distPos(0, n - 1);
    
    // Główna pętla wyżarzania.
    while (T > Tmin) {
        for (int iter = 0; iter < iterPerTemp; ++iter) {
            std::vector<int> neighbor = currentSequence;
            // Losowo wybieramy dwa indeksy i zamieniamy miejscami zadania – to prosta operacja sąsiedztwa.
            int i = distPos(rng);
            int j = distPos(rng);
            std::swap(neighbor[i], neighbor[j]);
            
            // Obliczamy C_max dla rozwiązania sąsiedniego.
            int neighborCost = computeCmax(neighbor, instanceIndex);
            int delta = neighborCost - currentCost;
            
            // Jeśli sąsiad jest lepszy – lub zaakceptowany z pewnym prawdopodobieństwem, przyjmujemy.
            if (delta < 0) {
                currentSequence = neighbor;
                currentCost = neighborCost;
                
                if (currentCost < bestCost) {
                    bestCost = currentCost;
                    bestSequence = currentSequence;
                }
            } else {
                double acceptanceProb = std::exp(-delta / T);
                if (distReal(rng) < acceptanceProb) {
                    currentSequence = neighbor;
                    currentCost = neighborCost;
                }
            }
        }
        // Zmniejszamy temperaturę.
        T *= alpha;
    }
    
    return { bestSequence, bestCost };
}
