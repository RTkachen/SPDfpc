#ifndef ZADANIE_HPP
#define ZADANIE_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Struct.hpp"

class Zadanie {
    public:
        // Konstruktor wczytujący instancje z pliku
        explicit Zadanie(const std::string& filename);

        // Zwraca wektor wszystkich wczytanych instancji
        const std::vector<Instance>& getInstances() const;

        // Wypisuje zawartość wszystkich instancji na stdout
        void printInstances() const;
    
    private:
        // Pomocnicza funkcja do wczytywania pliku
        void readTailDat(const std::string& filename);
        
        // Przechowywane instancje
        std::vector<Instance> instances_;
    };

#endif // ZADANIE_HPP