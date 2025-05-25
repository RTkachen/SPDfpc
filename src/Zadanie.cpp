#include "Zadanie.hpp"

Zadanie::Zadanie(const std::string& filename) {
    readTailDat(filename);
}

const std::vector<Instance>& Zadanie::getInstances() const {
    return instances_;
}

void Zadanie::readTailDat(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) {
        std::cerr << "Blad: nie mozna otworzyc pliku " << filename << std::endl;
        return;
    }

    int total_instances;
    in >> total_instances;
    instances_.clear();
    instances_.reserve(total_instances);

    for (int inst = 0; inst < total_instances; ++inst) {
        Instance I;
        in >> I.n_jobs >> I.m_machines;
        I.times.assign(I.n_jobs, std::vector<int>(I.m_machines));

        for (int j = 0; j < I.n_jobs; ++j) {
            for (int k = 0; k < I.m_machines; ++k) {
                int machine_id, ptime;
                in >> machine_id >> ptime;
                I.times[j][machine_id] = ptime;
            }
        }
        instances_.push_back(std::move(I));
    }
}

void Zadanie::printInstances() const {
    for (size_t idx = 0; idx < instances_.size(); ++idx) {
        const auto& inst = instances_[idx];
        std::cout << "Instancja " << idx + 1
                  << ": n_jobs=" << inst.n_jobs
                  << ", m_machines=" << inst.m_machines << std::endl;
        for (int j = 0; j < inst.n_jobs; ++j) {
            std::cout << " Job " << j << ":";
            for (int k = 0; k < inst.m_machines; ++k) {
                std::cout << " [M" << k << ":" << inst.times[j][k] << "]";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}