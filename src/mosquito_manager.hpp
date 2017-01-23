#pragma once

#include "mosquito.hpp"
#include <vector>

class MosquitoManager
{
private:
    unsigned int numMosquitoes=0;
    std::vector<Mosquito>* mosquitoes;
    std::vector<unsigned int> inactiveMosquitoes;
    std::vector<unsigned int> activeMosquitoes;

public:
    void initialise(std::vector<Mosquito>* mosquitoesArray);
    unsigned int get_count() const { return numMosquitoes; }
    void remove_mosquito(unsigned int numToRemove = 1);
    void add_mosquito(unsigned int numToAdd = 1);
    void modify_population(int numToChange = 0); //Just selects remove_mosquito / add_mosquito as appropriate.
    unsigned int random_active_mos() const;
};
