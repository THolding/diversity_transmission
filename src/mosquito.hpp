#pragma once
#include "infection.hpp"
#include "host.hpp"

class Output;

BITE_FREQUENCY_TABLE initialise_cumulative_bite_frequency_distribution();

class Mosquito
{
public:
    unsigned int age = 0; //Age in days.
    Infection infection; //Infection::active = false, by default.

    void infect(const Strain& strain, bool allowRecombination);
    void age_mosquito(const PTABLE& pDeath);
    void kill();
    void update_infection();
    void feed(Host& host, Output* output = nullptr, bool allowRecombination = true);

    bool is_infected() const { return infection.infected; }
};
