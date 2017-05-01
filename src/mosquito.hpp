#pragma once
#include "infection.hpp"
#include "host.hpp"

class Output;

class Mosquito
{
public:
    unsigned int age = 0; //Age in days.
    Infection infection; //Infection::active = false, by default.
    bool active = true;

    void infect(const Strain& strain, bool allowRecombination, bool bypassGenerationRegister = false); //bypassGeneratioNRegister prevents antigens being registered as newly generated antigens
    void age_mosquito(const PTABLE& pDeath);
    void kill();
    void update_infection();
    void feed(Host& host, Output* output = nullptr, bool allowRecombination = true);

    bool is_infected() const { return infection.infected; }
    bool is_active() const { return active; }
};
