#pragma once
#include <array>
#include "infection.hpp"

class Host
{
public:
    unsigned int age = 0; //In days
    Infection infection1; //Infection::active = false, by default.
    Infection infection2;
    ImmuneState immuneState;

    void infect(const Strain& strain);
    void age_host(const PTABLE& pDeath);
    void kill();
    void update_infections();

    bool is_infected() const { return infection1.infected || infection2.infected; }

    Host() : immuneState(NUM_PHENOTYPES, 0.0) {  }
};
