#pragma once
#include "infection.hpp"
#include "param_manager.hpp"

class Host
{
public:
    unsigned int age = 0; //In days
    Infection infection1; //Infection::active = false, by default.
    Infection infection2;
    ImmuneState immuneState;

    Host() : immuneState(ParamManager::instance().get_int("num_phenotypes"), 0.0) {  }

    void infect(const Strain& strain);
    void age_host(const PTABLE& pDeath);
    void kill();
    void update_infections();

    bool is_infected() const { return infection1.infected || infection2.infected; }
};
