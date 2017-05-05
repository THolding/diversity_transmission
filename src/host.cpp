#include "host.hpp"
#include "utilities.hpp"
#include "diversity_monitor.hpp"
#include <cmath>

#include <iostream>

//Attempt to infect a host.
void Host::infect(const Strain& strain)
{
    #pragma omp critical (host_infection)
    {
        if (!infection1.infected)
        {
            unsigned int projectedDuration = duration_kernal(strain, immuneState);
            if (projectedDuration > 0) {
                infection1.infected = true;
                infection1.strain = strain;
                infection1.infectivity = infectivity_kernal(strain, immuneState);
                infection1.durationRemaining = duration_kernal(strain, immuneState);
                exposure_kernal(strain, immuneState); //not needed as duration_kernal does this now too...
                //if (infection1.durationRemaining > 0)
                    DiversityMonitor::register_new_strain(strain);
                //std::cout << "host infected#1\tduration:" << projectedDuration <<"\n";
            }

        } else if (!infection2.infected)
        {
            unsigned int projectedDuration = duration_kernal(strain, immuneState);
            if (projectedDuration > 0) {
                infection2.infected = true;
                infection2.strain = strain;
                infection2.infectivity = infectivity_kernal(strain, immuneState);
                infection2.durationRemaining = duration_kernal(strain, immuneState);
                exposure_kernal(strain, immuneState); //not needed as duration_kernal does this now too...
                //if (infection2.durationRemaining > 0)
                    DiversityMonitor::register_new_strain(strain);
                //std::cout << "host infected#2\tduration:" << projectedDuration <<"\n";
            }
        }
    }
}

//Age host and kill / replace it with newborn if necessary.
void Host::age_host(const PTABLE& pDeath)
{
    if (utilities::random_float01() < pDeath[std::floor(age / 365)]) { //If the host dies.
        kill();
        //++deathCount;
    }
    else
        ++age;
}

void Host::kill()
{
    age = 0;
    infection1.reset();
    infection2.reset();
    std::fill(immuneState.begin(), immuneState.end(), 0.0);
    //std::cout << "Host died\n";
}

void Host::update_infections()
{
    if (infection1.infected)
    {
        if (infection1.durationRemaining <= 0)
            infection1.reset();
        else
            infection1.durationRemaining--;
    }
    if (infection2.infected)
    {
        if (infection2.durationRemaining <= 0)
            infection2.reset();
        else
            infection2.durationRemaining--;
    }
}








