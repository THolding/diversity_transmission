#include "host.hpp"
#include "utilities.hpp"
#include <cmath>

#include <iostream>

//Attempt to infect a host.
void Host::infect(const Strain& strain)
{
    if (!infection1.infected)
    {
        infection1.infected = true;
        infection1.strain = strain;
        infection1.infectivity = infectivity_kernal(strain, immuneState);
        infection1.durationRemaining = duration_kernal(strain, immuneState);
        //if (infection1.durationRemaining != 60)
        //    std::cout << "Host::infect, duration is not 60: " << infection1.durationRemaining << "\n";
        //exposure_kernal(strain, immuneState); //not needed as duration_kernal does this now too...
    } else if (!infection2.infected)
    {
        infection2.infected = true;
        infection2.strain = strain;
        infection2.infectivity = infectivity_kernal(strain, immuneState);
        infection2.durationRemaining = duration_kernal(strain, immuneState);
        //exposure_kernal(strain, immuneState); //not needed as duration_kernal does this now too...
    }
}

//Age host and kill / replace it with newborn if necessary.
void Host::age_host(const PTABLE& pDeath)
{
    if (utilities::random_float01() < pDeath[std::floor(age / 365)]) //If the host dies.
        kill();
}

void Host::kill()
{
    age = 0;
    infection1.reset();
    infection2.reset();
    std::fill(immuneState.begin(), immuneState.end(), 0.0);
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








