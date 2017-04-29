#include "mosquito.hpp"
#include "output.hpp"
#include "utilities.hpp"
#include "diversity_monitor.hpp"
#include <cmath>

//Enacts infection event to a mosquito if possible. Assumes any probabilistic factors affecting infection chance have been accounted for and infection is still going ahead.
void Mosquito::infect(const Strain& strain, bool allowRecombination)
{
    if (infection.infected == false) //Can only be infected once.
    {
        infection.infected = true;
        if (allowRecombination) {
            infection.strain = generate_recombinant_strain(strain);
            DiversityMonitor::register_new_strain(infection.strain);
        }
        else {
            infection.strain = strain;
            DiversityMonitor::register_new_strain(infection.strain);
        }
        infection.infectivity = 1.0;
        infection.durationRemaining = ParamManager::mosquito_eip;
    }
}

void Mosquito::age_mosquito(const PTABLE& pDeath)
{
    if (utilities::random_float01() < pDeath[age])
        kill();
    else
        ++age;
}

//killed and reborn
void Mosquito::kill()
{
    age = 0;
    infection.reset();
}

void Mosquito::update_infection()
{
    if (infection.infected && infection.durationRemaining > 0)
        --infection.durationRemaining;
}

void Mosquito::feed(Host& host, Output* output, bool allowRecombination)
{
    ///Host infecting mosquito (only if not already infected)
    if (infection.infected == false)
    {
        //If host has two infections then intergenic recombination occurs.
        if (host.infection1.infected && host.infection2.infected && allowRecombination)
        {
            //Choose strain at random to be primary parent.
            if (utilities::urandom(0,2) == 0)
                infect(generate_recombinant_strain(host.infection1.strain, host.infection2.strain), allowRecombination);
            else
                infect(generate_recombinant_strain(host.infection2.strain, host.infection1.strain), allowRecombination);
        }
        else if (host.infection1.infected && utilities::random_float01() < host.infection1.infectivity) //Can only be one infection so no intergenic recombination.
            infect(host.infection1.strain, allowRecombination);
        else if (host.infection2.infected && utilities::random_float01() < host.infection2.infectivity) //Can still only be one infection so no intergenic recombination.
            infect(host.infection2.strain, allowRecombination);
    }
    ///Handle mosquito infecting host
    else if (infection.durationRemaining <= 0) //Mosquito was already infected, so transmit to host if infectious
    {
        host.infect(infection.strain);
        if (output != nullptr) //Count infectious bites (to calculate EIR)
            output->register_infectious_bite();
    }
}
