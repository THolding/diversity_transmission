#include "mosquito.hpp"
#include "output.hpp"
#include "utilities.hpp"
#include <cmath>

BITE_FREQUENCY_TABLE initialise_cumulative_bite_frequency_distribution()
{
    BITE_FREQUENCY_TABLE pdfPoisson;
    for (unsigned int k=0; k<pdfPoisson.size(); k++)
        pdfPoisson[k] = ((float)std::pow(BITE_RATE, k) * std::exp(-BITE_RATE)) / (float) utilities::factorial(k);
    utilities::arrayToFile(pdfPoisson, FILE_PATH+RUN_NAME+"_pdfBiteFrequency.csv");

    //Calculate cumulative of pdfPoisson.
    BITE_FREQUENCY_TABLE cumulativeBiteFrequency;
    cumulativeBiteFrequency[0] = pdfPoisson[0];
    for (unsigned int i=1; i<pdfPoisson.size(); i++)
        cumulativeBiteFrequency[i] = cumulativeBiteFrequency[i-1]+pdfPoisson[i];

    return cumulativeBiteFrequency;
}

void Mosquito::infect(const Strain& strain, bool allowRecombination)
{
    if (!infection.infected) //Can only be infected once.
    {
        //float p = utilities::random_float01();
        //if (p < infection.infectivity)
        //if (utilities::random_float01() < infection.infectivity) {
            infection.infected = true;
            if (allowRecombination)
                infection.strain = generate_recombinant_strain(strain);
            else
                infection.strain = strain;
            infection.infectivity = 1.0;
            infection.durationRemaining = MOSQUITO_EIP;
        //}
    }
}

void Mosquito::age_mosquito(const PTABLE& pDeath)
{
    if (utilities::random_float01() < pDeath[age])
        kill();
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
    if (!infection.infected)
    {
        //If host has two infections there is a chance of recombination.
        if (host.infection1.infected && host.infection2.infected && allowRecombination)
        {
            //Choose strain at random to be primary parent.
            if (utilities::urandom(0,2) == 0)
                infect(generate_recombinant_strain(host.infection1.strain, host.infection2.strain), allowRecombination);
            else
                infect(generate_recombinant_strain(host.infection2.strain, host.infection1.strain), allowRecombination);
        }
        else if (host.infection1.infected && utilities::random_float01() < host.infection1.infectivity) //Can only be one infection so clone it.
            infect(host.infection1.strain, allowRecombination);
        else if (host.infection2.infected && utilities::random_float01() < host.infection2.infectivity) //Can still only be one infection so clone it.
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
