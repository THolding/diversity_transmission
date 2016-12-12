#include "testing.hpp"
#include "global_parameters.hpp"
#include "host.hpp"
#include "strain.hpp"
#include "mosquito.hpp"
#include <iostream>

void testing::run_tests()
{
    test_immunity();
    test_host_infection();
}

void testing::test_immunity()
{
    //Does infection alter immune state?
    Host host;
    host.kill();
    Strain strain = strain_from_antigen_pool({0*128, 1*128, 2*128, 3*128, 4*128, 5*128});

    std::cout << "*****Testing immune state*****\n";
    std::cout << "strain phenotype:\n";
    for (const Antigen a : strain)
        std::cout << get_phenotype_id(a) << ", ";
    std::cout << "\n\n";

    std::cout << "old immune state:\n";
    for (unsigned int i=0; i<10; i++)
        std::cout << host.immuneState[i] << ", ";
    std::cout << "\n\n";

    host.infect(strain);
    std::cout << "new immune state:\n";
    for (unsigned int i=0; i<10; i++)
        std::cout << host.immuneState[i] << ", ";
    std::cout << "\n\n";

    std::cout << "induced infection length = " << host.infection1.durationRemaining << "\n";
    host.infection1.reset();
    std::cout << "after clearing, duration = " << host.infection1.durationRemaining << "\n";
    host.infect(strain);
    std::cout << "after reinfecting, duration = " << host.infection1.durationRemaining << "\n";

    std::cout << "\n\n\n";
}

void testing::test_host_infection()
{
    Host host;
    host.kill();
    Strain strain = strain_from_antigen_pool({0*128, 1*128, 2*128, 3*128, 4*128, 5*128});

    host.infect(strain);
    std::cout << "Infected strain. Duration = " << host.infection1.durationRemaining << " infection status: " << host.infection1.infected << "\n";
    host.update_infections();
    std::cout << "updated host infection once. Duration = " << host.infection1.durationRemaining << " infection status: " << host.infection1.infected << "\n";
    host.update_infections();
    std::cout << "updated host infection once. Duration = " << host.infection1.durationRemaining << " infection status: " << host.infection1.infected << "\n";
    host.update_infections();
    std::cout << "updated host infection once. Duration = " << host.infection1.durationRemaining << " infection status: " << host.infection1.infected << "\n";
    host.update_infections();
    std::cout << "updated host infection once. Duration = " << host.infection1.durationRemaining << " infection status: " << host.infection1.infected << "\n";
}
