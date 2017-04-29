#include "testing.hpp"
#include "host.hpp"
#include "strain.hpp"
#include "infection.hpp"
#include "mosquito.hpp"
#include "utilities.hpp"
#include "output.hpp"
#include "demographic_tools.hpp"
#include <iostream>
#include <vector>

void testing::new_tests()
{
    utilities::initialise_random();
    BITE_FREQUENCY_TABLE cumulativeBiteFrequencyDistribution;
    cumulativeBiteFrequencyDistribution = ParamManager::get_cumulative_bite_frequency_distribution();
    PTABLE pDeathMosquitoes;
    pDeathMosquitoes = generate_mosquito_ptable();
    Output output(nullptr);
    output.preinitialise_output_storage();

    Strain strain = strain_from_antigen_pool({0*128, 1*128, 2*128, 3*128, 4*128, 5*128});
    Host host;
    host.kill();
    host.age = 2;
    host.infect(strain);
    std::cout << "host infection1: " << strain_phenotype_str(host.infection1.strain) << "\n";
    std::cout << "host infection2: " << strain_phenotype_str(host.infection2.strain) << "\n";

    std::vector<Mosquito> mosquitoes;
    for (unsigned int i=0; i<20; i++) {
        Mosquito newMos;
        newMos.kill();
        newMos.age = 1;
        newMos.active = true;

        mosquitoes.push_back(newMos);
    }

    std::cout << "Mosquito infection status\n";
    for (unsigned int i=0; i<mosquitoes.size(); ++i) {
        if (mosquitoes[i].is_infected())
            std::cout << i << ", infected\n";
        else
            std::cout << i << ", uninfected\n";
    }

    //feed once
    for (unsigned int i=0; i<mosquitoes.size(); ++i)
    {
        float p = utilities::random_float01()*cumulativeBiteFrequencyDistribution.back();
        unsigned int numBites = 0;
        while (p > cumulativeBiteFrequencyDistribution[numBites])
            ++numBites;

        //For each feed, select a host and feed.
        std::cout << "mosquito " << i << " feeds " << numBites << " times...\n";
        for (unsigned int b=0; b<numBites; ++b)
        {
            mosquitoes[i].feed(host, &output, false);
        }
    }

    std::cout << "Mosquito infection status\n";
    std::vector<unsigned int> infectedMossys;
    for (unsigned int i=0; i<mosquitoes.size(); ++i) {
        if (mosquitoes[i].is_infected()) {
            std::cout << i << ", infected\n";
            infectedMossys.push_back(i);
        }
        else
            std::cout << i << ", uninfected\n";
    }
    std::cout << "\n\nThere are " << infectedMossys.size() << " infected mosquitoes.\n";

    if (infectedMossys.size() >= 1) {
        Mosquito infMos = mosquitoes[infectedMossys[0]];
        if (infMos.is_infected())
            std::cout << "At least one infected mosquito aged: " << infMos.age << "\n" << "Infecting strain:\n" << strain_phenotype_str(infMos.infection.strain) << "\n\n\n";

        std::cout << "\nAging first infected mosquito...\n";
        while (infMos.age != 0)
        {
            infMos.age_mosquito(pDeathMosquitoes);
            std::cout << "New infected mosquito age: " << infMos.age << "\n";
        }

        std::cout << "Infected mosquito has died!\n";
        std::cout << "Replacement individual specifications:\n";
        std::cout << "\tage :" << infMos.age << "\n";
        std::cout << "\tinfection status : ";
        if (infMos.is_infected()) {
            std::cout << "infected\n";
            std::cout << "strain:\n";
            std::cout << strain_phenotype_str(infMos.infection.strain) << "\n\n";
        }
        else std::cout << "uninfected\n";
    }
}






















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
