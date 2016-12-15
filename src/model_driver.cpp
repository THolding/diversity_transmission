#include "model_driver.hpp"
#include "demographic_tools.hpp"
#include "strain.hpp"
#include "utilities.hpp"
#include <cmath>
#include <vector>

void ModelDriver::initialise_model()
{
    std::cout << "initialising random number generator" << std::endl;
    utilities::initialise_random();

    //Initialise PTABLEs.
    std::cout << "initialising PTABLEs" << std::endl;
    pDeathHosts =  generate_host_ptable();
    pDeathMosquitoes = generate_mosquito_ptable();
    cdfHosts = calculate_host_cdf(pDeathHosts);
    cdfMosquitoes = calculate_mosquito_cdf(pDeathMosquitoes);
    cumulativeBiteFrequencyDistribution = initialise_cumulative_bite_frequency_distribution();

    //Initialise hosts
    std::cout << "initialising host demographics" << std::endl;
    hosts.reserve(NUM_HOSTS);
    for (unsigned int h=0; h<NUM_HOSTS; ++h)
    {
        Host host;
        host.kill();
        host.age = random_host_equilibrum_age(cdfHosts);
        hosts.push_back(host);
    }

    //Initialise mosquitoes
    std::cout << "initialising mosquito demographics" << std::endl;
    mosquitoes.reserve(NUM_MOSQUITOES_MAX);
    for (unsigned int m=0; m<NUM_MOSQUITOES; ++m)
    {
        Mosquito mosquito;
        mosquito.kill();
        mosquito.age = random_mosquito_equilibrium_age(cdfMosquitoes);
        mosquito.active = true;
        mosquitoes.push_back(mosquito);
    }

    mManager.initialise(&mosquitoes);

    //Generate initial pool of antigen diversity
    std::cout << "initialising antigen pool" << std::endl;
    std::vector<Antigen> initialAntigenPool;
    initialAntigenPool.reserve(INITIAL_ANTIGEN_DIVERSITY);
    for (unsigned int a=0; a<INITIAL_ANTIGEN_DIVERSITY; ++a)
        initialAntigenPool.push_back(random_antigen());

    //Generate initial pool of strains from available initial antigen diversity
    std::cout << "initialising strain pool" << std::endl;
    std::vector<Strain> initialStrainPool;
    initialStrainPool.reserve(INITIAL_NUM_STRAINS);
    for (unsigned int s=0; s<INITIAL_NUM_STRAINS; ++s)
        initialStrainPool.push_back(strain_from_antigen_pool(initialAntigenPool));

    //Initial infections (mosquitoes)
    std::cout << "initialising mosquito infections" << std::endl;
    for (unsigned int i=0; i<INITIAL_NUM_MOSQUITO_INFECTIONS; ++i)
    {
        unsigned int iS = utilities::urandom(0, initialStrainPool.size());
        unsigned int iM = mManager.random_active_mos();
        mosquitoes[iM].infect(initialStrainPool[iS], false);
    }
}

void ModelDriver::run_model()
{
    initialise_model();
    output.preinitialise_output_storage(RUN_TIME, OUTPUT_INTERVAL);

    std::cout << "cumulative bite frequency distribution:\n";
    for (float f : cumulativeBiteFrequencyDistribution)
        std::cout << f << ", ";
    std::cout << "\n";

    bool finished = false;
    unsigned int time = 0;
    burnInPeriod = BURN_IN_PERIOD;
    output.append_output(time, hosts, mosquitoes);
    while (!finished)
    {
        if (VERBOSE) {
            std::cout << "t=" << time << "\n";
            if (burnInPeriod > 0)
                std::cout << "burnIn left: " << burnInPeriod << "\n";
        }
        else if (time % OUTPUT_INTERVAL == 0) {
            std::cout << "t=" << time << ". ";
            if (burnInPeriod > 0)
                std::cout << "burnIn left: " << burnInPeriod << ". ";
        }

        //Dynamic parameters
        update_parameters(time);

        //Host demographics
        //std::cout << "aging hosts...\n";
        age_hosts();

        //Mosquito demographics
        //std::cout << "aging mosquitoes...\n";
        age_mosquitoes();

        //Update infections in hosts
        //std::cout << "updating host infections...\n";
        update_host_infections();

        //Update infections in mosquitoes
        //std::cout << "updating mosquito infections...\n";
        update_mosquito_infections();

        //mosquitoes feed
        //std::cout << "feeding mosquitoes...\n";
        feed_mosquitoes();

        //Update logging / data collection.
        if (time % OUTPUT_INTERVAL== 0)
            output.append_output(time, hosts, mosquitoes);

        //Update time and check stop condition.
        //std::cout << "incrementing time...\n";
        if (burnInPeriod > 0)
            --burnInPeriod;

        ++time;
        if (time > RUN_TIME)
            finished = true;
    }

    output.export_output();
}

void ModelDriver::age_hosts()
{
    for (Host& host : hosts)
        host.age_host(pDeathHosts);
}

void ModelDriver::age_mosquitoes()
{
    for (Mosquito& mosquito : mosquitoes)
        if (mosquito.is_active())
            mosquito.age_mosquito(pDeathMosquitoes);
}

void ModelDriver::update_host_infections()
{
    for (Host& host : hosts)
        host.update_infections();
}

void ModelDriver::update_mosquito_infections()
{
    for (Mosquito& mosquito : mosquitoes)
        if (mosquito.is_active())
            mosquito.update_infection();
}

void ModelDriver::feed_mosquitoes()
{
    bool allowRecombination;
    if (burnInPeriod <= 0)
    {
        allowRecombination = true;
    }
    else
    {
        allowRecombination = false;
    }

    for (Mosquito& mosquito : mosquitoes)
    {
        if (mosquito.is_active()) {
            //Calculate number of times the mosquito feeds.
            float p = utilities::random_float01()*cumulativeBiteFrequencyDistribution.back();
            unsigned int numBites = 0;
            while (p > cumulativeBiteFrequencyDistribution[numBites])
                ++numBites;

            //For each feed, select a host and feed.
            for (unsigned int b=0; b<numBites; ++b)
            {
                unsigned int iH = utilities::urandom(0, hosts.size());
                mosquito.feed(hosts[iH], &output, allowRecombination);
            }
        }
    }
}

//Updates dynamic parameters
void ModelDriver::update_parameters(const unsigned int time)
{
    if (DYN_NUM_MOSQUITOES == true)
    {
        if (time >= START_MOS_INC && time < STOP_MOS_INC)
        {
            int numLeftToAdd = NUM_MOSQUITOES_MAX - mManager.get_count();
            int timeLeftToAdd = STOP_MOS_INC - time;
            float fractionToAdd = (float)numLeftToAdd / (float)timeLeftToAdd;
            mosChangeRemainder += fractionToAdd;
            int toAdd = std::floor(mosChangeRemainder);
            mosChangeRemainder -= toAdd;
            mManager.add_mosquito(toAdd);
            if (VERBOSE || (!VERBOSE && time % OUTPUT_INTERVAL == 0))
                std::cout << "Mosquitoes added: " << toAdd << ", leftToAdd: " << numLeftToAdd << ", timeLeftToAdd: " << timeLeftToAdd << ", total: " << mManager.get_count() << "\n";
        }

        if (time >= START_MOS_DCR && time < STOP_MOS_DCR)
        {
            int numLeftToRemove = mManager.get_count() - NUM_MOSQUITOES;
            int timeLeftToRemove = STOP_MOS_DCR - time;
            float fractionToRemove = (float)numLeftToRemove / (float)timeLeftToRemove;
            mosChangeRemainder -= fractionToRemove;
            int toRemove = std::abs(std::ceil(mosChangeRemainder));
            mosChangeRemainder += toRemove;
            mManager.remove_mosquito(toRemove);
            if (VERBOSE || (!VERBOSE && time % OUTPUT_INTERVAL == 0))
                std::cout << "Mosquitoes removed: " << toRemove << ", leftToRemove: " << numLeftToRemove << ", timeLeftToRemove: " << timeLeftToRemove << ", total: " << mManager.get_count() << "\n";
        }
    }
}






void MosquitoManager::initialise(std::vector<Mosquito>* mosquitoesArray)
{
    numMosquitoes = 0;
    mosquitoes = mosquitoesArray;
    activeMosquitoes.reserve(NUM_MOSQUITOES_MAX);
    inactiveMosquitoes.reserve(NUM_MOSQUITOES_MAX);
    for (unsigned int i=0; i<mosquitoes->size(); ++i)
    {
        if (mosquitoes->at(i).is_active())
        {
            activeMosquitoes.push_back(i);
            numMosquitoes++;
        } else
            inactiveMosquitoes.push_back(i);
    }
}

void MosquitoManager::remove_mosquito(unsigned int numToRemove)
{
    for (unsigned int i=0; i<numToRemove; ++i)
    {
        if (!activeMosquitoes.empty())
        {
            unsigned int iM = activeMosquitoes.back();
            mosquitoes->at(iM).kill();
            mosquitoes->at(iM).active = false;
            activeMosquitoes.pop_back();
            inactiveMosquitoes.push_back(iM);
            --numMosquitoes;
        }
    }
}
void MosquitoManager::add_mosquito(unsigned int numToAdd)
{
    for (unsigned int i=0; i<numToAdd; ++i)
    {
        if (inactiveMosquitoes.empty()) //then we need to create a new mosquitoe...
        {
            Mosquito newMosquito;
            newMosquito.kill();
            newMosquito.active = true;
            mosquitoes->push_back(newMosquito);
            activeMosquitoes.push_back(mosquitoes->size()-1);
            numMosquitoes++;
        }
        else //inactive mosquito can be reactivated
        {
            unsigned int iM = inactiveMosquitoes.back();
            mosquitoes->at(iM).kill();
            mosquitoes->at(iM).active = true;
            inactiveMosquitoes.pop_back();
            activeMosquitoes.push_back(iM);
            numMosquitoes++;
        }
    }
}

unsigned int MosquitoManager::random_active_mos() const
{
    return activeMosquitoes[utilities::random(0, activeMosquitoes.size())];
}
