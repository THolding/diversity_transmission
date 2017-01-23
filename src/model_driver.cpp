#include "model_driver.hpp"
#include "demographic_tools.hpp"
#include "strain.hpp"
#include "utilities.hpp"
#include "param_manager.hpp"
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

    //Initialise hosts
    std::cout << "initialising host demographics" << std::endl;
    hosts.reserve(ParamManager::instance().get_int("num_hosts"));
    for (unsigned int h=0; h<ParamManager::instance().get_int("num_hosts"); ++h)
    {
        Host host;
        host.kill();
        host.age = random_host_equilibrum_age(cdfHosts);
        hosts.push_back(host);
    }

    //Initialise mosquitoes
    std::cout << "initialising mosquito demographics" << std::endl;
    mosquitoes.reserve(ParamManager::instance().get_int("max_num_mosquitoes"));
    for (unsigned int m=0; m<ParamManager::instance().get_int("initial_num_mosquitoes"); ++m)
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
    initialAntigenPool.reserve(ParamManager::instance().get_int("initial_antigen_diversity"));
    for (unsigned int a=0; a<ParamManager::instance().get_int("initial_antigen_diversity"); ++a)
        initialAntigenPool.push_back(random_antigen());

    //Generate initial pool of strains from available initial antigen diversity
    std::cout << "initialising strain pool" << std::endl;
    std::vector<Strain> initialStrainPool;
    initialStrainPool.reserve(ParamManager::instance().get_int("initial_num_strains"));
    for (unsigned int s=0; s<ParamManager::instance().get_int("initial_num_strains"); ++s)
        initialStrainPool.push_back(strain_from_antigen_pool(initialAntigenPool));

    //Initial infections (mosquitoes)
    std::cout << "initialising mosquito infections" << std::endl;
    for (unsigned int i=0; i<ParamManager::instance().get_int("initial_num_mosquito_infections"); ++i)
    {
        unsigned int iS = utilities::urandom(0, initialStrainPool.size());
        unsigned int iM = mManager.random_active_mos();
        mosquitoes[iM].infect(initialStrainPool[iS], false);
    }
}

void ModelDriver::run_model()
{
    initialise_model();
    output.preinitialise_output_storage();

    bool finished = false;
    unsigned int timeElapsed = 0;
    unsigned int timeNextOutput = ParamManager::instance().get_int("output_interval");
    unsigned int lastOutputInterval = ParamManager::instance().get_int("output_interval");
    output.append_output(timeElapsed, hosts, mosquitoes);
    burnInPeriod = ParamManager::instance().get_int("burn_in_period");
    while (!finished)
    {
        //Dynamic parameters
        ParamManager::instance().update_adaptors(timeElapsed);
        if (lastOutputInterval != ParamManager::instance().get_int("output_interval")) {
            lastOutputInterval = ParamManager::instance().get_int("output_interval");
            timeNextOutput = timeElapsed;
        }

        //Update output times
        if (ParamManager::instance().get_bool("verbose")) {
            std::cout << "t=" << timeElapsed << "\n";
            if (burnInPeriod > 0)
                std::cout << "burnIn left: " << burnInPeriod << "\n";
        }
        else if (timeElapsed == timeNextOutput) {
            std::cout << "t=" << timeElapsed << ". ";
            if (burnInPeriod > 0)
                std::cout << "burnIn left: " << burnInPeriod << ". ";
        }

        //Host demographics
        //std::cout << "aging hosts...\n";
        age_hosts();
        //

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
        if (timeElapsed == timeNextOutput)
        {
            timeNextOutput = timeElapsed + ParamManager::instance().get_int("output_interval");
            output.append_output(timeElapsed, hosts, mosquitoes);
        }

        //Update time and check stop condition.
        //std::cout << "incrementing time...\n";
        if (burnInPeriod > 0)
            --burnInPeriod;

        ++timeElapsed;
        if (timeElapsed > ParamManager::instance().get_int("run_time"))
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
        allowRecombination = true;
    else
        allowRecombination = false;

    for (Mosquito& mosquito : mosquitoes)
    {
        if (mosquito.is_active()) {
            //Calculate number of times the mosquito feeds.
            float p = utilities::random_float01()*ParamManager::instance().get_cumulative_bite_frequency_distribution().back();
            unsigned int numBites = 0;
            while (p > ParamManager::instance().get_cumulative_bite_frequency_distribution()[numBites])
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
