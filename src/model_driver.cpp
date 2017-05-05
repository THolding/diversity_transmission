#include "model_driver.hpp"
#include "demographic_tools.hpp"
#include "strain.hpp"
#include "utilities.hpp"
#include "param_manager.hpp"
#include <algorithm>
#include <cmath>
#include <vector>
#include <omp.h>

//toremove:
#include "diversity_monitor.hpp"
#include "testing.hpp"

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
    hosts.reserve(ParamManager::num_hosts);
    for (unsigned int h=0; h<ParamManager::num_hosts; ++h)
    {
        Host host;
        host.kill();
        host.age = random_host_equilibrum_age(cdfHosts);
        hosts.push_back(host);
    }

    //Initialise mosquitoes
    std::cout << "initialising mosquito demographics" << std::endl;
    mosquitoes.reserve(ParamManager::max_num_mosquitoes);
    for (unsigned int m=0; m<ParamManager::initial_num_mosquitoes; ++m)
    {
        Mosquito mosquito;
        mosquito.kill();
        mosquito.age = random_mosquito_equilibrium_age(cdfMosquitoes);
        mosquito.active = true;
        mosquitoes.push_back(mosquito);
    }

    mManager.initialise(&mosquitoes);

    //Create initial pool of strains
    std::vector<Strain> initialStrainPool;
    if (ParamManager::unique_initial_strains == true)
        create_unique_initial_strains(initialStrainPool);
    else //randomly select initial strains and infections (default).
        create_random_initial_strains(initialStrainPool);

    //Initial infections (mosquitoes)
    std::cout << "initialising mosquito infections" << std::endl;
    for (unsigned int i=0; i<ParamManager::initial_num_mosquito_infections; ++i)
    {
        unsigned int iM = mManager.random_active_mos();
        mosquitoes[iM].infect(initialStrainPool[i % initialStrainPool.size()], false);
    }
}

void ModelDriver::run_model()
{
    initialise_model();
    output.preinitialise_output_storage();

    bool finished = false;
    unsigned int timeElapsed = 0;
    unsigned int timeNextOutput = ParamManager::output_interval;
    unsigned int lastOutputInterval = ParamManager::output_interval;
    output.append_output(timeElapsed, hosts, mosquitoes);
    burnInPeriod = ParamManager::burn_in_period;
    while (!finished)
    {
        //Dynamic parameters
        ParamManager::update_adaptors(timeElapsed);
        if (lastOutputInterval != ParamManager::output_interval) {
            lastOutputInterval = ParamManager::output_interval;
            timeNextOutput = timeElapsed;
        }

        //Update output times
        if (ParamManager::verbose) {
            std::cout << "t=" << timeElapsed << "\n";
            if (burnInPeriod > 0)
                std::cout << "burnIn left: " << burnInPeriod << "\n";
        }
        else if (timeElapsed == timeNextOutput) {
            std::cout << "t=" << timeElapsed << ". ";
            if (burnInPeriod > 0)
                std::cout << "burnIn left: " << burnInPeriod << ". ";
        }

        //tmp debug
        //std::cout << timeElapsed << ": UniqueAntigens: " << DiversityMonitor::get_num_unique_antigens() << "\n";

        //Host demographics
        //std::cout << "aging hosts...\n";
        #pragma omp barrier
        age_hosts();

        //Mosquito demographics
        //std::cout << "aging mosquitoes...\n";
        #pragma omp barrier
        age_mosquitoes();

        //Update infections in hosts
        //std::cout << "updating host infections...\n";
        #pragma omp barrier
        update_host_infections();

        //Update infections in mosquitoes
        //std::cout << "updating mosquito infections...\n";
        #pragma omp barrier
        update_mosquito_infections();

        //mosquitoes feed
        //std::cout << "feeding mosquitoes...\n";
        #pragma omp barrier
        feed_mosquitoes();

        //if appropriate, reintroduce an extinct initial strain.
        #pragma omp barrier
        attempt_reintroduction(timeElapsed);

        //Update logging / data collection.
        if (timeElapsed == timeNextOutput)
        {
            timeNextOutput = timeElapsed + ParamManager::output_interval;
            output.append_output(timeElapsed, hosts, mosquitoes);
        }

        //Update time and check stop condition.
        //std::cout << "incrementing time...\n";
        if (burnInPeriod > 0)
            --burnInPeriod;

        ++timeElapsed;
        if (timeElapsed > ParamManager::run_time)
            finished = true;

        #pragma omp barrier
    }

    output.export_output();
}

void ModelDriver::age_hosts()
{
    #pragma omp parallel for
    for (unsigned int i=0; i<hosts.size(); ++i)
    {
        hosts[i].age_host(pDeathHosts);
    }
}

void ModelDriver::age_mosquitoes()
{
    #pragma omp parallel for
    for (unsigned int i=0; i<mosquitoes.size(); ++i)
    {
        if (mosquitoes[i].is_active())
            mosquitoes[i].age_mosquito(pDeathMosquitoes);
    }
}

void ModelDriver::update_host_infections()
{
    #pragma omp parallel for
    for (unsigned int i=0; i<hosts.size(); ++i)
    {
        hosts[i].update_infections();
    }
}

void ModelDriver::update_mosquito_infections()
{
    #pragma omp parallel for
    for (unsigned int i=0; i<mosquitoes.size(); ++i)
    {
        if (mosquitoes[i].is_active())
            mosquitoes[i].update_infection();
    }
}

void ModelDriver::feed_mosquitoes()
{
    bool allowRecombination;
    if (burnInPeriod <= 0)
        allowRecombination = true;
    else
        allowRecombination = false;

    #pragma omp parallel for
    for (unsigned int i=0; i<mosquitoes.size(); ++i)
    {
        if (mosquitoes[i].is_active())
        {
            float p = utilities::random_float01()*ParamManager::get_cumulative_bite_frequency_distribution().back();
            unsigned int numBites = 0;
            while (p > ParamManager::get_cumulative_bite_frequency_distribution()[numBites])
                ++numBites;

            for (unsigned int b=0; b<numBites; ++b)
            {
                unsigned int iH = utilities::urandom(0, hosts.size());
                mosquitoes[i].feed(hosts[iH], &output, allowRecombination);
            }
        }
    }
}

//Attempts to reintroduce a strain IF and only if it is time to do so
//unique_initial_strains == true, static diversity is in use (i.e. intra and intergenic recombination = 0.0)
void ModelDriver::attempt_reintroduction(const unsigned int time)
{
    if (ParamManager::reintroduction_interval != 0 && time % ParamManager::reintroduction_interval == 0)
    {
        std::cout << "Attempting reintroduction at t=" << time << "\n";
        if (ParamManager::unique_initial_strains && ParamManager::intragenic_recombination_p == 0.0f && ParamManager::intergenic_recombination_p == 0.0f)
        {
            //Choose a random initial strain and if it is extinct try to infect a random mosquito (only works if mosquito is uninfected).
            unsigned int iS = utilities::random(0, cachedInitialStrainPool.size());
            Strain strain = cachedInitialStrainPool[iS];

            //If not extinct then ignore this...
            if (DiversityMonitor::get_antigen_count(get_phenotype_id(strain[0])) != 0)
                return; //Don't reintroduce strains which aren't extinct!

            unsigned int iM = mManager.random_active_mos();
            if (mosquitoes[iM].is_infected() == false) {
                mosquitoes[iM].infect(cachedInitialStrainPool[iS], false, true);
                std::cout << "Reintroduction successful!\n";
            }
        }
        else
            throw std::runtime_error("CANNOT REINTRODUCE INITIAL STRAINS: unique_initial_strains may be false, or intra/intergenic recombination may be non-zero.");
    }
}

//Guarentees that each strain generated is unique and non-overlapping, provided there are still unused antigens left during strain creation.
//If not enough antigens are available they will be reused such that
void ModelDriver::create_unique_initial_strains(std::vector<Strain>& _initialStrainPool)
{
    //Produce a random of all possible antigens
    std::cout << "initialising antigen pool" << std::endl;
    std::vector<Antigen> allAntigens;
    for (unsigned int i=0; i<ParamManager::num_phenotypes; ++i)
        allAntigens.push_back(i << ParamManager::num_genotype_only_bits);
    std::random_shuffle(allAntigens.begin(), allAntigens.end());

    //Create (unique) strain list
    std::cout << "initialising strain pool" << std::endl;
    _initialStrainPool.reserve(ParamManager::initial_num_strains);
    cachedInitialStrainPool.clear();
    cachedInitialStrainPool.reserve(ParamManager::initial_num_strains);
    unsigned int currentIndex=0;
    for (unsigned int s=0; s<ParamManager::initial_num_strains; ++s)
    {
        Strain newStrain;
        for (unsigned int i=0; i<ParamManager::repertoire_size; ++i)
        {
            if (currentIndex >= ParamManager::initial_antigen_diversity)
            {
                std::random_shuffle(allAntigens.begin(), allAntigens.begin()+ParamManager::initial_antigen_diversity);
                currentIndex = 0;
            }

            newStrain.push_back(allAntigens[currentIndex++]);
        }
        _initialStrainPool.push_back(newStrain);
        cachedInitialStrainPool.push_back(newStrain);
    }
}

//Generates strains by randomly sampling antigens to form a antigen pool, then randomly sampling the antigen pool to create each strain.
void ModelDriver::create_random_initial_strains(std::vector<Strain>& _initialStrainPool)
{
    cachedInitialStrainPool.clear();

    //Generate initial pool of antigen diversity
    std::cout << "initialising antigen pool" << std::endl;
    std::vector<Antigen> initialAntigenPool;
    initialAntigenPool.reserve(ParamManager::initial_antigen_diversity);
    for (unsigned int a=0; a<ParamManager::initial_antigen_diversity; ++a)
        initialAntigenPool.push_back(random_antigen() << ParamManager::num_genotype_only_bits); //turn antigen into gene

    //Generate initial pool of strains from available initial antigen diversity
    std::cout << "initialising strain pool" << std::endl;
    _initialStrainPool.reserve(ParamManager::initial_num_strains);
    for (unsigned int s=0; s<ParamManager::initial_num_strains; ++s)
        _initialStrainPool.push_back(strain_from_antigen_pool(initialAntigenPool));
}



//temp todelete
void ModelDriver::test()
{
    ParamManager::num_hosts = 10000;
    ParamManager::initial_num_mosquitoes = 10000;
    ParamManager::initial_num_mosquito_infections = 0;//10000;
    ParamManager::unique_initial_strains = true;
    ParamManager::intragenic_recombination_p = 0.0f;
    ParamManager::intergenic_recombination_p = 0.0f;
    ParamManager::num_phenotypes = 10000;
    ParamManager::initial_antigen_diversity = 10000;
    ParamManager::initial_num_strains = 1000;
    ParamManager::recalculate_derived_parameters();
    initialise_model();

    unsigned int uniqueCount;
    unsigned int totalCount;
    testing::long_diversity_count(uniqueCount, totalCount, hosts, mosquitoes);
    std::cout << "Mosquitoes infected:\n";
    std::cout << "Long method: " << uniqueCount << "\t" << totalCount << "\n";
    std::cout << "Quik method: " << DiversityMonitor::get_num_unique_antigens() << "\t" << DiversityMonitor::get_total_antigens() << "\n";

    for (unsigned int i=0; i<hosts.size(); ++i)
        hosts[i].infect(cachedInitialStrainPool[utilities::random(0, cachedInitialStrainPool.size())]);
    testing::long_diversity_count(uniqueCount, totalCount, hosts, mosquitoes);
    std::cout << "Mosquitoes+hosts infected:\n";
    std::cout << "Long method: " << uniqueCount << "\t" << totalCount << "\n";
    std::cout << "Quik method: " << DiversityMonitor::get_num_unique_antigens() << "\t" << DiversityMonitor::get_total_antigens() << "\n";


    for (unsigned int i=0; i<10000; ++i)
    {
        age_hosts();
        //age_mosquitoes();
    }
    testing::long_diversity_count(uniqueCount, totalCount, hosts, mosquitoes);
    std::cout << "After loop host+mosquito:\n";
    std::cout << "Long method: " << uniqueCount << "\t" << totalCount << "\n";
    std::cout << "Quik method: " << DiversityMonitor::get_num_unique_antigens() << "\t" << DiversityMonitor::get_total_antigens() << "\n";

    //age_hosts();
    //age_mosquitoes();
    //update_host_infections();
    //update_mosquito_infections();
    //feed_mosquitoes();


    //std::cout << "Counted unique antigens: " << uniqueAntigenCount << "\tIncremented antigens: " << DiversityMonitor::get_num_unique_antigens() << "\n";
    //std::cout << "Counted total antigens:  " << antigenTotal << "\tIncremented antigens: " << DiversityMonitor::get_total_antigens() << "\n";
}
