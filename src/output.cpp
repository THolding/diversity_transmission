#include "output.hpp"
#include "model_driver.hpp"
#include "strain.hpp"
#include "utilities.hpp"
#include <cmath>
#include <unordered_map>


void Output::preinitialise_output_storage()
{
    //unsigned int sizeNeeded = (numTimeSteps / outputInterval)+1; //+1 for initial conditions
    unsigned int sizeNeeded = ParamManager::instance().get_int("output_size_needed");

    timeLog.reserve(sizeNeeded);

    hPrevalence.reserve(sizeNeeded);
    hMultipliticyOfInfection.reserve(sizeNeeded);
    hImmunityMean.reserve(sizeNeeded);
    hImmunityVariance.reserve(sizeNeeded);

    mPrevalence.reserve(sizeNeeded);

    circulatingAntigenDiversity.reserve(sizeNeeded);
    shannonEntropy.reserve(sizeNeeded);
    //antigenicHostDiversity.reserve(sizeNeeded);

    eir.reserve(sizeNeeded);

    if (ParamManager::instance().get_bool("output_antigen_frequency"))
        antigenFrequency.reserve(sizeNeeded);

    if (ParamManager::instance().get_bool("dyn_num_mosquitoes"))
        numMosquitoesList.reserve(sizeNeeded);

    if (ParamManager::instance().get_bool("dyn_bite_rate"))
        biteRateList.reserve(sizeNeeded);

    curNumInfectiousBites = 0;
    lastUpdateTime = -1;
}

void Output::append_output(const unsigned int timestep, const Hosts& hosts, const Mosquitoes& mosquitoes)
{
    timeLog.push_back(timestep);

    float hPrev, hMOI;
    calc_host_infection_metrics(hosts, hPrev, hMOI);
    hPrevalence.push_back(hPrev);
    hMultipliticyOfInfection.push_back(hMOI);

    float immunityMean, immunityVariance;
    calc_host_immunity_metrics(hosts, immunityMean, immunityVariance);
    hImmunityMean.push_back(immunityMean);
    hImmunityVariance.push_back(immunityVariance);

    float mPrev;
    calc_mosquito_infection_metrics(mosquitoes, mPrev);
    mPrevalence.push_back(mPrev);

    float curShannonEntropy;
    circulatingAntigenDiversity.push_back(calc_circulating_antigen_diversity(hosts, mosquitoes, curShannonEntropy));
    shannonEntropy.push_back(curShannonEntropy);
    //antigenicHostDiversity.push_back(calc_antigen_host_diversity(hosts));

    eir.push_back(calc_eir(timestep)); //Used callback counter and self-resets.

    log_antigen_frequency(hosts, mosquitoes);

    log_dyn_params();

    lastUpdateTime = timestep;
}

void Output::export_output(const std::string runName, const std::string filePath)
{
    utilities::arrayToFile(timeLog, filePath+runName+"_timesteps.csv");

    utilities::arrayToFile(hPrevalence, filePath+runName+"_host_prevalences.csv");
    utilities::arrayToFile(hMultipliticyOfInfection, filePath+runName+"_host_moi.csv");
    utilities::arrayToFile(hImmunityMean, filePath+runName+"_host_immunity_mean.csv");
    utilities::arrayToFile(hImmunityVariance, filePath+runName+"_host_immunity_variance.csv");

    utilities::arrayToFile(mPrevalence, filePath+runName+"_mosquito_prevalences.csv");

    utilities::arrayToFile(circulatingAntigenDiversity, filePath+runName+"_circulating_antigen_diversity.csv");
    utilities::arrayToFile(shannonEntropy, filePath+runName+"_shannon_entropy_diversity.csv");

    //utilities::arrayToFile(hostDiversity, filePath+runName+"_host_diversity.csv");
    //utilities::arrayToFile(antigenicHostDiversity, filePath+runName+"_antigenic_host_diversity.csv");

    utilities::arrayToFile(eir, filePath+runName+"_eir.csv");

    if (ParamManager::instance().get_bool("output_antigen_frequency"))
        utilities::matrixToFile(antigenFrequency, filePath+runName+"_circulating_antigen_frequency.csv", ", ");

    if (ParamManager::instance().get_bool("dyn_num_mosquitoes"))
        utilities::arrayToFile(numMosquitoesList, filePath+runName+"_num_mosquitoes.csv");

    if (ParamManager::instance().get_bool("dyn_bite_rate"))
        utilities::arrayToFile(biteRateList, filePath+runName+"_bite_rate.csv");
}

void Output::register_infectious_bite()
{
    ++curNumInfectiousBites;
}

//
void Output::calc_host_infection_metrics(const Hosts& hosts, float& prevalence, float& multiplicityOfInfection)
{
    prevalence = 0.0f;
    multiplicityOfInfection = 0.0f;

    for (const Host& host : hosts)
    {
        if (host.is_infected())
        {
            ++prevalence;
            multiplicityOfInfection += host.infection1.infected;
            multiplicityOfInfection += host.infection2.infected;
        }
    }

    prevalence = prevalence / ParamManager::instance().get_int("num_hosts");
    multiplicityOfInfection = multiplicityOfInfection / ParamManager::instance().get_int("num_hosts");

    std::cout << "Host prevalence: " << prevalence << std::endl;
}

void Output::calc_mosquito_infection_metrics(const Mosquitoes& mosquitoes, float& prevalence)
{
    prevalence=0.0f;

    //Loop through mosquitoes
    for (const Mosquito& mosquito : mosquitoes)
    {
        //Count number of mosquitoes that are infected
        if (mosquito.is_active())
            prevalence += mosquito.is_infected();
    }
    prevalence = prevalence / model->get_mos_manager()->get_count();
}

void Output::calc_host_immunity_metrics(const Hosts& hosts, float& immunityMean, float& immunityVariance)
{
    immunityMean = 0.0f;
    immunityVariance = 0.0f;

    std::vector<float> individualTotals(ParamManager::instance().get_int("num_hosts"));

    //Calculate mean
    for (unsigned int iH=0; iH<hosts.size(); ++iH)
    {
        float curTotalImmunity = 0;
        for (auto immunity : hosts[iH].immuneState) //Sum immunity
        {
            curTotalImmunity += immunity;
        }
        curTotalImmunity = curTotalImmunity / hosts[iH].immuneState.size(); // Total immunity
        immunityMean += curTotalImmunity;
        individualTotals[iH] = curTotalImmunity; //Used to calculate variance
    }
    immunityMean = immunityMean / hosts.size();

    //Calculate variance
    for (float val : individualTotals)
        immunityVariance += std::pow(immunityMean - val, 2);
    immunityVariance = immunityVariance / ParamManager::instance().get_int("num_hosts");
}

/*float Output::calc_host_diversity(const Hosts& hosts)
{
    std::unordered_map<std::string, unsigned int> diversityPool;

    for (const Host& host : hosts)
    {
        if (host.infection1.infected)
        {
            std::string strainStr = strain_phenotype_str(host.infection1.strain);
            if (diversityPool.count(strainStr) == 0)
                diversityPool[strainStr] = 1;
        }
        if (host.infection2.infected)
        {
            std::string strainStr = strain_phenotype_str(host.infection2.strain);
            if (diversityPool.count(strainStr) == 0)
                diversityPool[strainStr] = 1;
        }
    }

    return diversityPool.size();
}*/

/*float Output::calc_antigen_host_diversity(const Hosts& hosts)
{
    std::unordered_map<Antigen, unsigned int> diversityPool;

    for (const Host& host : hosts)
    {
        if (host.infection1.infected)
        {
            for (const Antigen antigen : host.infection1.strain)
            {
                if (diversityPool.count(get_phenotype_id(antigen)) == 0)
                {
                    diversityPool[get_phenotype_id(antigen)] = 1;
                }
            }
        }
        if (host.infection2.infected)
        {
            for (const Antigen antigen : host.infection2.strain)
            {
                if (diversityPool.count(get_phenotype_id(antigen)) == 0)
                {
                    diversityPool[get_phenotype_id(antigen)] = 1;
                }
            }
        }
    }

    float antigenDiversity = (float)diversityPool.size() / (float)ParamManager::instance().get_int("num_phenotypes");
    return antigenDiversity;
}*/

float Output::calc_circulating_antigen_diversity(const Hosts& hosts, const Mosquitoes& mosquitoes, float& shannonEntropy)
{
    std::unordered_map<Antigen, unsigned int> diversityPool;

    for (const Host& host : hosts)
    {
        if (host.infection1.infected)
        {
            for (const Antigen antigen : host.infection1.strain)
            {
                if (diversityPool.count(get_phenotype_id(antigen)) == 0)
                {
                    diversityPool[get_phenotype_id(antigen)] = 1;
                }
                else
                    diversityPool[get_phenotype_id(antigen)] = diversityPool[get_phenotype_id(antigen)]+1;
            }
        }
        if (host.infection2.infected)
        {
            for (const Antigen antigen : host.infection2.strain)
            {
                if (diversityPool.count(get_phenotype_id(antigen)) == 0)
                {
                    diversityPool[get_phenotype_id(antigen)] = 1;
                }
                else
                    diversityPool[get_phenotype_id(antigen)] = diversityPool[get_phenotype_id(antigen)]+1;
            }
        }
    }

    for (const Mosquito& mosquito : mosquitoes)
    {
        if (mosquito.infection.infected)
        {
            for (const Antigen antigen : mosquito.infection.strain)
            {
                if (diversityPool.count(get_phenotype_id(antigen)) == 0)
                    diversityPool[get_phenotype_id(antigen)] = 1;
                else
                    diversityPool[get_phenotype_id(antigen)] = diversityPool[get_phenotype_id(antigen)]+1;
            }
        }
    }

    float antigenDiversity = (float)diversityPool.size() / (float)ParamManager::instance().get_int("num_phenotypes");

    shannonEntropy = calc_shannon_entropy(diversityPool);
    return antigenDiversity;
}

float Output::calc_shannon_entropy(const std::unordered_map<Antigen, unsigned int>& diversityPool)
{
    std::vector<float> proportions;
    float total = 0.0;
    for (const auto& item : diversityPool)
    {
        proportions.push_back(item.second);
        total += item.second;
    }
    if (total == 0.0)
        return 0.0;
    for (unsigned int i=0; i<proportions.size(); ++i)
        proportions[i] = proportions[i] / total;

    float entropy = 0.0;
    for (float proportion : proportions)
        entropy -= proportion * std::log(proportion);

    return entropy;
}

//returns daily EIR.
float Output::calc_eir(const unsigned int currentTime)
{
    //Calculates this because output interval can change over the course of a simulation.
    unsigned int timeSinceLastUpdate = currentTime - lastUpdateTime;

    float eir = (float)curNumInfectiousBites / (float)ParamManager::instance().get_int("num_hosts");
    eir = eir / (float)timeSinceLastUpdate;

    //Reset any counters.
    curNumInfectiousBites = 0;

    return eir;
}

void Output::antigen_counter_helper(std::vector<unsigned int>& antigenFreqs, const Strain& strain)
{
    for (const Antigen& antigen : strain)
        antigenFreqs[get_phenotype_id(antigen)] += 1;
}

void Output::log_antigen_frequency(const Hosts& hosts, const Mosquitoes& mosquitoes)
{
    if (ParamManager::instance().get_bool("output_antigen_frequency") == false)
        return;

    std::vector<unsigned int> antigenFreqs(ParamManager::instance().get_int("num_phenotypes"));

    //All host infections
    for (const Host& host : hosts)
    {
        if (host.infection1.infected)
            antigen_counter_helper(antigenFreqs, host.infection1.strain);
        if (host.infection2.infected)
            antigen_counter_helper(antigenFreqs, host.infection2.strain);
    }

    //All mosquito infections
    for (const Mosquito& mosquito : mosquitoes)
    {
        if (mosquito.is_active() && mosquito.infection.infected)
            antigen_counter_helper(antigenFreqs, mosquito.infection.strain);
    }

    antigenFrequency.push_back(antigenFreqs);
}

void Output::log_dyn_params()
{
    if (ParamManager::instance().get_bool("dyn_num_mosquitoes"))
        numMosquitoesList.push_back(model->get_mos_manager()->get_count());

    if (ParamManager::instance().get_bool("dyn_bite_rate"))
        biteRateList.push_back(ParamManager::instance().get_float("bite_rate"));
}
