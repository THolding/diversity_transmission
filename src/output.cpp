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

    hostDiversity.reserve(sizeNeeded);
    antigenicHostDiversity.reserve(sizeNeeded);

    eir.reserve(sizeNeeded);

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

    hostDiversity.push_back(calc_host_diversity(hosts));
    antigenicHostDiversity.push_back(calc_antigen_host_diversity(hosts));

    eir.push_back(calc_eir(timestep)); //Used callback counter and self-resets.

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

    utilities::arrayToFile(hostDiversity, filePath+runName+"_host_diversity.csv");
    utilities::arrayToFile(antigenicHostDiversity, filePath+runName+"_antigenic_host_diversity.csv");

    utilities::arrayToFile(eir, filePath+runName+"_eir.csv");

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

float Output::calc_host_diversity(const Hosts& hosts)
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
}

float Output::calc_antigen_host_diversity(const Hosts& hosts)
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
}

//returns daily EIR.
float Output::calc_eir(const unsigned int currentTime)
{
    unsigned int timeSinceLastUpdate = currentTime - lastUpdateTime;
    ////////THIS MIGHT BE BROKEN! Must cache previous call time and calculate EIR accordingly.
    float eir = (float)curNumInfectiousBites / (float)ParamManager::instance().get_int("num_hosts");
    eir = eir / (float)timeSinceLastUpdate;

    //Reset any counters.
    curNumInfectiousBites = 0;

    return eir;
}

void Output::log_dyn_params()
{
    if (ParamManager::instance().get_bool("dyn_num_mosquitoes"))
        numMosquitoesList.push_back(model->get_mos_manager()->get_count());

    if (ParamManager::instance().get_bool("dyn_bite_rate"))
        biteRateList.push_back(ParamManager::instance().get_float("bite_rate"));
}
