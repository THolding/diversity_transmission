#include "output.hpp"
#include "model_driver.hpp"
#include "strain.hpp"
#include "utilities.hpp"
#include <cmath>
#include <unordered_map>


void Output::preinitialise_output_storage()
{
    curNumInfectiousBites = 0;
    lastUpdateTime = -1;

    //unsigned int sizeNeeded = (numTimeSteps / outputInterval)+1; //+1 for initial conditions
    unsigned int sizeNeeded = ParamManager::instance().get_int("output_size_needed");


    timeLog.reserve(sizeNeeded);
    hPrevalence.reserve(sizeNeeded);
    mPrevalence.reserve(sizeNeeded);
    moi.reserve(sizeNeeded);
    eir.reserve(sizeNeeded);
    proportionCirculatingAntigens.reserve(sizeNeeded);
    shannonEntropy.reserve(sizeNeeded);
    absoluteImmunity.reserve(sizeNeeded);

    if (ParamManager::instance().get_bool("output_antigen_frequency"))
        antigenFrequency.reserve(sizeNeeded);

    if (ParamManager::instance().get_bool("output_parasite_adaptedness"))
        parasiteAdaptedness.reserve(sizeNeeded);

    if (ParamManager::instance().get_bool("dyn_num_mosquitoes"))
        numMosquitoesList.reserve(sizeNeeded);

    if (ParamManager::instance().get_bool("dyn_bite_rate"))
        biteRateList.reserve(sizeNeeded);

    if (ParamManager::instance().get_bool("dyn_intragenic_recombination_p"))
        intragenicRecombinationPList.reserve(sizeNeeded);
}

void Output::append_output(const unsigned int timestep, const Hosts& hosts, const Mosquitoes& mosquitoes)
{
    calc_host_dependent_metrics(hosts);
    calc_mosquito_dependent_metrics(mosquitoes);
    calc_host_mosquito_dependent_metrics(hosts, mosquitoes);
    calc_time_dependent_metrics(timestep);
    calc_dyn_metrics();

    lastUpdateTime = timestep;
}

void Output::export_output(const std::string runName, const std::string filePath)
{
    utilities::arrayToFile(timeLog, filePath+runName+"_timesteps.csv");
    utilities::arrayToFile(hPrevalence, filePath+runName+"_host_prevalences.csv");
    utilities::arrayToFile(mPrevalence, filePath+runName+"_mosquito_prevalences.csv");
    utilities::arrayToFile(moi, filePath+runName+"_moi.csv");
    utilities::arrayToFile(eir, filePath+runName+"_eir.csv");
    utilities::arrayToFile(proportionCirculatingAntigens, filePath+runName+"_num_circulating_antigens.csv");
    utilities::arrayToFile(shannonEntropy, filePath+runName+"_shannon_entropy_diversity.csv");
    utilities::arrayToFile(absoluteImmunity, filePath+runName+"_absolute_immunity.csv");

    if (ParamManager::instance().get_bool("output_antigen_frequency"))
        utilities::matrixToFile(antigenFrequency, filePath+runName+"_circulating_antigen_frequency.csv", ", ");

    if (ParamManager::instance().get_bool("output_parasite_adaptedness"))
        utilities::arrayToFile(parasiteAdaptedness, filePath+runName+"_parasite_adaptedness.csv");

    if (ParamManager::instance().get_bool("dyn_num_mosquitoes"))
        utilities::arrayToFile(numMosquitoesList, filePath+runName+"_num_mosquitoes.csv");

    if (ParamManager::instance().get_bool("dyn_bite_rate"))
        utilities::arrayToFile(biteRateList, filePath+runName+"_bite_rate.csv");

    if (ParamManager::instance().get_bool("dyn_intragenic_recombination_p"))
        utilities::arrayToFile(intragenicRecombinationPList, filePath+runName+"_intragenic_recombination_p.csv");
}

void Output::register_infectious_bite()
{
    ++curNumInfectiousBites;
}

//responsible for: host prevalence, host immunity, moi
void Output::calc_host_dependent_metrics(const Hosts& hosts)
{
    //Calculate prevalence and multiplicity of infection
    float prevalence = 0.0f;
    float multiplicityOfInfection = 0.0f;
    float absImmunity = 0.0f;

    for (unsigned int i=0; i<hosts.size(); ++i)
    {
        if (hosts[i].is_infected())
        {
            //Update counts for prevalence and moi
            ++prevalence;
            multiplicityOfInfection += hosts[i].infection1.infected;
            multiplicityOfInfection += hosts[i].infection2.infected;
        }

        //Calculate absolute immunity
        float curTotalImmunity = 0;
        for (const float immunity : hosts[i].immuneState) //Sum immunity
            curTotalImmunity += immunity;
        curTotalImmunity = curTotalImmunity / hosts[i].immuneState.size(); // Total immunity
        absImmunity += curTotalImmunity;
    }

    prevalence = prevalence / ParamManager::instance().get_int("num_hosts");
    multiplicityOfInfection = multiplicityOfInfection / ParamManager::instance().get_int("num_hosts");

    absImmunity = absImmunity / ParamManager::instance().get_int("num_hosts");
    hPrevalence.push_back(prevalence);

    moi.push_back(multiplicityOfInfection);
    absoluteImmunity.push_back(absImmunity);

    std::cout << "Host prevalence: " << prevalence << "\tabsImmunity: " << absImmunity << std::endl;
}

//mosquito prevalence
void Output::calc_mosquito_dependent_metrics(const Mosquitoes& mosquitoes)
{
    float prevalence=0.0f; //Mosquito prevalence

    //Loop through mosquitoes
    for (const Mosquito& mosquito : mosquitoes)
    {
        //Count number of mosquitoes that are infected
        if (mosquito.is_active())
            prevalence += mosquito.is_infected();
    }
    prevalence = prevalence / model->get_mos_manager()->get_count();
    mPrevalence.push_back(prevalence);
}


//time, daily EIR
void Output::calc_time_dependent_metrics(const unsigned int currentTime)
{
    timeLog.push_back(currentTime);

    //Calculate eir
    //Tracks time since last update because output interval can change over the course of a simulation.
    unsigned int timeSinceLastUpdate = currentTime - lastUpdateTime;
    float curEir = (float)curNumInfectiousBites / (float)ParamManager::instance().get_int("num_hosts");
    curEir = curEir / (float)timeSinceLastUpdate;
    //Reset infectious bite counteer
    curNumInfectiousBites = 0;
    eir.push_back(curEir);
}

//Track any dynamic (time dependent) parameters over time.
void Output::calc_dyn_metrics()
{
    if (ParamManager::instance().get_bool("dyn_num_mosquitoes"))
        numMosquitoesList.push_back(model->get_mos_manager()->get_count());

    if (ParamManager::instance().get_bool("dyn_bite_rate"))
        biteRateList.push_back(ParamManager::instance().get_float("bite_rate"));

    if (ParamManager::instance().get_bool("dyn_intragenic_recombination_p"))
        intragenicRecombinationPList.push_back(ParamManager::instance().get_float("intragenic_recombination_p"));
}

//numCirculatingAntigens, shannon entropy, antigen frequency, parasite adaptedness
//Optional: antigen frequency, parasite adaptedness
void Output::calc_host_mosquito_dependent_metrics(const Hosts& hosts, const Mosquitoes& mosquitoes)
{
    std::vector<unsigned int> curAntigenFrequencies;
    curAntigenFrequencies = std::vector<unsigned int>(ParamManager::instance().get_int("num_phenotypes"));

    std::unordered_map<Antigen, unsigned int> antigenCounter;

    //All host infections
    unsigned int antigenTotal = 0;
    for (const Host& host : hosts)
    {
        if (host.infection1.infected) {
            count_individual_antigens(curAntigenFrequencies, antigenCounter, host.infection1.strain);
            antigenTotal += ParamManager::instance().get_int("repertoire_size");
        }
        if (host.infection2.infected) {
            count_individual_antigens(curAntigenFrequencies, antigenCounter, host.infection2.strain);
            antigenTotal += ParamManager::instance().get_int("repertoire_size");
        }
    }


    //All mosquito infections
    for (const Mosquito& mosquito : mosquitoes)
    {
        if (mosquito.is_active() && mosquito.infection.infected)
            count_individual_antigens(curAntigenFrequencies, antigenCounter, mosquito.infection.strain);
    }

    //Calculate antigen frequencies by normalising by total
    for (unsigned int i=0; i<curAntigenFrequencies.size(); ++i)
        curAntigenFrequencies[i] = (float) curAntigenFrequencies[i] / (float) antigenTotal;
    //if (ParamManager::instance().get_bool("output_antigen_frequency")) //Turns out we need to do this for shannon entropy anyway!
    antigenFrequency.push_back(curAntigenFrequencies);

    proportionCirculatingAntigens.push_back(((float)antigenCounter.size()) / ((float)ParamManager::instance().get_int("num_phenotypes")));
    shannonEntropy.push_back(calc_shannon_entropy(curAntigenFrequencies));

    //Calculate parasite adaptability
    if (ParamManager::instance().get_bool("output_parasite_adaptedness"))// && ParamManager::instance().get_bool("output_antigen_frequency"))
        parasiteAdaptedness.push_back(calc_parasite_adaptedness(curAntigenFrequencies, antigenTotal, hosts));
}

float Output::calc_parasite_adaptedness(const std::vector<unsigned int> curAntigenFrequencies, const unsigned int antigenTotal, const Hosts& hosts)
{
    //Calculate susceptibility
    std::vector<float> susceptibility;
    for (unsigned int a=0; a<ParamManager::instance().get_int("num_phenotypes"); ++a)
    {
        susceptibility.push_back(0.0);
        for (unsigned int h=0; h<hosts.size(); ++h)
            susceptibility[a] = susceptibility[a] + hosts[h].immuneState[a];
        susceptibility[a] = susceptibility[a] / (float)hosts.size();
    }

    float parasiteAdaptedness = 0.0f;
    //For each antigen apply the forumlat and sum
    for (unsigned int a=0; a<ParamManager::instance().get_int("num_phenotypes"); ++a)
        parasiteAdaptedness = (curAntigenFrequencies[a] / (float) antigenTotal) * susceptibility[a];

    return parasiteAdaptedness;
}

void Output::count_individual_antigens(std::vector<unsigned int>& antigenFreqs, std::unordered_map<Antigen, unsigned int>& antigenCounter, const Strain& strain)
{
    for (const Antigen& antigen : strain)
    {
        if (ParamManager::instance().get_bool("output_antigen_frequency"))
            antigenFreqs[get_phenotype_id(antigen)] += 1;
        if (antigenCounter.find(get_phenotype_id(antigen)) == antigenCounter.end()) //if antigen type hasn't already been counted...
        //if (antigenCounter[get_phenotype_id(antigen)].count == 0)
            antigenCounter[get_phenotype_id(antigen)] = 1;
    }
}

float Output::calc_shannon_entropy(const std::vector<unsigned int>& curAntigenFrequency)//const std::unordered_map<Antigen, const unsigned int>& curAntigenFrequency)
{
    std::vector<float> proportions;
    float total = 0.0;
    for (const auto& item : curAntigenFrequency)
    {
        proportions.push_back(item);
        total += item;
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

