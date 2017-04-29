#include "output.hpp"
#include "model_driver.hpp"
#include "strain.hpp"
#include "utilities.hpp"
#include <cmath>
#include <sstream>
#include <numeric>
#include <unordered_map>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>


void Output::preinitialise_output_storage()
{
    curNumInfectiousBites = 0;
    cumulativeOutputCount = 0;
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

    if (ParamManager::instance().get_bool("output_host_susceptibility"))
        hostSusceptibility.reserve(sizeNeeded);

//    if (ParamManager::instance().get_bool("output_parasite_adaptedness"))
//        parasiteAdaptedness.reserve(sizeNeeded);

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

    process_strain_structure_output(hosts, mosquitoes);

    lastUpdateTime = timestep;
    ++cumulativeOutputCount;
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

    if (ParamManager::instance().get_bool("output_host_susceptibility"))
        utilities::arrayToFile(hostSusceptibility, filePath+runName+"_host_susceptibility.csv");

//    if (ParamManager::instance().get_bool("output_parasite_adaptedness"))
//        utilities::arrayToFile(parasiteAdaptedness, filePath+runName+"_parasite_adaptedness.csv");

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

    #pragma omp parallel for reduction(+:prevalence, multiplicityOfInfection, absImmunity)
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
    #pragma omp parallel for reduction(+:prevalence)
    for (unsigned int i=0; i<mosquitoes.size(); ++i)
    {
        //Count number of mosquitoes that are infected
        if (mosquitoes[i].is_active())
            prevalence += mosquitoes[i].is_infected();
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

    //All host infections
    unsigned int uniqueAntigenCount = 0;
    unsigned int antigenTotal = 0;
    //#pragma omp parallel for reduction (+:antigenTotal, uniqueAntigenCount) //can't really be parallelised effectively because each must read and modify curAntigenFrequencies
    //Parallelised count_individual_antigens instead
    for (unsigned int i=0; i<hosts.size(); ++i)
    {
        if (hosts[i].infection1.infected) {
            count_individual_antigens(curAntigenFrequencies, uniqueAntigenCount, hosts[i].infection1.strain);
            antigenTotal += ParamManager::instance().get_int("repertoire_size");
        }
        if (hosts[i].infection2.infected) {
            count_individual_antigens(curAntigenFrequencies, uniqueAntigenCount, hosts[i].infection2.strain);
            antigenTotal += ParamManager::instance().get_int("repertoire_size");
        }
    }


    //All mosquito infections
    unsigned int mosUniqueAntigenCount = 0;
    unsigned int mosAntigenTotal = 0;
    //#pragma omp parallel for reduction (+:mosUniqueAntigenCount, mosAntigenTotal) //can't really be parallelised effectively because each must read and modify curAntigenFrequencies
    //Parallelised count_individual_antigens instead
    for (unsigned int i=0; i<mosquitoes.size(); ++i)
    {
        if (mosquitoes[i].is_active() && mosquitoes[i].infection.infected) {
            count_individual_antigens(curAntigenFrequencies, mosUniqueAntigenCount, mosquitoes[i].infection.strain);
            mosAntigenTotal += 60;
        }
    }

    uniqueAntigenCount += mosUniqueAntigenCount;
    antigenTotal += mosAntigenTotal;

    //Use frequency to calculate some outputs
    //std::cout << uniqueAntigenCount << ", " << ((float)uniqueAntigenCount) / ((float)ParamManager::instance().get_int("num_phenotypes")) << "\n";
    proportionCirculatingAntigens.push_back(((float)uniqueAntigenCount) / ((float)ParamManager::instance().get_int("num_phenotypes")));
    shannonEntropy.push_back(Output::calc_shannon_entropy(curAntigenFrequencies, antigenTotal));

    //Calculate parasite adaptability
//    if (ParamManager::instance().get_bool("output_parasite_adaptedness"))// && ParamManager::instance().get_bool("output_antigen_frequency"))
//        parasiteAdaptedness.push_back(calc_parasite_adaptedness(curAntigenFrequencies, antigenTotal, hosts));

    //Calculate antigen proportions by normalising by total
    if (ParamManager::instance().get_bool("output_antigen_frequency")) //Turns out we need to do this for shannon entropy anyway!
        antigenFrequency.push_back(curAntigenFrequencies);

    //No need to calculate antigen proportions from antigen frequencies?
    if (ParamManager::instance().get_bool("output_host_susceptibility"))
    {
        //std::vector<float> curAntigenProportions;
        //curAntigenProportions.reserve(curAntigenFrequencies.size());
        //for (unsigned int i=0; i<curAntigenFrequencies.size(); ++i)
        //    curAntigenProportions.push_back((float) curAntigenFrequencies[i] / (float) antigenTotal);
        hostSusceptibility.push_back(calc_host_susceptibility(curAntigenFrequencies, antigenTotal, hosts));
    }
}

//TODO: remove antigenFreqs, and calculate it (if neccessary - we can probably just assume 0 if antigenCounter.find(a) == antigenCounter.end() - recal the vector once at the end).
//TODO: reconsider parallelisation
void Output::count_individual_antigens(std::vector<unsigned int>& antigenFreqs, unsigned int& antigenCounter, const Strain& strain)
{
    unsigned int antigenCounterChange = 0;

    #pragma omp parallel for reduction(+:antigenCounterChange) //probably not very efficient because of read/writes to shared antigenFrequs
    for (unsigned int i=0; i<strain.size(); ++i)
    //for (const Antigen& antigen : strain)
    {
        //if antigen type hasn't already been counted... No need to increment because we're just marking whether or not it exists
        if (antigenFreqs[get_phenotype_id(strain[i])] == 0)
            ++antigenCounterChange;
        //if (antigenFreqs[antigen] == 0)
        //    ++antigenCounter;

        //if (ParamManager::instance().get_bool("output_antigen_frequency"))
        antigenFreqs[get_phenotype_id(strain[i])] += 1;
    }

    antigenCounter += antigenCounterChange;
}


float Output::calc_shannon_entropy(const std::vector<unsigned int>& curAntigenFrequency, const unsigned int totalAntigens)//const std::unordered_map<Antigen, const unsigned int>& curAntigenFrequency)
{
    if (totalAntigens == 0)
        return 0.0;

    std::vector<float> proportions;
    for (unsigned int i=0; i<curAntigenFrequency.size(); ++i) {
        proportions.push_back((float)curAntigenFrequency[i]/(float)totalAntigens);
    }

    float entropy = 0.0;
    #pragma omp parallel for reduction(+:entropy)
    for (unsigned int i=0; i<proportions.size(); ++i)
    {
        if (proportions[i] != 0)
            entropy -= proportions[i] * std::log(proportions[i]);
    }
    //std::cout << "Entropy = " << entropy << "\n";

    return entropy;
}

//Measure of how susceptible the host population is (ranges between 0 and 1). I.e. take away from 1 to give host adaptedness.
float Output::calc_host_susceptibility(const std::vector<unsigned int>& curAntigenFrequencies, const unsigned int antigenTotal, const Hosts& hosts)
{
    if (antigenTotal == 0)
        return 0;

    //Calculate immunity to each antigen
    std::vector<float> immunity(ParamManager::instance().get_int("num_phenotypes"), 0.0f);

    #pragma omp parallel for
    for (unsigned int a=0; a<ParamManager::instance().get_int("num_phenotypes"); ++a)
    {
        for (unsigned int h=0; h<hosts.size(); ++h)
            immunity[a] += 1.0 - hosts[h].immuneState[a];
        immunity[a] = immunity[a] / (float)hosts.size();
        //if (immunity[a] > 1.0)
        //    std::cout << "WARNING: immunity[a] = " << immunity[a] << "\n";
    }


    //Calculate host susceptibility
    float hostSusceptibility = 0.0f;
    #pragma omp parallel for reduction (+:hostSusceptibility)
    for (unsigned int a=0; a<ParamManager::instance().get_int("num_phenotypes"); ++a) {
        //if (((float)curAntigenFrequencies[a] / (float) antigenTotal) > 1.0)
        //    std::cout << "WARNING: p(a) = " << ((float)curAntigenFrequencies[a] / (float) antigenTotal) << "\n";
        hostSusceptibility += ((float)curAntigenFrequencies[a] / (float) antigenTotal) * immunity[a];
    }

    //std::cout << "Host susceptibility = " << hostSusceptibility << "\n";
    return hostSusceptibility;
}


//Counts and outputs the frequency of strain repertoires.
//TODO: optimise with omp. http://stackoverflow.com/questions/15855609/openmpc-c-efficient-way-of-sharing-an-unordered-mapstring-vectorint-an for hints on sharing the unordered_map
void Output::process_strain_structure_output(const Hosts& hosts, const Mosquitoes& mosquitoes)
{
    if (ParamManager::instance().get_bool("output_strain_structure") == false)
        return;

    std::unordered_map<std::string, unsigned int> strainFrequencies;

    //Need to count strains in host and mosquito infections
    for (unsigned int iH=0; iH<hosts.size(); ++iH)
    {
        if (hosts[iH].infection1.infected)
            strainFrequencies[strain_phenotype_str_ordered(hosts[iH].infection1.strain)] += 1;

        if (hosts[iH].infection2.infected)
            strainFrequencies[strain_phenotype_str_ordered(hosts[iH].infection2.strain)] += 1;
    }

    for (unsigned int iM=0; iM<mosquitoes.size(); ++iM)
    {
        if (mosquitoes[iM].is_active() && mosquitoes[iM].infection.infected)
            strainFrequencies[strain_phenotype_str_ordered(mosquitoes[iM].infection.strain)] += 1;
    }

    //Output to file
    int status = mkdir("strain_structure", 0755); //https://linux.die.net/man/2/mkdir

    //create file name
    std::ostringstream outputFilename;
    outputFilename << ParamManager::instance().file_path() << "strain_structure/" << ParamManager::instance().run_name() << "_strain_structure" << cumulativeOutputCount << ".csv";


    std::ofstream file;
    file.open(outputFilename.str(), std::ofstream::out | std::ofstream::trunc);
    for (auto entry : strainFrequencies)
        file << entry.first << ", " << entry.second << "\n";

    file.flush();
    file.close();
}

