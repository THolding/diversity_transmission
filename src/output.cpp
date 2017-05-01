#include "output.hpp"
#include "diversity_monitor.hpp"
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

//tmp
#include "testing.hpp"


void Output::preinitialise_output_storage()
{
    curNumInfectiousBites = 0;
    cumulativeOutputCount = 0;
    lastUpdateTime = -1;

    //unsigned int sizeNeeded = (numTimeSteps / outputInterval)+1; //+1 for initial conditions
    unsigned int sizeNeeded = ParamManager::output_size_needed;


    timeLog.reserve(sizeNeeded);
    hPrevalence.reserve(sizeNeeded);
    mPrevalence.reserve(sizeNeeded);
    moi.reserve(sizeNeeded);
    eir.reserve(sizeNeeded);
    proportionCirculatingAntigens.reserve(sizeNeeded);
    shannonEntropy.reserve(sizeNeeded);
    absoluteImmunity.reserve(sizeNeeded);
    antigenGenerationRate.reserve(sizeNeeded);
    antigenLossRate.reserve(sizeNeeded);

    if (ParamManager::output_antigen_frequency)
        antigenFrequency.reserve(sizeNeeded);

    if (ParamManager::output_host_susceptibility)
        hostSusceptibility.reserve(sizeNeeded);

//    if (ParamManager::instance().get_bool("output_parasite_adaptedness"))
//        parasiteAdaptedness.reserve(sizeNeeded);

    if (ParamManager::dyn_num_mosquitoes)
        numMosquitoesList.reserve(sizeNeeded);

    if (ParamManager::dyn_bite_rate)
        biteRateList.reserve(sizeNeeded);

    if (ParamManager::dyn_intragenic_recombination_p)
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

    DiversityMonitor::reset_loss_gen_count();
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
    utilities::arrayToFile(antigenGenerationRate, filePath+runName+"_antigen_generation_rate.csv");
    utilities::arrayToFile(antigenLossRate, filePath+runName+"_antigen_loss_rate.csv");

    if (ParamManager::output_antigen_frequency)
        utilities::matrixToFile(antigenFrequency, filePath+runName+"_circulating_antigen_frequency.csv", ", ");

    if (ParamManager::output_host_susceptibility)
        utilities::arrayToFile(hostSusceptibility, filePath+runName+"_host_susceptibility.csv");

//    if (ParamManager::instance().get_bool("output_parasite_adaptedness"))
//        utilities::arrayToFile(parasiteAdaptedness, filePath+runName+"_parasite_adaptedness.csv");

    if (ParamManager::dyn_num_mosquitoes)
        utilities::arrayToFile(numMosquitoesList, filePath+runName+"_num_mosquitoes.csv");

    if (ParamManager::dyn_bite_rate)
        utilities::arrayToFile(biteRateList, filePath+runName+"_bite_rate.csv");

    if (ParamManager::dyn_intragenic_recombination_p)
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

    prevalence = prevalence / (float) ParamManager::num_hosts;
    multiplicityOfInfection = multiplicityOfInfection / (float) ParamManager::num_hosts;
    absImmunity = absImmunity / (float) ParamManager::num_hosts;

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
    float curEir = (float)curNumInfectiousBites / (float)ParamManager::num_hosts;
    curEir = curEir / (float)timeSinceLastUpdate;
    //Reset infectious bite counteer
    curNumInfectiousBites = 0;
    eir.push_back(curEir);

    //Calculate rate of new antigen generation and loss
    antigenGenerationRate.push_back(((float)DiversityMonitor::get_current_generation_count()) / (float)timeSinceLastUpdate);
    antigenLossRate.push_back(((float)DiversityMonitor::get_current_loss_count()) / (float)timeSinceLastUpdate);

    //std::cout << "\ngen: " << antigenGenerationRate.back() << "\n";
    //std::cout << "loss: " << antigenLossRate.back() << "\n";
}

//Track any dynamic (time dependent) parameters over time.
void Output::calc_dyn_metrics()
{
    if (ParamManager::dyn_num_mosquitoes)
        numMosquitoesList.push_back(model->get_mos_manager()->get_count());

    if (ParamManager::dyn_bite_rate)
        biteRateList.push_back(ParamManager::bite_rate);

    if (ParamManager::dyn_intragenic_recombination_p)
        intragenicRecombinationPList.push_back(ParamManager::intragenic_recombination_p);
}


//numCirculatingAntigens, shannon entropy, antigen frequency, parasite adaptedness
//Optional: antigen frequency, parasite adaptedness
void Output::calc_host_mosquito_dependent_metrics(const Hosts& hosts, const Mosquitoes& mosquitoes)
{
    //Use frequency to calculate some outputs
    //std::cout << uniqueAntigenCount << ", " << ((float)uniqueAntigenCount) / ((float)ParamManager::instance().get_int("num_phenotypes")) << "\n";
    proportionCirculatingAntigens.push_back(((float)DiversityMonitor::get_num_unique_antigens()) / ((float)ParamManager::num_phenotypes));
    shannonEntropy.push_back(Output::calc_shannon_entropy(DiversityMonitor::get_antigen_counts(), DiversityMonitor::get_total_antigens()));

    //Calculate antigen proportions by normalising by total
    if (ParamManager::output_antigen_frequency) //Turns out we need to do this for shannon entropy anyway!
        antigenFrequency.push_back(DiversityMonitor::get_antigen_counts());

    //No need to calculate antigen proportions from antigen frequencies?
    if (ParamManager::output_host_susceptibility)
    {
        hostSusceptibility.push_back(calc_host_susceptibility(DiversityMonitor::get_antigen_counts(), DiversityMonitor::get_total_antigens(), hosts));
    }

    //unsigned int uniqueCount;
    //unsigned int totalCount;
    //testing::long_diversity_count(uniqueCount, totalCount, hosts, mosquitoes);
    //std::cout << "Long method: " << uniqueCount << "\t" << totalCount << "\n";
    //std::cout << "Quik method: " << DiversityMonitor::get_num_unique_antigens() << "\t" << DiversityMonitor::get_total_antigens() << "\n";
}





float Output::calc_shannon_entropy(const std::vector<unsigned int>& curAntigenFrequency, const unsigned int totalAntigens)//const std::unordered_map<Antigen, const unsigned int>& curAntigenFrequency)
{
    if (totalAntigens == 0)
        return 0.0;

    float entropy = 0.0;
    #pragma omp parallel for reduction(+:entropy)
    for (unsigned int i=0; i<curAntigenFrequency.size(); ++i)
    {
        float proportion = (float)curAntigenFrequency[i]/(float)totalAntigens;
        if (proportion != 0)
            entropy -= proportion * std::log(proportion);
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
    std::vector<float> immunity(ParamManager::num_phenotypes, 0.0f);

    #pragma omp parallel for
    for (unsigned int a=0; a<ParamManager::num_phenotypes; ++a)
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
    for (unsigned int a=0; a<ParamManager::num_phenotypes; ++a) {
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
    if (ParamManager::output_strain_structure == false)
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
    outputFilename << ParamManager::file_path() << "strain_structure/" << ParamManager::run_name() << "_strain_structure" << cumulativeOutputCount << ".csv";


    std::ofstream file;
    file.open(outputFilename.str(), std::ofstream::out | std::ofstream::trunc);
    for (auto entry : strainFrequencies)
        file << entry.first << ", " << entry.second << "\n";

    file.flush();
    file.close();
}

