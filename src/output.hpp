#pragma once
#include "host.hpp"
#include "mosquito.hpp"

class ModelDriver;

class Output
{
private:
    typedef std::vector<Host> Hosts;
    typedef std::vector<Mosquito> Mosquitoes;

    ModelDriver* model;

    //Counters
    int lastUpdateTime = -1;
    unsigned int curNumInfectiousBites;
    unsigned int cumulativeOutputCount;

    std::vector<unsigned int> timeLog;
    std::vector<float> hPrevalence; //Proportion of hosts which are infected
    std::vector<float> mPrevalence; //Proportion of mosquitoes which are infected
    std::vector<float> moi; //Mean number of infections of those infected
    std::vector<float> eir; //EIR (daily)
    std::vector<float> proportionCirculatingAntigens; //Proportion of antigen space represented by the antigens in circulation.
    std::vector<float> shannonEntropy; //shannon entropy of antigens in circulation
    std::vector<float> absoluteImmunity; //Total number of antigens to which
    std::vector<float> antigenGenerationRate; //Rate of new antigen generation (per day)
    std::vector<float> antigenLossRate; //Rate of antigen loss (per day)

    //Optional output
    std::vector<float> hostSusceptibility; //Outputs a number (ranging between 0 and 1) indicating the mean susceptibility of the host popualtion to currently circulating parasite population.
    //std::vector<float> parasiteAdaptedness; //Measure of how adapted the parasite population is to the current host immunity
    std::vector<std::vector<unsigned int>> antigenFrequency; //frequency of each antigen type, at each output time interval

    //Dynamic parameters
    std::vector<unsigned int> numMosquitoesList; //Tracks number of mosquitoes over time
    std::vector<float> biteRateList; //Tracks bite rate over time
    std::vector<float> intragenicRecombinationPList; //Tracks intragenic recombination rate over time

    void calc_host_dependent_metrics(const Hosts& hosts); //host prevalence, host immunity, moi
    void calc_mosquito_dependent_metrics(const Mosquitoes& mosquitoes); //mosquito prevalence
    void calc_host_mosquito_dependent_metrics(const Hosts& hosts, const Mosquitoes& mosquitoes); //antigen diversity, shannon entropy, antigen frequency, parasite adaptedness
    void calc_time_dependent_metrics(const unsigned int currentTime); //time, daily EIR
    void calc_dyn_metrics();

    //Calculates and outputs repertoire frequencies
    void process_strain_structure_output(const Hosts& hosts, const Mosquitoes& mosquitoes);

    void count_individual_antigens(std::vector<unsigned int>& antigenFreqs, unsigned int& antigenCounter, const Strain& strain);
    float calc_host_susceptibility(const std::vector<unsigned int>& curAntigenFrequencies, const unsigned int totalAntigens, const Hosts& hosts);
    float calc_parasite_adaptedness(const std::vector<unsigned int>& curAntigenFrequencies, const unsigned int antigenTotal, const Hosts& hosts);
    static float calc_shannon_entropy(const std::vector<unsigned int>& curAntigenFrequency, const unsigned int totalAntigens);
    //float calc_shannon_entropy(const std::unordered_map<Antigen, unsigned int>& diversityPool);
    //float calc_shannon_entropy(const std::vector<unsigned int>& curAntigenFrequency);
    //void antigen_counter_helper(std::vector<unsigned int>& antigenFreqs, const Strain& strain);
    void log_dyn_params();

public:

    Output(ModelDriver* _model) : model(_model) {  }
    void preinitialise_output_storage();
    void append_output(const unsigned int timestep, const Hosts& hosts, const Mosquitoes& mosquitoes);
    void export_output(const std::string runName=ParamManager::run_name(), const std::string filePath=ParamManager::file_path());
    void register_infectious_bite();
};


