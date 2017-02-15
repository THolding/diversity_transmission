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
    std::vector<unsigned int> timeLog;

    std::vector<float> hPrevalence;
    std::vector<float> hImmunityMean;
    std::vector<float> hImmunityVariance;
    std::vector<float> hMultipliticyOfInfection;

    std::vector<float> mPrevalence;

    std::vector<float> circulatingAntigenDiversity;
    std::vector<float> shannonEntropy; //shannon entropy of antigens in circulation
    //std::vector<float> antigenicHostDiversity; //antigen count in circulation

    std::vector<float> eir;

    std::vector<std::vector<unsigned int>> antigenFrequency;

    std::vector<unsigned int> numMosquitoesList;
    std::vector<float> biteRateList;

    void calc_host_infection_metrics(const Hosts& hosts, float& prevalence, float& multiplicityOfInfection);
    void calc_mosquito_infection_metrics(const Mosquitoes& mosquitoes, float& prevalence);

    void calc_host_immunity_metrics(const Hosts& hosts, float& immunityMean, float& immunityVariance);

    float calc_circulating_antigen_diversity(const Hosts& hosts, const Mosquitoes& mosquitoes, float& shannonEntropy);
    //float calc_antigen_host_diversity(const Hosts& hosts);

    float calc_eir(const unsigned int currentTime);
    float calc_shannon_entropy(const std::unordered_map<Antigen, unsigned int>& diversityPool);

    void antigen_counter_helper(std::vector<unsigned int>& antigenFreqs, const Strain& strain);
    void log_antigen_frequency(const Hosts& hosts, const Mosquitoes& mosquitos);
    void log_dyn_params();

public:
    Output(ModelDriver* _model) : model(_model) {  }
    void preinitialise_output_storage();
    void append_output(const unsigned int timestep, const Hosts& hosts, const Mosquitoes& mosquitoes);
    void export_output(const std::string runName=ParamManager::instance().run_name(), const std::string filePath=ParamManager::instance().file_path());

    void register_infectious_bite();
};
