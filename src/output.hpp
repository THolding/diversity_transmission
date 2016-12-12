#pragma once
#include "host.hpp"
#include "mosquito.hpp"

class Output
{
private:
    typedef std::vector<Host> Hosts;
    typedef std::vector<Mosquito> Mosquitoes;

    std::vector<unsigned int> timeLog;

    std::vector<float> hPrevalence;
    std::vector<float> hImmunityMean;
    std::vector<float> hImmunityVariance;
    std::vector<float> hMultipliticyOfInfection;

    std::vector<float> mPrevalence;

    std::vector<unsigned int> hostDiversity;
    std::vector<float> antigenicHostDiversity;

    std::vector<float> eir;

    //Counters
    unsigned int curNumInfectiousBites;

    void calc_host_infection_metrics(const Hosts& hosts, float& prevalence, float& multiplicityOfInfection);
    void calc_mosquito_infection_metrics(const Mosquitoes& mosquitoes, float& prevalence);

    void calc_host_immunity_metrics(const Hosts& hosts, float& immunityMean, float& immunityVariance);

    float calc_host_diversity(const Hosts& hosts);
    float calc_antigen_host_diversity(const Hosts& hosts);

    float calc_eir();

public:
    void preinitialise_output_storage(const unsigned int numTimeSteps, const unsigned int outputInterval);
    void append_output(const unsigned int timestep, const Hosts& hosts, const Mosquitoes& mosquitoes);
    void export_output(const std::string runName=RUN_NAME, const std::string filePath=FILE_PATH);

    void register_infectious_bite();
};
