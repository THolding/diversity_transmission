#pragma once
#include "global_parameters.hpp"
#include "host.hpp"
#include "output.hpp"
#include "mosquito.hpp"



class ModelDriver
{
private:
    PTABLE pDeathHosts;
    PTABLE pDeathMosquitoes;
    PTABLE cdfHosts;
    PTABLE cdfMosquitoes;
    BITE_FREQUENCY_TABLE cumulativeBiteFrequencyDistribution;

    std::vector<Host> hosts;
    std::vector<Mosquito> mosquitoes;

    Output output;
    int burnInPeriod;

    void age_hosts();
    void age_mosquitoes();
    void update_host_infections();
    void update_mosquito_infections();
    void feed_mosquitoes();
public:
    void initialise_model();
    void run_model();
};
