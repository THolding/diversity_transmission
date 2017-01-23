#pragma once
#include "host.hpp"
#include "output.hpp"
#include "mosquito_manager.hpp"

class ModelDriver
{
private:
    PTABLE pDeathHosts;
    PTABLE pDeathMosquitoes;
    PTABLE cdfHosts;
    PTABLE cdfMosquitoes;

    std::vector<Host> hosts;
    std::vector<Mosquito> mosquitoes;
    MosquitoManager mManager;
    float mosChangeRemainder = 0.0;

    Output output;
    int burnInPeriod;

    void age_hosts();
    void age_mosquitoes();
    void update_host_infections();
    void update_mosquito_infections();
    void feed_mosquitoes();
    void update_parameters(const unsigned int time);
public:
    ModelDriver() : output(Output(this)) {  }
    void initialise_model();
    void run_model();
    MosquitoManager* get_mos_manager() {  return &mManager; }
};
